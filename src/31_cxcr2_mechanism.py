"""CXCR2 mechanism MR: expression -> blood-cell traits -> IBD.

The goal is not to force a neat mediation story. It tests whether the robust
CXCR2 protective cis-MR signal can be explained by circulating neutrophil
count/percentage using GWAS Catalog blood-cell trait GWAS.
"""
from __future__ import annotations

import gzip
import math
import re
import shlex
import subprocess
import time
import urllib.parse
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from pydantic import BaseModel, ConfigDict
from scipy import stats
from scipy.special import logsumexp

from paths import P


PMIDS = {"Vuckovic_2020": "32888494", "Astle_2016": "27863252"}
TRAIT_TARGETS = {
    "neutrophil_count": ("Neutrophil count",),
    "neutrophil_percent": ("Neutrophil percentage of white cells",),
    "lymphocyte_count": ("Lymphocyte count",),
    "lymphocyte_percent": ("Lymphocyte percentage of white cells",),
    "monocyte_count": ("Monocyte count",),
    "monocyte_percent": ("Monocyte percentage of white cells",),
    "eosinophil_count": ("Eosinophil count",),
    "eosinophil_percent": ("Eosinophil percentage of white cells",),
    "white_blood_cell_count": ("White blood cell count",),
}
CORE_TRAITS = {"neutrophil_count", "neutrophil_percent"}
DOWNLOAD_TRAITS = CORE_TRAITS | {"white_blood_cell_count"}
OUTCOMES = {"IBD": "IBD.h.tsv.gz", "CD": "CD.h.tsv.gz", "UC": "UC.h.tsv.gz"}
W_EQTL = 0.15**2
W_GWAS = 0.2**2
P1, P2, P12 = 1e-4, 1e-4, 1e-5


class Config(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    out_dir: Path = P.out("cxcr2_mechanism")
    bloodtrait_dir: Path = P.raw / "gwas" / "bloodtraits"
    gwas_dir: Path = P.raw / "gwas"
    n_jobs: int = 32
    p_instrument: float = 5e-8
    clump_bp: int = 500_000


CFG = Config()


def run_text(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def fetch_json(url: str, tries: int = 8) -> dict:
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return __import__("json").load(r)
        except Exception:
            time.sleep(min(180, 10 * 2**attempt))
    raise RuntimeError(f"failed JSON request after retries: {url}")


def url_exists(url: str) -> bool:
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=45) as r:
            return 200 <= r.status < 400
    except Exception:
        return False


def ftp_group(accession: str) -> str:
    n = int(accession[4:])
    lo = (n - 1) // 1000 * 1000 + 1
    hi = lo + 999
    width = len(accession) - 4
    return f"GCST{lo:0{width}d}-GCST{hi:0{width}d}"


def discover_studies() -> pd.DataFrame:
    cache = CFG.out_dir / "gwas_catalog_bloodtrait_accessions.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    rows = []
    for publication, pmid in PMIDS.items():
        qs = urllib.parse.urlencode({"pubmedId": pmid, "size": 300})
        url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByPublicationIdPubmedId?{qs}"
        data = fetch_json(url)
        for s in data.get("_embedded", {}).get("studies", []):
            trait = s.get("diseaseTrait", {}).get("trait", "")
            for key, labels in TRAIT_TARGETS.items():
                if trait in labels:
                    acc = s["accessionId"]
                    rows.append(
                        {
                            "publication": publication,
                            "pmid": pmid,
                            "trait_key": key,
                            "trait": trait,
                            "accession": acc,
                            "sample_size": s.get("initialSampleSize", ""),
                            "full_pvalue_set": s.get("fullPvalueSet", ""),
                            "source": "GWAS Catalog REST",
                        }
                    )
    df = pd.DataFrame(rows).sort_values(["trait_key", "publication"])
    df.to_csv(cache, sep="\t", index=False)
    return df


def list_harmonised_url(accession: str, pmid: str) -> str | None:
    group = ftp_group(accession)
    base = f"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{group}/{accession}/harmonised/"
    try:
        with urllib.request.urlopen(base, timeout=60) as r:
            text = r.read().decode("utf-8", "replace")
    except Exception:
        return None
    names = re.findall(r'href="([^"]+\.h\.tsv\.gz)"', text)
    names = [n for n in names if n.startswith(f"{pmid}-{accession}-")]
    if not names:
        return None
    return base + sorted(names, key=len)[0]


def download(url: str, path: Path, tries: int = 10) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return "already_present"
    tmp = path.with_suffix(path.suffix + ".part")
    last = ""
    for attempt in range(tries):
        try:
            subprocess.run(
                [
                    "curl",
                    "-L",
                    "--fail",
                    "--retry",
                    "10",
                    "--retry-all-errors",
                    "--connect-timeout",
                    "60",
                    "--speed-time",
                    "900",
                    "--speed-limit",
                    "1024",
                    "--continue-at",
                    "-",
                    "--output",
                    str(tmp),
                    url,
                ],
                check=True,
            )
            tmp.replace(path)
            subprocess.run(["gzip", "-t", str(path)], check=True)
            return "downloaded"
        except subprocess.CalledProcessError as exc:
            last = str(exc)
            time.sleep(min(600, 30 * 2**attempt))
    raise RuntimeError(f"download failed: {url}: {last}")


def ensure_trait_files(studies: pd.DataFrame) -> pd.DataFrame:
    preferred = []
    for key in TRAIT_TARGETS:
        sub = studies[studies["trait_key"].eq(key)].copy()
        if sub.empty:
            continue
        sub["rank"] = sub["publication"].map({"Vuckovic_2020": 0, "Astle_2016": 1}).fillna(9)
        preferred.append(sub.sort_values("rank").iloc[0])
    rows = []
    for r in preferred:
        url = list_harmonised_url(r.accession, str(r.pmid))
        path = CFG.bloodtrait_dir / f"{r.trait_key}_{r.accession}.h.tsv.gz"
        status = "optional_not_downloaded"
        if r.trait_key not in DOWNLOAD_TRAITS:
            rows.append({**r.to_dict(), "url": url or "", "path": str(path), "download_status": status, "usable": path.exists() and path.stat().st_size > 0})
            continue
        status = "no_harmonised_url"
        if url:
            try:
                status = download(url, path)
            except Exception as exc:
                status = f"download_failed: {exc}"
        rows.append({**r.to_dict(), "url": url or "", "path": str(path), "download_status": status, "usable": path.exists() and path.stat().st_size > 0})
    df = pd.DataFrame(rows)
    df.to_csv(CFG.out_dir / "bloodtrait_gwas_downloads.tsv", sep="\t", index=False)
    return df


def clean_allele(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.upper().replace({"NA": "", "NAN": ""})


def rsid_lookup(path: Path, rsids: set[str], tag: str) -> pd.DataFrame:
    rsfile = CFG.out_dir / f"_{tag}_rsids.txt"
    rsfile.write_text("\n".join(sorted(rsids)) + "\n")
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    ix = {c: i + 1 for i, c in enumerate(hdr)}
    rs_col = "hm_rsid" if "hm_rsid" in ix else "variant_id"
    beta_col = "hm_beta" if "hm_beta" in ix else "beta"
    ea_col = "hm_effect_allele" if "hm_effect_allele" in ix else "effect_allele"
    oa_col = "hm_other_allele" if "hm_other_allele" in ix else "other_allele"
    chr_col = "hm_chrom" if "hm_chrom" in ix else "chromosome"
    pos_col = "hm_pos" if "hm_pos" in ix else "base_pair_location"
    se_col = "standard_error"
    p_col = "p_value"
    cmd = (
        f"gunzip -c {shlex.quote(str(path))} | awk -F'\\t' -v OFS='\\t' "
        f"-v rs={ix[rs_col]} -v b={ix[beta_col]} -v ea={ix[ea_col]} -v oa={ix[oa_col]} "
        f"-v se={ix[se_col]} -v pv={ix[p_col]} -v c={ix[chr_col]} -v pos={ix[pos_col]} "
        "'NR==FNR{want[$1]=1; next} FNR==1{next} "
        "{n=split($rs,a,\",\"); for(i=1;i<=n;i++){if(a[i] in want){print a[i],$c,$pos,$b,$se,$pv,$ea,$oa; break}}}' "
        f"{shlex.quote(str(rsfile))} -"
    )
    rows = []
    for line in run_text(cmd).splitlines():
        rsid, chrom, pos, beta, se, p, ea, oa = line.split("\t")
        try:
            rows.append(
                {
                    "rsid": rsid,
                    "chrom": str(chrom).replace("chr", ""),
                    "pos": int(float(pos)),
                    "beta_out": float(beta),
                    "se_out": float(se),
                    "p_out": float(p),
                    "ea_out": ea.upper(),
                    "oa_out": oa.upper(),
                }
            )
        except ValueError:
            continue
    return pd.DataFrame(rows).drop_duplicates("rsid")


def align(beta: float, ea: str, oa: str, target_ea: str, target_oa: str) -> float | None:
    ea, oa, target_ea, target_oa = str(ea).upper(), str(oa).upper(), str(target_ea).upper(), str(target_oa).upper()
    if ea == target_ea and oa == target_oa:
        return float(beta)
    if ea == target_oa and oa == target_ea:
        return -float(beta)
    return None


def ivw(beta_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray) -> tuple[float, float, float, float, float]:
    w = 1 / np.square(se_y)
    denom = float(np.sum(w * np.square(beta_x)))
    if denom <= 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    theta = float(np.sum(w * beta_x * beta_y) / denom)
    se = float(np.sqrt(1 / denom))
    p = float(2 * stats.norm.sf(abs(theta / se)))
    q = float(np.sum(w * np.square(beta_y - theta * beta_x)))
    q_p = float(stats.chi2.sf(q, len(beta_x) - 1)) if len(beta_x) > 1 else np.nan
    return theta, se, p, q, q_p


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    v = values[order]
    w = weights[order] / weights.sum()
    return float(v[np.searchsorted(np.cumsum(w), 0.5)])


def mr_methods(d: pd.DataFrame, exposure: str, outcome: str, beta_x: str, beta_y: str, se_y: str, extra: dict) -> list[dict]:
    if d.empty:
        return [{**extra, "exposure": exposure, "outcome": outcome, "method": "IVW", "nsnp": 0, "theta": np.nan, "se": np.nan, "p": np.nan, "OR": np.nan, "q": np.nan, "q_p": np.nan, "status": "no harmonised instruments"}]
    bx, by, sy = d[beta_x].to_numpy(float), d[beta_y].to_numpy(float), d[se_y].to_numpy(float)
    rows = []
    th, se, p, q, q_p = ivw(bx, by, sy)
    rows.append({**extra, "exposure": exposure, "outcome": outcome, "method": "IVW", "nsnp": len(d), "theta": th, "se": se, "p": p, "OR": math.exp(th) if pd.notna(th) else np.nan, "q": q, "q_p": q_p, "status": "ok"})
    if len(d) >= 3:
        ratio = by / bx
        w = np.square(bx / sy)
        wm = weighted_median(ratio, w)
        rows.append({**extra, "exposure": exposure, "outcome": outcome, "method": "weighted_median", "nsnp": len(d), "theta": wm, "se": np.nan, "p": np.nan, "OR": math.exp(wm), "status": "ok"})
    else:
        rows.append({**extra, "exposure": exposure, "outcome": outcome, "method": "weighted_median", "nsnp": len(d), "theta": np.nan, "se": np.nan, "p": np.nan, "OR": np.nan, "status": "need >=3 SNPs"})
    return rows


def cxcr2_instruments() -> pd.DataFrame:
    blood = pd.read_csv(P.tables / "instruments.tsv", sep="\t")
    blood = blood[blood["gene"].eq("CXCR2")].copy()
    blood["context"] = "blood_eQTLGen"
    blood["effect_allele"] = blood["assessed"]
    blood["other_allele"] = blood["other"]
    pos = pd.read_csv(P.raw / "eqtlgen" / "cis_full_candidates.tsv", sep="\t")
    pos = pos[pos["GeneSymbol"].eq("CXCR2")][["SNP", "SNPChr", "SNPPos"]].drop_duplicates("SNP")
    blood = blood.merge(pos, left_on="rsid", right_on="SNP", how="left")
    blood["chrom"] = blood["SNPChr"]
    blood["pos"] = blood["SNPPos"]
    blood = blood[["gene", "context", "rsid", "chrom", "pos", "effect_allele", "other_allele", "eaf", "beta_eqtl", "se_eqtl", "p_eqtl"]]

    immune = pd.read_csv(P.out("causal_module") / "instruments_multicontext.tsv", sep="\t")
    immune = immune[(immune["gene"].eq("CXCR2")) & (immune["context"].eq("neutrophil"))].copy()
    immune["context"] = "neutrophil_BLUEPRINT"
    immune["eaf"] = immune["maf"]
    immune = immune[["gene", "context", "rsid", "chrom", "pos", "effect_allele", "other_allele", "eaf", "beta_eqtl", "se_eqtl", "p_eqtl"]]
    inst = pd.concat([blood, immune], ignore_index=True)
    inst.to_csv(CFG.out_dir / "cxcr2_expression_instruments.tsv", sep="\t", index=False)
    return inst


def scan_instruments(path: Path, trait_key: str) -> pd.DataFrame:
    usecols = [
        "hm_rsid",
        "hm_chrom",
        "hm_pos",
        "hm_other_allele",
        "hm_effect_allele",
        "hm_beta",
        "other_allele",
        "effect_allele",
        "beta",
        "standard_error",
        "p_value",
        "chromosome",
        "base_pair_location",
    ]
    chunks = []
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", usecols=lambda c: c in usecols, dtype=str, chunksize=750_000):
        p = pd.to_numeric(chunk["p_value"], errors="coerce")
        x = chunk.loc[p < CFG.p_instrument].copy()
        if x.empty:
            continue
        x["p_trait"] = p.loc[x.index].to_numpy()
        x["rsid"] = x["hm_rsid"].fillna("").replace({"NA": ""})
        x["beta_trait"] = pd.to_numeric(x["hm_beta"].where(x["hm_beta"] != "NA", x["beta"]), errors="coerce")
        x["se_trait"] = pd.to_numeric(x["standard_error"], errors="coerce")
        x["ea_trait"] = clean_allele(x["hm_effect_allele"].where(x["hm_effect_allele"] != "NA", x["effect_allele"]))
        x["oa_trait"] = clean_allele(x["hm_other_allele"].where(x["hm_other_allele"] != "NA", x["other_allele"]))
        x["chrom"] = x["hm_chrom"].where(x["hm_chrom"].notna() & (x["hm_chrom"] != "NA"), x.get("chromosome", ""))
        x["pos"] = pd.to_numeric(x["hm_pos"].where(x["hm_pos"].notna() & (x["hm_pos"] != "NA"), x.get("base_pair_location", "")), errors="coerce")
        chunks.append(x[["rsid", "chrom", "pos", "ea_trait", "oa_trait", "beta_trait", "se_trait", "p_trait"]])
    if not chunks:
        return pd.DataFrame()
    hits = pd.concat(chunks, ignore_index=True).dropna()
    hits = hits[(hits["rsid"] != "") & (hits["se_trait"] > 0)].sort_values("p_trait").drop_duplicates("rsid")
    keep, chosen = [], {}
    for r in hits.itertuples(index=False):
        chrom, pos = str(r.chrom).replace("chr", ""), int(r.pos)
        if any(abs(pos - old) <= CFG.clump_bp for old in chosen.get(chrom, [])):
            continue
        chosen.setdefault(chrom, []).append(pos)
        keep.append(r._asdict())
    return pd.DataFrame(keep).assign(trait_key=trait_key)


def cxcr2_to_traits(inst: pd.DataFrame, downloads: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, detail = [], []
    for r in downloads[downloads["usable"]].itertuples(index=False):
        path = Path(r.path)
        out = rsid_lookup(path, set(inst["rsid"]), f"cxcr2_{r.trait_key}")
        h = inst.merge(out, on="rsid", how="inner")
        aligned = []
        for x in h.itertuples(index=False):
            b = align(x.beta_out, x.ea_out, x.oa_out, x.effect_allele, x.other_allele)
            if b is None or x.beta_eqtl == 0 or x.se_out <= 0:
                continue
            aligned.append({**x._asdict(), "beta_trait_aligned": b, "trait_key": r.trait_key, "trait": r.trait, "accession": r.accession})
        d = pd.DataFrame(aligned)
        if not d.empty:
            detail.append(d)
        for context, g in (d.groupby("context") if not d.empty else []):
            rows.extend(mr_methods(g, f"CXCR2_expression_{context}", r.trait_key, "beta_eqtl", "beta_trait_aligned", "se_out", {"trait": r.trait, "accession": r.accession, "context": context}))
        if d.empty:
            for context in inst["context"].unique():
                rows.extend(mr_methods(d, f"CXCR2_expression_{context}", r.trait_key, "beta_eqtl", "beta_trait_aligned", "se_out", {"trait": r.trait, "accession": r.accession, "context": context}))
    res = pd.DataFrame(rows)
    det = pd.concat(detail, ignore_index=True) if detail else pd.DataFrame()
    res.to_csv(CFG.out_dir / "cxcr2_to_bloodtraits_mr.tsv", sep="\t", index=False)
    det.to_csv(CFG.out_dir / "cxcr2_to_bloodtraits_instruments.tsv", sep="\t", index=False)
    return res, det


def trait_to_ibd(downloads: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    usable = downloads[downloads["usable"] & downloads["trait_key"].isin(CORE_TRAITS)].copy()
    insts = Parallel(n_jobs=min(CFG.n_jobs, max(1, len(usable))))(
        delayed(scan_instruments)(Path(r.path), r.trait_key) for r in usable.itertuples(index=False)
    )
    inst = pd.concat([x for x in insts if not x.empty], ignore_index=True) if insts else pd.DataFrame()
    inst = inst.merge(usable[["trait_key", "trait", "accession"]], on="trait_key", how="left")
    inst.to_csv(CFG.out_dir / "bloodtrait_instruments.tsv", sep="\t", index=False)
    rows, details = [], []
    for trait_key, x0 in (inst.groupby("trait_key") if not inst.empty else []):
        for outcome, fname in OUTCOMES.items():
            out = rsid_lookup(CFG.gwas_dir / fname, set(x0["rsid"]), f"{trait_key}_{outcome}")
            h = x0.merge(out, on="rsid", how="inner")
            aligned = []
            for r in h.itertuples(index=False):
                b = align(r.beta_out, r.ea_out, r.oa_out, r.ea_trait, r.oa_trait)
                if b is None or r.beta_trait == 0 or r.se_out <= 0:
                    continue
                aligned.append({**r._asdict(), "beta_ibd_aligned": b, "outcome": outcome})
            d = pd.DataFrame(aligned)
            if not d.empty:
                details.append(d)
            extra = {"trait": x0["trait"].iloc[0], "accession": x0["accession"].iloc[0]}
            rows.extend(mr_methods(d, trait_key, outcome, "beta_trait", "beta_ibd_aligned", "se_out", extra))
    res = pd.DataFrame(rows)
    det = pd.concat(details, ignore_index=True) if details else pd.DataFrame()
    res.to_csv(CFG.out_dir / "bloodtrait_to_ibd_mr.tsv", sep="\t", index=False)
    det.to_csv(CFG.out_dir / "bloodtrait_to_ibd_instruments_harmonised.tsv", sep="\t", index=False)
    return res, det


def approx_coloc(exposure: pd.DataFrame, outcome_path: Path, trait_key: str, context: str) -> dict:
    eq = exposure.copy()
    if eq.empty:
        return {"trait_key": trait_key, "context": context, "nsnp_coloc": 0, "PP4": np.nan, "status": "no eqtl variants"}
    out = rsid_lookup(outcome_path, set(eq["rsid"].dropna()), f"coloc_{context}_{trait_key}")
    h = eq.merge(out, on="rsid", how="inner")
    rows = []
    for r in h.itertuples(index=False):
        b = align(r.beta_out, r.ea_out, r.oa_out, r.effect_allele, r.other_allele)
        if b is not None and r.se_eqtl > 0 and r.se_out > 0:
            rows.append({**r._asdict(), "beta_trait_aligned": b})
    d = pd.DataFrame(rows)
    if len(d) < 20:
        return {"trait_key": trait_key, "context": context, "nsnp_coloc": len(d), "PP4": np.nan, "status": "too few overlapping SNPs"}
    z1 = d["beta_eqtl"].to_numpy(float) / d["se_eqtl"].to_numpy(float)
    z2 = d["beta_trait_aligned"].to_numpy(float) / d["se_out"].to_numpy(float)
    v1 = np.square(d["se_eqtl"].to_numpy(float))
    v2 = np.square(d["se_out"].to_numpy(float))
    l1 = 0.5 * (np.log(v1 / (v1 + W_EQTL)) + z1 * z1 * W_EQTL / (v1 + W_EQTL))
    l2 = 0.5 * (np.log(v2 / (v2 + W_GWAS)) + z2 * z2 * W_GWAS / (v2 + W_GWAS))
    n = len(d)
    lh0 = 0.0
    lh1 = math.log(P1) + logsumexp(l1)
    lh2 = math.log(P2) + logsumexp(l2)
    pairs = logsumexp(l1) + logsumexp(l2)
    same = logsumexp(l1 + l2)
    diff = np.log(max(np.exp(pairs - same) - 1.0, 1e-300)) + math.log(P1 * P2)
    lh4 = math.log(P12) + same - math.log(n)
    den = logsumexp([lh0, lh1, lh2, diff, lh4])
    return {"trait_key": trait_key, "context": context, "nsnp_coloc": n, "PP4": float(np.exp(lh4 - den)), "status": "ok"}


def coloc_tables(downloads: pd.DataFrame) -> pd.DataFrame:
    blood = pd.read_csv(P.raw / "eqtlgen" / "cis_full_candidates.tsv", sep="\t")
    blood = blood[blood["GeneSymbol"].eq("CXCR2")].copy()
    blood["rsid"] = blood["SNP"]
    blood["effect_allele"] = blood["AssessedAllele"]
    blood["other_allele"] = blood["OtherAllele"]
    blood["beta_eqtl"] = pd.read_csv(P.tables / "instruments.tsv", sep="\t").query("gene == 'CXCR2'")["beta_eqtl"].iloc[0]
    blood["se_eqtl"] = np.nan
    # Reconstruct approximate beta/se from z for coloc ABF; sign/scale is enough for locus-level support.
    blood["se_eqtl"] = abs(blood["beta_eqtl"].iloc[0] / blood["Zscore"].replace(0, np.nan).astype(float))
    blood["beta_eqtl"] = blood["Zscore"].astype(float) * blood["se_eqtl"]
    blood = blood[["rsid", "effect_allele", "other_allele", "beta_eqtl", "se_eqtl"]].dropna()
    neut = pd.read_csv(P.out("causal_module") / "multicontext_eqtl_pairs.tsv", sep="\t")
    neut = neut[(neut["gene"].eq("CXCR2")) & (neut["context"].eq("neutrophil"))].copy()
    neut = neut[["rsid", "effect_allele", "other_allele", "beta_eqtl", "se_eqtl"]].dropna()
    rows = []
    for r in downloads[downloads["usable"] & downloads["trait_key"].isin(CORE_TRAITS)].itertuples(index=False):
        rows.append(approx_coloc(blood, Path(r.path), r.trait_key, "blood_eQTLGen"))
        rows.append(approx_coloc(neut, Path(r.path), r.trait_key, "neutrophil_BLUEPRINT"))
    df = pd.DataFrame(rows)
    df.to_csv(CFG.out_dir / "cxcr2_to_bloodtraits_coloc.tsv", sep="\t", index=False)
    return df


def mediation_summary(cx: pd.DataFrame, bt: pd.DataFrame, coloc: pd.DataFrame) -> pd.DataFrame:
    direct = pd.read_csv(P.tables / "module_causal_map_multicontext.tsv", sep="\t")
    direct = direct[(direct["gene"].eq("CXCR2")) & (direct["context"].isin(["blood", "neutrophil"]))]
    rows = []
    for trait in sorted(CORE_TRAITS):
        a = cx[(cx["outcome"].eq(trait)) & (cx["method"].eq("IVW"))].copy()
        b = bt[(bt["exposure"].eq(trait)) & (bt["outcome"].eq("IBD")) & (bt["method"].eq("IVW"))].copy()
        for ar in a.itertuples(index=False):
            br = b.iloc[0] if len(b) else pd.Series(dtype=object)
            prod = ar.theta * br.get("theta", np.nan) if pd.notna(ar.theta) and len(b) else np.nan
            pp4 = coloc[(coloc["trait_key"].eq(trait)) & (coloc["context"].eq(ar.context))]["PP4"]
            rows.append(
                {
                    "trait_key": trait,
                    "cxcr2_context": ar.context,
                    "cxcr2_to_trait_theta": ar.theta,
                    "cxcr2_to_trait_p": ar.p,
                    "cxcr2_to_trait_direction": "higher_CXCR2_increases_trait" if pd.notna(ar.theta) and ar.theta > 0 else "higher_CXCR2_decreases_trait" if pd.notna(ar.theta) else "unknown",
                    "trait_to_IBD_theta": br.get("theta", np.nan),
                    "trait_to_IBD_p": br.get("p", np.nan),
                    "trait_to_IBD_direction": "trait_increases_IBD_risk" if len(b) and br.get("theta", np.nan) > 0 else "trait_protective" if len(b) and br.get("theta", np.nan) < 0 else "unknown",
                    "product_theta": prod,
                    "product_direction_vs_cxcr2_protection": "supports_protection" if pd.notna(prod) and prod < 0 else "opposes_protection" if pd.notna(prod) and prod > 0 else "unknown",
                    "cxcr2_trait_coloc_PP4": pp4.iloc[0] if len(pp4) else np.nan,
                    "direct_CXCR2_IBD_OR_blood": direct[direct["context"].eq("blood")]["MR_OR"].iloc[0],
                    "direct_CXCR2_IBD_OR_neutrophil": direct[direct["context"].eq("neutrophil")]["MR_OR"].iloc[0],
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(CFG.out_dir / "mediation_summary.tsv", sep="\t", index=False)
    return df


def plot_summary(med: pd.DataFrame, bt: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    x = med.copy()
    x["label"] = x["trait_key"].str.replace("_", " ") + "\n" + x["cxcr2_context"].str.replace("_", "\n")
    colors = np.where(x["product_theta"] < 0, "#2b8cbe", "#d95f0e")
    axes[0].barh(np.arange(len(x)), x["product_theta"], color=colors)
    axes[0].axvline(0, color="black", lw=0.8)
    axes[0].set_yticks(np.arange(len(x)), x["label"])
    axes[0].set_xlabel("Network product theta")
    axes[0].set_title("CXCR2 -> trait -> IBD")
    y = bt[(bt["method"].eq("IVW")) & (bt["outcome"].eq("IBD"))].copy()
    y["label"] = y["exposure"].str.replace("_", " ")
    axes[1].barh(np.arange(len(y)), y["theta"], color="#756bb1")
    axes[1].axvline(0, color="black", lw=0.8)
    axes[1].set_yticks(np.arange(len(y)), y["label"])
    axes[1].set_xlabel("Trait -> IBD theta")
    axes[1].set_title("Blood-cell trait causal leg")
    fig.tight_layout()
    path = CFG.out_dir / "Fig_cxcr2_mechanism.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def write_summary(studies: pd.DataFrame, downloads: pd.DataFrame, cx: pd.DataFrame, bt: pd.DataFrame, coloc: pd.DataFrame, med: pd.DataFrame, fig: Path) -> None:
    neut = med[med["trait_key"].isin(CORE_TRAITS)].copy()
    resolved = bool((neut["product_direction_vs_cxcr2_protection"].eq("supports_protection") & (neut["cxcr2_trait_coloc_PP4"].fillna(0) >= 0.5)).any())
    lines = [
        "# CXCR2 neutrophil mechanism summary",
        "",
        "## GWAS Catalog accessions",
        "",
    ]
    for r in studies[studies["trait_key"].isin(CORE_TRAITS)].sort_values(["trait_key", "publication"]).itertuples(index=False):
        lines.append(f"- {r.publication}: {r.trait} = `{r.accession}` ({r.sample_size}).")
    key_block = med.to_csv(sep="\t", index=False)
    lines.extend(
        [
            "",
            "## Mechanism read",
            "",
            f"- Direct prior result: higher CXCR2 expression is protective for IBD (blood OR 0.753; neutrophil OR 0.851).",
            f"- Network verdict: {'the circulating neutrophil trait chain is directionally compatible with CXCR2 protection' if resolved else 'the circulating neutrophil trait chain does not cleanly resolve the CXCR2 protective paradox'}.",
            "- Interpret cautiously: these are circulating blood-cell traits; they do not directly measure intestinal neutrophil trafficking, clearance, or resolution programs.",
            "",
            "## Key rows",
            "",
            "```tsv",
            key_block.rstrip("\n"),
            "```",
            "",
            "## Outputs",
            "",
            "- `cxcr2_to_bloodtraits_mr.tsv`",
            "- `bloodtrait_to_ibd_mr.tsv`",
            "- `mediation_summary.tsv`",
            f"- `{fig.name}`",
        ]
    )
    (CFG.out_dir / "SUMMARY.md").write_text("\n".join(lines) + "\n")


def promote_and_journal(fig: Path) -> None:
    P.promote_figure(fig)
    for name in ["cxcr2_to_bloodtraits_mr.tsv", "bloodtrait_to_ibd_mr.tsv", "mediation_summary.tsv"]:
        P.promote_table(CFG.out_dir / name)
    note = P.journal / "docs" / "cxcr2_neutrophil_mechanism_2026-06-23.md"
    note.write_text(
        "# CXCR2 neutrophil mediation mechanism\n\n"
        "Promoted `outputs/cxcr2_mechanism/` artifacts as the current final mechanism read because they directly test the requested chain: "
        "CXCR2 expression -> circulating neutrophil/blood-cell traits -> de Lange IBD/CD/UC. "
        "The analysis is intentionally honest about non-resolution if the direction, significance, or coloc support is weak.\n"
    )
    status = P.journal / "status" / "overnight_autorun_log.md"
    with status.open("a") as fh:
        fh.write(
            "\n- **2026-06-23 — CXCR2 mechanism B2 complete** (`src/31`, `outputs/cxcr2_mechanism/`). "
            "GWAS Catalog accessions confirmed for Vuckovic 2020 and Astle 2016 neutrophil traits; ran CXCR2 expression -> blood-cell trait MR, "
            "blood-cell trait -> de Lange IBD/CD/UC MR, approximate CXCR2-trait coloc, and network mediation read. See `outputs/cxcr2_mechanism/SUMMARY.md`.\n"
        )


def main() -> None:
    CFG.out_dir.mkdir(parents=True, exist_ok=True)
    studies = discover_studies()
    downloads = ensure_trait_files(studies)
    inst = cxcr2_instruments()
    cx, _ = cxcr2_to_traits(inst, downloads)
    bt, _ = trait_to_ibd(downloads)
    coloc = coloc_tables(downloads)
    med = mediation_summary(cx, bt, coloc)
    fig = plot_summary(med, bt)
    write_summary(studies, downloads, cx, bt, coloc, med, fig)
    promote_and_journal(fig)


if __name__ == "__main__":
    main()
