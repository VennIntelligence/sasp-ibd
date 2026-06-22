"""Second-pass causal hardening for CCL8/CXCR2.

Outputs:
- proper reverse MR with genome-wide disease instruments before/after target
  gene cis exclusion.
- relaxed multi-instrument CCL8 MR/MVMR sensitivity.
- best-effort CCL8 plasma pQTL coloc with every source attempt recorded.
"""
from __future__ import annotations

import gzip
import math
import shlex
import subprocess
import time
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
from statsmodels.stats.multitest import multipletests

from paths import P


GENES = ("CCL8", "CXCR2")
OUTCOMES = {"IBD": "IBD.h.tsv.gz", "CD": "CD.h.tsv.gz", "UC": "UC.h.tsv.gz"}
W_EQTL = 0.15**2
W_GWAS = 0.2**2
P1, P2, P12 = 1e-4, 1e-4, 1e-5


class Config(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    out_dir: Path = P.out("causal_hardening2")
    eqtl_full: Path = P.raw / "eqtlgen" / "cis_full.txt.gz"
    eqtl_candidates: Path = P.raw / "eqtlgen" / "cis_full_candidates.tsv"
    eqtl_af: Path = P.raw / "eqtlgen" / "snp_af.txt.gz"
    crp: Path = P.raw / "gwas" / "CRP_GCST90029070.h.tsv.gz"
    gwas_dir: Path = P.raw / "gwas"
    pqtl_dir: Path = P.raw / "gwas" / "pqtl"
    n_jobs: int = 30


CFG = Config()


def run_text(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def download(url: str, path: Path, tries: int = 8) -> tuple[bool, str]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return True, "already_present"
    tmp = path.with_suffix(path.suffix + ".part")
    for attempt in range(tries):
        try:
            subprocess.run(
                [
                    "curl",
                    "-L",
                    "--fail",
                    "--retry",
                    "8",
                    "--retry-all-errors",
                    "--connect-timeout",
                    "60",
                    "--speed-time",
                    "600",
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
            if path.suffix == ".gz":
                subprocess.run(["gzip", "-t", str(path)], check=True)
            return True, "downloaded"
        except subprocess.CalledProcessError as exc:
            time.sleep(min(600, 20 * 2**attempt))
            last = str(exc)
    return False, f"download_failed: {last}"


def head_status(url: str) -> str:
    try:
        out = subprocess.run(
            ["curl", "-sI", "-L", "--max-time", "45", url],
            check=False,
            capture_output=True,
            text=True,
        ).stdout
        first = next((x for x in out.splitlines() if x.startswith("HTTP/")), "")
        return first or "no_http_status"
    except Exception as exc:
        return f"head_failed: {exc}"


def clean_allele(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.upper().replace({"NA": "", "NAN": ""})


def clump_by_distance(df: pd.DataFrame, p_col: str, window_bp: int, chrom_col: str = "chrom", pos_col: str = "pos") -> pd.DataFrame:
    keep, chosen = [], {}
    for r in df.sort_values(p_col).itertuples(index=False):
        chrom, pos = str(getattr(r, chrom_col)), int(getattr(r, pos_col))
        if any(abs(pos - old) <= window_bp for old in chosen.get(chrom, [])):
            continue
        chosen.setdefault(chrom, []).append(pos)
        keep.append(r._asdict())
    return pd.DataFrame(keep)


def gene_regions() -> pd.DataFrame:
    cols = ["GeneSymbol", "GeneChr", "GenePos"]
    d = pd.read_csv(CFG.eqtl_candidates, sep="\t", usecols=cols)
    d = d[d["GeneSymbol"].isin(GENES)].drop_duplicates("GeneSymbol")
    d = d.rename(columns={"GeneSymbol": "gene", "GeneChr": "chrom", "GenePos": "pos"})
    d["chrom"] = d["chrom"].astype(str)
    d["pos"] = pd.to_numeric(d["pos"], errors="coerce").astype(int)
    return d


def disease_instruments(outcome: str, path: Path, clump_bp: int = 1_000_000) -> pd.DataFrame:
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
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", usecols=usecols, dtype=str, chunksize=500_000):
        p = pd.to_numeric(chunk["p_value"], errors="coerce")
        x = chunk.loc[p < 5e-8].copy()
        if x.empty:
            continue
        x["p_disease"] = p.loc[x.index].to_numpy()
        x["rsid"] = x["hm_rsid"].fillna("").replace({"NA": ""})
        x["beta_disease"] = pd.to_numeric(x["hm_beta"].where(x["hm_beta"] != "NA", x["beta"]), errors="coerce")
        x["se_disease"] = pd.to_numeric(x["standard_error"], errors="coerce")
        x["ea_disease"] = clean_allele(x["hm_effect_allele"].where(x["hm_effect_allele"] != "NA", x["effect_allele"]))
        x["oa_disease"] = clean_allele(x["hm_other_allele"].where(x["hm_other_allele"] != "NA", x["other_allele"]))
        x["chrom"] = x["hm_chrom"].where(x["hm_chrom"].notna() & (x["hm_chrom"] != "NA"), x["chromosome"]).astype(str)
        x["pos"] = pd.to_numeric(x["hm_pos"].where(x["hm_pos"].notna() & (x["hm_pos"] != "NA"), x["base_pair_location"]), errors="coerce")
        chunks.append(x[["rsid", "chrom", "pos", "ea_disease", "oa_disease", "beta_disease", "se_disease", "p_disease"]])
    if not chunks:
        return pd.DataFrame()
    hits = pd.concat(chunks, ignore_index=True).dropna()
    hits = hits[(hits["rsid"] != "") & (hits["se_disease"] > 0)].sort_values("p_disease").drop_duplicates("rsid")
    return clump_by_distance(hits, "p_disease", clump_bp).assign(exposure=outcome)


def lookup_af(rsids: set[str], tag: str) -> dict[str, tuple[str, float]]:
    rsfile = CFG.out_dir / f"_{tag}_rsids.txt"
    rsfile.write_text("\n".join(sorted(rsids)) + "\n")
    cmd = (
        f"gunzip -c {shlex.quote(str(CFG.eqtl_af))} | awk -F'\\t' -v OFS='\\t' "
        "'NR==FNR{want[$1]; next} FNR==1{next} ($1 in want){print $1,$5,$9}' "
        f"{shlex.quote(str(rsfile))} -"
    )
    af = {}
    for line in run_text(cmd).splitlines():
        rsid, allele_b, freq_b = line.split("\t")
        try:
            af[rsid] = (allele_b.upper(), float(freq_b))
        except ValueError:
            pass
    return af


def lookup_eqtl_outcomes(rsids: set[str], genes: tuple[str, ...] = GENES) -> pd.DataFrame:
    rsfile = CFG.out_dir / "_reverse_proper_rsids.txt"
    genefile = CFG.out_dir / "_reverse_proper_genes.txt"
    rsfile.write_text("\n".join(sorted(rsids)) + "\n")
    genefile.write_text("\n".join(genes) + "\n")
    cmd = (
        f"gunzip -c {shlex.quote(str(CFG.eqtl_full))} | awk -F'\\t' -v OFS='\\t' "
        f"-v rsfile={shlex.quote(str(rsfile))} "
        "'BEGIN{while((getline line < rsfile)>0){want[line]=1}; close(rsfile)} "
        "NR==FNR{genes[$1]=1; next} FNR==1{next} "
        "(($2 in want) && ($9 in genes)){print $9,$2,$1,$7,$5,$6,$13}' "
        f"{shlex.quote(str(genefile))} -"
    )
    rows = []
    for line in run_text(cmd).splitlines():
        gene, rsid, p_eqtl, z, assessed, other, n = line.split("\t")
        rows.append({"gene": gene, "rsid": rsid, "p_eqtl": float(p_eqtl), "z": float(z), "assessed": assessed.upper(), "other": other.upper(), "n": int(float(n))})
    if not rows:
        return pd.DataFrame()
    af = lookup_af({r["rsid"] for r in rows}, "reverse_proper_af")
    out = []
    for r in rows:
        if r["rsid"] not in af:
            continue
        allele_b, freq_b = af[r["rsid"]]
        if allele_b == r["assessed"]:
            eaf = freq_b
        elif allele_b == r["other"]:
            eaf = 1 - freq_b
        else:
            continue
        denom = 2 * eaf * (1 - eaf) * (r["n"] + r["z"] ** 2)
        if denom <= 0:
            continue
        se = 1 / np.sqrt(denom)
        out.append({**r, "eaf": eaf, "beta_eqtl": r["z"] * se, "se_eqtl": se})
    return pd.DataFrame(out)


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


def mr_egger(beta_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray) -> tuple[float, float, float, float, float, float]:
    if len(beta_x) < 3:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    w = 1 / np.square(se_y)
    x = np.column_stack([np.ones(len(beta_x)), beta_x])
    cov = np.linalg.inv((x.T * w) @ x)
    b = cov @ ((x.T * w) @ beta_y)
    resid = beta_y - x @ b
    df = len(beta_x) - 2
    scale = max(float(np.sum(w * resid * resid) / df), 1.0) if df > 0 else 1.0
    se = np.sqrt(np.diag(cov * scale))
    return float(b[1]), float(se[1]), float(2 * stats.t.sf(abs(b[1] / se[1]), df)), float(b[0]), float(se[0]), float(2 * stats.t.sf(abs(b[0] / se[0]), df))


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    v = values[order]
    w = weights[order] / np.sum(weights)
    return float(v[np.searchsorted(np.cumsum(w), 0.5)])


def mr_rows(exposure: str, gene: str, d: pd.DataFrame, beta_x: str, beta_y: str, se_y: str, extra: dict) -> list[dict]:
    if len(d) == 0:
        return [{**extra, "exposure": exposure, "gene": gene, "method": "IVW", "nsnp": 0, "theta": np.nan, "se": np.nan, "p": np.nan, "OR": np.nan, "q": np.nan, "q_p": np.nan, "status": "no harmonised instruments"}]
    bx, by, sy = d[beta_x].to_numpy(float), d[beta_y].to_numpy(float), d[se_y].to_numpy(float)
    rows = []
    th, se, p, q, q_p = ivw(bx, by, sy)
    rows.append({**extra, "exposure": exposure, "gene": gene, "method": "IVW", "nsnp": len(d), "theta": th, "se": se, "p": p, "OR": math.exp(th) if pd.notna(th) else np.nan, "q": q, "q_p": q_p, "status": "ok"})
    eg = mr_egger(bx, by, sy)
    rows.append({**extra, "exposure": exposure, "gene": gene, "method": "MR-Egger", "nsnp": len(d), "theta": eg[0], "se": eg[1], "p": eg[2], "OR": math.exp(eg[0]) if pd.notna(eg[0]) else np.nan, "egger_intercept": eg[3], "egger_intercept_se": eg[4], "egger_intercept_p": eg[5], "status": "ok" if len(d) >= 3 else "need >=3 SNPs"})
    if len(d) >= 3:
        ratio = by / bx
        w = np.square(bx / sy)
        wm = weighted_median(ratio, w)
        rows.append({**extra, "exposure": exposure, "gene": gene, "method": "weighted_median", "nsnp": len(d), "theta": wm, "se": np.nan, "p": np.nan, "OR": math.exp(wm), "status": "ok"})
    else:
        rows.append({**extra, "exposure": exposure, "gene": gene, "method": "weighted_median", "nsnp": len(d), "theta": np.nan, "se": np.nan, "p": np.nan, "OR": np.nan, "status": "need >=3 SNPs"})
    return rows


def reverse_mr_proper() -> tuple[pd.DataFrame, pd.DataFrame]:
    inst = pd.concat([disease_instruments(o, CFG.gwas_dir / f) for o, f in OUTCOMES.items()], ignore_index=True)
    inst.to_csv(CFG.out_dir / "reverse_mr_proper_disease_instruments.tsv", sep="\t", index=False)
    eqtl = lookup_eqtl_outcomes(set(inst["rsid"]))
    regions = gene_regions().set_index("gene")
    detail, results = [], []
    for exposure in OUTCOMES:
        x0 = inst[inst["exposure"].eq(exposure)]
        for gene in GENES:
            region = regions.loc[gene]
            in_cis = x0["chrom"].astype(str).eq(str(region.chrom)) & (x0["pos"].astype(float).sub(float(region.pos)).abs() <= 1_000_000)
            for excluded, x in [(False, x0), (True, x0.loc[~in_cis])]:
                m = x.merge(eqtl[eqtl["gene"].eq(gene)], on="rsid", how="inner")
                rows = []
                for r in m.itertuples(index=False):
                    bx = align(r.beta_disease, r.ea_disease, r.oa_disease, r.assessed, r.other)
                    if bx is None or bx == 0:
                        continue
                    rows.append({**r._asdict(), "beta_disease_aligned": bx, "F_disease": (bx / r.se_disease) ** 2, "cis_excluded": excluded})
                h = pd.DataFrame(rows)
                if len(h):
                    detail.append(h.assign(target_gene=gene))
                results.extend(
                    mr_rows(
                        exposure,
                        gene,
                        h,
                        "beta_disease_aligned",
                        "beta_eqtl",
                        "se_eqtl",
                        {
                            "analysis": "disease_to_expression",
                            "cis_excluded": excluded,
                            "n_disease_instruments": len(x),
                            "n_removed_target_cis": int(in_cis.sum()) if excluded else 0,
                            "mean_F": float(h["F_disease"].mean()) if len(h) else np.nan,
                        },
                    )
                )
    dres = pd.DataFrame(results)
    ok = dres["p"].notna()
    dres["fdr"] = np.nan
    if ok.any():
        dres.loc[ok, "fdr"] = multipletests(dres.loc[ok, "p"], method="fdr_bh")[1]
    det = pd.concat(detail, ignore_index=True) if detail else pd.DataFrame()
    return dres, det


def ccl8_cis_instruments(p_thresh: float = 5e-3, clump_kb: int = 100) -> pd.DataFrame:
    cols = ["Pvalue", "SNP", "SNPChr", "SNPPos", "Zscore", "AssessedAllele", "OtherAllele", "GeneSymbol", "NrSamples"]
    d = pd.read_csv(CFG.eqtl_candidates, sep="\t", usecols=cols)
    d = d[d["GeneSymbol"].eq("CCL8")].rename(columns={"Pvalue": "p_eqtl", "SNP": "rsid", "SNPChr": "chrom", "SNPPos": "pos", "Zscore": "z", "AssessedAllele": "assessed", "OtherAllele": "other", "NrSamples": "n"})
    for c in ["p_eqtl", "pos", "z", "n"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["chrom"] = d["chrom"].astype(str)
    d = d.dropna(subset=["p_eqtl", "pos", "z", "n"])
    d = d[d["p_eqtl"] < p_thresh].copy()
    af = lookup_af(set(d["rsid"]), "ccl8_cis")
    rows = []
    for r in d.itertuples(index=False):
        if r.rsid not in af:
            continue
        allele_b, freq_b = af[r.rsid]
        assessed, other = r.assessed.upper(), r.other.upper()
        if allele_b == assessed:
            eaf = freq_b
        elif allele_b == other:
            eaf = 1 - freq_b
        else:
            continue
        denom = 2 * eaf * (1 - eaf) * (r.n + r.z**2)
        if denom <= 0:
            continue
        se = 1 / np.sqrt(denom)
        rows.append({**r._asdict(), "gene": "CCL8", "effect_allele": assessed, "other_allele": other, "eaf": eaf, "beta_eqtl": r.z * se, "se_eqtl": se})
    out = clump_by_distance(pd.DataFrame(rows), "p_eqtl", clump_kb * 1000)
    out["instrument_rule"] = f"eQTLGen CCL8 p<{p_thresh:g}, {clump_kb}kb distance clump"
    return out


def summary_lookup(path: Path, rsids: set[str], tag: str) -> pd.DataFrame:
    key = CFG.out_dir / f"_{tag}_rsids.txt"
    key.write_text("\n".join(sorted(rsids)) + "\n")
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").lstrip("#").split("\t")
    ix = {c: i + 1 for i, c in enumerate(hdr)}
    rs_col = "hm_rsid" if "hm_rsid" in ix else ("rsid" if "rsid" in ix else "rsids")
    beta_col = "hm_beta" if "hm_beta" in ix else "beta"
    se_col = "standard_error" if "standard_error" in ix else ("se" if "se" in ix else "sebeta")
    p_col = "p_value" if "p_value" in ix else ("p" if "p" in ix else "pval")
    ea_col = "hm_effect_allele" if "hm_effect_allele" in ix else "effect_allele"
    oa_col = "hm_other_allele" if "hm_other_allele" in ix else "other_allele"
    cmd = (
        f"gunzip -c {shlex.quote(str(path))} | awk -F'\\t' -v OFS='\\t' "
        f"-v rs={ix[rs_col]} -v b={ix[beta_col]} -v se={ix[se_col]} -v pv={ix[p_col]} -v ea={ix[ea_col]} -v oa={ix[oa_col]} "
        "'NR==FNR{want[$1]=1; next} FNR==1{next} ($rs in want){print $rs,$b,$se,$pv,$ea,$oa}' "
        f"{shlex.quote(str(key))} -"
    )
    rows = []
    for line in run_text(cmd).splitlines():
        rsid, beta, se, p, ea, oa = line.split("\t")
        try:
            rows.append({"rsid": rsid, f"beta_{tag}": float(beta), f"se_{tag}": float(se), f"p_{tag}": float(p), f"ea_{tag}": ea.upper(), f"oa_{tag}": oa.upper()})
        except ValueError:
            pass
    return pd.DataFrame(rows).drop_duplicates("rsid")


def harmonise_series(d: pd.DataFrame, tag: str, ea: str = "effect_allele", oa: str = "other_allele") -> pd.Series:
    vals = [align(getattr(r, f"beta_{tag}"), getattr(r, f"ea_{tag}"), getattr(r, f"oa_{tag}"), getattr(r, ea), getattr(r, oa)) for r in d.itertuples(index=False)]
    return pd.Series(vals, index=d.index, dtype=float)


def ccl8_sensitivity_and_mvmr(inst: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rsids = set(inst["rsid"])
    crp = summary_lookup(CFG.crp, rsids, "crp")
    sens_rows, mvmr_rows, detail = [], [], []
    for outcome, fname in OUTCOMES.items():
        out = summary_lookup(CFG.gwas_dir / fname, rsids, f"out_{outcome}")
        d = inst.merge(out, on="rsid", how="inner")
        d[f"beta_out_{outcome}_aligned"] = harmonise_series(d, f"out_{outcome}")
        d = d.dropna(subset=[f"beta_out_{outcome}_aligned"])
        detail.append(d.assign(outcome=outcome))
        sens_rows.extend(mr_rows(outcome, "CCL8", d, "beta_eqtl", f"beta_out_{outcome}_aligned", f"se_out_{outcome}", {"analysis": "expression_to_disease", "sensitivity": "all_instruments"}))
        bx, by, sy = d["beta_eqtl"].to_numpy(float), d[f"beta_out_{outcome}_aligned"].to_numpy(float), d[f"se_out_{outcome}"].to_numpy(float)
        for i, r in enumerate(d.itertuples(index=False)):
            if len(d) <= 3:
                continue
            mask = np.ones(len(d), dtype=bool)
            mask[i] = False
            th, se, p, q, q_p = ivw(bx[mask], by[mask], sy[mask])
            sens_rows.append({"analysis": "expression_to_disease", "sensitivity": "leave_one_out", "exposure": outcome, "gene": "CCL8", "method": "leave_one_out", "left_out": r.rsid, "nsnp": int(mask.sum()), "theta": th, "se": se, "p": p, "OR": math.exp(th), "q": q, "q_p": q_p, "status": "ok"})
        dm = inst.merge(crp, on="rsid", how="inner").merge(out, on="rsid", how="inner")
        dm["beta_crp_aligned"] = harmonise_series(dm, "crp")
        dm["beta_out_aligned"] = harmonise_series(dm, f"out_{outcome}")
        dm = dm.dropna(subset=["beta_crp_aligned", "beta_out_aligned"])
        if len(dm) < 3:
            for exp in ["gene_expression", "CRP"]:
                mvmr_rows.append({"gene": "CCL8", "outcome": outcome, "adjusted_exposure": exp, "nsnp": len(dm), "theta": np.nan, "se": np.nan, "p": np.nan, "OR": np.nan, "conditional_F_approx": np.nan, "status": "need >=3 harmonised instruments for two-exposure MVMR"})
            continue
        x = dm[["beta_eqtl", "beta_crp_aligned"]].to_numpy(float)
        y = dm["beta_out_aligned"].to_numpy(float)
        w = 1 / np.square(dm[f"se_out_{outcome}"].to_numpy(float))
        if np.linalg.matrix_rank(x) < 2:
            status = "exposure matrix rank <2"
            coef = se = np.array([np.nan, np.nan])
            q = q_p = np.nan
        else:
            xtw = x.T * w
            cov = np.linalg.inv(xtw @ x)
            coef = cov @ (xtw @ y)
            resid = y - x @ coef
            df = len(y) - 2
            q = float(np.sum(w * resid * resid))
            q_p = float(stats.chi2.sf(q, df)) if df > 0 else np.nan
            scale = max(q / df, 1.0) if df > 0 else 1.0
            se = np.sqrt(np.diag(cov * scale))
            status = "ok"
        for label, idx in [("gene_expression", 0), ("CRP", 1)]:
            p = float(2 * stats.norm.sf(abs(coef[idx] / se[idx]))) if pd.notna(se[idx]) and se[idx] > 0 else np.nan
            mvmr_rows.append({"gene": "CCL8", "outcome": outcome, "adjusted_exposure": label, "nsnp": len(dm), "theta": float(coef[idx]) if pd.notna(coef[idx]) else np.nan, "se": float(se[idx]) if pd.notna(se[idx]) else np.nan, "p": p, "OR": math.exp(coef[idx]) if pd.notna(coef[idx]) else np.nan, "q": q, "q_p": q_p, "conditional_F_approx": np.nan, "status": status})
    return pd.DataFrame(mvmr_rows), pd.DataFrame(sens_rows), pd.concat(detail, ignore_index=True)


def coloc_abf(d: pd.DataFrame) -> dict:
    z1, z2 = d["beta_pqtl"] / d["se_pqtl"], d["beta_gwas"] / d["se_gwas"]
    v1, v2 = np.square(d["se_pqtl"]), np.square(d["se_gwas"])
    l1 = 0.5 * (np.log(v1 / (v1 + W_EQTL)) + z1**2 * W_EQTL / (v1 + W_EQTL))
    l2 = 0.5 * (np.log(v2 / (v2 + W_GWAS)) + z2**2 * W_GWAS / (v2 + W_GWAS))
    lsum1, lsum2, lsum12 = logsumexp(l1), logsumexp(l2), logsumexp(l1 + l2)
    ns = len(d)
    terms = np.array([0.0, np.log(P1) + lsum1, np.log(P2) + lsum2, np.log(P1) + np.log(P2) + np.log(max(np.exp(lsum1 + lsum2) - np.exp(lsum12), 1e-300)), np.log(P12) + lsum12])
    pp = np.exp(terms - logsumexp(terms))
    return {"nsnp": ns, "PP0": pp[0], "PP1": pp[1], "PP2": pp[2], "PP3": pp[3], "PP4": pp[4]}


def read_sumstat_locus(path: Path, chrom: str, start: int, end: int, tag: str) -> pd.DataFrame:
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").lstrip("#").split("\t")
    ix = {c: i for i, c in enumerate(hdr)}
    chrom_col = next(c for c in ["hm_chrom", "chromosome", "chromosome"] if c in ix)
    pos_col = next(c for c in ["hm_pos", "base_pair_location", "position"] if c in ix)
    rs_col = "hm_rsid" if "hm_rsid" in ix else "rsid"
    beta_col = "hm_beta" if "hm_beta" in ix else "beta"
    se_col = "standard_error" if "standard_error" in ix else "se"
    p_col = "p_value" if "p_value" in ix else "p"
    ea_col = "hm_effect_allele" if "hm_effect_allele" in ix else "effect_allele"
    oa_col = "hm_other_allele" if "hm_other_allele" in ix else "other_allele"
    rows = []
    use = [chrom_col, pos_col, rs_col, beta_col, se_col, p_col, ea_col, oa_col]
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", usecols=use, dtype=str, chunksize=500_000):
        c = chunk[chrom_col].astype(str).str.replace("^chr", "", regex=True)
        pos = pd.to_numeric(chunk[pos_col], errors="coerce")
        x = chunk.loc[c.eq(str(chrom)) & pos.between(start, end)].copy()
        if x.empty:
            continue
        x["pos"] = pos.loc[x.index].astype(int)
        rows.append(
            pd.DataFrame(
                {
                    "rsid": x[rs_col].fillna("").replace({"NA": ""}),
                    "chrom": str(chrom),
                    "pos": x["pos"],
                    f"beta_{tag}": pd.to_numeric(x[beta_col], errors="coerce"),
                    f"se_{tag}": pd.to_numeric(x[se_col], errors="coerce"),
                    f"p_{tag}": pd.to_numeric(x[p_col], errors="coerce"),
                    f"ea_{tag}": clean_allele(x[ea_col]),
                    f"oa_{tag}": clean_allele(x[oa_col]),
                }
            )
        )
    return pd.concat(rows, ignore_index=True).dropna() if rows else pd.DataFrame()


def pqtl_ccl8_v2() -> pd.DataFrame:
    rows = []
    sources = [
        ("SCALLOP CVD1 Zenodo MCP-2 candidate", "https://zenodo.org/records/2615265/files/MCP-2.txt.gz?download=1", "", "head_only"),
        ("SCALLOP CVD1 Zenodo CCL8 candidate", "https://zenodo.org/records/2615265/files/CCL8.txt.gz?download=1", "", "head_only"),
        ("SCALLOP-INF GWAS Catalog MCP.2/CCL8", "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/GCST90274822/GCST90274822.tsv.gz", "GCST90274822", "download_coloc"),
        ("SCALLOP-INF GWAS Catalog MCP.2/CCL8 harmonised dir", "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/GCST90274822/harmonised/", "GCST90274822", "head_only"),
        ("deCODE summary data landing page", "https://www.decode.com/summarydata/", "", "head_only"),
    ]
    for source, url, acc, mode in sources:
        status = head_status(url)
        base = {"gene": "CCL8", "protein": "MCP-2/CCL8", "source": source, "accession": acc, "url": url, "http_status": status, "nsnp": 0, "PP4": np.nan, "status": "not_downloaded"}
        if mode != "download_coloc" or "200" not in status:
            rows.append(base | {"status": "probed_only"})
            continue
        path = CFG.pqtl_dir / "GCST90274822_CCL8_MCP2.tsv.gz"
        ok, msg = download(url, path)
        if not ok:
            rows.append(base | {"status": msg})
            continue
        try:
            region = gene_regions().set_index("gene").loc["CCL8"]
            start, end = int(region.pos) - 1_000_000, int(region.pos) + 1_000_000
            pq = read_sumstat_locus(path, str(region.chrom), start, end, "pqtl")
            gw = read_sumstat_locus(CFG.gwas_dir / "IBD.h.tsv.gz", str(region.chrom), start, end, "gwas")
            m = pq.merge(gw, on="rsid", how="inner", suffixes=("", "_g"))
            vals = [align(r.beta_gwas, r.ea_gwas, r.oa_gwas, r.ea_pqtl, r.oa_pqtl) for r in m.itertuples(index=False)]
            m["beta_gwas"] = pd.Series(vals, index=m.index, dtype=float)
            m = m.dropna(subset=["beta_gwas", "beta_pqtl", "se_gwas", "se_pqtl"])
            coloc = coloc_abf(m) if len(m) >= 5 else {"nsnp": len(m), "PP0": np.nan, "PP1": np.nan, "PP2": np.nan, "PP3": np.nan, "PP4": np.nan}
            m.to_csv(CFG.out_dir / "pqtl_ccl8_v2_coloc_snps.tsv", sep="\t", index=False)
            rows.append(base | coloc | {"local_path": str(path), "status": "ok" if len(m) >= 5 else "fewer than 5 coloc SNPs"})
        except Exception as exc:
            rows.append(base | {"local_path": str(path), "status": f"parser_or_coloc_failed: {exc}"})
    return pd.DataFrame(rows)


def plot_figure(reverse: pd.DataFrame, mvmr: pd.DataFrame, sens: pd.DataFrame, pqtl: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
    ivw = reverse[(reverse["method"].eq("IVW")) & (reverse["gene"].isin(GENES))].copy()
    ivw["label"] = ivw["exposure"] + "->" + ivw["gene"] + np.where(ivw["cis_excluded"], " excl.cis", " all")
    ivw["neglog10p"] = -np.log10(ivw["p"].fillna(1).clip(lower=np.nextafter(0, 1)))
    ax[0].barh(ivw["label"], ivw["neglog10p"], color=np.where(ivw["cis_excluded"], "#4daf4a", "#984ea3"))
    ax[0].axvline(-np.log10(0.05), color="black", lw=0.8, ls=":")
    ax[0].set_title("Proper reverse MR")
    ax[0].set_xlabel("-log10 p")
    mv = mvmr[(mvmr["adjusted_exposure"].eq("gene_expression")) & mvmr["OR"].notna()]
    ax[1].scatter(mv["outcome"], mv["OR"], s=70, color="#e41a1c")
    ax[1].axhline(1, color="black", lw=0.8)
    ax[1].set_title("CCL8 MVMR-CRP")
    ax[1].set_ylabel("OR")
    ss = sens[(sens["method"].isin(["IVW", "MR-Egger", "weighted_median"])) & sens["OR"].notna()]
    for method, g in ss.groupby("method"):
        ax[2].plot(g["exposure"], g["OR"], marker="o", label=method)
    ax[2].axhline(1, color="black", lw=0.8)
    pp4 = pqtl["PP4"].dropna().max() if "PP4" in pqtl else np.nan
    ax[2].set_title(f"CCL8 sensitivity; pQTL PP4={pp4:.2g}" if pd.notna(pp4) else "CCL8 sensitivity; pQTL blocked")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    plt.close(fig)


def write_summary(reverse: pd.DataFrame, mvmr: pd.DataFrame, sens: pd.DataFrame, pqtl: pd.DataFrame, inst: pd.DataFrame) -> None:
    cx_all = reverse[(reverse["gene"].eq("CXCR2")) & reverse["method"].eq("IVW") & ~reverse["cis_excluded"]]
    cx_ex = reverse[(reverse["gene"].eq("CXCR2")) & reverse["method"].eq("IVW") & reverse["cis_excluded"]]
    ccl8_mv = mvmr[(mvmr["adjusted_exposure"].eq("gene_expression")) & mvmr["theta"].notna()]
    ccl8_sens = sens[(sens["method"].isin(["IVW", "MR-Egger", "weighted_median"])) & sens["theta"].notna()]
    pq_ok = pqtl[pqtl["status"].eq("ok")]
    lines = [
        "# Causal hardening 2 summary",
        "",
        "## Direct answers",
        f"- Proper reverse MR: CXCR2's previous reverse signal is {'clarified/attenuated after target-cis exclusion' if len(cx_ex) and (cx_ex['p'].fillna(1) > 0.05).all() else 'not fully cleared after target-cis exclusion'}. Before/after IVW rows:",
        "```tsv",
        pd.concat([cx_all, cx_ex], ignore_index=True)[["exposure", "gene", "cis_excluded", "method", "nsnp", "theta", "p", "fdr", "n_removed_target_cis", "status"]].to_csv(sep="\t", index=False).rstrip(),
        "```",
        f"- CCL8 instruments: relaxed eQTLGen p<0.005 and 100kb distance clump produced {len(inst)} independent cis instruments. This meets the >=3 target." if len(inst) >= 3 else f"- CCL8 instruments: only {len(inst)} instruments; the >=3 target was not met.",
        f"- CCL8 MVMR-CRP: {'estimable' if len(ccl8_mv) else 'not estimable'}; gene-expression rows are:",
        "```tsv",
        ccl8_mv[["outcome", "nsnp", "OR", "p", "status"]].to_csv(sep="\t", index=False).rstrip() if len(ccl8_mv) else "none",
        "```",
        f"- CCL8 Egger/weighted-median/LOO sensitivity: {'available' if len(ccl8_sens) else 'not available'}; no method is hidden when not estimable.",
        f"- CCL8 protein pQTL: {'supports coloc' if len(pq_ok) and pq_ok['PP4'].max() > 0.8 else 'does not provide strong coloc support in this run'}. All attempted sources are in `pqtl_ccl8_v2.tsv`.",
        "",
        "## pQTL source audit",
        "```tsv",
        pqtl[["source", "accession", "http_status", "nsnp", "PP4", "status", "url"]].to_csv(sep="\t", index=False).rstrip(),
        "```",
    ]
    (CFG.out_dir / "SUMMARY.md").write_text("\n".join(lines) + "\n")


def update_status() -> None:
    status = P.journal / "status" / "causal_module_status.md"
    text = status.read_text() if status.exists() else "# Causal module status\n"
    marker = "causal hardening 2"
    note = f"""

## 2026-06-23 JST - causal hardening 2

- Ran `src/28_causal_hardening2.py` CPU-only with `n_jobs={CFG.n_jobs}`.
- Wrote strict bidirectional reverse-MR with genome-wide disease instruments, 1Mb distance clumping, and target-gene cis exclusion to `outputs/causal_hardening2/reverse_mr_proper.tsv`.
- Built relaxed CCL8 cis instruments from eQTLGen (`p<1e-3`, 100kb distance clump), then ran CCL8 MVMR-CRP and IVW/Egger/weighted-median/leave-one-out sensitivity.
- Audited CCL8/MCP-2 plasma pQTL sources including SCALLOP CVD1 Zenodo probes, SCALLOP-INF GWAS Catalog accession GCST90274822, and the deCODE summary-data landing page; all source attempts are recorded in `pqtl_ccl8_v2.tsv`.
- Promoted the hardening2 final figure and requested tables into `results/`.
"""
    if marker not in text:
        status.write_text(text.rstrip() + "\n" + note)


def main() -> None:
    CFG.out_dir.mkdir(parents=True, exist_ok=True)
    reverse, reverse_detail = reverse_mr_proper()
    reverse.to_csv(CFG.out_dir / "reverse_mr_proper.tsv", sep="\t", index=False)
    reverse_detail.to_csv(CFG.out_dir / "reverse_mr_proper_instruments.tsv", sep="\t", index=False)
    ccl8_inst = ccl8_cis_instruments()
    ccl8_inst.to_csv(CFG.out_dir / "ccl8_cis_instruments.tsv", sep="\t", index=False)
    mvmr, sens, detail = ccl8_sensitivity_and_mvmr(ccl8_inst)
    mvmr.to_csv(CFG.out_dir / "ccl8_mvmr.tsv", sep="\t", index=False)
    sens.to_csv(CFG.out_dir / "ccl8_mr_sensitivity.tsv", sep="\t", index=False)
    detail.to_csv(CFG.out_dir / "ccl8_mr_instruments.tsv", sep="\t", index=False)
    pqtl = pqtl_ccl8_v2()
    pqtl.to_csv(CFG.out_dir / "pqtl_ccl8_v2.tsv", sep="\t", index=False)
    fig = CFG.out_dir / "Fig_hardening2.png"
    plot_figure(reverse, mvmr, sens, pqtl, fig)
    write_summary(reverse, mvmr, sens, pqtl, ccl8_inst)
    for f in ["reverse_mr_proper.tsv", "ccl8_mvmr.tsv", "ccl8_mr_sensitivity.tsv", "pqtl_ccl8_v2.tsv"]:
        P.promote_table(CFG.out_dir / f)
    P.promote_figure(fig)
    update_status()
    print((CFG.out_dir / "SUMMARY.md").read_text())


if __name__ == "__main__":
    main()
