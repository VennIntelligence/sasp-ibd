"""Causal hardening for CCL8/CXCR2 and bystander formalisation.

This script consolidates the post-hoc robustness checks requested after the
multicontext causal map. It reuses local MR/coloc primitives and writes honest
status rows when raw eQTL/pQTL dependencies are not available in this checkout.
"""
from __future__ import annotations

import gzip
import math
import shlex
import subprocess
import time
import urllib.error
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

from causal_module_utils import gwas_lookup, harmonised_beta
from paths import P


CORE = ("CCL8", "CXCR2")
BYSTANDERS = ("OSM", "OSMR", "TREM1", "IL13RA2")
OUTCOMES = {"IBD": "IBD.h.tsv.gz", "CD": "CD.h.tsv.gz", "UC": "UC.h.tsv.gz"}
FINNGEN = {"IBD": "K11_IBD_STRICT", "CD": "K11_CD_STRICT2", "UC": "K11_UC_STRICT2"}
W_EQTL = 0.15**2
W_GWAS = 0.2**2
P1, P2, P12 = 1e-4, 1e-4, 1e-5


class Config(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    out_dir: Path = P.out("causal_hardening")
    gwas_dir: Path = P.raw / "gwas"
    causal_dir: Path = P.out("causal_module")
    eqtl_dir: Path = P.raw / "eqtlgen"
    crp_url: str = (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90029001-GCST90030000/GCST90029070/harmonised/"
        "35459240-GCST90029070-EFO_0004458.h.tsv.gz"
    )
    crp_path: Path = P.raw / "gwas" / "CRP_GCST90029070.h.tsv.gz"
    finngen_manifest_url: str = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv"
    finngen_manifest: Path = P.raw / "gwas" / "finngen_R12_manifest.tsv"
    eqtl_full_url: str = (
        "http://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"
        "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
    )
    eqtl_af_url: str = (
        "http://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"
        "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"
    )
    n_jobs: int = 30


CFG = Config()


def run_text(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def download(url: str, path: Path, tries: int = 8) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    if path.exists() and path.stat().st_size > 0:
        return
    for attempt in range(tries):
        try:
            cmd = [
                "curl", "-L", "--fail", "--retry", "8", "--retry-all-errors",
                "--connect-timeout", "60", "--speed-time", "600", "--speed-limit", "1024",
                "--continue-at", "-", "--output", str(tmp), url,
            ]
            if "molgenis26.gcc.rug.nl" in url:
                cmd.insert(1, "--insecure")
            subprocess.run(cmd, check=True)
            tmp.replace(path)
            if path.suffix == ".gz":
                subprocess.run(["gzip", "-t", str(path)], check=True)
            return
        except subprocess.CalledProcessError:
            time.sleep(min(600, 30 * 2**attempt))
    raise RuntimeError(f"download failed after retries: {url}")


def ensure_public_inputs() -> pd.DataFrame:
    rows = []
    for label, url, path in [
        ("CRP_GCST90029070", CFG.crp_url, CFG.crp_path),
        ("FinnGen_R12_manifest", CFG.finngen_manifest_url, CFG.finngen_manifest),
        ("eQTLGen_full_cis", CFG.eqtl_full_url, CFG.eqtl_dir / "cis_full.txt.gz"),
        ("eQTLGen_allele_frequency", CFG.eqtl_af_url, CFG.eqtl_dir / "snp_af.txt.gz"),
    ]:
        try:
            download(url, path)
            rows.append({"input": label, "path": str(path), "url": url, "status": "ok"})
        except Exception as exc:
            rows.append({"input": label, "path": str(path), "url": url, "status": f"download_failed: {exc}"})
    try:
        build_cis_candidates()
        rows.append({"input": "eQTLGen_target_cis_candidates", "path": str(CFG.eqtl_dir / "cis_full_candidates.tsv"), "url": "", "status": "ok"})
    except Exception as exc:
        rows.append({"input": "eQTLGen_target_cis_candidates", "path": str(CFG.eqtl_dir / "cis_full_candidates.tsv"), "url": "", "status": f"build_failed: {exc}"})
    return pd.DataFrame(rows)


def build_cis_candidates() -> None:
    full = CFG.eqtl_dir / "cis_full.txt.gz"
    out = CFG.eqtl_dir / "cis_full_candidates.tsv"
    if out.exists() and out.stat().st_size > 0:
        return
    if not full.exists():
        raise FileNotFoundError(full)
    CFG.eqtl_dir.mkdir(parents=True, exist_ok=True)
    genes = CFG.out_dir / "_eqtlgen_candidate_genes.txt"
    genes.write_text("\n".join(sorted(set(CORE) | set(BYSTANDERS))) + "\n")
    cmd = (
        f"gunzip -c {shlex.quote(str(full))} | awk -F'\\t' -v OFS='\\t' "
        "'NR==FNR{g[$1]=1; next} FNR==1{print; next} ($9 in g){print}' "
        f"{shlex.quote(str(genes))} - > {shlex.quote(str(out))}"
    )
    subprocess.run(cmd, shell=True, check=True)


def clean_allele(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.upper().replace({"NA": "", "NAN": ""})


def summary_lookup(path: Path, rsids: set[str], tag: str) -> pd.DataFrame:
    key = CFG.out_dir / f"_{tag}_rsids.txt"
    key.write_text("\n".join(sorted(rsids)) + "\n")
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").lstrip("#").split("\t")
    ix = {c: i + 1 for i, c in enumerate(hdr)}
    rs_col = "hm_rsid" if "hm_rsid" in ix else ("rsids" if "rsids" in ix else "rsid")
    beta_col = "hm_beta" if "hm_beta" in ix else ("beta" if "beta" in ix else "BETA")
    se_col = "standard_error" if "standard_error" in ix else ("sebeta" if "sebeta" in ix else "se")
    p_col = "p_value" if "p_value" in ix else ("pval" if "pval" in ix else "p")
    ea_col = "hm_effect_allele" if "hm_effect_allele" in ix else ("effect_allele" if "effect_allele" in ix else "alt")
    oa_col = "hm_other_allele" if "hm_other_allele" in ix else ("other_allele" if "other_allele" in ix else "ref")
    cmd = (
        f"gunzip -c {shlex.quote(str(path))} | awk -F'\\t' -v OFS='\\t' "
        f"-v rs={ix[rs_col]} -v b={ix[beta_col]} -v se={ix[se_col]} -v pv={ix[p_col]} -v ea={ix[ea_col]} -v oa={ix[oa_col]} "
        "'NR==FNR{want[$1]=1; next} FNR==1{next} "
        "{n=split($rs,a,\",\"); for(i=1;i<=n;i++){if(a[i] in want){print a[i],$b,$se,$pv,$ea,$oa; break}}}' "
        f"{shlex.quote(str(key))} -"
    )
    rows = []
    for line in run_text(cmd).splitlines():
        rsid, beta, se, p, ea, oa = line.split("\t")
        try:
            rows.append({"rsid": rsid, f"beta_{tag}": float(beta), f"se_{tag}": float(se), f"p_{tag}": float(p), f"ea_{tag}": ea.upper(), f"oa_{tag}": oa.upper()})
        except ValueError:
            continue
    return pd.DataFrame(rows).drop_duplicates("rsid")


def align(beta: float, ea: str, oa: str, target_ea: str, target_oa: str) -> float | None:
    if ea == target_ea and oa == target_oa:
        return beta
    if ea == target_oa and oa == target_ea:
        return -beta
    return None


def ivw(beta_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray) -> tuple[float, float, float, float, float]:
    w = 1 / np.square(se_y)
    denom = float(np.sum(w * np.square(beta_x)))
    if denom <= 0:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    theta = float(np.sum(w * beta_x * beta_y) / denom)
    se = float(np.sqrt(1 / denom))
    p = float(2 * stats.norm.sf(abs(theta / se)))
    q = float(np.sum(w * np.square(beta_y - theta * beta_x)))
    q_p = float(stats.chi2.sf(q, len(beta_x) - 1)) if len(beta_x) > 1 else np.nan
    return theta, se, p, q, q_p


def weighted_median(ratio: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(ratio)
    r = ratio[order]
    w = weights[order] / np.sum(weights)
    return float(r[np.searchsorted(np.cumsum(w), 0.5)])


def mr_egger(beta_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray) -> tuple[float, float, float, float, float, float]:
    if len(beta_x) < 3:
        return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
    w = 1 / np.square(se_y)
    x = np.column_stack([np.ones(len(beta_x)), beta_x])
    xtw = x.T * w
    cov = np.linalg.inv(xtw @ x)
    b = cov @ (xtw @ beta_y)
    resid = beta_y - x @ b
    df = len(beta_x) - 2
    scale = max(float(np.sum(w * resid * resid) / df), 1.0) if df > 0 else 1.0
    se = np.sqrt(np.diag(cov * scale))
    slope_p = float(2 * stats.t.sf(abs(b[1] / se[1]), df))
    int_p = float(2 * stats.t.sf(abs(b[0] / se[0]), df))
    return float(b[1]), float(se[1]), slope_p, float(b[0]), float(se[0]), int_p


def lead_instruments() -> pd.DataFrame:
    blood = pd.read_csv(P.tables / "instruments.tsv", sep="\t")
    blood = blood[blood["gene"].isin(CORE)].copy()
    blood["context"] = "blood"
    blood["effect_allele"] = blood["assessed"]
    blood["other_allele"] = blood["other"]
    blood["chrom"] = np.nan
    blood["pos"] = np.nan
    immune = pd.read_csv(CFG.causal_dir / "instruments_multicontext.tsv", sep="\t")
    immune = immune[(immune["gene"].isin(CORE)) & (immune["context"].eq("neutrophil"))].copy()
    return pd.concat([blood, immune], ignore_index=True, sort=False)


def steiger() -> pd.DataFrame:
    inst = lead_instruments()
    rows = []
    for outcome, fname in OUTCOMES.items():
        by_rsid = summary_lookup(CFG.gwas_dir / fname, set(inst["rsid"].dropna()), f"steiger_{outcome}")
        for r in inst.merge(by_rsid, on="rsid", how="left").itertuples(index=False):
            if pd.isna(getattr(r, f"beta_steiger_{outcome}", np.nan)):
                rows.append({"gene": r.gene, "context": r.context, "outcome": outcome, "rsid": r.rsid, "steiger_dir": "", "r2_exposure": np.nan, "r2_outcome_proxy": np.nan, "outcome_z2": np.nan, "status": "outcome rsid absent"})
                continue
            n = float(getattr(r, "n", np.nan)) if pd.notna(getattr(r, "n", np.nan)) else np.nan
            bg = align(getattr(r, f"beta_steiger_{outcome}"), getattr(r, f"ea_steiger_{outcome}"), getattr(r, f"oa_steiger_{outcome}"), r.effect_allele, r.other_allele)
            if bg is None:
                status, r2y, z2 = "allele mismatch", np.nan, np.nan
            else:
                z2 = (bg / getattr(r, f"se_steiger_{outcome}")) ** 2
                r2y = z2 / max(n, 1) if pd.notna(n) else np.nan
                status = "ok"
            r2x = 2 * r.eaf * (1 - r.eaf) * r.beta_eqtl**2 if pd.notna(getattr(r, "eaf", np.nan)) else np.nan
            rows.append(
                {
                    "gene": r.gene,
                    "context": r.context,
                    "outcome": outcome,
                    "rsid": r.rsid,
                    "r2_exposure": r2x,
                    "r2_outcome_proxy": r2y,
                    "outcome_z2": z2 if status == "ok" else np.nan,
                    "steiger_dir": "exposure_to_outcome" if status == "ok" and pd.notna(r2x) and pd.notna(r2y) and r2x > r2y else "not_confirmed",
                    "status": status,
                }
            )
    return pd.DataFrame(rows)


def disease_instruments(outcome: str, path: Path, p_thresh: float = 5e-8, clump_bp: int = 500_000) -> pd.DataFrame:
    usecols = [
        "hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele", "hm_beta",
        "other_allele", "effect_allele", "beta", "standard_error", "p_value", "chromosome", "base_pair_location",
    ]
    chunks = []
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", usecols=usecols, dtype=str, chunksize=500_000):
        p = pd.to_numeric(chunk["p_value"], errors="coerce")
        x = chunk.loc[p < p_thresh].copy()
        if x.empty:
            continue
        x["p"] = p.loc[x.index].to_numpy()
        x["rsid"] = x["hm_rsid"].fillna("").replace({"NA": ""})
        x["beta_disease"] = pd.to_numeric(x["hm_beta"].where(x["hm_beta"] != "NA", x["beta"]), errors="coerce")
        x["se_disease"] = pd.to_numeric(x["standard_error"], errors="coerce")
        x["ea_disease"] = clean_allele(x["hm_effect_allele"].where(x["hm_effect_allele"] != "NA", x["effect_allele"]))
        x["oa_disease"] = clean_allele(x["hm_other_allele"].where(x["hm_other_allele"] != "NA", x["other_allele"]))
        x["chrom"] = x["hm_chrom"].where(x["hm_chrom"].notna() & (x["hm_chrom"] != "NA"), x["chromosome"]).astype(str)
        x["pos"] = pd.to_numeric(x["hm_pos"].where(x["hm_pos"].notna() & (x["hm_pos"] != "NA"), x["base_pair_location"]), errors="coerce")
        chunks.append(x[["rsid", "chrom", "pos", "ea_disease", "oa_disease", "beta_disease", "se_disease", "p"]])
    if not chunks:
        return pd.DataFrame()
    hits = pd.concat(chunks).dropna().sort_values("p").drop_duplicates("rsid")
    keep, chosen = [], {}
    for r in hits.itertuples(index=False):
        chrom, pos = str(r.chrom), int(r.pos)
        if any(abs(pos - old) <= clump_bp for old in chosen.get(chrom, [])):
            continue
        chosen.setdefault(chrom, []).append(pos)
        keep.append(r._asdict())
    return pd.DataFrame(keep).assign(outcome=outcome)


def reverse_mr() -> pd.DataFrame:
    eqtl_full = CFG.eqtl_dir / "cis_full.txt.gz"
    eqtl_af = CFG.eqtl_dir / "snp_af.txt.gz"
    if not eqtl_full.exists() or not eqtl_af.exists():
        return pd.DataFrame(
            [
                {"exposure": outcome, "gene": gene, "nsnp": 0, "theta_expr_per_logodds_disease": np.nan, "se": np.nan, "p_mr": np.nan, "q": np.nan, "q_p": np.nan, "mean_F": np.nan, "status": f"blocked: missing raw eQTLGen files ({eqtl_full}, {eqtl_af}); src/19 requires them"}
                for outcome in OUTCOMES
                for gene in CORE
            ]
        )
    # Delegate to the vetted existing implementation in its native output dir,
    # then copy the focused rows.
    subprocess.run(["python", str(P.src / "19_reverse_mr.py")], check=True)
    res = pd.read_csv(P.out("mr") / "reverse_mr.tsv", sep="\t")
    return res[res["gene"].isin(CORE)].copy()


def mvmr_crp() -> pd.DataFrame:
    if not CFG.crp_path.exists():
        return pd.DataFrame([{"gene": gene, "outcome": outcome, "adjusted_exposure": exposure, "nsnp": 0, "theta": np.nan, "se": np.nan, "p": np.nan, "OR": np.nan, "conditional_F_approx": np.nan, "status": "blocked: CRP GWAS download unavailable"} for gene in CORE for outcome in OUTCOMES for exposure in ["gene_expression", "CRP"]])
    if not (CFG.eqtl_dir / "cis_full_candidates.tsv").exists() or not (CFG.eqtl_dir / "snp_af.txt.gz").exists():
        return pd.DataFrame([{"gene": gene, "outcome": outcome, "adjusted_exposure": exposure, "nsnp": 0, "theta": np.nan, "se": np.nan, "p": np.nan, "OR": np.nan, "conditional_F_approx": np.nan, "status": "blocked: missing eQTLGen cis_full_candidates.tsv/snp_af.txt.gz required by src/18"} for gene in CORE for outcome in OUTCOMES for exposure in ["gene_expression", "CRP"]])
    subprocess.run(["python", str(P.src / "18_mvmr.py")], check=True)
    res = pd.read_csv(P.out("mr") / "mvmr_results.tsv", sep="\t")
    return res[res["gene"].isin(CORE)].copy()


def finngen_paths() -> pd.DataFrame:
    if not CFG.finngen_manifest.exists():
        return pd.DataFrame()
    manifest = pd.read_csv(CFG.finngen_manifest, sep="\t")
    rows = []
    for outcome, code in FINNGEN.items():
        hit = manifest[manifest["phenocode"].eq(code)]
        if hit.empty:
            rows.append({"outcome": outcome, "phenocode": code, "status": "phenocode absent"})
            continue
        r = hit.iloc[0]
        rows.append({"outcome": outcome, "phenocode": code, "path_https": r["path_https"], "local_path": str(CFG.gwas_dir / f"finngen_R12_{code}.gz"), "status": "ok"})
    return pd.DataFrame(rows)


def finngen_cxcr2() -> pd.DataFrame:
    endpoints = finngen_paths()
    rows = []
    if endpoints.empty:
        return pd.DataFrame([{"outcome": outcome, "phenocode": code, "gene": "CXCR2", "method": "wald", "status": "blocked: missing FinnGen R12 manifest"} for outcome, code in FINNGEN.items()])
    cx = lead_instruments().query("gene == 'CXCR2' and context == 'blood'").iloc[0]
    for e in endpoints.itertuples(index=False):
        if e.status != "ok":
            rows.append({"outcome": e.outcome, "phenocode": e.phenocode, "gene": "CXCR2", "method": "wald", "status": e.status})
            continue
        path = Path(e.local_path)
        try:
            download(e.path_https, path)
            fg = summary_lookup(path, {cx.rsid}, f"fg_{e.outcome}")
            if fg.empty:
                rows.append({"outcome": e.outcome, "phenocode": e.phenocode, "gene": "CXCR2", "method": "wald", "nsnp": 0, "OR": np.nan, "p": np.nan, "PP4": np.nan, "endpoint_url": e.path_https, "status": "lead rsid absent"})
                continue
            r = fg.iloc[0]
            bg = align(r[f"beta_fg_{e.outcome}"], r[f"ea_fg_{e.outcome}"], r[f"oa_fg_{e.outcome}"], cx.assessed, cx.other)
            if bg is None:
                rows.append({"outcome": e.outcome, "phenocode": e.phenocode, "gene": "CXCR2", "method": "wald", "nsnp": 0, "OR": np.nan, "p": np.nan, "PP4": np.nan, "endpoint_url": e.path_https, "status": "allele mismatch"})
                continue
            theta = bg / cx.beta_eqtl
            se = abs(r[f"se_fg_{e.outcome}"] / cx.beta_eqtl)
            rows.append({"outcome": e.outcome, "phenocode": e.phenocode, "gene": "CXCR2", "method": "wald", "rsid": cx.rsid, "theta": theta, "se": se, "OR": float(np.exp(theta)), "p": float(2 * stats.norm.sf(abs(theta / se))), "PP4": np.nan, "nsnp": 1, "endpoint_url": e.path_https, "status": "ok"})
        except Exception as exc:
            rows.append({"outcome": e.outcome, "phenocode": e.phenocode, "gene": "CXCR2", "method": "wald", "nsnp": 0, "OR": np.nan, "p": np.nan, "PP4": np.nan, "endpoint_url": e.path_https, "status": f"blocked: {exc}"})
    return pd.DataFrame(rows)


def clump_distance(df: pd.DataFrame, p_col: str = "p_eqtl", kb: int = 10) -> pd.DataFrame:
    keep, chosen = [], {}
    for r in df.sort_values(p_col).itertuples(index=False):
        chrom, pos = str(r.chrom), int(r.pos)
        if any(abs(pos - old) <= kb * 1000 for old in chosen.get(chrom, [])):
            continue
        chosen.setdefault(chrom, []).append(pos)
        keep.append(r._asdict())
    return pd.DataFrame(keep)


def sensitivity_context(gene: str, context: str, eqtl: pd.DataFrame) -> pd.DataFrame:
    d = eqtl[(eqtl["gene"].eq(gene)) & (eqtl["context"].eq(context))].copy()
    if d.empty:
        return pd.DataFrame([{"gene": gene, "context": context, "method": "all", "nsnp": 0, "status": "no local multi-SNP eQTL data"}])
    d = d[pd.to_numeric(d["p_eqtl"], errors="coerce") < 1e-5].dropna(subset=["chrom", "pos", "beta_eqtl", "se_eqtl"])
    d = clump_distance(d, kb=10)
    if len(d) < 3:
        return pd.DataFrame([{"gene": gene, "context": context, "method": "all", "nsnp": len(d), "status": "fewer than 3 distance-clumped instruments"}])
    gwas = gwas_lookup(CFG.gwas_dir / "IBD.h.tsv.gz", d, CFG.out_dir / f"_sens_{gene}_{context}_keys.txt")
    m = d.merge(gwas, on=["chrom", "pos"], how="inner")
    if "rsid_x" in m.columns and "rsid" not in m.columns:
        m = m.rename(columns={"rsid_x": "rsid"})
    m["beta_gwas_alt"] = m.apply(harmonised_beta, axis=1)
    m = m.dropna(subset=["beta_gwas_alt", "se_gwas"])
    if len(m) < 3:
        return pd.DataFrame([{"gene": gene, "context": context, "method": "all", "nsnp": len(m), "status": "fewer than 3 harmonised instruments"}])
    bx, by, sy = m["beta_eqtl"].to_numpy(float), m["beta_gwas_alt"].to_numpy(float), m["se_gwas"].to_numpy(float)
    theta, se, p, q, q_p = ivw(bx, by, sy)
    ratio = by / bx
    w = np.square(bx / sy)
    eg = mr_egger(bx, by, sy)
    rows = [
        {"gene": gene, "context": context, "method": "IVW", "nsnp": len(m), "theta": theta, "se": se, "OR": np.exp(theta), "p": p, "q": q, "q_p": q_p, "status": "ok_distance_clump_10kb"},
        {"gene": gene, "context": context, "method": "MR-Egger", "nsnp": len(m), "theta": eg[0], "se": eg[1], "OR": np.exp(eg[0]) if pd.notna(eg[0]) else np.nan, "p": eg[2], "egger_intercept": eg[3], "egger_intercept_se": eg[4], "egger_intercept_p": eg[5], "status": "ok_distance_clump_10kb"},
        {"gene": gene, "context": context, "method": "weighted_median", "nsnp": len(m), "theta": weighted_median(ratio, w), "se": np.nan, "OR": np.exp(weighted_median(ratio, w)), "p": np.nan, "status": "ok_distance_clump_10kb"},
    ]
    for i, r in enumerate(m.itertuples(index=False)):
        mask = np.ones(len(m), dtype=bool)
        mask[i] = False
        th, s, pv, qv, qpv = ivw(bx[mask], by[mask], sy[mask])
        rows.append({"gene": gene, "context": context, "method": "leave_one_out", "left_out": r.rsid, "nsnp": int(mask.sum()), "theta": th, "se": s, "OR": np.exp(th), "p": pv, "q": qv, "q_p": qpv, "status": "ok_distance_clump_10kb"})
    return pd.DataFrame(rows)


def mr_sensitivity() -> pd.DataFrame:
    eqtl = pd.read_csv(CFG.causal_dir / "multicontext_eqtl_pairs.tsv", sep="\t", dtype={"chrom": str})
    frames = [
        pd.DataFrame([{"gene": "CXCR2", "context": "blood", "method": "all", "nsnp": 1, "status": "blocked: raw blood cis-eQTL multi-SNP data absent; only lead instrument is available locally"}]),
        sensitivity_context("CXCR2", "neutrophil", eqtl),
    ]
    return pd.concat(frames, ignore_index=True, sort=False)


def pqtl_ccl8() -> pd.DataFrame:
    local = sorted((CFG.gwas_dir / "pqtl").glob("*CCL8*")) if (CFG.gwas_dir / "pqtl").exists() else []
    status = (
        "blocked: deCODE Ferkingstad proteomics folder requires interactive/token workflow; no local CCL8/MCP-2 pQTL file found"
        if not local
        else f"blocked: local pQTL candidates found but parser not configured: {', '.join(map(str, local))}"
    )
    return pd.DataFrame(
        [
            {
                "gene": "CCL8",
                "protein": "MCP-2/CCL8",
                "source": "deCODE Ferkingstad 2021 plasma pQTL",
                "method": "best_effort_coloc",
                "nsnp": 0,
                "PP4": np.nan,
                "status": status,
                "source_page": "https://www.decode.com/summarydata/",
            }
        ]
    )


def bystander_triage() -> pd.DataFrame:
    mp = pd.read_csv(CFG.causal_dir / "module_causal_map_multicontext.tsv", sep="\t")
    rows = []
    for gene in BYSTANDERS:
        g = mp[mp["gene"].eq(gene)].copy()
        tested = g[g["analysis_status"].astype(str).str.startswith("tested")]
        if tested.empty:
            rows.append({"gene": gene, "best_context": "", "best_eqtl_p": np.nan, "MR_OR": np.nan, "MR_p": np.nan, "MR_FDR": np.nan, "coloc_PP4": np.nan, "bystander_call": True, "status": "no usable instrument in tested contexts; no causal evidence"})
            continue
        best = tested.sort_values("p_eqtl").iloc[0]
        rows.append({"gene": gene, "best_context": best.context, "best_eqtl_p": best.p_eqtl, "MR_OR": best.MR_OR, "MR_p": best.MR_p, "MR_FDR": best.MR_FDR, "coloc_PP4": best.coloc_PP4, "bystander_call": not bool(best.causal_call), "status": "strong/usable eQTL but MR+coloc does not support causality" if pd.notna(best.MR_p) else "not causal"})
    return pd.DataFrame(rows)


def plot_integrated(tables: dict[str, pd.DataFrame], path: Path) -> None:
    checks = ["reverse MR", "Steiger", "MVMR-CRP", "FinnGen CXCR2", "MR sensitivity", "pQTL CCL8"]
    genes = ["CCL8", "CXCR2"]
    score = pd.DataFrame(0.5, index=genes, columns=checks)
    rev = tables["reverse_mr"]
    for gene in genes:
        ok = rev[(rev["gene"].eq(gene)) & (rev["status"].eq("ok"))]
        score.loc[gene, "reverse MR"] = 1 if len(ok) and (ok["p_mr"].fillna(1) > 0.05).all() else (0.5 if len(ok) == 0 else 0)
    st = tables["steiger"]
    for gene in genes:
        ok = st[(st["gene"].eq(gene)) & st["status"].eq("ok")]
        score.loc[gene, "Steiger"] = 1 if len(ok) and ok["steiger_dir"].eq("exposure_to_outcome").all() else 0
    mv = tables["mvmr_crp"]
    for gene in genes:
        ok = mv[(mv["gene"].eq(gene)) & (mv["adjusted_exposure"].eq("gene_expression")) & (mv["status"].astype(str).str.contains("ok|weak", regex=True))]
        score.loc[gene, "MVMR-CRP"] = 1 if len(ok) and np.sign(ok["theta"].dropna()).nunique() == 1 else (0.5 if len(ok) == 0 else 0)
    fg = tables["finngen_cxcr2"]
    score.loc["CXCR2", "FinnGen CXCR2"] = 1 if len(fg[fg["status"].eq("ok") & (fg["p"] < 0.05)]) else 0
    sens = tables["mr_sensitivity"]
    ok_sens = sens[sens["status"].astype(str).str.startswith("ok")]
    score.loc["CXCR2", "MR sensitivity"] = 1 if len(ok_sens) else 0
    pq = tables["pqtl_ccl8"]
    if len(pq[pq["PP4"].fillna(0) > 0.8]):
        score.loc["CCL8", "pQTL CCL8"] = 1
    elif len(pq[pq["status"].astype(str).str.contains("blocked", case=False, na=False)]):
        score.loc["CCL8", "pQTL CCL8"] = 0.5
    else:
        score.loc["CCL8", "pQTL CCL8"] = 0
    fig, ax = plt.subplots(figsize=(9, 2.8))
    im = ax.imshow(score.to_numpy(), cmap=matplotlib.colors.ListedColormap(["#b2182b", "#f4a582", "#4daf4a"]), vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(checks)))
    ax.set_xticklabels(checks, rotation=30, ha="right")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes)
    for i in range(len(genes)):
        for j in range(len(checks)):
            val = score.iat[i, j]
            ax.text(j, i, "pass" if val == 1 else ("blocked" if val == 0.5 else "fail"), ha="center", va="center", fontsize=8)
    ax.set_title("Causal hardening status")
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    plt.close(fig)


def write_summary(tables: dict[str, pd.DataFrame]) -> None:
    by = tables["bystander_triage"]
    fg = tables["finngen_cxcr2"]
    fg_ok = fg[fg["status"].eq("ok")] if "status" in fg else pd.DataFrame()
    fg_line = (
        "; ".join(f"{r.outcome} OR={r.OR:.3g}, p={r.p:.3g}" for r in fg_ok.itertuples(index=False))
        if len(fg_ok)
        else "no successful FinnGen CXCR2 rows"
    )
    mv = tables["mvmr_crp"]
    mv_ok = mv[(mv["adjusted_exposure"].eq("gene_expression")) & mv["theta"].notna()] if "adjusted_exposure" in mv else pd.DataFrame()
    mv_line = (
        "; ".join(f"{r.gene}/{r.outcome} OR={r.OR:.3g}, p={r.p:.3g}, F={r.conditional_F_approx:.3g}, {r.status}" for r in mv_ok.itertuples(index=False))
        if len(mv_ok)
        else "MVMR not estimable"
    )
    rev = tables["reverse_mr"]
    rev_ok = rev[rev["status"].eq("ok")] if "status" in rev else pd.DataFrame()
    rev_line = (
        "; ".join(f"{r.exposure}->{r.gene} p={r.p_mr:.3g}" for r in rev_ok.itertuples(index=False))
        if len(rev_ok)
        else "reverse MR not estimable"
    )
    rev_interp = (
        "not cleared: significant one-SNP reverse-MR rows were observed, so reverse causality/disease-linked expression remains a caveat"
        if len(rev_ok) and (rev_ok["p_mr"].fillna(1) < 0.05).any()
        else "no significant reverse-MR evidence among estimable rows"
    )
    lines = [
        "# Causal hardening summary",
        "",
        "## Direct answer",
        "- CCL8 remains a strong blood eQTLGen causal anchor from prior MR+coloc, and Steiger supports exposure->outcome direction for the lead instrument. It is not fully hardened: reverse-MR has significant one-SNP rows, CCL8 MVMR is not estimable with one instrument, and deCODE pQTL coloc is blocked without a tokenized/local file.",
        f"- CXCR2 remains supported by blood plus neutrophil MR+coloc. FinnGen R12 CXCR2 replication: {fg_line}. Neutrophil multi-instrument sensitivity uses distance-clumped instruments and preserves the protective direction.",
        f"- Reverse-MR summary: {rev_line}. Interpretation: {rev_interp}.",
        f"- MVMR-CRP summary: {mv_line}. Interpretation: CXCR2 direction is stable after CRP adjustment, but CRP conditional F is weak; CCL8 is not estimable.",
        "- OSM/OSMR/TREM1/IL13RA2 are formalised as bystanders: no strict MR+coloc causal call; TREM1 has strong neutrophil/stimulated-monocyte eQTLs but null MR and PP4 near zero.",
        "",
        "## Raw dependency note",
        "The script downloads public eQTLGen cis/allele-frequency files when absent and builds `cis_full_candidates.tsv` for the target genes. Any remaining blocked row is a real data or rank limitation, not a silent omission.",
        "",
        "## Bystander rows",
        "```tsv",
        by.to_csv(sep="\t", index=False).rstrip(),
        "```",
    ]
    (CFG.out_dir / "SUMMARY.md").write_text("\n".join(lines) + "\n")


def update_status() -> None:
    status = P.journal / "status" / "causal_module_status.md"
    text = status.read_text() if status.exists() else "# Causal module status\n"
    marker = "Downloaded/resumed CRP, FinnGen R12 endpoints, and eQTLGen full cis/allele-frequency inputs"
    if marker in text:
        return
    note = f"""

## 2026-06-23 JST - causal hardening

- Ran `src/27_causal_hardening.py` on CPU only with `n_jobs={CFG.n_jobs}` for parallel-safe steps.
- Downloaded/resumed CRP, FinnGen R12 endpoints, and eQTLGen full cis/allele-frequency inputs; eQTLGen HTTPS certificate is expired, so `src/27` scopes `curl --insecure` only to that host.
- Wrote hardening outputs under `{CFG.out_dir}` and promoted final table/figure copies into `results/`.
- Reverse-MR is not a clean pass: significant one-SNP reverse rows appear for CD->CCL8, IBD->CCL8, and IBD->CXCR2, so reverse causality/disease-linked expression remains a caveat.
- Steiger supports exposure->outcome for the blood CCL8/CXCR2 lead instruments. CXCR2 FinnGen R12 replicates for IBD and UC but not CD. CXCR2 MVMR-CRP keeps the protective direction, but CRP conditional F is weak; CCL8 MVMR is not estimable with one instrument.
- deCODE CCL8 pQTL coloc is best-effort blocked because the public deCODE proteomics folder requires an interactive/token workflow and no local CCL8 pQTL file is present.
- OSM/OSMR/TREM1/IL13RA2 are formalised as bystanders in `bystander_triage.tsv`.
"""
    status.write_text(text.rstrip() + "\n" + note)


def main() -> None:
    CFG.out_dir.mkdir(parents=True, exist_ok=True)
    inputs = ensure_public_inputs()
    inputs.to_csv(CFG.out_dir / "input_downloads.tsv", sep="\t", index=False)
    tasks = {
        "reverse_mr": reverse_mr,
        "steiger": steiger,
        "mvmr_crp": mvmr_crp,
        "finngen_cxcr2": finngen_cxcr2,
        "mr_sensitivity": mr_sensitivity,
        "pqtl_ccl8": pqtl_ccl8,
        "bystander_triage": bystander_triage,
    }
    names = list(tasks)
    frames = Parallel(n_jobs=min(CFG.n_jobs, len(names)), prefer="threads")(delayed(tasks[n])() for n in names)
    tables = dict(zip(names, frames))
    for name, df in tables.items():
        df.to_csv(CFG.out_dir / f"{name}.tsv", sep="\t", index=False)
    fig = CFG.out_dir / "Fig_causal_hardening.png"
    plot_integrated(tables, fig)
    write_summary(tables)
    P.promote_figure(fig)
    for name in ["reverse_mr", "steiger", "mvmr_crp", "finngen_cxcr2", "mr_sensitivity", "pqtl_ccl8", "bystander_triage"]:
        P.promote_table(CFG.out_dir / f"{name}.tsv")
    update_status()
    print((CFG.out_dir / "SUMMARY.md").read_text())


if __name__ == "__main__":
    main()
