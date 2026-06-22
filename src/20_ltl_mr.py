"""MR: Codd 2021 leukocyte telomere length (LTL) -> IBD/CD/UC.

The exposure file is the GWAS Catalog harmonised summary statistics for
GCST90002398. Instruments are genome-wide significant and distance-clumped
because no local LD panel/R stack is available.
"""
from __future__ import annotations

import gzip
import shlex
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict
from scipy import stats

from paths import P


class Config(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    out_dir: Path = P.out("mr")
    exposure: Path = P.raw / "gwas" / "LTL_Codd2021_GCST90002398.h.tsv.gz"
    exposure_url: str = (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90002001-GCST90003000/GCST90002398/harmonised/"
        "32888494-GCST90002398-EFO_0004833.h.tsv.gz"
    )
    gwas: dict[str, Path] = {
        "IBD": P.raw / "gwas" / "IBD.h.tsv.gz",
        "CD": P.raw / "gwas" / "CD.h.tsv.gz",
        "UC": P.raw / "gwas" / "UC.h.tsv.gz",
    }
    gws_p: float = 5e-8
    clump_kb: int = 500
    bootstrap_n: int = 1000
    seed: int = 20260622


CFG = Config()


def run_text(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def blocked_table(reason: str) -> None:
    rows = [
        {
            "outcome": outcome,
            "method": method,
            "nsnp": 0,
            "theta_per_higher_LTL": np.nan,
            "se": np.nan,
            "p": np.nan,
            "OR_per_higher_LTL": np.nan,
            "OR_per_shorter_LTL": np.nan,
            "q": np.nan,
            "q_p": np.nan,
            "egger_intercept": np.nan,
            "egger_intercept_p": np.nan,
            "status": reason,
        }
        for outcome in CFG.gwas
        for method in ["IVW", "MR_Egger", "weighted_median"]
    ]
    pd.DataFrame(rows).to_csv(CFG.out_dir / "mr_LTL.tsv", sep="\t", index=False)


def clean_allele(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.upper().replace({"NA": "", "NAN": ""})


def clump_by_distance(df: pd.DataFrame, window_bp: int) -> pd.DataFrame:
    keep = []
    chosen: dict[str, list[int]] = {}
    for r in df.sort_values("p").itertuples(index=False):
        chrom = str(r.chrom)
        pos = int(r.pos)
        if any(abs(pos - old) <= window_bp for old in chosen.get(chrom, [])):
            continue
        chosen.setdefault(chrom, []).append(pos)
        keep.append(r._asdict())
    return pd.DataFrame(keep)


def read_ltl_instruments() -> pd.DataFrame:
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
    for chunk in pd.read_csv(CFG.exposure, sep="\t", compression="gzip", usecols=usecols, dtype=str, chunksize=750_000):
        p = pd.to_numeric(chunk["p_value"], errors="coerce")
        m = p < CFG.gws_p
        if not m.any():
            continue
        x = chunk.loc[m].copy()
        x["p"] = p.loc[m].to_numpy()
        x["rsid"] = x["hm_rsid"].fillna("").replace({"NA": ""})
        x["beta_ltl"] = pd.to_numeric(x["hm_beta"].where(x["hm_beta"] != "NA", x["beta"]), errors="coerce")
        x["se_ltl"] = pd.to_numeric(x["standard_error"], errors="coerce")
        x["ea_ltl"] = clean_allele(x["hm_effect_allele"].where(x["hm_effect_allele"] != "NA", x["effect_allele"]))
        x["oa_ltl"] = clean_allele(x["hm_other_allele"].where(x["hm_other_allele"] != "NA", x["other_allele"]))
        x["chrom"] = x["hm_chrom"].where(x["hm_chrom"].notna() & (x["hm_chrom"] != "NA"), x["chromosome"])
        x["pos"] = pd.to_numeric(
            x["hm_pos"].where(x["hm_pos"].notna() & (x["hm_pos"] != "NA"), x["base_pair_location"]),
            errors="coerce",
        )
        chunks.append(x[["rsid", "chrom", "pos", "ea_ltl", "oa_ltl", "beta_ltl", "se_ltl", "p"]])
    if not chunks:
        return pd.DataFrame(columns=["rsid", "chrom", "pos", "ea_ltl", "oa_ltl", "beta_ltl", "se_ltl", "p"])
    hits = pd.concat(chunks, ignore_index=True)
    hits = hits.dropna(subset=["rsid", "chrom", "pos", "beta_ltl", "se_ltl", "p"])
    hits = hits[(hits["rsid"] != "") & (hits["se_ltl"] > 0)]
    hits = hits.sort_values("p").drop_duplicates("rsid")
    return clump_by_distance(hits, CFG.clump_kb * 1000)


def write_rsids(rsids: list[str]) -> Path:
    rsfile = CFG.out_dir / "_ltl_rsids.txt"
    rsfile.write_text("\n".join(sorted(rsids)) + "\n")
    return rsfile


def gwas_lookup(path: Path, rsfile: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    idx = {c: i + 1 for i, c in enumerate(hdr)}
    required = [
        "hm_rsid",
        "hm_beta",
        "hm_effect_allele",
        "hm_other_allele",
        "beta",
        "effect_allele",
        "other_allele",
        "standard_error",
        "p_value",
    ]
    missing = [c for c in required if c not in idx]
    if missing:
        raise ValueError(f"{path} missing required columns: {missing}")
    cmd = (
        f"gunzip -c {shlex.quote(str(path))} | "
        "awk -F'\\t' -v OFS='\\t' "
        f"-v rs={idx['hm_rsid']} -v hb={idx['hm_beta']} "
        f"-v hea={idx['hm_effect_allele']} -v hoa={idx['hm_other_allele']} "
        f"-v b={idx['beta']} -v ea={idx['effect_allele']} -v oa={idx['other_allele']} "
        f"-v se={idx['standard_error']} -v pv={idx['p_value']} "
        "'NR==FNR{want[$1]=1; next} FNR==1{next} "
        "($rs in want){print $rs,$hb,$hea,$hoa,$b,$ea,$oa,$se,$pv}' "
        f"{shlex.quote(str(rsfile))} -"
    )
    rows = []
    for line in run_text(cmd).splitlines():
        rsid, hbeta, hea, hoa, beta, ea, oa, se, p = line.split("\t")
        b = hbeta if hbeta != "NA" else beta
        effect = hea if hbeta != "NA" else ea
        other = hoa if hbeta != "NA" else oa
        try:
            rows.append(
                {
                    "rsid": rsid,
                    "beta_out": float(b),
                    "ea_out": effect.upper(),
                    "oa_out": other.upper(),
                    "se_out": float(se),
                    "p_out": float(p),
                }
            )
        except ValueError:
            continue
    return pd.DataFrame(rows)


def align_outcome(row: pd.Series) -> float | None:
    if row.ea_out == row.ea_ltl and row.oa_out == row.oa_ltl:
        return row.beta_out
    if row.ea_out == row.oa_ltl and row.oa_out == row.ea_ltl:
        return -row.beta_out
    return None


def ivw(beta_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray) -> tuple[float, float, float, float, float]:
    w = 1 / np.square(se_y)
    denom = np.sum(w * np.square(beta_x))
    theta = float(np.sum(w * beta_x * beta_y) / denom)
    se = float(np.sqrt(1 / denom))
    p = float(2 * stats.norm.sf(abs(theta / se)))
    q = float(np.sum(w * np.square(beta_y - theta * beta_x)))
    q_p = float(stats.chi2.sf(q, len(beta_x) - 1)) if len(beta_x) > 1 else np.nan
    return theta, se, p, q, q_p


def egger(beta_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray) -> tuple[float, float, float, float, float, float]:
    w = 1 / np.square(se_y)
    x = np.column_stack([np.ones_like(beta_x), beta_x])
    xtw = x.T * w
    cov = np.linalg.inv(xtw @ x)
    coef = cov @ (xtw @ beta_y)
    resid = beta_y - x @ coef
    if len(beta_x) > 2:
        sigma2 = max(float(np.sum(w * np.square(resid)) / (len(beta_x) - 2)), 1.0)
    else:
        sigma2 = 1.0
    cov = cov * sigma2
    se = np.sqrt(np.diag(cov))
    intercept, slope = float(coef[0]), float(coef[1])
    intercept_se, slope_se = float(se[0]), float(se[1])
    slope_p = float(2 * stats.norm.sf(abs(slope / slope_se)))
    intercept_p = float(2 * stats.norm.sf(abs(intercept / intercept_se)))
    return slope, slope_se, slope_p, intercept, intercept_se, intercept_p


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cutoff = weights.sum() / 2
    return float(values[np.searchsorted(np.cumsum(weights), cutoff)])


def weighted_median_with_bootstrap(beta_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray) -> tuple[float, float, float]:
    ratios = beta_y / beta_x
    weights = np.square(beta_x) / np.square(se_y)
    theta = weighted_median(ratios, weights)
    rng = np.random.default_rng(CFG.seed)
    probs = weights / weights.sum()
    boot = []
    n = len(ratios)
    for _ in range(CFG.bootstrap_n):
        idx = rng.choice(n, size=n, replace=True, p=probs)
        boot.append(weighted_median(ratios[idx], weights[idx]))
    se = float(np.std(boot, ddof=1))
    p = float(2 * stats.norm.sf(abs(theta / se))) if se > 0 else np.nan
    return theta, se, p


def mr_rows(outcome: str, d: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    rows = []
    loo_rows = []
    bx, by, sy = d["beta_ltl"].to_numpy(), d["beta_out_aligned"].to_numpy(), d["se_out"].to_numpy()
    theta, se, p, q, q_p = ivw(bx, by, sy)
    rows.append(result_row(outcome, "IVW", len(d), theta, se, p, q, q_p))
    if len(d) >= 3:
        slope, slope_se, slope_p, intercept, intercept_se, intercept_p = egger(bx, by, sy)
        rows.append(result_row(outcome, "MR_Egger", len(d), slope, slope_se, slope_p, np.nan, np.nan, intercept, intercept_p))
    else:
        rows.append(result_row(outcome, "MR_Egger", len(d), np.nan, np.nan, np.nan, np.nan, np.nan, status="need >=3 SNPs"))
    if len(d) >= 3:
        wm, wm_se, wm_p = weighted_median_with_bootstrap(bx, by, sy)
        rows.append(result_row(outcome, "weighted_median", len(d), wm, wm_se, wm_p, np.nan, np.nan))
    else:
        rows.append(result_row(outcome, "weighted_median", len(d), np.nan, np.nan, np.nan, np.nan, np.nan, status="need >=3 SNPs"))
    if len(d) > 2:
        for rsid in d["rsid"]:
            sub = d[d["rsid"] != rsid]
            t, s, p_loo, _, _ = ivw(sub["beta_ltl"].to_numpy(), sub["beta_out_aligned"].to_numpy(), sub["se_out"].to_numpy())
            loo_rows.append(
                {
                    "outcome": outcome,
                    "left_out_rsid": rsid,
                    "nsnp": len(sub),
                    "theta_per_higher_LTL": t,
                    "se": s,
                    "p": p_loo,
                    "OR_per_higher_LTL": float(np.exp(t)),
                    "OR_per_shorter_LTL": float(np.exp(-t)),
                }
            )
    return rows, loo_rows


def result_row(
    outcome: str,
    method: str,
    nsnp: int,
    theta: float,
    se: float,
    p: float,
    q: float,
    q_p: float,
    egger_intercept: float = np.nan,
    egger_intercept_p: float = np.nan,
    status: str = "ok",
) -> dict:
    return {
        "outcome": outcome,
        "method": method,
        "nsnp": nsnp,
        "theta_per_higher_LTL": theta,
        "se": se,
        "p": p,
        "OR_per_higher_LTL": float(np.exp(theta)) if np.isfinite(theta) else np.nan,
        "OR_per_shorter_LTL": float(np.exp(-theta)) if np.isfinite(theta) else np.nan,
        "q": q,
        "q_p": q_p,
        "egger_intercept": egger_intercept,
        "egger_intercept_p": egger_intercept_p,
        "status": status,
    }


def main() -> None:
    CFG.out_dir.mkdir(parents=True, exist_ok=True)
    if not CFG.exposure.exists():
        reason = f"missing LTL exposure file; download {CFG.exposure_url} to {CFG.exposure}"
        blocked_table(reason)
        print(reason)
        return

    try:
        ltl = read_ltl_instruments()
    except (OSError, EOFError, pd.errors.ParserError) as exc:
        reason = f"LTL exposure unreadable or partial: {exc}"
        blocked_table(reason)
        print(reason)
        return

    ltl.to_csv(CFG.out_dir / "ltl_instruments.tsv", sep="\t", index=False)
    print(f"LTL instruments after distance clump: {len(ltl)}")
    if ltl.empty:
        blocked_table("no genome-wide significant LTL instruments after clumping")
        return
    rsfile = write_rsids(list(ltl["rsid"]))

    all_details = []
    all_rows = []
    all_loo = []
    for outcome, path in CFG.gwas.items():
        out = gwas_lookup(path, rsfile)
        d = ltl.merge(out, on="rsid", how="inner")
        if d.empty:
            all_rows.extend(
                [
                    result_row(outcome, "IVW", 0, np.nan, np.nan, np.nan, np.nan, np.nan, status="no outcome overlap"),
                    result_row(outcome, "MR_Egger", 0, np.nan, np.nan, np.nan, np.nan, np.nan, status="no outcome overlap"),
                    result_row(outcome, "weighted_median", 0, np.nan, np.nan, np.nan, np.nan, np.nan, status="no outcome overlap"),
                ]
            )
            continue
        d["beta_out_aligned"] = d.apply(align_outcome, axis=1)
        d = d.dropna(subset=["beta_out_aligned"])
        d = d[(d["se_out"] > 0) & (d["beta_ltl"] != 0)]
        print(f"{outcome}: {len(d)} harmonised LTL instruments")
        if d.empty:
            all_rows.extend(
                [
                    result_row(outcome, "IVW", 0, np.nan, np.nan, np.nan, np.nan, np.nan, status="allele harmonisation removed all SNPs"),
                    result_row(outcome, "MR_Egger", 0, np.nan, np.nan, np.nan, np.nan, np.nan, status="allele harmonisation removed all SNPs"),
                    result_row(outcome, "weighted_median", 0, np.nan, np.nan, np.nan, np.nan, np.nan, status="allele harmonisation removed all SNPs"),
                ]
            )
            continue
        all_details.append(d.assign(outcome=outcome))
        rows, loo = mr_rows(outcome, d)
        all_rows.extend(rows)
        all_loo.extend(loo)

    if all_details:
        pd.concat(all_details, ignore_index=True).to_csv(CFG.out_dir / "mr_LTL_instruments.tsv", sep="\t", index=False)
    if all_loo:
        pd.DataFrame(all_loo).to_csv(CFG.out_dir / "mr_LTL_leave_one_out.tsv", sep="\t", index=False)
    res = pd.DataFrame(all_rows)
    res.to_csv(CFG.out_dir / "mr_LTL.tsv", sep="\t", index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
