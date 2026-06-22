"""Minimal local MVMR: gene expression vs CRP for IBD/CD/UC outcomes.

This is intentionally conservative and local: it uses target-gene cis-eQTL
variants already prefetched for coloc, then looks up CRP and disease effects by
rsid with awk hash joins. If CRP data is absent/partial or the local cis region
does not provide enough rank for two-exposure MVMR, the output records that
status instead of fabricating an adjusted estimate.
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
    cis_candidates: Path = P.raw / "eqtlgen" / "cis_full_candidates.tsv"
    eqtl_af: Path = P.raw / "eqtlgen" / "snp_af.txt.gz"
    crp: Path = P.raw / "gwas" / "CRP_GCST90029070.h.tsv.gz"
    crp_url: str = (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90029001-GCST90030000/GCST90029070/harmonised/"
        "35459240-GCST90029070-EFO_0004458.h.tsv.gz"
    )
    gwas: dict[str, Path] = {
        "IBD": P.raw / "gwas" / "IBD.h.tsv.gz",
        "CD": P.raw / "gwas" / "CD.h.tsv.gz",
        "UC": P.raw / "gwas" / "UC.h.tsv.gz",
    }
    genes: tuple[str, ...] = ("CCL8", "CXCR2")
    eqtl_p: float = 5e-8
    clump_kb: int = 100


CFG = Config()


def run_text(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def blocked(reason: str) -> None:
    rows = []
    for gene in CFG.genes:
        for outcome in CFG.gwas:
            for exposure in ["gene_expression", "CRP"]:
                rows.append(
                    {
                        "gene": gene,
                        "outcome": outcome,
                        "adjusted_exposure": exposure,
                        "nsnp": 0,
                        "theta": np.nan,
                        "se": np.nan,
                        "p": np.nan,
                        "OR": np.nan,
                        "q": np.nan,
                        "q_p": np.nan,
                        "conditional_F_approx": np.nan,
                        "status": reason,
                    }
                )
    pd.DataFrame(rows).to_csv(CFG.out_dir / "mvmr_results.tsv", sep="\t", index=False)


def clump_by_distance(df: pd.DataFrame, window_bp: int) -> pd.DataFrame:
    keep = []
    chosen: dict[str, list[int]] = {}
    for r in df.sort_values("p_eqtl").itertuples(index=False):
        chrom = str(r.snp_chr)
        pos = int(r.snp_pos)
        if any(abs(pos - old) <= window_bp for old in chosen.get(chrom, [])):
            continue
        chosen.setdefault(chrom, []).append(pos)
        keep.append(r._asdict())
    return pd.DataFrame(keep)


def lookup_af(rsids: set[str]) -> dict[str, tuple[str, float]]:
    rsfile = CFG.out_dir / "_mvmr_rsids.txt"
    rsfile.write_text("\n".join(sorted(rsids)) + "\n")
    cmd = (
        f"gunzip -c {shlex.quote(str(CFG.eqtl_af))} | "
        "awk -F'\\t' -v OFS='\\t' "
        "'NR==FNR{want[$1]; next} FNR==1{next} ($1 in want){print $1,$5,$9}' "
        f"{shlex.quote(str(rsfile))} -"
    )
    af = {}
    for line in run_text(cmd).splitlines():
        rsid, allele_b, freq_b = line.split("\t")
        try:
            af[rsid] = (allele_b.upper(), float(freq_b))
        except ValueError:
            continue
    return af


def read_gene_cis() -> pd.DataFrame:
    cols = ["Pvalue", "SNP", "SNPChr", "SNPPos", "Zscore", "AssessedAllele", "OtherAllele", "GeneSymbol", "NrSamples"]
    df = pd.read_csv(CFG.cis_candidates, sep="\t", usecols=cols)
    df = df[df["GeneSymbol"].isin(CFG.genes)].copy()
    df = df.rename(
        columns={
            "Pvalue": "p_eqtl",
            "SNP": "rsid",
            "SNPChr": "snp_chr",
            "SNPPos": "snp_pos",
            "Zscore": "z",
            "AssessedAllele": "assessed",
            "OtherAllele": "other",
            "GeneSymbol": "gene",
            "NrSamples": "n",
        }
    )
    for c in ["p_eqtl", "snp_chr", "snp_pos", "z", "n"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["p_eqtl", "snp_chr", "snp_pos", "z", "n"])
    df = df[df["p_eqtl"] < CFG.eqtl_p]
    af = lookup_af(set(df["rsid"]))
    rows = []
    for r in df.itertuples(index=False):
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
        rows.append({**r._asdict(), "assessed": assessed, "other": other, "eaf": eaf, "beta_eqtl": r.z * se, "se_eqtl": se})
    out = []
    for gene, g in pd.DataFrame(rows).groupby("gene"):
        c = clump_by_distance(g, CFG.clump_kb * 1000)
        out.append(c)
        print(f"{gene}: {len(g)} eQTL p<{CFG.eqtl_p:g}, {len(c)} after {CFG.clump_kb}kb distance clump")
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def write_rsids(rsids: list[str]) -> Path:
    rsfile = CFG.out_dir / "_mvmr_rsids.txt"
    rsfile.write_text("\n".join(sorted(rsids)) + "\n")
    return rsfile


def summary_lookup(path: Path, rsfile: Path, prefix: str) -> pd.DataFrame:
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
                    f"beta_{prefix}": float(b),
                    f"se_{prefix}": float(se),
                    f"p_{prefix}": float(p),
                    f"ea_{prefix}": effect.upper(),
                    f"oa_{prefix}": other.upper(),
                }
            )
        except ValueError:
            continue
    return pd.DataFrame(rows)


def align_beta(beta: float, ea: str, oa: str, assessed: str, other: str) -> float | None:
    if ea == assessed and oa == other:
        return beta
    if ea == other and oa == assessed:
        return -beta
    return None


def harmonise(d: pd.DataFrame, prefix: str) -> pd.Series:
    vals = []
    for r in d.itertuples(index=False):
        vals.append(align_beta(getattr(r, f"beta_{prefix}"), getattr(r, f"ea_{prefix}"), getattr(r, f"oa_{prefix}"), r.assessed, r.other))
    return pd.Series(vals, index=d.index, dtype=float)


def wls_no_intercept(x: np.ndarray, y: np.ndarray, se_y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    w = 1 / np.square(se_y)
    xtw = x.T * w
    xtwx = xtw @ x
    beta = np.linalg.inv(xtwx) @ (xtw @ y)
    resid = y - x @ beta
    q = float(np.sum(w * np.square(resid)))
    df = len(y) - x.shape[1]
    q_p = float(stats.chi2.sf(q, df)) if df > 0 else np.nan
    scale = max(q / df, 1.0) if df > 0 else 1.0
    se = np.sqrt(np.diag(np.linalg.inv(xtwx) * scale))
    return beta, se, q, q_p


def conditional_f(xj: np.ndarray, xother: np.ndarray, sej: np.ndarray) -> float:
    w = 1 / np.square(sej)
    xo = xother.reshape(-1, 1)
    denom = float(((xo.T * w) @ xo)[0, 0])
    if denom <= 0:
        return np.nan
    b = float(((xo.T * w) @ xj)[0] / denom)
    resid = xj - b * xother
    df = len(xj) - 2
    if df <= 0:
        return np.nan
    return float(np.sum(w * np.square(resid)) / df)


def mvmr_for_gene_outcome(gene: str, outcome: str, d: pd.DataFrame) -> list[dict]:
    if len(d) < 3:
        return status_rows(gene, outcome, len(d), "need >=3 harmonised instruments for two-exposure MVMR")
    x = d[["beta_eqtl", "beta_crp_aligned"]].to_numpy()
    if np.linalg.matrix_rank(x) < 2:
        return status_rows(gene, outcome, len(d), "exposure matrix rank <2 after harmonisation")
    coef, se, q, q_p = wls_no_intercept(x, d["beta_out_aligned"].to_numpy(), d["se_out"].to_numpy())
    f_gene = conditional_f(d["beta_eqtl"].to_numpy(), d["beta_crp_aligned"].to_numpy(), d["se_eqtl"].to_numpy())
    f_crp = conditional_f(d["beta_crp_aligned"].to_numpy(), d["beta_eqtl"].to_numpy(), d["se_crp"].to_numpy())
    model_status = "ok" if np.nanmin([f_gene, f_crp]) >= 10 else "weak conditional F for at least one exposure"
    rows = []
    for label, i, fval in [("gene_expression", 0, f_gene), ("CRP", 1, f_crp)]:
        p = float(2 * stats.norm.sf(abs(coef[i] / se[i])))
        rows.append(
            {
                "gene": gene,
                "outcome": outcome,
                "adjusted_exposure": label,
                "nsnp": len(d),
                "theta": float(coef[i]),
                "se": float(se[i]),
                "p": p,
                "OR": float(np.exp(coef[i])),
                "q": q,
                "q_p": q_p,
                "conditional_F_approx": fval,
                "status": model_status,
            }
        )
    return rows


def status_rows(gene: str, outcome: str, nsnp: int, status: str) -> list[dict]:
    return [
        {
            "gene": gene,
            "outcome": outcome,
            "adjusted_exposure": exposure,
            "nsnp": nsnp,
            "theta": np.nan,
            "se": np.nan,
            "p": np.nan,
            "OR": np.nan,
            "q": np.nan,
            "q_p": np.nan,
            "conditional_F_approx": np.nan,
            "status": status,
        }
        for exposure in ["gene_expression", "CRP"]
    ]


def main() -> None:
    CFG.out_dir.mkdir(parents=True, exist_ok=True)
    if not CFG.crp.exists():
        reason = f"missing CRP exposure file; download {CFG.crp_url} to {CFG.crp}"
        blocked(reason)
        print(reason)
        return
    try:
        gene_cis = read_gene_cis()
    except (OSError, EOFError, pd.errors.ParserError) as exc:
        reason = f"cannot read local cis candidate/eQTL data: {exc}"
        blocked(reason)
        print(reason)
        return
    if gene_cis.empty:
        blocked("no target-gene cis instruments")
        print("no target-gene cis instruments")
        return
    rsfile = write_rsids(list(gene_cis["rsid"]))
    try:
        crp = summary_lookup(CFG.crp, rsfile, "crp")
    except (OSError, EOFError, pd.errors.ParserError, subprocess.CalledProcessError) as exc:
        reason = f"CRP exposure unreadable or partial: {exc}"
        blocked(reason)
        print(reason)
        return

    result_rows = []
    detail_rows = []
    for outcome, path in CFG.gwas.items():
        out = summary_lookup(path, rsfile, "out")
        for gene in CFG.genes:
            d = gene_cis[gene_cis["gene"] == gene].merge(crp, on="rsid", how="inner").merge(out, on="rsid", how="inner")
            if d.empty:
                result_rows.extend(status_rows(gene, outcome, 0, "no CRP/outcome overlap for clumped cis instruments"))
                continue
            d["beta_crp_aligned"] = harmonise(d, "crp")
            d["beta_out_aligned"] = harmonise(d, "out")
            d = d.dropna(subset=["beta_crp_aligned", "beta_out_aligned"])
            d = d[(d["se_crp"] > 0) & (d["se_out"] > 0) & (d["se_eqtl"] > 0)]
            detail_rows.append(d.assign(outcome=outcome))
            result_rows.extend(mvmr_for_gene_outcome(gene, outcome, d))

    if detail_rows:
        pd.concat(detail_rows, ignore_index=True).to_csv(CFG.out_dir / "mvmr_instruments.tsv", sep="\t", index=False)
    res = pd.DataFrame(result_rows).sort_values(["outcome", "gene", "adjusted_exposure"])
    res.to_csv(CFG.out_dir / "mvmr_results.tsv", sep="\t", index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
