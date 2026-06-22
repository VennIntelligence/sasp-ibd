"""Reverse MR: IBD/CD/UC liability -> SenMayo candidate gene expression.

Disease instruments come from genome-wide significant de Lange GWAS variants,
distance-clumped without an LD panel. eQTLGen full cis statistics are queried
with an awk hash join against only the needed rsids and genes.
"""
from __future__ import annotations

import shlex
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict
from scipy import stats
from statsmodels.stats.multitest import multipletests

from paths import P


class Config(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    out_dir: Path = P.out("mr")
    eqtl_full: Path = P.raw / "eqtlgen" / "cis_full.txt.gz"
    eqtl_af: Path = P.raw / "eqtlgen" / "snp_af.txt.gz"
    gwas: dict[str, Path] = {
        "IBD": P.raw / "gwas" / "IBD.h.tsv.gz",
        "CD": P.raw / "gwas" / "CD.h.tsv.gz",
        "UC": P.raw / "gwas" / "UC.h.tsv.gz",
    }
    gws_p: float = 5e-8
    clump_kb: int = 500


CFG = Config()


def run_text(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def _clean_allele(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.upper().replace({"NA": "", "NAN": ""})


def target_genes() -> list[str]:
    genes: set[str] = {"CCL8", "CXCR2"}
    for outcome in CFG.gwas:
        path = CFG.out_dir / f"mr_{outcome}.tsv"
        if path.exists():
            mr = pd.read_csv(path, sep="\t")
            genes |= set(mr.loc[mr["fdr"] < 0.05, "gene"])
    coloc_path = CFG.out_dir / "coloc_IBD.tsv"
    if coloc_path.exists():
        coloc = pd.read_csv(coloc_path, sep="\t")
        genes |= set(coloc.loc[coloc["PP4"] > 0.5, "gene"])
    return sorted(genes)


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


def read_disease_instruments(outcome: str, path: Path) -> pd.DataFrame:
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
        m = p < CFG.gws_p
        if not m.any():
            continue
        x = chunk.loc[m].copy()
        x["p"] = p.loc[m].to_numpy()
        x["rsid"] = x["hm_rsid"].fillna("").replace({"NA": ""})
        x["beta_exp"] = pd.to_numeric(x["hm_beta"].where(x["hm_beta"] != "NA", x["beta"]), errors="coerce")
        x["se_exp"] = pd.to_numeric(x["standard_error"], errors="coerce")
        x["ea"] = _clean_allele(x["hm_effect_allele"].where(x["hm_effect_allele"] != "NA", x["effect_allele"]))
        x["oa"] = _clean_allele(x["hm_other_allele"].where(x["hm_other_allele"] != "NA", x["other_allele"]))
        x["chrom"] = x["hm_chrom"].where(x["hm_chrom"].notna() & (x["hm_chrom"] != "NA"), x["chromosome"])
        x["pos"] = pd.to_numeric(
            x["hm_pos"].where(x["hm_pos"].notna() & (x["hm_pos"] != "NA"), x["base_pair_location"]),
            errors="coerce",
        )
        chunks.append(x[["rsid", "chrom", "pos", "ea", "oa", "beta_exp", "se_exp", "p"]])
    if not chunks:
        return pd.DataFrame(columns=["outcome", "rsid", "chrom", "pos", "ea", "oa", "beta_exp", "se_exp", "p"])
    hits = pd.concat(chunks, ignore_index=True)
    hits = hits.dropna(subset=["rsid", "chrom", "pos", "beta_exp", "se_exp", "p"])
    hits = hits[(hits["rsid"] != "") & (hits["se_exp"] > 0)]
    hits = hits.sort_values("p").drop_duplicates("rsid")
    clumped = clump_by_distance(hits, CFG.clump_kb * 1000)
    clumped.insert(0, "outcome", outcome)
    return clumped


def write_values(path: Path, values: list[str]) -> None:
    path.write_text("\n".join(values) + "\n")


def lookup_af(rsids: set[str]) -> dict[str, tuple[str, float]]:
    rsfile = CFG.out_dir / "_reverse_rsids.txt"
    write_values(rsfile, sorted(rsids))
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


def lookup_eqtl(rsids: set[str], genes: list[str]) -> pd.DataFrame:
    rsfile = CFG.out_dir / "_reverse_rsids.txt"
    genefile = CFG.out_dir / "_reverse_genes.txt"
    write_values(rsfile, sorted(rsids))
    write_values(genefile, genes)
    cmd = (
        f"gunzip -c {shlex.quote(str(CFG.eqtl_full))} | "
        "awk -F'\\t' -v OFS='\\t' "
        f"-v rsfile={shlex.quote(str(rsfile))} "
        "'BEGIN{while((getline line < rsfile)>0){want[line]=1}; close(rsfile)} "
        "NR==FNR{genes[$1]=1; next} FNR==1{next} "
        "(($2 in want) && ($9 in genes)){print $9,$2,$1,$5,$6,$7,$13}' "
        f"{shlex.quote(str(genefile))} -"
    )
    rows = []
    for line in run_text(cmd).splitlines():
        gene, rsid, p_eqtl, z, assessed, other, n = line.split("\t")
        rows.append(
            {
                "gene": gene,
                "rsid": rsid,
                "p_eqtl": float(p_eqtl),
                "z": float(z),
                "assessed": assessed.upper(),
                "other": other.upper(),
                "n": int(float(n)),
            }
        )
    return pd.DataFrame(rows)


def add_eqtl_beta(eqtl: pd.DataFrame, af: dict[str, tuple[str, float]]) -> pd.DataFrame:
    rows = []
    for r in eqtl.itertuples(index=False):
        if r.rsid not in af:
            continue
        allele_b, freq_b = af[r.rsid]
        if allele_b == r.assessed:
            eaf = freq_b
        elif allele_b == r.other:
            eaf = 1 - freq_b
        else:
            continue
        denom = 2 * eaf * (1 - eaf) * (r.n + r.z**2)
        if denom <= 0:
            continue
        se = 1 / np.sqrt(denom)
        beta = r.z * se
        rows.append({**r._asdict(), "eaf": eaf, "beta_eqtl": beta, "se_eqtl": se})
    return pd.DataFrame(rows)


def align_beta(beta: float, ea: str, oa: str, assessed: str, other: str) -> float | None:
    ea, oa, assessed, other = ea.upper(), oa.upper(), assessed.upper(), other.upper()
    if ea == assessed and oa == other:
        return beta
    if ea == other and oa == assessed:
        return -beta
    return None


def ivw(beta_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray) -> tuple[float, float, float, float, float]:
    w = 1 / np.square(se_y)
    denom = np.sum(w * np.square(beta_x))
    theta = float(np.sum(w * beta_x * beta_y) / denom)
    se = float(np.sqrt(1 / denom))
    z = theta / se
    p = float(2 * stats.norm.sf(abs(z)))
    q = float(np.sum(w * np.square(beta_y - theta * beta_x)))
    q_p = float(stats.chi2.sf(q, len(beta_x) - 1)) if len(beta_x) > 1 else np.nan
    return theta, se, p, q, q_p


def main() -> None:
    CFG.out_dir.mkdir(parents=True, exist_ok=True)
    genes = target_genes()
    print(f"reverse MR target genes ({len(genes)}): {', '.join(genes)}")

    inst = []
    for outcome, path in CFG.gwas.items():
        d = read_disease_instruments(outcome, path)
        print(f"{outcome}: {len(d)} distance-clumped disease instruments")
        inst.append(d)
    disease_inst = pd.concat(inst, ignore_index=True)
    disease_inst.to_csv(CFG.out_dir / "reverse_mr_disease_instruments.tsv", sep="\t", index=False)

    rsids = set(disease_inst["rsid"])
    af = lookup_af(rsids)
    eqtl = add_eqtl_beta(lookup_eqtl(rsids, genes), af)
    eqtl.to_csv(CFG.out_dir / "reverse_mr_eqtl_outcomes.tsv", sep="\t", index=False)
    print(f"eQTL outcome rows after AF alignment: {len(eqtl)}")

    detail_rows = []
    result_rows = []
    for outcome in CFG.gwas:
        x = disease_inst[disease_inst["outcome"] == outcome]
        for gene in genes:
            m = x.merge(eqtl[eqtl["gene"] == gene], on="rsid", how="inner")
            rows = []
            for r in m.itertuples(index=False):
                bx = align_beta(r.beta_exp, r.ea, r.oa, r.assessed, r.other)
                if bx is None or bx == 0 or r.se_eqtl <= 0:
                    continue
                rows.append(
                    {
                        "exposure": outcome,
                        "gene": gene,
                        "rsid": r.rsid,
                        "beta_disease": bx,
                        "se_disease": r.se_exp,
                        "p_disease": r.p,
                        "beta_eqtl": r.beta_eqtl,
                        "se_eqtl": r.se_eqtl,
                        "p_eqtl": r.p_eqtl,
                        "assessed": r.assessed,
                        "other": r.other,
                        "F_disease": (bx / r.se_exp) ** 2,
                    }
                )
            detail_rows.extend(rows)
            if not rows:
                result_rows.append(
                    {
                        "exposure": outcome,
                        "gene": gene,
                        "nsnp": 0,
                        "theta_expr_per_logodds_disease": np.nan,
                        "se": np.nan,
                        "p_mr": np.nan,
                        "q": np.nan,
                        "q_p": np.nan,
                        "mean_F": np.nan,
                        "status": "no overlapping cis-eQTL outcome after allele harmonisation",
                    }
                )
                continue
            d = pd.DataFrame(rows)
            theta, se, p, q, q_p = ivw(d["beta_disease"].to_numpy(), d["beta_eqtl"].to_numpy(), d["se_eqtl"].to_numpy())
            result_rows.append(
                {
                    "exposure": outcome,
                    "gene": gene,
                    "nsnp": len(d),
                    "theta_expr_per_logodds_disease": theta,
                    "se": se,
                    "p_mr": p,
                    "q": q,
                    "q_p": q_p,
                    "mean_F": float(d["F_disease"].mean()),
                    "status": "ok",
                }
            )

    details = pd.DataFrame(detail_rows)
    if not details.empty:
        details.to_csv(CFG.out_dir / "reverse_mr_instruments.tsv", sep="\t", index=False)

    res = pd.DataFrame(result_rows)
    ok = res["p_mr"].notna()
    res["fdr"] = np.nan
    if ok.any():
        res.loc[ok, "fdr"] = multipletests(res.loc[ok, "p_mr"], method="fdr_bh")[1]
    res = res.sort_values(["exposure", "p_mr", "gene"], na_position="last")
    res.to_csv(CFG.out_dir / "reverse_mr.tsv", sep="\t", index=False)
    print(res[["exposure", "gene", "nsnp", "theta_expr_per_logodds_disease", "p_mr", "fdr", "status"]].to_string(index=False))


if __name__ == "__main__":
    main()
