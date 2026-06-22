"""Approximate coloc for gut eQTL module pairs vs de Lange IBD GWAS.

The intended full cis allpairs GTEx v8 files are not anonymously accessible in
this environment. This script therefore performs an ABF coloc calculation on
GTEx significant variant-gene pairs only and labels the method accordingly.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, field_validator
from scipy.special import logsumexp

from causal_module_utils import gwas_lookup, harmonised_beta
from paths import P


W_EQTL = 0.15**2
W_GWAS = 0.2**2
P1, P2, P12 = 1e-4, 1e-4, 1e-5
METHOD = "significant_pairs_approx"


class Inputs(BaseModel):
    model_config = ConfigDict(frozen=True)

    out_dir: Path = P.out("causal_module")
    gwas_ibd: Path = P.raw / "gwas" / "IBD.h.tsv.gz"

    @field_validator("out_dir")
    @classmethod
    def have_pairs(cls, v: Path) -> Path:
        if not (v / "gtex_colon_module_sig_pairs.tsv").exists():
            raise FileNotFoundError("run src/13_build_instruments.py first")
        return v


def labf(beta: np.ndarray, varbeta: np.ndarray, prior_var: float) -> np.ndarray:
    z2 = beta * beta / varbeta
    r = prior_var / (varbeta + prior_var)
    return 0.5 * (np.log1p(-r) + r * z2)


def coloc_one(df: pd.DataFrame) -> dict[str, float | int | str]:
    if len(df) < 2:
        return {"nsnps": len(df), "PP0": np.nan, "PP1": np.nan, "PP2": np.nan, "PP3": np.nan, "PP4": np.nan}
    l1 = labf(df["beta_eqtl"].to_numpy(), np.square(df["se_eqtl"].to_numpy()), W_EQTL)
    l2 = labf(df["beta_gwas_alt"].to_numpy(), np.square(df["se_gwas"].to_numpy()), W_GWAS)
    h1, h2 = logsumexp(l1), logsumexp(l2)
    h4 = logsumexp(l1 + l2)
    s = h1 + h2
    h3 = s + np.log1p(-np.exp(min(h4 - s, -1e-12)))
    logs = np.array([0, np.log(P1) + h1, np.log(P2) + h2, np.log(P1) + np.log(P2) + h3, np.log(P12) + h4])
    pp = np.exp(logs - logsumexp(logs))
    return {"nsnps": len(df), "PP0": pp[0], "PP1": pp[1], "PP2": pp[2], "PP3": pp[3], "PP4": pp[4]}


def main() -> None:
    cfg = Inputs()
    pairs = pd.read_csv(cfg.out_dir / "gtex_colon_module_sig_pairs.tsv", sep="\t", dtype={"chrom": str})
    if not len(pairs):
        raise RuntimeError("no GTEx significant pairs for coloc")

    gwas = gwas_lookup(cfg.gwas_ibd, pairs, cfg.out_dir / "_coloc_gwas_keys_IBD.txt")
    merged = pairs.merge(gwas, on=["chrom", "pos"], how="left")
    merged["beta_gwas_alt"] = merged.apply(harmonised_beta, axis=1)
    merged = merged.dropna(subset=["beta_gwas_alt", "se_gwas"]).copy()
    merged = merged[merged["se_gwas"] > 0]

    rows = []
    for (gene, tissue), g in merged.groupby(["gene", "tissue"], sort=True):
        rec = {"gene": gene, "tissue": tissue, "method": METHOD}
        rec.update(coloc_one(g))
        rows.append(rec)
    res = pd.DataFrame(rows).sort_values(["PP4", "gene", "tissue"], ascending=[False, True, True])
    res.to_csv(cfg.out_dir / "coloc_gut_IBD.tsv", sep="\t", index=False)
    print("=== gut coloc IBD (significant-pairs restricted) ===")
    print(res.head(30).round(4).to_string(index=False))
    print(f"PP4>0.8: {(res['PP4'] > 0.8).sum() if len(res) else 0}")


if __name__ == "__main__":
    main()
