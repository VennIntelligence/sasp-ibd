"""MR for GTEx colon eQTL instruments against de Lange IBD/CD/UC GWAS."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, field_validator
from scipy import stats
from statsmodels.stats.multitest import multipletests

from causal_module_utils import gwas_lookup, harmonised_beta
from paths import P


OUTCOMES = {
    "IBD": "IBD.h.tsv.gz",
    "CD": "CD.h.tsv.gz",
    "UC": "UC.h.tsv.gz",
}


class Inputs(BaseModel):
    model_config = ConfigDict(frozen=True)

    out_dir: Path = P.out("causal_module")
    gwas_dir: Path = P.raw / "gwas"

    @field_validator("out_dir")
    @classmethod
    def have_instruments(cls, v: Path) -> Path:
        if not (v / "instruments_gut.tsv").exists():
            raise FileNotFoundError("run src/13_build_instruments.py first")
        return v


def run_outcome(outcome: str, gwas_path: Path, inst: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    gwas = gwas_lookup(gwas_path, inst, out_dir / f"_gwas_keys_{outcome}.txt")
    if not len(gwas):
        raise RuntimeError(f"no GWAS hits for {outcome}")
    merged = inst.merge(gwas, on=["chrom", "pos"], how="left")
    merged["beta_gwas_alt"] = merged.apply(harmonised_beta, axis=1)
    merged = merged.dropna(subset=["beta_gwas_alt", "se_gwas"]).copy()
    merged = merged[(merged["beta_eqtl"] != 0) & (merged["se_gwas"] > 0)]
    theta = merged["beta_gwas_alt"] / merged["beta_eqtl"]
    se = (merged["se_gwas"] / merged["beta_eqtl"]).abs()
    z = theta / se
    res = merged[
        [
            "gene", "tissue", "variant_id", "rsid", "chrom", "pos", "ref", "alt",
            "beta_eqtl", "se_eqtl", "p_eqtl", "beta_gwas_alt", "se_gwas", "p_gwas",
        ]
    ].copy()
    res["theta"] = theta
    res["se"] = se
    res["OR"] = np.exp(theta.clip(-20, 20))
    res["p_mr"] = 2 * stats.norm.sf(np.abs(z))
    res["fdr"] = multipletests(res["p_mr"], method="fdr_bh")[1] if len(res) else []
    res = res.sort_values(["p_mr", "gene", "tissue"])
    res.to_csv(out_dir / f"mr_gut_{outcome}.tsv", sep="\t", index=False)
    print(f"{outcome}: tested {res['gene'].nunique()} genes / {len(res)} gene-tissue rows")
    if len(res):
        print(res.head(12)[["gene", "tissue", "OR", "p_mr", "fdr", "p_gwas"]].to_string(index=False))
    return res


def main() -> None:
    cfg = Inputs()
    inst = pd.read_csv(cfg.out_dir / "instruments_gut.tsv", sep="\t", dtype={"chrom": str})
    if not len(inst):
        raise RuntimeError("no GTEx gut instruments")
    for outcome, name in OUTCOMES.items():
        path = cfg.gwas_dir / name
        if not path.exists():
            raise FileNotFoundError(path)
        run_outcome(outcome, path, inst, cfg.out_dir)


if __name__ == "__main__":
    main()
