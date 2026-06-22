"""Shared helpers for the refractory-module causal pipeline."""
from __future__ import annotations

import shlex
import subprocess
from pathlib import Path

import pandas as pd


def gwas_lookup(path: Path, variants: pd.DataFrame, key_file: Path) -> pd.DataFrame:
    keys = variants[["chrom", "pos"]].drop_duplicates().copy()
    keys["key"] = keys["chrom"].astype(str) + ":" + keys["pos"].astype(str)
    key_file.write_text("\n".join(keys["key"].sort_values()) + "\n")

    hdr = subprocess.run(
        f"gunzip -c {shlex.quote(str(path))} | head -1",
        shell=True,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.rstrip("\n").split("\t")
    ix = {c: i + 1 for i, c in enumerate(hdr)}
    required = [
        "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele",
        "hm_beta", "hm_rsid", "standard_error", "p_value",
    ]
    missing = [c for c in required if c not in ix]
    if missing:
        raise ValueError(f"{path} missing columns {missing}")
    cmd = (
        f"gunzip -c {shlex.quote(str(path))} | "
        "awk -F'\\t' -v OFS='\\t' "
        f"-v c={ix['hm_chrom']} -v p={ix['hm_pos']} -v oa={ix['hm_other_allele']} "
        f"-v ea={ix['hm_effect_allele']} -v b={ix['hm_beta']} -v rs={ix['hm_rsid']} "
        f"-v se={ix['standard_error']} -v pv={ix['p_value']} "
        "'NR==FNR{want[$1]; next} FNR==1{next} "
        "{key=$c \":\" $p; if(key in want) print $c,$p,$oa,$ea,$b,$rs,$se,$pv}' "
        f"{shlex.quote(str(key_file))} -"
    )
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    rows = []
    for line in out.splitlines():
        p = line.split("\t")
        if len(p) != 8:
            continue
        chrom, pos, other, effect, beta, rsid, se, pv = p
        try:
            rows.append(
                {
                    "chrom": chrom.replace("chr", ""),
                    "pos": int(pos),
                    "gwas_other": other.upper(),
                    "gwas_effect": effect.upper(),
                    "beta_gwas": float(beta),
                    "rsid": rsid,
                    "se_gwas": float(se),
                    "p_gwas": float(pv),
                }
            )
        except ValueError:
            continue
    return pd.DataFrame(rows)


def harmonised_beta(row: pd.Series) -> float | None:
    if row["gwas_effect"] == row["effect_allele"] and row["gwas_other"] == row["other_allele"]:
        return row["beta_gwas"]
    if row["gwas_effect"] == row["other_allele"] and row["gwas_other"] == row["effect_allele"]:
        return -row["beta_gwas"]
    return None
