"""Build GTEx colon cis-eQTL instruments for the refractory IBD module.

GTEx v8 significant variant-gene pairs encode variants as
chr_pos_ref_alt_b38. TensorQTL slopes are ALT-dosage effects on expression, so
MR harmonisation later aligns GWAS beta to ALT.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, field_validator

from paths import P


MODULE_GENES = [
    "OSM", "OSMR", "TREM1", "IL13RA2", "CXCR2", "CCL8", "IL11", "IL24",
    "CCL2", "CXCL10", "TNFRSF1A", "FPR1", "FPR2", "G0S2", "MNDA", "SELE",
    "ADGRE2", "FCGR2A", "FCGR3B", "CSF3R", "C5AR1", "TLR2", "S100A8",
    "S100A9", "S100A12", "SRGN", "BCL2A1", "TFPI2", "PROK2", "AQP9",
    "PTGS2", "FFAR2", "SERPINE1", "MMP1", "MMP3", "TNC", "FAP", "IL1B",
]
TISSUES = ("Colon_Sigmoid", "Colon_Transverse")


class Inputs(BaseModel):
    model_config = ConfigDict(frozen=True)

    gtex_dir: Path = P.raw / "gtex_colon"
    out_dir: Path = P.out("causal_module")

    @field_validator("gtex_dir")
    @classmethod
    def have_gtex(cls, v: Path) -> Path:
        missing = [
            v / f"{t}.v8.signif_variant_gene_pairs.txt.gz"
            for t in TISSUES
            if not (v / f"{t}.v8.signif_variant_gene_pairs.txt.gz").exists()
        ]
        if missing:
            raise FileNotFoundError(f"missing GTEx files: {missing}")
        return v


def parse_variant_id(variant_id: str) -> tuple[str, int, str, str]:
    chrom, pos, ref, alt, build = variant_id.split("_", 4)
    if build != "b38":
        raise ValueError(f"unexpected GTEx build in {variant_id}")
    return chrom.replace("chr", ""), int(pos), ref.upper(), alt.upper()


def gene_map(gtex_dir: Path) -> pd.DataFrame:
    rows = []
    for tissue in TISSUES:
        path = gtex_dir / f"{tissue}.v8.egenes.txt.gz"
        with gzip.open(path, "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(head)}
            for line in fh:
                p = line.rstrip("\n").split("\t")
                rows.append((p[ix["gene_id"]], p[ix["gene_name"]], tissue))
    gm = pd.DataFrame(rows, columns=["gene_id", "gene", "tissue"])
    gm["gene_base"] = gm["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    return gm


def module_genes() -> list[str]:
    genes = list(dict.fromkeys(MODULE_GENES))
    nr = P.outputs / "nr_unbiased" / "NR_up_robust.tsv"
    if nr.exists():
        try:
            df = pd.read_csv(nr, sep="\t")
            col = "gene" if "gene" in df.columns else df.columns[0]
            for g in df[col].dropna().astype(str).head(50):
                genes.append(g)
        except Exception as exc:
            print(f"warning: could not read {nr}: {exc}")
    else:
        print(f"note: optional unbiased NR-up table absent: {nr}")
    return sorted(dict.fromkeys(genes))


def read_tissue_pairs(path: Path, gene_ids: set[str], id_to_gene: dict[str, str], tissue: str) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(head)}
        need = [
            "variant_id", "gene_id", "tss_distance", "ma_samples", "ma_count",
            "maf", "pval_nominal", "slope", "slope_se", "pval_beta",
        ]
        missing = [c for c in need if c not in ix]
        if missing:
            raise ValueError(f"{path} missing columns {missing}")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            gid = p[ix["gene_id"]]
            if gid not in gene_ids:
                continue
            chrom, pos, ref, alt = parse_variant_id(p[ix["variant_id"]])
            beta = float(p[ix["slope"]])
            se = float(p[ix["slope_se"]])
            if not np.isfinite(beta) or not np.isfinite(se) or se <= 0:
                continue
            rows.append(
                {
                    "tissue": tissue,
                    "gene": id_to_gene[gid],
                    "gene_id": gid,
                    "variant_id": p[ix["variant_id"]],
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "effect_allele": alt,
                    "other_allele": ref,
                    "eaf": float(p[ix["maf"]]),
                    "tss_distance": int(p[ix["tss_distance"]]),
                    "ma_samples": int(p[ix["ma_samples"]]),
                    "ma_count": int(p[ix["ma_count"]]),
                    "beta_eqtl": beta,
                    "se_eqtl": se,
                    "z_eqtl": beta / se,
                    "p_eqtl": float(p[ix["pval_nominal"]]),
                    "p_beta": float(p[ix["pval_beta"]]),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    cfg = Inputs()
    genes = module_genes()
    gm = gene_map(cfg.gtex_dir)
    target = gm[gm["gene"].isin(genes)].copy()
    id_to_gene = dict(zip(target["gene_id"], target["gene"]))
    gene_ids = set(id_to_gene)
    print(f"module genes requested: {len(genes)}")
    print(f"GTEx-resolvable genes: {target['gene'].nunique()}/{len(genes)}")

    all_pairs = []
    for tissue in TISSUES:
        path = cfg.gtex_dir / f"{tissue}.v8.signif_variant_gene_pairs.txt.gz"
        df = read_tissue_pairs(path, gene_ids, id_to_gene, tissue)
        print(f"{tissue}: {len(df)} significant pairs across {df['gene'].nunique() if len(df) else 0} genes")
        all_pairs.append(df)
    pairs = pd.concat(all_pairs, ignore_index=True) if all_pairs else pd.DataFrame()
    pairs.to_csv(cfg.out_dir / "gtex_colon_module_sig_pairs.tsv", sep="\t", index=False)

    if len(pairs):
        lead = (
            pairs.assign(abs_z=lambda x: x["z_eqtl"].abs())
            .sort_values(["tissue", "gene", "p_eqtl", "abs_z"], ascending=[True, True, True, False])
            .drop_duplicates(["tissue", "gene"])
            .drop(columns=["abs_z"])
            .sort_values(["gene", "tissue"])
        )
    else:
        lead = pd.DataFrame()
    lead.to_csv(cfg.out_dir / "instruments_gut.tsv", sep="\t", index=False)

    have = set(lead["gene"]) if len(lead) else set()
    no_inst = pd.DataFrame({"gene": genes})
    no_inst["gtex_gene_found"] = no_inst["gene"].isin(set(target["gene"]))
    no_inst["has_gut_instrument"] = no_inst["gene"].isin(have)
    no_inst = no_inst[~no_inst["has_gut_instrument"]].sort_values("gene")
    no_inst.to_csv(cfg.out_dir / "genes_without_gut_instruments.tsv", sep="\t", index=False)

    print(f"lead gut instruments: {len(lead)} gene-tissue rows for {len(have)} genes")
    if len(no_inst):
        print("no gut instruments:", ", ".join(no_inst["gene"].tolist()))


if __name__ == "__main__":
    main()
