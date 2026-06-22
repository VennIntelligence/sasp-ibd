"""Triangulate gut-eQTL MR/coloc results and compare with blood eQTL results."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, field_validator

from paths import P


OUTCOMES = ("IBD", "CD", "UC")


class Inputs(BaseModel):
    model_config = ConfigDict(frozen=True)

    out_dir: Path = P.out("causal_module")

    @field_validator("out_dir")
    @classmethod
    def have_results(cls, v: Path) -> Path:
        needed = [v / f"mr_gut_{o}.tsv" for o in OUTCOMES] + [v / "coloc_gut_IBD.tsv"]
        missing = [p for p in needed if not p.exists()]
        if missing:
            raise FileNotFoundError(f"missing pipeline outputs: {missing}")
        return v


def load_mr(out_dir: Path) -> pd.DataFrame:
    base = pd.read_csv(out_dir / "mr_gut_IBD.tsv", sep="\t")
    base = base.rename(columns={"OR": "MR_OR", "p_mr": "MR_p", "fdr": "MR_FDR"})
    keep = [
        "gene", "tissue", "variant_id", "rsid", "chrom", "pos", "ref", "alt",
        "MR_OR", "MR_p", "MR_FDR", "p_gwas", "p_eqtl",
    ]
    base = base[[c for c in keep if c in base.columns]]
    for outcome in ("CD", "UC"):
        df = pd.read_csv(out_dir / f"mr_gut_{outcome}.tsv", sep="\t")
        df = df[["gene", "tissue", "OR", "p_mr", "fdr"]].rename(
            columns={"OR": f"{outcome}_MR_OR", "p_mr": f"{outcome}_MR_p", "fdr": f"{outcome}_MR_FDR"}
        )
        base = base.merge(df, on=["gene", "tissue"], how="left")
    return base


def blood_status(genes: pd.Series) -> pd.DataFrame:
    blood = pd.DataFrame({"gene": sorted(set(genes))})
    mr_path = P.tables / "mr_IBD.tsv"
    coloc_path = P.tables / "coloc_IBD.tsv"
    if mr_path.exists():
        mr = pd.read_csv(mr_path, sep="\t")
        mr = mr[["gene", "OR", "p_mr", "fdr"]].rename(
            columns={"OR": "blood_MR_OR", "p_mr": "blood_MR_p", "fdr": "blood_MR_FDR"}
        )
        blood = blood.merge(mr, on="gene", how="left")
    if coloc_path.exists():
        coloc = pd.read_csv(coloc_path, sep="\t")
        coloc = coloc[["gene", "PP4"]].rename(columns={"PP4": "blood_coloc_PP4"})
        blood = blood.merge(coloc, on="gene", how="left")
    blood["blood_causal"] = (blood.get("blood_MR_FDR", np.nan) < 0.05) & (blood.get("blood_coloc_PP4", np.nan) > 0.8)
    blood["blood_tested"] = blood.get("blood_MR_p", pd.Series(np.nan, index=blood.index)).notna()
    return blood


def compare_label(gut: bool, blood: bool, blood_tested: bool) -> str:
    if gut and blood:
        return "both_consistent"
    if gut and not blood:
        return "gut_only" if blood_tested else "gut_only_blood_not_tested"
    if blood and not gut:
        return "blood_only"
    return "neither" if blood_tested else "gut_noncausal_blood_not_tested"


def plot_map(df: pd.DataFrame, path: Path) -> None:
    d = df[df["analysis_status"].eq("tested")].sort_values(["causal_call", "MR_p"], ascending=[False, True]).copy()
    if not len(d):
        fig, ax = plt.subplots(figsize=(7, 3))
        ax.text(0.5, 0.5, "No harmonised gut-eQTL MR rows", ha="center", va="center")
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(path, dpi=220)
        plt.close(fig)
        return
    d["neglog10_MR_p"] = -np.log10(d["MR_p"].clip(lower=np.nextafter(0, 1)))
    tissues = {t: i for i, t in enumerate(sorted(d["tissue"].unique()))}
    colors = d["causal_call"].map({True: "#b2182b", False: "#2166ac"})
    fig, ax = plt.subplots(figsize=(10, max(5, 0.28 * d["gene"].nunique() + 1.5)))
    y_labels = [f"{r.gene} | {r.tissue.replace('Colon_', '')}" for r in d.itertuples()]
    y = np.arange(len(d))
    sizes = 35 + 220 * d["coloc_PP4"].fillna(0).clip(0, 1)
    ax.scatter(d["neglog10_MR_p"], y, s=sizes, c=colors, alpha=0.85, edgecolor="black", linewidth=0.35)
    for i, r in enumerate(d.itertuples()):
        if r.MR_OR < 1:
            ax.text(d["neglog10_MR_p"].iloc[i] + 0.04, i, "protect", va="center", fontsize=7, color="#2166ac")
    ax.axvline(-np.log10(0.05), color="grey", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("-log10 IBD MR p-value")
    ax.set_title("Refractory inflammatory module: gut eQTL MR + coloc")
    causal = plt.Line2D([], [], marker="o", linestyle="", color="#b2182b", label="MR FDR<0.05 + PP4>0.8")
    non = plt.Line2D([], [], marker="o", linestyle="", color="#2166ac", label="not causal by rule")
    ax.legend(handles=[causal, non], frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_summary(df: pd.DataFrame, out_dir: Path) -> None:
    tested = df[df["analysis_status"].eq("tested")]
    causal = tested[tested["causal_call"]].sort_values("MR_p")
    gut_only = df[df["blood_compare"].str.startswith("gut_only")]
    no_genes = sorted(df.loc[df["analysis_status"].eq("no_gut_instrument"), "gene"].unique())
    no_harm = sorted(df.loc[df["analysis_status"].eq("no_gwas_harmonized"), "gene"].unique())
    lines = [
        "# Gut eQTL causal-module MR+coloc summary",
        "",
        "Method: GTEx v8 Colon_Sigmoid and Colon_Transverse significant variant-gene pairs were used for lead cis-eQTL MR instruments. GWAS harmonisation used GRCh38 chr/pos/ref/alt from GTEx variant_id against de Lange harmonised hm_* fields. Coloc is ABF restricted to GTEx significant pairs because full allpairs were not anonymously accessible in this environment.",
        "",
        f"Genes with harmonised gut MR rows: {tested['gene'].nunique()}",
        f"Genes with gut instruments but no harmonised GWAS row: {len(no_harm)}" + (f" ({', '.join(no_harm)})" if no_harm else ""),
        f"Genes without gut instruments: {len(no_genes)}" + (f" ({', '.join(no_genes)})" if no_genes else ""),
        f"Causal gene-tissue calls (IBD MR FDR<0.05 and coloc PP4>0.8): {len(causal)}",
    ]
    if len(causal):
        lines.append("")
        lines.append("Causal calls:")
        for r in causal.itertuples():
            lines.append(
                f"- {r.gene} / {r.tissue}: OR={r.MR_OR:.3g}, MR p={r.MR_p:.3g}, FDR={r.MR_FDR:.3g}, PP4={r.coloc_PP4:.3g}, blood_compare={r.blood_compare}"
            )
    else:
        lines.append("")
        lines.append("No module gene passed the strict causal rule under significant-pairs restricted gut coloc.")
    lines.extend(
        [
            "",
            "Blood comparison:",
            f"- gut_only rows: {len(gut_only)}",
            f"- blood_only rows: {(df['blood_compare'] == 'blood_only').sum()}",
            f"- both_consistent rows: {(df['blood_compare'] == 'both_consistent').sum()}",
            f"- neither rows: {(df['blood_compare'] == 'neither').sum()}",
            "",
            "Caveat: PP4 values here are not full-locus coloc PP4; they are restricted to significant GTEx pairs and should be treated as a conservative/approximate screen until allpairs are available.",
        ]
    )
    (out_dir / "SUMMARY.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    cfg = Inputs()
    mr = load_mr(cfg.out_dir)
    mr["analysis_status"] = "tested"
    coloc = pd.read_csv(cfg.out_dir / "coloc_gut_IBD.tsv", sep="\t")
    coloc = coloc[["gene", "tissue", "PP3", "PP4", "method"]].rename(columns={"PP4": "coloc_PP4", "PP3": "coloc_PP3", "method": "coloc_method"})
    df = mr.merge(coloc, on=["gene", "tissue"], how="left")

    inst = pd.read_csv(cfg.out_dir / "instruments_gut.tsv", sep="\t")
    tested_pairs = set(zip(df["gene"], df["tissue"]))
    missing_harm = inst[~inst.set_index(["gene", "tissue"]).index.isin(tested_pairs)].copy()
    if len(missing_harm):
        extra = missing_harm[["gene", "tissue", "variant_id", "chrom", "pos", "ref", "alt", "p_eqtl"]].copy()
        extra["analysis_status"] = "no_gwas_harmonized"
        df = pd.concat([df, extra], ignore_index=True, sort=False)

    no_inst_path = cfg.out_dir / "genes_without_gut_instruments.tsv"
    no_inst = pd.read_csv(no_inst_path, sep="\t") if no_inst_path.exists() else pd.DataFrame(columns=["gene"])
    if len(no_inst):
        extra = no_inst[["gene"]].copy()
        extra["tissue"] = "none"
        extra["analysis_status"] = "no_gut_instrument"
        df = pd.concat([df, extra], ignore_index=True, sort=False)

    blood = blood_status(df["gene"])
    df = df.merge(blood, on="gene", how="left")
    df["causal_call"] = (df["MR_FDR"] < 0.05) & (df["coloc_PP4"] > 0.8)
    df["blood_compare"] = [
        compare_label(bool(g), bool(b), bool(t))
        for g, b, t in zip(df["causal_call"], df["blood_causal"], df["blood_tested"])
    ]
    status_order = {"tested": 0, "no_gwas_harmonized": 1, "no_gut_instrument": 2}
    df["_status_order"] = df["analysis_status"].map(status_order).fillna(9)
    df = df.sort_values(["_status_order", "causal_call", "MR_p", "gene", "tissue"], ascending=[True, False, True, True, True])
    df = df.drop(columns=["_status_order"])
    df.to_csv(cfg.out_dir / "module_causal_map.tsv", sep="\t", index=False)

    plot_map(df, cfg.out_dir / "Fig_module_causal_map.png")
    write_summary(df, cfg.out_dir)

    P.promote_table(cfg.out_dir / "module_causal_map.tsv")
    P.promote_table(cfg.out_dir / "mr_gut_IBD.tsv")
    P.promote_table(cfg.out_dir / "mr_gut_CD.tsv")
    P.promote_table(cfg.out_dir / "mr_gut_UC.tsv")
    P.promote_table(cfg.out_dir / "coloc_gut_IBD.tsv")
    P.promote_figure(cfg.out_dir / "Fig_module_causal_map.png")

    print("=== module causal map ===")
    cols = ["gene", "tissue", "analysis_status", "MR_OR", "MR_p", "MR_FDR", "coloc_PP4", "causal_call", "blood_compare"]
    print(df[cols].head(40).to_string(index=False))


if __name__ == "__main__":
    main()
