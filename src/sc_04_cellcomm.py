"""Task 2 focused SASP ligand-receptor communication and integrated figure."""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from pydantic import BaseModel, ConfigDict, Field

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paths import P
from sc_task2_utils import HEALTH_ORDER, append_status, group_expression, normalize_copy, write_tsv


OUT = P.out("sc_04")
SC01 = P.out("sc_01")
SC02 = P.out("sc_02")
SC03 = P.out("sc_03")

LR_PAIRS = [
    ("CCL8", "CCR1", "monocyte/macrophage recruitment"),
    ("CCL8", "CCR2", "monocyte/macrophage recruitment"),
    ("CCL8", "CCR3", "chemokine recruitment"),
    ("CCL8", "CCR5", "chemokine recruitment"),
    ("CXCL8", "CXCR1", "neutrophil axis"),
    ("CXCL8", "CXCR2", "neutrophil axis"),
    ("CXCL1", "CXCR2", "neutrophil axis"),
    ("CXCL2", "CXCR2", "neutrophil axis"),
    ("CXCL3", "CXCR2", "neutrophil axis"),
    ("IL1B", "IL1R1", "inflammatory SASP"),
    ("IL6", "IL6R", "inflammatory SASP"),
    ("TNF", "TNFRSF1A", "inflammatory SASP"),
    ("TNF", "TNFRSF1B", "inflammatory SASP"),
    ("AREG", "EGFR", "epithelial repair SASP"),
    ("HGF", "MET", "stromal epithelial crosstalk"),
]


class CommConfig(BaseModel):
    model_config = ConfigDict(frozen=True)

    min_sender_cells: int = Field(default=20, ge=1)
    min_receiver_cells: int = Field(default=20, ge=1)
    max_celltypes_per_dataset: int = Field(default=24, ge=4)


def input_paths(use_subsample: bool) -> list[Path]:
    names = (
        ("smillie_uc_subsample.h5ad", "martin_cd_subsample.h5ad")
        if use_subsample
        else ("smillie_uc_qc.h5ad", "martin_cd_qc.h5ad")
    )
    return [SC01 / x for x in names if (SC01 / x).exists()]


def communication_for_dataset(path: Path, cfg: CommConfig) -> pd.DataFrame:
    adata = normalize_copy(ad.read_h5ad(path))
    top_ct = adata.obs["cell_type"].value_counts().head(cfg.max_celltypes_per_dataset).index
    adata = adata[adata.obs["cell_type"].isin(top_ct)].copy()
    genes = sorted({x for pair in LR_PAIRS for x in pair[:2]})
    ge = group_expression(adata, genes, ["dataset", "disease", "Health", "cell_type"])
    rows = []
    for (dataset, disease, health), block in ge.groupby(["dataset", "disease", "Health"], observed=True):
        by_gene = {(r.cell_type, r.gene): r for r in block.itertuples(index=False)}
        celltypes = sorted(block["cell_type"].unique())
        for ligand, receptor, axis in LR_PAIRS:
            for sender in celltypes:
                l = by_gene.get((sender, ligand))
                if l is None or l.n_cells < cfg.min_sender_cells:
                    continue
                for receiver in celltypes:
                    r = by_gene.get((receiver, receptor))
                    if r is None or r.n_cells < cfg.min_receiver_cells:
                        continue
                    score = l.mean_expr * r.mean_expr * l.frac_expr * r.frac_expr
                    rows.append(
                        {
                            "dataset": dataset,
                            "disease": disease,
                            "Health": health,
                            "axis": axis,
                            "ligand": ligand,
                            "receptor": receptor,
                            "sender_cell_type": sender,
                            "receiver_cell_type": receiver,
                            "sender_n_cells": int(l.n_cells),
                            "receiver_n_cells": int(r.n_cells),
                            "ligand_mean": float(l.mean_expr),
                            "receptor_mean": float(r.mean_expr),
                            "ligand_fraction": float(l.frac_expr),
                            "receptor_fraction": float(r.frac_expr),
                            "communication_score": float(score),
                        }
                    )
    return pd.DataFrame(rows)


def plot_cellcomm(comm: pd.DataFrame) -> Path:
    inflamed = comm[comm["Health"].eq("Inflamed")].copy()
    top = inflamed.sort_values("communication_score", ascending=False).head(30)
    top["interaction"] = (
        top["ligand"]
        + "->"
        + top["receptor"]
        + " | "
        + top["sender_cell_type"].str.slice(0, 24)
        + " -> "
        + top["receiver_cell_type"].str.slice(0, 24)
    )
    fig, ax = plt.subplots(figsize=(9.5, 7.5))
    sns.barplot(data=top, y="interaction", x="communication_score", hue="dataset", dodge=False, ax=ax)
    ax.set_title("Top inflamed SASP ligand-receptor axes")
    ax.set_xlabel("focused communication score")
    ax.set_ylabel("")
    fig.tight_layout()
    out = OUT / "cellcomm_SASP_axis.png"
    fig.savefig(out, dpi=190)
    plt.close(fig)
    return out


def plot_integrated() -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10.5))
    sen_path = SC02 / "senescence_per_celltype.tsv"
    health_path = SC02 / "senescence_by_health_celltype.tsv"
    pert_path = SC03 / "insilico_perturbation.tsv"
    comm_path = OUT / "cellcomm_SASP_axis.tsv"

    ax = axes[0, 0]
    if sen_path.exists():
        sen = pd.read_csv(sen_path, sep="\t")
        top = sen.sort_values("true_senescent_fraction", ascending=False).head(15)
        ax.barh(top["dataset"] + " | " + top["cell_type"], top["true_senescent_fraction"], color="#4c78a8")
        ax.invert_yaxis()
        ax.set_xlabel("candidate fraction")
        ax.set_title("Bona fide senescence candidates")
    else:
        ax.text(0.5, 0.5, "sc_02 missing", ha="center")

    ax = axes[0, 1]
    if health_path.exists():
        h = pd.read_csv(health_path, sep="\t")
        h = h[h["Health"].isin(HEALTH_ORDER)]
        h["Health"] = pd.Categorical(h["Health"], HEALTH_ORDER, ordered=True)
        top_ct = h.groupby("cell_type", observed=True)["n_cells"].sum().sort_values(ascending=False).head(16).index
        mat = h[h["cell_type"].isin(top_ct)].pivot_table(
            index="cell_type", columns="Health", values="sasp_score_mean", aggfunc="mean"
        )
        sns.heatmap(mat.reindex(columns=HEALTH_ORDER), cmap="Reds", ax=ax, cbar_kws={"shrink": 0.7})
        ax.set_title("SASP score by inflammation state")
        ax.set_xlabel("")
        ax.set_ylabel("")
    else:
        ax.text(0.5, 0.5, "sc_02 missing", ha="center")

    ax = axes[1, 0]
    if pert_path.exists():
        pert = pd.read_csv(pert_path, sep="\t")
        p = pert.groupby(["dataset", "gene", "perturbation"], observed=True)[
            "projection_toward_inflamed_mean"
        ].mean().reset_index()
        p["label"] = p["dataset"] + " | " + p["gene"] + " " + p["perturbation"]
        ax.barh(np.arange(len(p)), p["projection_toward_inflamed_mean"], color=np.where(p["projection_toward_inflamed_mean"] >= 0, "#d55e00", "#0072b2"))
        ax.set_yticks(np.arange(len(p)))
        ax.set_yticklabels(p["label"], fontsize=8)
        ax.axvline(0, color="#333333", lw=0.8)
        ax.set_title("Geneformer perturbation")
        ax.set_xlabel("toward inflamed centroid")
    else:
        ax.text(0.5, 0.5, "sc_03 missing", ha="center")

    ax = axes[1, 1]
    if comm_path.exists():
        comm = pd.read_csv(comm_path, sep="\t")
        top = comm[comm["Health"].eq("Inflamed")].sort_values("communication_score", ascending=False).head(15)
        labels = top["dataset"] + " | " + top["ligand"] + "->" + top["receptor"]
        ax.barh(labels, top["communication_score"], color="#59a14f")
        ax.invert_yaxis()
        ax.set_title("Inflamed SASP communication")
        ax.set_xlabel("score")
    else:
        ax.text(0.5, 0.5, "sc_04 missing", ha="center")

    fig.suptitle("Task 2 single-cell SASP/senescence and causal-gene perturbation", fontweight="bold")
    fig.tight_layout()
    out = OUT / "Fig_task2_singlecell.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def write_summary(comm: pd.DataFrame) -> Path:
    def md_tsv(df: pd.DataFrame) -> str:
        return "```tsv\n" + df.to_csv(sep="\t", index=False).strip() + "\n```"

    sen = pd.read_csv(SC02 / "senescence_per_celltype.tsv", sep="\t") if (SC02 / "senescence_per_celltype.tsv").exists() else pd.DataFrame()
    sen_health = pd.read_csv(SC02 / "senescence_by_health_celltype.tsv", sep="\t") if (SC02 / "senescence_by_health_celltype.tsv").exists() else pd.DataFrame()
    pert = pd.read_csv(SC03 / "insilico_perturbation.tsv", sep="\t") if (SC03 / "insilico_perturbation.tsv").exists() else pd.DataFrame()
    auc_path = P.tables / "multicohort_auc.tsv"
    tri_path = P.tables / "triangulation.tsv"
    auc = pd.read_csv(auc_path, sep="\t") if auc_path.exists() else pd.DataFrame()
    tri = pd.read_csv(tri_path, sep="\t") if tri_path.exists() else pd.DataFrame()
    lines = [
        "# Task2 single-cell honest summary",
        "",
        "## Bottom line",
        "- IBD gut contains small-to-moderate candidate non-proliferative SASP/arrest-high cells, not a dominant global senescent compartment.",
        "- The signal is most interpretable in epithelial/stromal/endothelial axes: Smillie UC highlights M cells, Goblet/immature enterocytes and fibroblast-like WNT2B/WNT5B groups; Martin CD highlights inflamed fibroblast, endothelial, myeloid and epithelial groups.",
        "- This supports the bulk conclusion that most mucosal SenMayo/SASP signal is inflammatory-secretory, while a minority of cells meet a stricter arrest-plus-SASP definition.",
    ]
    if not sen.empty:
        top = sen.sort_values("true_senescent_fraction", ascending=False).head(8)
        if not sen_health.empty:
            inflamed_top = (
                sen_health[sen_health["Health"].eq("Inflamed")]
                .sort_values("true_senescent_fraction", ascending=False)
                .head(8)
            )
        else:
            inflamed_top = pd.DataFrame()
        lines += [
            "",
            "## Bona fide senescence",
            "Candidate cells require SASP-high, arrest-high, positive absolute SASP/arrest scores, non-positive proliferation score, and no MKI67 expression.",
            "",
            md_tsv(top[["dataset", "cell_type", "n_cells", "true_senescent_fraction", "sasp_score_mean", "arrest_score_mean"]]),
            "",
        ]
        if not inflamed_top.empty:
            lines += [
                "Top inflamed-state groups:",
                "",
                md_tsv(
                    inflamed_top[
                        [
                            "dataset",
                            "Health",
                            "cell_type",
                            "n_cells",
                            "true_senescent_fraction",
                            "sasp_score_mean",
                            "arrest_score_mean",
                            "proliferation_score_mean",
                            "MKI67_mean",
                        ]
                    ]
                ),
                "",
            ]
    if not pert.empty:
        g = pert.groupby(["gene", "perturbation"], observed=True)["projection_toward_inflamed_mean"].mean().reset_index()
        gd = (
            pert.groupby(["dataset", "gene", "perturbation"], observed=True)[
                "projection_toward_inflamed_mean"
            ]
            .mean()
            .reset_index()
        )
        def effect(gene: str, mode: str) -> float | None:
            vals = g[(g["gene"].eq(gene)) & (g["perturbation"].eq(mode))][
                "projection_toward_inflamed_mean"
            ]
            return None if vals.empty else float(vals.iloc[0])

        ccl8_del = effect("CCL8", "delete")
        ccl8_oe = effect("CCL8", "overexpress")
        cxcr2_del = effect("CXCR2", "delete")
        cxcr2_oe = effect("CXCR2", "overexpress")
        lines += [
            "## Geneformer perturbation",
            "Positive projection means movement toward the inflamed centroid; negative means away from inflammation.",
            "",
            f"- CCL8 overexpression moves cells toward inflammation ({ccl8_oe:.4g}); CCL8 deletion is near-null overall ({ccl8_del:.4g}) and differs by dataset, so the perturbation gives partial rather than clean support for CCL8 risk biology.",
            f"- CXCR2 deletion moves toward inflammation ({cxcr2_del:.4g}) while CXCR2 overexpression moves away ({cxcr2_oe:.4g}), matching the genetic-protective interpretation and arguing against naive CXCR2 blockade.",
            "",
            md_tsv(g),
            "",
            "Dataset-specific perturbation means:",
            "",
            md_tsv(gd),
            "",
        ]
    if not auc.empty or not tri.empty:
        lines += ["## Triangulation", ""]
        if not auc.empty:
            analysis_col = "analysis" if "analysis" in auc.columns else "model"
            pooled = auc[(auc[analysis_col].eq("pooled_random_effects")) & (auc["score"].eq("senmayo"))]
            if not pooled.empty:
                r = pooled.iloc[0]
                lines.append(
                    f"- Bulk mucosal SenMayo predicts biologic non-response across cohorts "
                    f"(pooled AUC {r['auc']:.3f}; direction={r['direction']})."
                )
        if not tri.empty:
            keep = tri[tri["gene"].isin(["CCL8", "CXCR2"])].copy()
            if not keep.empty:
                tri_cols = [
                    c
                    for c in [
                        "gene",
                        "MR_OR",
                        "MR_p",
                        "coloc_PP4",
                        "GSE16879_FC",
                        "GSE73661_FC",
                        "triangulated",
                        "convergent",
                    ]
                    if c in keep.columns
                ]
                lines += [
                    "- Genetics/transcriptomics prior used for interpretation:",
                    "",
                    md_tsv(keep[tri_cols]),
                    "",
                ]
        lines.append(
            "- The single-cell datasets used here do not contain biologic response labels, so response triangulation is borrowed from the existing bulk-response results rather than estimated de novo at single-cell level."
        )
        lines.append("")
    if not comm.empty:
        topc = comm[comm["Health"].eq("Inflamed")].sort_values("communication_score", ascending=False).head(10)
        ccl8 = (
            comm[comm["Health"].eq("Inflamed") & comm["ligand"].eq("CCL8")]
            .sort_values("communication_score", ascending=False)
            .head(8)
        )
        cxcr2 = (
            comm[comm["Health"].eq("Inflamed") & comm["receptor"].eq("CXCR2")]
            .sort_values("communication_score", ascending=False)
            .head(8)
        )
        lines += [
            "## SASP communication",
            "The focused ligand-receptor score is dominated by IL1B/TNF inflammatory SASP axes, with CCL8->CCR and CXCL chemokine->CXCR2 axes present but not always top-ranked.",
            "",
            md_tsv(topc[["dataset", "ligand", "receptor", "sender_cell_type", "receiver_cell_type", "communication_score"]]),
            "",
        ]
        if not ccl8.empty:
            lines += ["Top inflamed CCL8->CCR rows:", "", md_tsv(ccl8[["dataset", "ligand", "receptor", "sender_cell_type", "receiver_cell_type", "communication_score"]]), ""]
        if not cxcr2.empty:
            lines += ["Top inflamed CXCR2-receptor rows:", "", md_tsv(cxcr2[["dataset", "ligand", "receptor", "sender_cell_type", "receiver_cell_type", "communication_score"]]), ""]
    lines += [
        "## Caveats",
        "- Martin CD cell types are broad marker-derived labels because the local tar contains raw 10x matrices but no author per-cell annotation.",
        "- sc_04 uses a focused, auditable ligand-receptor score for SASP axes. CellPhoneDB was kept out of the main venv as required; this focused score avoids polluting the analysis environment while preserving the CCL8/CXCR2 biology needed here.",
    ]
    out = OUT / "task2_honest_summary.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--subsample", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = CommConfig()
    paths = input_paths(args.subsample)
    if not paths:
        raise FileNotFoundError("No sc_01 h5ad files found. Run src/sc_01_qc.py first.")
    comm = pd.concat([communication_for_dataset(p, cfg) for p in paths], ignore_index=True)
    write_tsv(comm, OUT / "cellcomm_SASP_axis.tsv")
    fig_comm = plot_cellcomm(comm)
    fig_main = plot_integrated()
    summary = write_summary(comm)
    P.promote_table(OUT / "cellcomm_SASP_axis.tsv")
    P.promote_figure(fig_main)
    append_status(
        "## sc_04 cell communication and integrated figure\n"
        f"- inputs: {', '.join(str(p.relative_to(P.root)) for p in paths)}\n"
        f"- outputs: {fig_comm.relative_to(P.root)}, {fig_main.relative_to(P.root)}, {summary.relative_to(P.root)}\n"
        "- communication_score = ligand mean * receptor mean * ligand fraction * receptor fraction "
        "within Health x sender x receiver groups for curated SASP ligand-receptor pairs."
    )
    print("wrote", OUT / "cellcomm_SASP_axis.tsv", fig_main)


if __name__ == "__main__":
    main()
