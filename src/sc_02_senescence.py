"""Task 2 senescence scoring: arrest arm vs SASP arm by cell type."""
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
from sc_task2_utils import (
    HEALTH_ORDER,
    add_gene_scores,
    append_status,
    load_senmayo_symbols,
    normalize_copy,
    sparse_col_frame,
    write_tsv,
)


OUT = P.out("sc_02")
SC01 = P.out("sc_01")

ARREST_UP = ["CDKN2A", "CDKN1A", "CDKN2B", "GLB1", "SERPINE1"]
ARREST_DOWN = ["LMNB1", "MKI67"]
PROLIF = ["MKI67", "TOP2A", "PCNA", "MCM5", "MCM6", "TYMS"]


class ScoreConfig(BaseModel):
    model_config = ConfigDict(frozen=True)

    min_cells_per_group: int = Field(default=20, ge=1)
    true_z: float = 0.5


def input_paths(use_subsample: bool) -> list[Path]:
    names = (
        ("smillie_uc_subsample.h5ad", "martin_cd_subsample.h5ad")
        if use_subsample
        else ("smillie_uc_qc.h5ad", "martin_cd_qc.h5ad")
    )
    return [SC01 / x for x in names if (SC01 / x).exists()]


def add_scores(adata: ad.AnnData) -> tuple[ad.AnnData, dict[str, int]]:
    work = normalize_copy(adata)
    counts = {
        "senmayo": add_gene_scores(work, load_senmayo_symbols(), "sasp_score"),
        "arrest_up": add_gene_scores(work, ARREST_UP, "arrest_up_score"),
        "arrest_down": add_gene_scores(work, ARREST_DOWN, "arrest_down_score"),
        "proliferation": add_gene_scores(work, PROLIF, "proliferation_score"),
    }
    work.obs["arrest_score"] = work.obs["arrest_up_score"] - work.obs["arrest_down_score"]
    expr = sparse_col_frame(work, ["MKI67", "LMNB1", "CDKN1A", "CDKN2A", "CCL8", "CXCR2"])
    for col in expr:
        work.obs[f"{col}_logexpr"] = expr[col].values
    return work, counts


def z_by_celltype(df: pd.DataFrame, value: str) -> pd.Series:
    g = df.groupby(["dataset", "cell_type"], observed=True)[value]
    mu = g.transform("mean")
    sd = g.transform("std").replace(0, np.nan)
    return ((df[value] - mu) / sd).fillna(0.0)


def summarize(work: ad.AnnData, cfg: ScoreConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    obs = work.obs.copy()
    for col in ["sasp_score", "arrest_score", "proliferation_score", "MKI67_logexpr"]:
        obs[f"{col}_z_celltype"] = z_by_celltype(obs, col)
    obs["true_senescent_candidate"] = (
        (obs["sasp_score_z_celltype"] >= cfg.true_z)
        & (obs["arrest_score_z_celltype"] >= cfg.true_z)
        & (obs["sasp_score"] > 0.0)
        & (obs["arrest_score"] > 0.0)
        & (obs["proliferation_score"] <= 0.0)
        & (obs["MKI67_logexpr"] <= 0.0)
    )
    by_ct = (
        obs.groupby(["dataset", "disease", "cell_type"], observed=True)
        .agg(
            n_cells=("sasp_score", "size"),
            sasp_score_mean=("sasp_score", "mean"),
            arrest_score_mean=("arrest_score", "mean"),
            proliferation_score_mean=("proliferation_score", "mean"),
            true_senescent_fraction=("true_senescent_candidate", "mean"),
            inflamed_fraction=("Health", lambda s: float((s == "Inflamed").mean())),
        )
        .reset_index()
    )
    by_health = (
        obs.groupby(["dataset", "disease", "Health", "cell_type"], observed=True)
        .agg(
            n_cells=("sasp_score", "size"),
            sasp_score_mean=("sasp_score", "mean"),
            arrest_score_mean=("arrest_score", "mean"),
            proliferation_score_mean=("proliferation_score", "mean"),
            MKI67_mean=("MKI67_logexpr", "mean"),
            CCL8_mean=("CCL8_logexpr", "mean"),
            CXCR2_mean=("CXCR2_logexpr", "mean"),
            true_senescent_fraction=("true_senescent_candidate", "mean"),
        )
        .reset_index()
    )
    by_ct = by_ct[by_ct["n_cells"] >= cfg.min_cells_per_group]
    by_health = by_health[by_health["n_cells"] >= cfg.min_cells_per_group]
    return by_ct, by_health


def plot_scores(by_ct: pd.DataFrame, by_health: pd.DataFrame) -> list[Path]:
    paths: list[Path] = []
    top = (
        by_ct.sort_values("n_cells", ascending=False)
        .groupby("dataset", observed=True)
        .head(18)
        .copy()
    )
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    sns.scatterplot(
        data=top,
        x="arrest_score_mean",
        y="sasp_score_mean",
        size="true_senescent_fraction",
        hue="dataset",
        sizes=(25, 260),
        ax=ax,
    )
    for _, r in top.iterrows():
        ax.text(r["arrest_score_mean"], r["sasp_score_mean"], r["cell_type"], fontsize=6)
    ax.axhline(0, color="#999999", lw=0.7)
    ax.axvline(0, color="#999999", lw=0.7)
    ax.set_title("Cell-type senescence arms")
    fig.tight_layout()
    p = OUT / "senescence_celltype_scatter.png"
    fig.savefig(p, dpi=190)
    plt.close(fig)
    paths.append(p)

    h = by_health[by_health["Health"].isin(HEALTH_ORDER)].copy()
    h["Health"] = pd.Categorical(h["Health"], HEALTH_ORDER, ordered=True)
    top_ct = (
        h.groupby("cell_type", observed=True)["n_cells"].sum().sort_values(ascending=False).head(22).index
    )
    h = h[h["cell_type"].isin(top_ct)]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), sharey=True)
    for ax, val, title in [
        (axes[0], "sasp_score_mean", "SASP arm"),
        (axes[1], "arrest_score_mean", "Arrest arm"),
    ]:
        mat = h.pivot_table(index="cell_type", columns="Health", values=val, aggfunc="mean")
        mat = mat.reindex(columns=HEALTH_ORDER)
        sns.heatmap(mat, cmap="vlag", center=0, ax=ax, cbar_kws={"shrink": 0.7})
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.tight_layout()
    p = OUT / "senescence_health_heatmap.png"
    fig.savefig(p, dpi=190)
    plt.close(fig)
    paths.append(p)
    return paths


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--subsample", action="store_true")
    p.add_argument("--min-cells-per-group", type=int, default=20)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = ScoreConfig(min_cells_per_group=args.min_cells_per_group)
    paths = input_paths(args.subsample)
    if not paths:
        raise FileNotFoundError("No sc_01 h5ad files found. Run src/sc_01_qc.py first.")
    scored = []
    gene_counts = {}
    for path in paths:
        a = ad.read_h5ad(path)
        work, counts = add_scores(a)
        scored.append(work)
        gene_counts[path.stem] = counts
        per_cell = work.obs[
            [
                "dataset",
                "disease",
                "Health",
                "Subject",
                "Sample",
                "cell_type",
                "sasp_score",
                "arrest_score",
                "proliferation_score",
                "MKI67_logexpr",
                "LMNB1_logexpr",
                "CDKN1A_logexpr",
                "CDKN2A_logexpr",
                "CCL8_logexpr",
                "CXCR2_logexpr",
            ]
        ].reset_index(names="cell")
        write_tsv(per_cell, OUT / f"{path.stem}_senescence_per_cell.tsv")
    combo = ad.concat(scored, axis=0, join="inner", merge="same")
    by_ct, by_health = summarize(combo, cfg)
    write_tsv(by_ct, OUT / "senescence_per_celltype.tsv")
    write_tsv(by_health, OUT / "senescence_by_health_celltype.tsv")
    figs = plot_scores(by_ct, by_health)
    P.promote_table(OUT / "senescence_per_celltype.tsv")
    append_status(
        "## sc_02 senescence scoring\n"
        f"- inputs: {', '.join(str(p.relative_to(P.root)) for p in paths)}\n"
        f"- gene coverage: {gene_counts}\n"
        f"- outputs: outputs/sc_02/senescence_per_celltype.tsv, "
        f"outputs/sc_02/senescence_by_health_celltype.tsv, {', '.join(p.name for p in figs)}\n"
        "- true_senescent_candidate is defined within each dataset/cell type as SASP z>=0.5, "
        "arrest z>=0.5, positive absolute SASP/arrest scores, non-positive proliferation score, "
        "and no MKI67 expression."
    )
    print("wrote", OUT / "senescence_per_celltype.tsv")


if __name__ == "__main__":
    main()
