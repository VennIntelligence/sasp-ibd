"""Single-cell CXCR1/CXCR2/IL8 mechanism localization.

Loads SCP259 Smillie UC atlas (Imm compartment, MTX format).
Computes per-cell-type dotplot statistics for the CXCL8 axis:
  CXCR1, CXCR2 (receptors), IL8/CXCL8, CXCL1, CXCL5 (CXCR2 ligands)
Stratified by Health (Healthy / Non-inflamed / Inflamed).

Run on remote:  cd ~/mycode/sasp-ibd && .venv/bin/python src/51_scrnaseq_cxcr2_mechanism.py
CPU only. Output written to outputs/51_scrnaseq_cxcr2/.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io

# ── paths ──────────────────────────────────────────────────────────────────
REPO = Path(__file__).parent.parent
SCP  = REPO / "data/scrna/SCP259/SCP259"
IMM  = SCP / "expression/5cdc540d328cee7a2efc234a"
META = SCP / "metadata/all.meta2.txt"
OUT  = REPO / "outputs/51_scrnaseq_cxcr2"
OUT.mkdir(parents=True, exist_ok=True)

# ── genes of interest ──────────────────────────────────────────────────────
RECEPTORS = ["CXCR1", "CXCR2"]
LIGANDS   = ["IL8", "CXCL1", "CXCL5"]   # IL8 = CXCL8 (old HGNC name)
GENES_OI  = RECEPTORS + LIGANDS

HEALTH_ORDER = ["Healthy", "Non-inflamed", "Inflamed"]
HEALTH_LABEL = {"Healthy": "Healthy", "Non-inflamed": "Non-inflamed", "Inflamed": "Inflamed (UC)"}

# cell-type order for the Imm compartment (roughly myeloid → lymphoid)
CELLTYPE_ORDER = [
    "Inflammatory Monocytes", "Cycling Monocytes", "Macrophages",
    "DC1", "DC2",
    "CD69+ Mast", "CD69- Mast",
    "ILCs", "NKs",
    "CD4+ Memory", "CD4+ Activated Fos-hi", "CD4+ Activated Fos-lo",
    "CD4+ PD1+", "Tregs",
    "CD8+ LP", "CD8+ IELs", "CD8+ IL17+",
    "Cycling T", "Cycling B",
    "Follicular", "GC", "Plasma",
    "MT-hi",
]


# ── helpers ────────────────────────────────────────────────────────────────

def load_imm() -> tuple[np.ndarray, list[str], list[str]]:
    """Return (dense_submatrix cells×genes, barcodes, gene_names_kept).

    The SCP MTX is genes × cells (standard 10x orientation).
    We read it as sparse, subset to GENES_OI rows, then densify only that slice.
    """
    mtx_path = IMM / "gene_sorted-Imm.matrix.mtx"
    bc_path  = IMM / "Imm.barcodes2.tsv"
    g_path   = IMM / "Imm.genes.tsv"

    print("reading sparse MTX …", flush=True)
    mat = scipy.io.mmread(mtx_path).tocsr()   # genes × cells
    genes   = pd.read_csv(g_path, header=None)[0].tolist()
    barcodes = pd.read_csv(bc_path, header=None)[0].tolist()

    # If mat.shape is (cells × genes) the barcodes length should match rows
    # Determine orientation by comparing sizes
    if mat.shape[0] == len(barcodes) and mat.shape[1] == len(genes):
        # already cells × genes — transpose to genes × cells for consistent slicing
        mat = mat.T.tocsr()
    elif mat.shape[0] == len(genes) and mat.shape[1] == len(barcodes):
        pass  # already genes × cells — correct
    else:
        warnings.warn(f"MTX shape {mat.shape} vs genes {len(genes)} / barcodes {len(barcodes)}")

    gene_map = {g: i for i, g in enumerate(genes)}
    kept_genes = [g for g in GENES_OI if g in gene_map]
    missing    = [g for g in GENES_OI if g not in gene_map]
    if missing:
        print(f"  genes absent from Imm compartment: {missing}", flush=True)

    idx = [gene_map[g] for g in kept_genes]
    # col_sums before subsetting (for library-size normalisation)
    col_sums = np.asarray(mat.sum(axis=0)).flatten()   # length = n_cells
    sub = mat[idx, :].toarray().T.astype(np.float32)   # cells × kept_genes

    # log1p-normalise (library size → 10k)
    col_sums_safe = np.where(col_sums == 0, 1, col_sums).astype(np.float32)
    sub = sub / col_sums_safe[:, None] * 1e4
    sub = np.log1p(sub)

    print(f"  loaded: {sub.shape[0]} cells × {sub.shape[1]} genes", flush=True)
    return sub, barcodes, kept_genes


def load_meta() -> pd.DataFrame:
    meta = pd.read_csv(META, sep="\t", skiprows=[1])  # drop TYPE row
    return meta


def merge_meta(sub: np.ndarray, barcodes: list[str], genes: list[str],
               meta: pd.DataFrame) -> pd.DataFrame:
    df = pd.DataFrame(sub, index=barcodes, columns=genes)
    df.index.name = "NAME"
    merged = df.join(meta.set_index("NAME")[["Cluster", "Health"]])
    return merged.dropna(subset=["Cluster", "Health"])


def compute_stats(merged: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    for (cluster, health), grp in merged.groupby(["Cluster", "Health"]):
        if health not in HEALTH_ORDER:
            continue
        row = {"Cluster": cluster, "Health": health, "n_cells": len(grp)}
        for g in genes:
            row[f"{g}_mean"] = float(grp[g].mean())
            row[f"{g}_pct"]  = float((grp[g] > 0).mean() * 100)
        rows.append(row)
    return pd.DataFrame(rows)


def dotplot(stats: pd.DataFrame, genes: list[str], title: str, path: Path) -> None:
    ct_present = [c for c in CELLTYPE_ORDER if c in stats["Cluster"].values]
    ct_other   = sorted(set(stats["Cluster"]) - set(CELLTYPE_ORDER))
    celltypes  = ct_present + ct_other

    n_ct = len(celltypes)
    n_h  = len(HEALTH_ORDER)
    n_g  = len(genes)

    # figure: rows = celltypes, col-groups = health, within-group = genes
    fig_w = max(10, n_g * n_h * 0.9 + 3)
    fig_h = max(6, n_ct * 0.45 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap = plt.get_cmap("Reds")
    all_means = [stats.get(f"{g}_mean", pd.Series([])).dropna().values for g in genes]
    vmax = max(np.percentile(v, 95) if len(v) else 1 for v in all_means)
    vmax = max(vmax, 0.5)

    x_tick_labels = []
    x_ticks = []
    col = 0
    for h in HEALTH_ORDER:
        hs = stats[stats["Health"] == h].set_index("Cluster")
        for g in genes:
            for row_i, ct in enumerate(celltypes):
                y = n_ct - 1 - row_i
                if ct in hs.index:
                    mean_val = hs.loc[ct, f"{g}_mean"]
                    pct_val  = hs.loc[ct, f"{g}_pct"]
                    size = (pct_val / 100) * 220
                    color = cmap(min(mean_val / vmax, 1.0))
                    ax.scatter(col, y, s=size, c=[color], edgecolors="none", zorder=3)
            x_ticks.append(col)
            x_tick_labels.append(f"{g}\n{HEALTH_LABEL[h]}")
            col += 1
        col += 0.6  # gap between health groups

    ax.set_xlim(-0.7, col - 0.3)
    ax.set_ylim(-0.8, n_ct - 0.2)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels, fontsize=7, rotation=45, ha="right")
    ax.set_yticks(range(n_ct))
    ax.set_yticklabels(list(reversed(celltypes)), fontsize=8)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.grid(axis="x", color="#eeeeee", zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

    # colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.35, pad=0.02, aspect=15)
    cbar.set_label("mean log-norm expr", fontsize=8)

    # size legend
    for pct, lbl in [(25, "25%"), (50, "50%"), (75, "75%")]:
        ax.scatter([], [], s=(pct / 100) * 220, c="gray", edgecolors="none",
                   label=f"{lbl} expressing")
    ax.legend(title="% cells > 0", loc="lower right", frameon=False,
              fontsize=7, title_fontsize=7)

    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {path}", flush=True)


def summary_bar(stats: pd.DataFrame, genes: list[str], path: Path) -> None:
    """Mean expression of each gene across Health categories (all immune cells)."""
    rows = []
    all_merged = None
    # recompute per health (mean over all cells in that health category)
    for h in HEALTH_ORDER:
        hs = stats[stats["Health"] == h]
        for g in genes:
            mean_overall = np.average(hs[f"{g}_mean"], weights=hs["n_cells"])
            pct_overall  = np.average(hs[f"{g}_pct"], weights=hs["n_cells"])
            rows.append({"Health": h, "Gene": g, "mean": mean_overall, "pct": pct_overall})
    df = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, col, ylabel in zip(axes, ["mean", "pct"], ["mean log-norm expr", "% cells expressing"]):
        pivot = df.pivot(index="Gene", columns="Health", values=col)[HEALTH_ORDER]
        pivot.plot(kind="bar", ax=ax, color=["#4dac26", "#b8e186", "#d01c8b"],
                   edgecolor="none", width=0.7)
        ax.set_title(ylabel)
        ax.set_xlabel("")
        ax.set_xticklabels(pivot.index, rotation=30, ha="right", fontsize=9)
        ax.legend(fontsize=8, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("CXCL8 axis: expression across IBD health states (Smillie UC atlas immune cells)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {path}", flush=True)


def main() -> None:
    print("=== CXCR1/CXCR2 single-cell mechanism (SCP259 UC atlas) ===", flush=True)

    sub, barcodes, kept_genes = load_imm()
    meta = load_meta()
    merged = merge_meta(sub, barcodes, kept_genes, meta)
    print(f"  merged: {len(merged)} cells with metadata", flush=True)

    print("  computing dotplot stats …", flush=True)
    stats = compute_stats(merged, kept_genes)
    stats.to_csv(OUT / "dotplot_stats.tsv", sep="\t", index=False)

    print("  plotting …", flush=True)
    dotplot(
        stats, kept_genes,
        title="CXCL8 axis receptor/ligand expression — SCP259 Smillie UC atlas (immune compartment)",
        path=OUT / "Fig_cxcr2_dotplot.png",
    )
    summary_bar(stats, kept_genes, OUT / "Fig_cxcr2_health_bar.png")

    # write SUMMARY.md
    receptor_stats = stats[stats["Health"] == "Inflamed"].copy()
    if "CXCR2_mean" in receptor_stats.columns:
        top_cxcr2 = receptor_stats.nlargest(5, "CXCR2_mean")[["Cluster", "CXCR2_mean", "CXCR2_pct", "n_cells"]]
        top_il8   = receptor_stats.nlargest(5, "IL8_mean")[["Cluster", "IL8_mean", "IL8_pct", "n_cells"]] if "IL8_mean" in receptor_stats.columns else None
    else:
        top_cxcr2 = top_il8 = None

    lines = [
        "# CXCR1/CXCR2 axis single-cell mechanism (SCP259 Smillie UC atlas)",
        "",
        "## Source",
        "Smillie et al. 2019 UC atlas (SCP259). Immune (Imm) compartment. 365,493 cells; 3 health states: Healthy, Non-inflamed, Inflamed (UC).",
        "",
        "## Genes assayed",
        f"Receptors: {RECEPTORS}. Ligands (CXCR2-activating): {LIGANDS} (IL8 = CXCL8 by old HGNC nomenclature).",
        "",
        "## Key findings",
    ]
    if top_cxcr2 is not None:
        lines.append("Top CXCR2-expressing cell types (inflamed tissue):")
        for _, r in top_cxcr2.iterrows():
            lines.append(f"  - {r['Cluster']}: mean={r['CXCR2_mean']:.3f}, pct={r['CXCR2_pct']:.1f}%, n={r['n_cells']}")
    if top_il8 is not None:
        lines.append("Top IL8/CXCL8-expressing cell types (inflamed tissue):")
        for _, r in top_il8.iterrows():
            lines.append(f"  - {r['Cluster']}: mean={r['IL8_mean']:.3f}, pct={r['IL8_pct']:.1f}%, n={r['n_cells']}")
    lines += [
        "",
        "## Caveats",
        "- Neutrophils are absent from SCP259 (expected for scRNA-seq of IBD tissue — fragility/preparation).",
        "  The principal CXCR2-expressing cells captured are myeloid (monocytes/macrophages).",
        "- CXCL8/IL8 is named 'IL8' in SCP259 gene matrix (old HGNC; identical gene product).",
        "- Analyses is restricted to Imm compartment. Epithelial IL8 (major CXCL8 source in IBD) not shown here.",
    ]
    (OUT / "SUMMARY.md").write_text("\n".join(lines) + "\n")
    print("done — outputs in", OUT, flush=True)


if __name__ == "__main__":
    main()
