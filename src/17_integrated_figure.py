"""Integrated causal-genetics figure: MR forest (IBD), coloc PP4, cross-outcome
(IBD/CD/UC) map, and genetics x transcriptomics triangulation.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MR = "/Users/ujs/Downloads/lzy/outputs/mr"
mr = {oc: pd.read_csv(f"{MR}/mr_{oc}.tsv", sep="\t") for oc in ["IBD", "CD", "UC"]
      if os.path.exists(f"{MR}/mr_{oc}.tsv")}
coloc = pd.read_csv(f"{MR}/coloc_IBD.tsv", sep="\t") if os.path.exists(f"{MR}/coloc_IBD.tsv") else pd.DataFrame()
tri = pd.read_csv(f"{MR}/triangulation.tsv", sep="\t") if os.path.exists(f"{MR}/triangulation.tsv") else pd.DataFrame()

fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# (A) MR forest IBD, FDR<0.1
ax = axes[0, 0]
r = mr["IBD"][mr["IBD"].fdr < 0.10].sort_values("OR")
y = np.arange(len(r))
lo = np.exp(np.log(r.OR) - 1.96 * r.se); hi = np.exp(np.log(r.OR) + 1.96 * r.se)
cols = ["#c0392b" if o > 1 else "#2471a3" for o in r.OR]
ax.errorbar(r.OR, y, xerr=[r.OR - lo, hi - r.OR], fmt="none", ecolor="grey", zorder=1)
ax.scatter(r.OR, y, c=cols, s=45, zorder=2)
ax.axvline(1, ls="--", color="k", lw=.8)
ax.set_yticks(y); ax.set_yticklabels(r.gene)
ax.set_xlabel("OR for IBD per +1 SD expression")
ax.set_title("(A) Cis-MR: causal senescence genes -> IBD (FDR<0.1)\nred=risk, blue=protective")

# (B) coloc PP4
ax = axes[0, 1]
if len(coloc):
    c = coloc.sort_values("PP4").tail(15)
    cc = ["#27ae60" if v > 0.75 else "#f39c12" if v > 0.5 else "#bbbbbb" for v in c.PP4]
    ax.barh(range(len(c)), c.PP4, color=cc)
    ax.set_yticks(range(len(c))); ax.set_yticklabels(c.gene)
    ax.axvline(0.75, ls="--", color="green", lw=.8)
    ax.set_xlabel("PP.H4 (shared causal variant)")
    ax.set_title("(B) Colocalization\ngreen=colocalized (PP4>0.75)")
else:
    ax.text(.5, .5, "coloc pending", ha="center")

# (C) cross-outcome map (signed -log10p, FDR<0.1 in any)
ax = axes[1, 0]
genes = set()
for oc, d in mr.items():
    genes |= set(d[d.fdr < 0.10].gene)
genes = sorted(genes)
M = np.full((len(genes), 3), np.nan)
for j, oc in enumerate(["IBD", "CD", "UC"]):
    if oc in mr:
        dd = mr[oc].set_index("gene")
        for i, g in enumerate(genes):
            if g in dd.index:
                M[i, j] = np.sign(np.log(dd.loc[g, "OR"])) * -np.log10(dd.loc[g, "p_mr"])
im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-6, vmax=6)
ax.set_xticks(range(3)); ax.set_xticklabels(["IBD", "CD", "UC"])
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=8)
ax.set_title("(C) Causal across subtypes\nsigned -log10 p (red=risk)")
fig.colorbar(im, ax=ax, fraction=0.046)

# (D) triangulation: MR direction vs mucosal fold-change
ax = axes[1, 1]
if len(tri) and "GSE16879_FC" in tri.columns:
    t = tri[tri.MR_p < 0.05].dropna(subset=["GSE16879_FC"]).copy()
    t["mr_dir"] = np.log(t.MR_OR)
    ax.scatter(t.mr_dir, t.GSE16879_FC, s=40, c=["#c0392b" if o > 1 else "#2471a3" for o in t.MR_OR])
    for _, rr in t.iterrows():
        ax.annotate(rr.gene, (rr.mr_dir, rr.GSE16879_FC), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    ax.axhline(0, ls="--", color="grey", lw=.6); ax.axvline(0, ls="--", color="grey", lw=.6)
    ax.set_xlabel("MR causal direction  ln(OR)  (>0 risk)")
    ax.set_ylabel("Mucosal log2FC (IBD vs control)")
    ax.set_title("(D) Genetics x transcriptomics\ntriangulation")
else:
    ax.text(.5, .5, "triangulation pending", ha="center")

fig.suptitle("Causal-genetic dissection of cellular senescence in IBD (MR + colocalization)",
             fontweight="bold", fontsize=13)
fig.tight_layout()
fig.savefig(f"{MR}/Fig_CAUSAL_integrated.png", dpi=170)
print("saved results/mr/Fig_CAUSAL_integrated.png")
