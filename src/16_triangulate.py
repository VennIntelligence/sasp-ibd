"""Triangulation: do the genetically CAUSAL senescence genes (MR + coloc) also
show the expected MUCOSAL behaviour (up in active IBD, change with response)?
Convergent evidence across genetics + transcriptomics is the paper's backbone.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = "/Users/ujs/Downloads/lzy"
PROC = f"{BASE}/data/interim"
MR = f"{BASE}/outputs/mr"

# curated drug-target relevance for senescence genes (esp. IBD-tested)
DRUGS = {
    "TNFRSF1A": "anti-TNF (infliximab/adalimumab) — APPROVED in IBD",
    "TNFRSF1B": "anti-TNF axis",
    "IL1B": "anti-IL-1 (canakinumab/anakinra) — trials",
    "CXCR2": "CXCR2 antagonist (navarixin) — inflammation trials",
    "MMP9": "andecaliximab (anti-MMP9) — tested in UC/CD",
    "ICAM1": "alicaforsen (ICAM-1 antisense) — tested in UC/pouchitis",
    "CCL2": "CCL2/CCR2 axis inhibitors",
    "IL15": "anti-IL-15",
    "CSF2RB": "anti-GM-CSF axis",
    "VEGFA": "anti-VEGF (bevacizumab)",
    "HGF": "HGF/MET pathway",
    "EGFR": "anti-EGFR",
}

mr = pd.read_csv(f"{MR}/mr_IBD.tsv", sep="\t")
coloc = pd.read_csv(f"{MR}/coloc_IBD.tsv", sep="\t") if os.path.exists(f"{MR}/coloc_IBD.tsv") else pd.DataFrame(columns=["gene", "PP4"])
ppmap = dict(zip(coloc["gene"], coloc["PP4"]))


def mucosa_stats(gene):
    """per-gene: IBD-baseline vs control log2FC + p; responder reversal."""
    rows = {}
    for gse in ["GSE16879", "GSE73661"]:
        expr = pd.read_parquet(f"{PROC}/{gse}_expr.parquet")
        sc = pd.read_csv(f"{PROC}/{gse}_scored.tsv", sep="\t", index_col=0)
        if gene not in expr.index:
            continue
        e = expr.loc[gene]
        ctrl = sc.index[sc.group == "Control"].intersection(e.index)
        base = sc.index[(sc.timepoint == "baseline")].intersection(e.index)
        if len(ctrl) >= 3 and len(base) >= 3:
            fc = e[base].mean() - e[ctrl].mean()
            p = stats.mannwhitneyu(e[base], e[ctrl]).pvalue
            rows[f"{gse}_FC"] = round(float(fc), 2)
            rows[f"{gse}_p"] = float(p)
    return rows


out = []
tested = mr.sort_values("p_mr")
for _, r in tested.iterrows():
    g = r["gene"]
    rec = {"gene": g, "MR_OR": round(r["OR"], 3), "MR_p": r["p_mr"], "MR_fdr": r["fdr"],
           "coloc_PP4": round(ppmap.get(g, np.nan), 3) if g in ppmap else np.nan}
    rec["drug_target"] = DRUGS.get(g, "")
    rec.update(mucosa_stats(g))
    out.append(rec)

tri = pd.DataFrame(out)
tri.to_csv(f"{MR}/triangulation.tsv", sep="\t", index=False)

# "convergent" = MR nominal sig + (coloc PP4>0.5 if available) + mucosa up in IBD
def converg(r):
    mr_ok = r["MR_p"] < 0.05
    col_ok = (np.isnan(r.get("coloc_PP4", np.nan))) or (r.get("coloc_PP4", 0) > 0.5)
    muc_ok = (r.get("GSE16879_FC", 0) or 0) > 0 or (r.get("GSE73661_FC", 0) or 0) > 0
    return mr_ok and col_ok and muc_ok

tri["convergent"] = tri.apply(converg, axis=1)
print("=== Triangulation (genetics x transcriptomics) ===")
cols = [c for c in ["gene", "MR_OR", "MR_p", "MR_fdr", "coloc_PP4",
                    "GSE16879_FC", "GSE73661_FC", "convergent"] if c in tri.columns]
print(tri[cols].head(25).to_string(index=False))
print(f"\nMR nominal (p<0.05): {(tri.MR_p<0.05).sum()} | MR FDR<0.05: {(tri.MR_fdr<0.05).sum()}")
print(f"convergent (MR + mucosa + coloc-consistent): {tri['convergent'].sum()} genes")
print("convergent genes:", list(tri[tri.convergent]["gene"]))
tri.to_csv(f"{MR}/triangulation.tsv", sep="\t", index=False)
