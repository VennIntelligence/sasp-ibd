"""Colocalization (coloc ABF, Giambartolomei 2014) for MR-candidate senescence
genes: does the eQTL signal and the IBD GWAS signal share the SAME causal
variant (PP.H4) rather than two linked variants (PP.H3)? This guards MR against
confounding by LD. Uses eQTLGen full cis (region) + AF + IBD GWAS region.
"""
import gzip, subprocess, os
import numpy as np
import pandas as pd
from scipy.special import logsumexp

BASE = "/Users/ujs/Downloads/lzy"
EG = f"{BASE}/data/raw/eqtlgen"
GWAS_IBD = f"{BASE}/data/raw/gwas/IBD.h.tsv.gz"
W_EQTL, W_GWAS = 0.15**2, 0.2**2
P1, P2, P12 = 1e-4, 1e-4, 1e-5

# ---- candidate genes: union of top MR hits across outcomes ----
cand = set()
for oc in ["IBD", "CD", "UC"]:
    f = f"{BASE}/outputs/mr/mr_{oc}.tsv"
    if os.path.exists(f):
        r = pd.read_csv(f, sep="\t")
        cand |= set(r.sort_values("p_mr").head(20)["gene"])
        cand |= set(r[r["fdr"] < 0.10]["gene"])
print(f"coloc candidate genes: {len(cand)} -> {sorted(cand)}")

# ---- 1. region cis-eQTL for candidates (from prefiltered TSV) ----
print("reading prefiltered candidate cis-eQTL...")
gene_snps = {g: [] for g in cand}
with open(f"{EG}/cis_full_candidates.tsv") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(hdr)}
    gi, si, zi, ni = ix["GeneSymbol"], ix["SNP"], ix["Zscore"], ix["NrSamples"]
    ai, oi = ix["AssessedAllele"], ix["OtherAllele"]
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if p[gi] in cand:
            gene_snps[p[gi]].append((p[si], float(p[zi]), int(p[ni]), p[ai], p[oi]))
for g in list(gene_snps):
    print(f"  {g}: {len(gene_snps[g])} region SNPs")

# ---- 2-3. AF + GWAS lookup via awk hash (single pass, O(1) lookup) ----
need = {s[0] for g in gene_snps for s in gene_snps[g]}
rsfile = f"{BASE}/outputs/mr/_coloc_rsids.txt"
open(rsfile, "w").write("\n".join(need) + "\n")
print(f"need {len(need)} unique region SNPs")

# AF file cols: SNP(1) hg19_chr hg19_pos AlleleA AlleleB(5) ... AlleleB_all(9)
awk_af = (rf"gunzip -c {EG}/snp_af.txt.gz | "
          rf"awk -F'\t' -v OFS='\t' 'NR==FNR{{a[$1];next}} ($1 in a){{print $1,$5,$9}}' {rsfile} -")
af = {}
for line in subprocess.run(awk_af, shell=True, capture_output=True, text=True).stdout.splitlines():
    try:
        rs, ab, fb = line.split("\t"); af[rs] = (ab, float(fb))
    except ValueError:
        pass
print(f"AF for {len(af)}/{len(need)} region SNPs")

# GWAS IBD: locate columns (1-based for awk)
hdr = subprocess.run(f"gunzip -c {GWAS_IBD} | head -1", shell=True,
                     capture_output=True, text=True).stdout.rstrip("\n").split("\t")
gi2 = {c: i + 1 for i, c in enumerate(hdr)}
rc, hb, bb, se = gi2["hm_rsid"], gi2["hm_beta"], gi2["beta"], gi2["standard_error"]
awk_g = (rf"gunzip -c {GWAS_IBD} | awk -F'\t' -v OFS='\t' -v rc={rc} -v hb={hb} -v bb={bb} -v se={se} "
         rf"'NR==FNR{{a[$1];next}} ($rc in a){{print $rc,$hb,$bb,$se}}' {rsfile} -")
gwas = {}
for line in subprocess.run(awk_g, shell=True, capture_output=True, text=True).stdout.splitlines():
    try:
        rs, b_h, b_o, s = line.split("\t")
        b = float(b_h) if b_h != "NA" else float(b_o)
        gwas[rs] = (b, float(s))
    except ValueError:
        pass
print(f"GWAS hits for {len(gwas)}/{len(need)} region SNPs")


def labf(z, V, W):
    r = W / (V + W)
    return 0.5 * (np.log(1 - r) + r * z * z)


# ---- 4. coloc per gene ----
res = []
for g, snps in gene_snps.items():
    l1, l2 = [], []
    for rs, z_e, n_e, aa, oa in snps:
        if rs not in af or rs not in gwas:
            continue
        alleleB, fB = af[rs]
        f = fB if alleleB == aa else (1 - fB) if alleleB == oa else None
        if f is None or f <= 0 or f >= 1:
            continue
        V_e = 1.0 / (2 * n_e * f * (1 - f))      # eQTLGen z is beta/se
        b_g, se_g = gwas[rs]
        if se_g <= 0:
            continue
        z_g, V_g = b_g / se_g, se_g**2
        l1.append(labf(z_e, V_e, W_EQTL))
        l2.append(labf(z_g, V_g, W_GWAS))
    if len(l1) < 5:
        continue
    l1, l2 = np.array(l1), np.array(l2)
    h1, h2 = logsumexp(l1), logsumexp(l2)
    h4 = logsumexp(l1 + l2)
    # h3 = log(exp(h1+h2) - exp(h4)) via stable log-diff-exp (avoids overflow)
    s = h1 + h2
    h3 = s + np.log1p(-np.exp(min(h4 - s, -1e-12)))
    logs = np.array([0, np.log(P1) + h1, np.log(P2) + h2,
                     np.log(P1) + np.log(P2) + h3, np.log(P12) + h4])
    pp = np.exp(logs - logsumexp(logs))
    res.append({"gene": g, "nsnps": len(l1),
                "PP0": pp[0], "PP1": pp[1], "PP2": pp[2], "PP3": pp[3], "PP4": pp[4]})

cr = pd.DataFrame(res).sort_values("PP4", ascending=False)
cr.to_csv(f"{BASE}/outputs/mr/coloc_IBD.tsv", sep="\t", index=False)
print("\n=== coloc results (PP.H4 = shared causal variant) ===")
print(cr.round(3).to_string(index=False))
print(f"\ncolocalized (PP4>0.75): {(cr['PP4']>0.75).sum()} genes")
