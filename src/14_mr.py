"""Mendelian randomization: senescence-gene expression (eQTLGen instruments)
-> IBD / CD / UC (de Lange 2017 GWAS). Single-instrument Wald ratio per gene,
allele-harmonised by rsid, FDR across genes. Outputs tables + forest plot.
"""
import os
import shlex
import subprocess
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/ujs/Downloads/lzy"
GWAS = {"IBD": f"{BASE}/data/raw/gwas/IBD.h.tsv.gz",
        "CD": f"{BASE}/data/raw/gwas/CD.h.tsv.gz",
        "UC": f"{BASE}/data/raw/gwas/UC.h.tsv.gz"}
inst = pd.read_csv(f"{BASE}/outputs/mr/instruments.tsv", sep="\t")
rsids = set(inst["rsid"])
rsfile = f"{BASE}/outputs/mr/_rsids.txt"
open(rsfile, "w").write("\n".join(sorted(rsids)) + "\n")
print(f"instruments: {len(inst)} genes / {len(rsids)} SNPs")


def gwas_lookup(path):
    """Look up instrument rsids in a harmonised GWAS with an awk hash join."""
    hdr = subprocess.run(f"gunzip -c {path} | head -1", shell=True,
                         capture_output=True, text=True).stdout.rstrip("\n").split("\t")
    i = {c: k for k, c in enumerate(hdr)}
    required = [
        "hm_rsid", "hm_beta", "hm_effect_allele", "hm_other_allele",
        "beta", "effect_allele", "other_allele", "standard_error", "p_value",
    ]
    missing = [c for c in required if c not in i]
    if missing:
        raise ValueError(f"{path} missing required GWAS columns: {missing}")
    c = {name: i[name] + 1 for name in required}
    cmd = (
        f"gunzip -c {shlex.quote(path)} | "
        "awk -F'\\t' -v OFS='\\t' "
        f"-v rs={c['hm_rsid']} -v hb={c['hm_beta']} "
        f"-v hea={c['hm_effect_allele']} -v hoa={c['hm_other_allele']} "
        f"-v b={c['beta']} -v ea={c['effect_allele']} -v oa={c['other_allele']} "
        f"-v se={c['standard_error']} -v pv={c['p_value']} "
        "'NR==FNR{want[$1]; next} FNR==1{next} "
        "($rs in want){print $rs,$hb,$hea,$hoa,$b,$ea,$oa,$se,$pv}' "
        f"{shlex.quote(rsfile)} -"
    )
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    d = {}
    for line in out.splitlines():
        p = line.split("\t")
        if len(p) != 9:
            continue
        rs, hm_beta, hm_ea, hm_oa, beta_raw, ea_raw, oa_raw, se_raw, p_raw = p
        if rs not in rsids:
            continue
        try:
            beta = hm_beta
            ea, oa = hm_ea, hm_oa
            if beta == "NA":
                beta, ea, oa = beta_raw, ea_raw.upper(), oa_raw.upper()
            d[rs] = {"ea": ea.upper(), "oa": oa.upper(),
                     "beta": float(beta), "se": float(se_raw),
                     "p": float(p_raw)}
        except (ValueError, KeyError, IndexError):
            continue
    return d


def harmonise(beta_g, ea, oa, assessed, other):
    """align GWAS beta to the eQTL assessed allele."""
    if ea == assessed and oa == other:
        return beta_g
    if ea == other and oa == assessed:
        return -beta_g
    return None  # allele mismatch / ambiguous


all_res = {}
for outcome, path in GWAS.items():
    if not os.path.exists(path):
        print(f"skip {outcome} (not downloaded)"); continue
    g = gwas_lookup(path)
    rows = []
    for _, r in inst.iterrows():
        hit = g.get(r["rsid"])
        if not hit:
            continue
        bg = harmonise(hit["beta"], hit["ea"], hit["oa"], r["assessed"], r["other"])
        if bg is None:
            continue
        theta = bg / r["beta_eqtl"]                 # effect of +1 expr unit on log-OR
        se = abs(hit["se"] / r["beta_eqtl"])        # Wald-ratio SE (1st order)
        z = theta / se
        rows.append({"gene": r["gene"], "rsid": r["rsid"],
                     "theta": theta, "se": se, "OR": np.exp(theta),
                     "p_mr": 2 * stats.norm.sf(abs(z)),
                     "gwas_p": hit["p"], "eqtl_p": r["p_eqtl"]})
    res = pd.DataFrame(rows)
    if len(res):
        res["fdr"] = multipletests(res["p_mr"], method="fdr_bh")[1]
        res = res.sort_values("p_mr")
        res.to_csv(f"{BASE}/outputs/mr/mr_{outcome}.tsv", sep="\t", index=False)
        all_res[outcome] = res
        nsig = (res.fdr < 0.05).sum()
        print(f"\n{outcome}: tested {len(res)} genes, FDR<0.05: {nsig}")
        print(res.head(12)[["gene", "OR", "p_mr", "fdr", "gwas_p"]].to_string(index=False))

# ---- forest plot of top IBD hits ----
if "IBD" in all_res:
    r = all_res["IBD"].head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(6, 6))
    y = np.arange(len(r))
    lo = np.exp(r["theta"] - 1.96 * r["se"]); hi = np.exp(r["theta"] + 1.96 * r["se"])
    ax.errorbar(r["OR"], y, xerr=[r["OR"] - lo, hi - r["OR"]], fmt="o", color="#c0392b")
    ax.axvline(1, ls="--", color="grey")
    ax.set_yticks(y); ax.set_yticklabels(r["gene"])
    ax.set_xlabel("OR for IBD per +1 SD senescence-gene expression")
    ax.set_title("Cis-MR: senescence genes -> IBD (top 15)")
    fig.tight_layout(); fig.savefig(f"{BASE}/outputs/mr/MR_forest_IBD.png", dpi=200)
    print("\nsaved results/mr/MR_forest_IBD.png")
