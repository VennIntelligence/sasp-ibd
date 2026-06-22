"""Fig 2: IBD mucosa shows accelerated transcriptomic aging, scaling with
disease activity.
- Baseline IBD vs control in both cohorts (UC/CD split for GSE16879).
- GSE73661: age acceleration vs endoscopic Mayo subscore (activity gradient).
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROC = "/Users/ujs/Downloads/lzy/data/interim"
RES = "/Users/ujs/Downloads/lzy/outputs"
out = {}

fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))

# (a) GSE16879 baseline by group
sc = pd.read_csv(f"{PROC}/GSE16879_scored.tsv", sep="\t", index_col=0)
base = sc[sc["timepoint"].isin(["baseline", "control"])]
order = ["Control", "UC", "CD"]
data = [base[base.group == g]["age_accel"].dropna() for g in order]
ax = axes[0]
ax.boxplot(data, tick_labels=order, showfliers=False)
for i, d in enumerate(data, 1):
    ax.scatter(np.random.normal(i, 0.06, len(d)), d, s=12, alpha=0.5)
ax.axhline(0, ls="--", color="grey", lw=0.8)
ax.set_ylabel("Age acceleration (z vs control)")
p_uc = stats.mannwhitneyu(data[1], data[0], alternative="greater").pvalue
p_cd = stats.mannwhitneyu(data[2], data[0], alternative="greater").pvalue
out["16879_UC_vs_ctrl_p"] = float(p_uc); out["16879_CD_vs_ctrl_p"] = float(p_cd)
ax.set_title(f"GSE16879 baseline\nUC p={p_uc:.1e}, CD p={p_cd:.1e}")

# (b) GSE73661 baseline by group
sc2 = pd.read_csv(f"{PROC}/GSE73661_scored.tsv", sep="\t", index_col=0)
base2 = sc2[sc2["timepoint"].isin(["baseline", "control"])]
order2 = ["Control", "UC"]
data2 = [base2[base2.group == g]["age_accel"].dropna() for g in order2]
ax = axes[1]
ax.boxplot(data2, tick_labels=order2, showfliers=False)
for i, d in enumerate(data2, 1):
    ax.scatter(np.random.normal(i, 0.06, len(d)), d, s=12, alpha=0.5)
ax.axhline(0, ls="--", color="grey", lw=0.8)
p2 = stats.mannwhitneyu(data2[1], data2[0], alternative="greater").pvalue
out["73661_UC_vs_ctrl_p"] = float(p2)
ax.set_ylabel("Age acceleration"); ax.set_title(f"GSE73661 baseline\nUC p={p2:.1e}")

# (c) GSE73661 activity gradient: age accel vs Mayo endoscopic subscore
b3 = sc2[sc2["timepoint"] == "baseline"].dropna(subset=["mayo", "age_accel"])
ax = axes[2]
if len(b3) > 5:
    rho, pr = stats.spearmanr(b3["mayo"], b3["age_accel"])
    out["73661_ageaccel_vs_mayo_rho"] = float(rho); out["73661_ageaccel_vs_mayo_p"] = float(pr)
    grades = sorted(b3["mayo"].dropna().unique())
    dd = [b3[b3.mayo == m]["age_accel"] for m in grades]
    ax.boxplot(dd, tick_labels=[int(m) for m in grades], showfliers=False)
    for i, d in enumerate(dd, 1):
        ax.scatter(np.random.normal(i, 0.06, len(d)), d, s=12, alpha=0.5)
    ax.set_xlabel("Mayo endoscopic subscore"); ax.set_ylabel("Age acceleration")
    ax.set_title(f"GSE73661 activity gradient\nSpearman rho={rho:.2f}, p={pr:.1e}")

fig.suptitle("IBD mucosa is transcriptomically 'older' and scales with disease activity",
             fontweight="bold")
fig.tight_layout()
fig.savefig(f"{RES}/Fig2_accelerated_aging.png", dpi=200)
json.dump(out, open(f"{RES}/activity_stats.json", "w"), indent=2)
print("saved results/Fig2_accelerated_aging.png and activity_stats.json")
print(json.dumps(out, indent=2))
