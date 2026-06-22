"""CORE NOVEL ANALYSIS: is transcriptomic age acceleration a REVERSIBLE,
treatment-coupled biomarker in IBD?

H1 (baseline): IBD mucosa shows higher age acceleration than controls.
H2 (reversal): in treatment RESPONDERS, age acceleration drops post-treatment
               toward control levels; in NON-responders it does not.
H3 (predict):  baseline age acceleration differs between future R and NR.
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
results = {}


def load(gse):
    return pd.read_csv(f"{PROC}/{gse}_scored.tsv", sep="\t", index_col=0)


def paired(df):
    """Return per-patient (baseline, post, response) for IBD patients."""
    rows = []
    for pid, sub in df[df["group"] != "Control"].groupby("patient"):
        b = sub[sub["timepoint"] == "baseline"]["age_accel"]
        p = sub[sub["timepoint"] == "post"]["age_accel"]
        resp = sub["response"].iloc[0]
        if len(b) and len(p):
            rows.append((pid, b.mean(), p.mean(), resp))
    return pd.DataFrame(rows, columns=["patient", "baseline", "post", "response"])

# ===================== GSE16879 (primary) =====================
d = load("GSE16879")
ctrl = d[d["group"] == "Control"]["age_accel"]
base_ibd = d[d["timepoint"] == "baseline"]["age_accel"]
u1 = stats.mannwhitneyu(base_ibd, ctrl, alternative="greater")
print(f"[16879] H1 baseline IBD vs control: median {base_ibd.median():.2f} vs "
      f"{ctrl.median():.2f}, MWU p={u1.pvalue:.2e}")
results["16879_H1_baseline_vs_control_p"] = float(u1.pvalue)

pr = paired(d)
pr["delta"] = pr["post"] - pr["baseline"]
for resp in ["R", "NR"]:
    g = pr[pr["response"] == resp]
    if len(g) >= 3:
        w = stats.wilcoxon(g["baseline"], g["post"])
        print(f"[16879] H2 {resp} (n={len(g)}): baseline {g['baseline'].median():.2f} "
              f"-> post {g['post'].median():.2f}, delta median {g['delta'].median():.2f}, "
              f"Wilcoxon p={w.pvalue:.3f}")
        results[f"16879_H2_{resp}_paired_p"] = float(w.pvalue)
        results[f"16879_H2_{resp}_delta_median"] = float(g["delta"].median())
# delta R vs NR
if (pr["response"] == "R").sum() >= 3 and (pr["response"] == "NR").sum() >= 3:
    md = stats.mannwhitneyu(pr[pr.response == "R"]["delta"],
                            pr[pr.response == "NR"]["delta"], alternative="less")
    print(f"[16879] H2 delta R vs NR: MWU p={md.pvalue:.3f}")
    results["16879_H2_deltaRvsNR_p"] = float(md.pvalue)
    # H3 baseline R vs NR
    mb = stats.mannwhitneyu(pr[pr.response == "R"]["baseline"],
                            pr[pr.response == "NR"]["baseline"])
    print(f"[16879] H3 baseline R vs NR: median {pr[pr.response=='R']['baseline'].median():.2f} "
          f"vs {pr[pr.response=='NR']['baseline'].median():.2f}, MWU p={mb.pvalue:.3f}")
    results["16879_H3_baselineRvsNR_p"] = float(mb.pvalue)

# ===================== GSE73661 (validation) =====================
d2 = load("GSE73661")
ctrl2 = d2[d2["group"] == "Control"]["age_accel"]
base2 = d2[d2["timepoint"] == "baseline"]["age_accel"]
u2 = stats.mannwhitneyu(base2, ctrl2, alternative="greater")
print(f"\n[73661] H1 baseline UC vs control: median {base2.median():.2f} vs "
      f"{ctrl2.median():.2f}, MWU p={u2.pvalue:.2e}")
results["73661_H1_baseline_vs_control_p"] = float(u2.pvalue)
pr2 = paired(d2)
pr2["delta"] = pr2["post"] - pr2["baseline"]
for resp in ["R", "NR"]:
    g = pr2[pr2["response"] == resp]
    if len(g) >= 3:
        w = stats.wilcoxon(g["baseline"], g["post"])
        print(f"[73661] H2 {resp} (n={len(g)}): baseline {g['baseline'].median():.2f} "
              f"-> post {g['post'].median():.2f}, Wilcoxon p={w.pvalue:.3f}")
        results[f"73661_H2_{resp}_paired_p"] = float(w.pvalue)

json.dump(results, open(f"{RES}/reversal_stats.json", "w"), indent=2)

# ===================== FIGURE 3 =====================
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))
# (a) baseline by group, 16879
ax = axes[0]
groups = [("Control", ctrl), ("IBD baseline", base_ibd)]
ax.boxplot([g[1] for g in groups], tick_labels=[g[0] for g in groups],
           showfliers=False)
for i, (_, vals) in enumerate(groups, 1):
    ax.scatter(np.random.normal(i, 0.06, len(vals)), vals, s=12, alpha=0.5)
ax.set_ylabel("Transcriptomic age acceleration (z vs control)")
ax.set_title(f"GSE16879 baseline\nIBD > control (p={u1.pvalue:.1e})")

# (b) paired reversal, 16879
ax = axes[1]
for resp, col, xoff in [("R", "#2ca02c", 0), ("NR", "#d62728", 1)]:
    g = pr[pr["response"] == resp]
    for _, r in g.iterrows():
        ax.plot([xoff*3+0, xoff*3+1], [r["baseline"], r["post"]],
                color=col, alpha=0.35, lw=0.8)
    ax.scatter([xoff*3]*len(g), g["baseline"], color=col, s=18)
    ax.scatter([xoff*3+1]*len(g), g["post"], color=col, s=18)
ax.set_xticks([0, 1, 3, 4])
ax.set_xticklabels(["R base", "R post", "NR base", "NR post"])
ax.axhline(0, ls="--", color="grey", lw=0.8)
ax.set_ylabel("Age acceleration")
ax.set_title("GSE16879 paired pre/post\nresponders rejuvenate")

# (c) delta by response, both cohorts
ax = axes[2]
data, labels = [], []
for tag, p in [("16879", pr), ("73661", pr2)]:
    for resp in ["R", "NR"]:
        g = p[p["response"] == resp]["delta"].dropna()
        if len(g):
            data.append(g); labels.append(f"{tag}\n{resp}")
ax.boxplot(data, tick_labels=labels, showfliers=False)
for i, g in enumerate(data, 1):
    ax.scatter(np.random.normal(i, 0.06, len(g)), g, s=12, alpha=0.5)
ax.axhline(0, ls="--", color="grey", lw=0.8)
ax.set_ylabel("Δ age acceleration (post − baseline)")
ax.set_title("Reversal magnitude by response")

fig.suptitle("Transcriptomic mucosal age is a reversible, treatment-coupled biomarker in IBD",
             fontweight="bold")
fig.tight_layout()
fig.savefig(f"{RES}/Fig3_reversal.png", dpi=200)
print("\nsaved results/Fig3_reversal.png and reversal_stats.json")
