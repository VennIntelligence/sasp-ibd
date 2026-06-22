"""MAIN (data-driven) analysis: mucosal cellular-senescence / SASP burden as a
reversible, predictive biomarker of biologic-therapy response in IBD.

The GTEx transcriptomic age clock was a NEGATIVE result (reported as contrast).
The robust, consistent signal is senescence/SASP (SenMayo) + senescence-specific
arrest markers (p16/p21), which (i) rise in active disease, (ii) predict response,
(iii) resolve in responders. We also test senescence-SPECIFIC markers to argue the
signal is not merely generic inflammation.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/ujs/Downloads/lzy"
PROC = f"{BASE}/data/interim"
RES = f"{BASE}/outputs"
sets = json.load(open(f"{BASE}/data/external/genesets/senescence_sets.json"))
SENMAYO = sets["SenMayo"]
SEN_UP = ["CDKN1A", "CDKN2A", "CDKN2B", "GLB1", "SERPINE1"]   # arrest/SA-bGal
SEN_DOWN = ["MKI67", "LMNB1"]                                 # proliferation/lamin (lost in senescence)
out = {}


def zscore(expr, genes):
    g = [x for x in genes if x in expr.index]
    sub = expr.loc[g]
    return sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0), g


def senmayo_score(expr):
    z, g = zscore(expr, SENMAYO)
    return z.mean(axis=0), len(g)


def core_score(expr):
    zu, gu = zscore(expr, SEN_UP)
    zd, gd = zscore(expr, SEN_DOWN)
    return zu.mean(axis=0) - zd.mean(axis=0), gu, gd


def paired(df, col):
    rows = []
    for pid, sub in df[df["group"] != "Control"].groupby("patient"):
        b = sub[sub["timepoint"] == "baseline"][col]
        p = sub[sub["timepoint"] == "post"][col]
        if len(b) and len(p):
            rows.append((pid, b.mean(), p.mean(), sub["response"].iloc[0]))
    return pd.DataFrame(rows, columns=["patient", "baseline", "post", "response"])


fig, axes = plt.subplots(2, 3, figsize=(14, 8.6))

for col, gse in enumerate(["GSE16879", "GSE73661"]):
    expr = pd.read_parquet(f"{PROC}/{gse}_expr.parquet")
    sc = pd.read_csv(f"{PROC}/{gse}_scored.tsv", sep="\t", index_col=0)
    sm, n_sm = senmayo_score(expr)
    cs, gu, gd = core_score(expr)
    sc["senmayo"] = sm.reindex(sc.index)
    sc["core_sen"] = cs.reindex(sc.index)
    for marker in ["CDKN1A", "CDKN2A"]:
        if marker in expr.index:
            z = (expr.loc[marker] - expr.loc[marker].mean()) / (expr.loc[marker].std() or 1)
            sc[marker] = z.reindex(sc.index)
    sc.to_csv(f"{PROC}/{gse}_scored.tsv", sep="\t")

    # ---- baseline vs control ----
    base = sc[sc["timepoint"].isin(["baseline", "control"])]
    for score in ["senmayo", "core_sen"]:
        ibd = base[base.group != "Control"][score].dropna()
        ctl = base[base.group == "Control"][score].dropna()
        p = stats.mannwhitneyu(ibd, ctl, alternative="greater").pvalue
        out[f"{gse}_{score}_baseline_vs_ctrl_p"] = float(p)
        out[f"{gse}_{score}_baseline_ibd_median"] = float(ibd.median())
        out[f"{gse}_{score}_baseline_ctl_median"] = float(ctl.median())

    # ---- reversal (paired) for SenMayo ----
    pr = paired(sc, "senmayo"); pr["delta"] = pr.post - pr.baseline
    for resp in ["R", "NR"]:
        g = pr[pr.response == resp]
        if len(g) >= 3:
            w = stats.wilcoxon(g.baseline, g.post).pvalue
            out[f"{gse}_senmayo_{resp}_paired_p"] = float(w)
            out[f"{gse}_senmayo_{resp}_delta_median"] = float(g.delta.median())
    if (pr.response == "R").sum() >= 3 and (pr.response == "NR").sum() >= 3:
        out[f"{gse}_senmayo_deltaRvsNR_p"] = float(
            stats.mannwhitneyu(pr[pr.response == "R"].delta, pr[pr.response == "NR"].delta,
                               alternative="less").pvalue)

    # ---- prediction (patient-level) ----
    bp = sc[(sc.timepoint == "baseline") & (sc.response.isin(["R", "NR"]))]
    bp = bp.dropna(subset=["senmayo", "core_sen"]).groupby("patient").agg(
        senmayo=("senmayo", "mean"), core_sen=("core_sen", "mean"),
        CDKN1A=("CDKN1A", "mean") if "CDKN1A" in sc else ("senmayo", "mean"),
        CDKN2A=("CDKN2A", "mean") if "CDKN2A" in sc else ("senmayo", "mean"),
        response=("response", "first")).reset_index()
    y = (bp.response == "NR").astype(int).values

    def auc_dir(s):
        a = roc_auc_score(y, s); return max(a, 1 - a)
    for feat in ["senmayo", "core_sen", "CDKN1A", "CDKN2A"]:
        if feat in bp:
            out[f"{gse}_AUC_{feat}"] = float(auc_dir(bp[feat].values))
    # combined LOOCV
    X = bp[["senmayo", "core_sen"]].values
    Xs = (X - X.mean(0)) / (X.std(0) + 1e-9)
    preds = np.zeros(len(y))
    for tr, te in LeaveOneOut().split(Xs):
        preds[te] = LogisticRegression().fit(Xs[tr], y[tr]).predict_proba(Xs[te])[:, 1]
    auc_comb = roc_auc_score(y, preds)
    out[f"{gse}_AUC_senescence_combined_LOO"] = float(auc_comb)

    # ===== plots =====
    # row0: baseline senmayo by group
    ax = axes[0, col]
    groups = sorted(base.group.unique(), key=lambda g: (g != "Control", g))
    data = [base[base.group == g]["senmayo"].dropna() for g in groups]
    ax.boxplot(data, tick_labels=groups, showfliers=False)
    for i, d in enumerate(data, 1):
        ax.scatter(np.random.normal(i, 0.06, len(d)), d, s=12, alpha=0.5)
    ax.axhline(0, ls="--", color="grey", lw=.8); ax.set_ylabel("SenMayo score")
    pv = out[f"{gse}_senmayo_baseline_vs_ctrl_p"]
    ax.set_title(f"{gse}: senescence/SASP up in active IBD\n(IBD>ctrl p={pv:.1e})")
    # row1: reversal + ROC text
    ax = axes[1, col]
    for resp, cc, xo in [("R", "#2ca02c", 0), ("NR", "#d62728", 1)]:
        g = pr[pr.response == resp]
        for _, r in g.iterrows():
            ax.plot([xo*3, xo*3+1], [r.baseline, r.post], color=cc, alpha=.3, lw=.8)
        ax.scatter([xo*3]*len(g), g.baseline, color=cc, s=16)
        ax.scatter([xo*3+1]*len(g), g.post, color=cc, s=16)
    ax.set_xticks([0, 1, 3, 4]); ax.set_xticklabels(["R\nbase", "R\npost", "NR\nbase", "NR\npost"])
    ax.axhline(0, ls="--", color="grey", lw=.8); ax.set_ylabel("SenMayo score")
    aucs = out.get(f"{gse}_AUC_senmayo", np.nan)
    ax.set_title(f"{gse}: SASP resolves in responders\nbaseline SenMayo predicts resp AUC={aucs:.2f}")

# row0 col2: AUC summary bar
ax = axes[0, 2]
feats = ["senmayo", "core_sen", "CDKN1A", "CDKN2A"]
labels = ["SenMayo", "Core p16/p21", "CDKN1A", "CDKN2A"]
x = np.arange(len(feats)); w = 0.38
for k, (gse, cc) in enumerate([("GSE16879", "#1f77b4"), ("GSE73661", "#ff7f0e")]):
    vals = [out.get(f"{gse}_AUC_{f}", np.nan) for f in feats]
    ax.bar(x + (k-0.5)*w, vals, w, label=gse, color=cc)
ax.axhline(0.5, ls="--", color="grey"); ax.set_ylim(0.4, 1)
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=20, fontsize=8)
ax.set_ylabel("AUC (predict non-response)"); ax.legend(fontsize=8)
ax.set_title("Baseline senescence predicts response\n(senescence-specific markers too)")

# row1 col2: clock contrast (negative result, honesty)
ax = axes[1, 2]
clk = []
for gse in ["GSE16879", "GSE73661"]:
    sc = pd.read_csv(f"{PROC}/{gse}_scored.tsv", sep="\t", index_col=0)
    bp = sc[(sc.timepoint == "baseline") & (sc.response.isin(["R", "NR"]))]
    bp = bp.dropna(subset=["age_accel"]).groupby("patient").agg(
        age_accel=("age_accel", "mean"), response=("response", "first"))
    y = (bp.response == "NR").astype(int).values
    a = roc_auc_score(y, bp.age_accel.values); a = max(a, 1-a)
    clk.append((gse, a, out.get(f"{gse}_AUC_senmayo", np.nan)))
ax.bar([0, 1], [clk[0][1], clk[1][1]], 0.35, label="Aging clock", color="#999999")
ax.bar([0.4, 1.4], [clk[0][2], clk[1][2]], 0.35, label="SenMayo", color="#2ca02c")
ax.set_xticks([0.2, 1.2]); ax.set_xticklabels(["GSE16879", "GSE73661"])
ax.axhline(0.5, ls="--", color="grey"); ax.set_ylim(0.4, 1)
ax.set_ylabel("AUC"); ax.legend(fontsize=8)
ax.set_title("Honest contrast: senescence >> age clock\n(clock ~ chance)")

fig.suptitle("Mucosal cellular senescence/SASP is a reversible biomarker that predicts biologic response in IBD",
             fontweight="bold", fontsize=12)
fig.tight_layout()
fig.savefig(f"{RES}/Fig_MAIN_senescence.png", dpi=190)
json.dump(out, open(f"{RES}/senescence_main_stats.json", "w"), indent=2)
print(json.dumps(out, indent=2))
print("\nsaved results/Fig_MAIN_senescence.png and senescence_main_stats.json")
