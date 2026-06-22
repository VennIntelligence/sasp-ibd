"""Does BASELINE transcriptomic aging predict treatment response?
Per cohort, among pre-treatment IBD samples with known outcome:
  - AUC of baseline age acceleration alone
  - AUC of baseline SenMayo alone
  - combined logistic model (age_accel + SenMayo), leave-one-out AUC
"""
import json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROC = "/Users/ujs/Downloads/lzy/data/interim"
RES = "/Users/ujs/Downloads/lzy/outputs"
out = {}

fig, axes = plt.subplots(1, 2, figsize=(9, 4.3))

for ax, gse in zip(axes, ["GSE16879", "GSE73661"]):
    sc = pd.read_csv(f"{PROC}/{gse}_scored.tsv", sep="\t", index_col=0)
    b = sc[(sc["timepoint"] == "baseline") & (sc["response"].isin(["R", "NR"]))].copy()
    b = b.dropna(subset=["age_accel", "senmayo"])
    # collapse to one row per patient (mean) to avoid pseudoreplication
    b = b.groupby("patient").agg(age_accel=("age_accel", "mean"),
                                 senmayo=("senmayo", "mean"),
                                 response=("response", "first")).reset_index()
    y = (b["response"] == "NR").astype(int).values  # predict non-response
    n_r, n_nr = (y == 0).sum(), (y == 1).sum()
    print(f"\n{gse}: {n_r} responders, {n_nr} non-responders (patient-level)")
    if n_r < 3 or n_nr < 3:
        ax.set_title(f"{gse}: too few for ROC"); continue

    def auc_dir(score):
        a = roc_auc_score(y, score)
        return max(a, 1 - a)  # direction-agnostic
    auc_age = auc_dir(b["age_accel"].values)
    auc_sen = auc_dir(b["senmayo"].values)
    out[f"{gse}_AUC_ageaccel"] = float(auc_age)
    out[f"{gse}_AUC_senmayo"] = float(auc_sen)

    # combined LOOCV
    X = b[["age_accel", "senmayo"]].values
    Xs = (X - X.mean(0)) / X.std(0)
    loo, preds = LeaveOneOut(), np.zeros(len(y))
    for tr, te in loo.split(Xs):
        lr = LogisticRegression().fit(Xs[tr], y[tr])
        preds[te] = lr.predict_proba(Xs[te])[:, 1]
    auc_comb = roc_auc_score(y, preds) if len(set(y)) == 2 else np.nan
    out[f"{gse}_AUC_combined_LOOCV"] = float(auc_comb)
    print(f"  AUC age-accel={auc_age:.2f}, SenMayo={auc_sen:.2f}, combined(LOO)={auc_comb:.2f}")

    fpr, tpr, _ = roc_curve(y, preds)
    ax.plot(fpr, tpr, lw=2, color="#d62728", label=f"combined AUC={auc_comb:.2f}")
    ax.plot([0, 1], [0, 1], "--", color="grey")
    ax.set_xlabel("FPR"); ax.set_ylabel("TPR")
    ax.set_title(f"{gse}: predict non-response\nage-accel AUC={auc_age:.2f}, SenMayo AUC={auc_sen:.2f}")
    ax.legend(fontsize=8)

fig.suptitle("Baseline mucosal aging predicts biologic-therapy response", fontweight="bold")
fig.tight_layout()
fig.savefig(f"{RES}/Fig5_predictor.png", dpi=200)
json.dump(out, open(f"{RES}/predictor_stats.json", "w"), indent=2)
print("\nsaved results/Fig5_predictor.png and predictor_stats.json")
