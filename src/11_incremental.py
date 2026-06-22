"""Is the SASP predictor just baseline severity?
GSE73661 has baseline Mayo endoscopic subscore (severity). Test whether SenMayo
adds predictive value beyond severity, report direction (which way SASP predicts),
and SASP-severity collinearity. Also report direction in GSE16879.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score

PROC = "/Users/ujs/Downloads/lzy/data/interim"
RES = "/Users/ujs/Downloads/lzy/outputs"
out = {}


def loo_auc(X, y):
    X = np.asarray(X, float)
    Xs = (X - X.mean(0)) / (X.std(0) + 1e-9)
    p = np.zeros(len(y))
    for tr, te in LeaveOneOut().split(Xs):
        p[te] = LogisticRegression().fit(Xs[tr], y[tr]).predict_proba(Xs[te])[:, 1]
    return roc_auc_score(y, p)


# ---- direction in both cohorts ----
for gse in ["GSE16879", "GSE73661"]:
    sc = pd.read_csv(f"{PROC}/{gse}_scored.tsv", sep="\t", index_col=0)
    bp = sc[(sc.timepoint == "baseline") & (sc.response.isin(["R", "NR"]))]
    bp = bp.dropna(subset=["senmayo"]).groupby("patient").agg(
        senmayo=("senmayo", "mean"), response=("response", "first"))
    mr = bp[bp.response == "R"].senmayo.median()
    mn = bp[bp.response == "NR"].senmayo.median()
    out[f"{gse}_baseline_SenMayo_R_median"] = float(mr)
    out[f"{gse}_baseline_SenMayo_NR_median"] = float(mn)
    out[f"{gse}_direction"] = "higher SASP -> responder" if mr > mn else "higher SASP -> non-responder"
    print(f"{gse}: baseline SenMayo  R={mr:.2f}  NR={mn:.2f}  -> {out[f'{gse}_direction']}")

# ---- incremental value over severity (GSE73661, has baseline Mayo) ----
sc = pd.read_csv(f"{PROC}/GSE73661_scored.tsv", sep="\t", index_col=0)
bp = sc[(sc.timepoint == "baseline") & (sc.response.isin(["R", "NR"]))]
bp = bp.dropna(subset=["senmayo", "mayo"]).groupby("patient").agg(
    senmayo=("senmayo", "mean"), mayo=("mayo", "mean"), response=("response", "first")).reset_index()
y = (bp.response == "NR").astype(int).values
print(f"\nGSE73661 incremental test: n={len(bp)} ({(y==0).sum()} R, {(y==1).sum()} NR)")

rho, pms = stats.spearmanr(bp.senmayo, bp.mayo)
out["GSE73661_senmayo_vs_mayo_rho"] = float(rho)
print(f"  SenMayo vs baseline Mayo collinearity: rho={rho:.2f} (p={pms:.2g})")

auc_mayo = loo_auc(bp[["mayo"]].values, y)
auc_sen = loo_auc(bp[["senmayo"]].values, y)
auc_both = loo_auc(bp[["mayo", "senmayo"]].values, y)
out["GSE73661_AUC_mayo_only"] = float(auc_mayo)
out["GSE73661_AUC_senmayo_only"] = float(auc_sen)
out["GSE73661_AUC_mayo_plus_senmayo"] = float(auc_both)
print(f"  LOO-AUC: Mayo only={auc_mayo:.2f} | SenMayo only={auc_sen:.2f} | Mayo+SenMayo={auc_both:.2f}")
print(f"  -> incremental ΔAUC over severity = {auc_both - auc_mayo:+.2f}")

json.dump(out, open(f"{RES}/incremental_stats.json", "w"), indent=2)
print("\nsaved results/incremental_stats.json")
