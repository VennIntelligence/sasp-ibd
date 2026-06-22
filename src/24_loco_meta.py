"""Task 3 multi-cohort baseline AUC and random-effects meta-analysis."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.metrics import roc_auc_score

from paths import P


OUT = P.out("24_loco_meta")
COHORTS = ["GSE16879", "GSE73661", "GSE12251", "GSE23597", "GSE92415"]
SCORES = ["senmayo", "core_sen", "inflammation_score"]
THERAPY_KEEP = {
    "GSE16879": {"IFX"},
    "GSE73661": {"IFX", "VDZ"},
    "GSE12251": {"IFX"},
    "GSE23597": {"IFX"},
    "GSE92415": {"GLM"},
}


def patient_baseline(cohort: str) -> pd.DataFrame:
    sc = pd.read_csv(P.interim / f"{cohort}_scored.tsv", sep="\t")
    b = sc[(sc["timepoint"] == "baseline") & (sc["response"].isin(["R", "NR"]))].copy()
    keep = THERAPY_KEEP.get(cohort)
    if keep is not None:
        b = b[b["therapy"].isin(keep)]
    cols = ["cohort", "patient", "response", "therapy"] + [c for c in SCORES if c in b]
    b = b[cols].dropna(subset=[c for c in SCORES if c in b])
    agg = {c: (c, "mean") for c in SCORES if c in b}
    agg.update({"response": ("response", "first"), "therapy": ("therapy", lambda x: ",".join(sorted(set(map(str, x)))))})
    out = b.groupby("patient", as_index=False).agg(**agg)
    out["cohort"] = cohort
    return out


def auc_ci(y: np.ndarray, score: np.ndarray, n_boot: int = 5000) -> tuple[float, float, float, float]:
    auc = roc_auc_score(y, score)
    rng = np.random.default_rng(20260622)
    vals = []
    idx = np.arange(len(y))
    for _ in range(n_boot):
        take = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[take])) < 2:
            continue
        vals.append(roc_auc_score(y[take], score[take]))
    if vals:
        lo, hi = np.quantile(vals, [0.025, 0.975])
        se = float(np.std(vals, ddof=1))
    else:
        lo = hi = se = np.nan
    return float(auc), float(lo), float(hi), se


def logit(x: float) -> float:
    x = min(max(x, 1e-4), 1 - 1e-4)
    return math.log(x / (1 - x))


def inv_logit(x: float) -> float:
    return 1 / (1 + math.exp(-x))


def random_effects(rows: pd.DataFrame) -> dict[str, float]:
    usable = rows.dropna(subset=["auc", "auc_se"]).copy()
    usable = usable[(usable["auc_se"] > 0) & (usable["auc"] > 0) & (usable["auc"] < 1)]
    if len(usable) < 2:
        return {"auc": np.nan, "ci_low": np.nan, "ci_high": np.nan, "i2": np.nan, "tau2": np.nan}
    yi = usable["auc"].map(logit).to_numpy()
    sei = (usable["auc_se"] / (usable["auc"] * (1 - usable["auc"]))).to_numpy()
    vi = sei**2
    wi = 1 / vi
    fixed = np.sum(wi * yi) / np.sum(wi)
    q = np.sum(wi * (yi - fixed) ** 2)
    df = len(yi) - 1
    c = np.sum(wi) - (np.sum(wi**2) / np.sum(wi))
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    w_re = 1 / (vi + tau2)
    mu = np.sum(w_re * yi) / np.sum(w_re)
    se_mu = math.sqrt(1 / np.sum(w_re))
    i2 = max(0.0, (q - df) / q) * 100 if q > 0 else 0.0
    return {
        "auc": inv_logit(mu),
        "ci_low": inv_logit(mu - 1.96 * se_mu),
        "ci_high": inv_logit(mu + 1.96 * se_mu),
        "i2": i2,
        "tau2": tau2,
    }


def main() -> None:
    base = pd.concat([patient_baseline(c) for c in COHORTS], ignore_index=True)
    rows = []
    for score in SCORES:
        for cohort, sub in base.dropna(subset=[score]).groupby("cohort"):
            y = (sub["response"] == "NR").astype(int).to_numpy()
            if len(np.unique(y)) < 2:
                continue
            auc, lo, hi, se = auc_ci(y, sub[score].to_numpy())
            rows.append(
                {
                    "analysis": "cohort",
                    "score": score,
                    "cohort": cohort,
                    "n_patients": len(sub),
                    "n_R": int((y == 0).sum()),
                    "n_NR": int((y == 1).sum()),
                    "therapy": ",".join(sorted(set(",".join(sub["therapy"]).split(",")))),
                    "auc": auc,
                    "ci_low": lo,
                    "ci_high": hi,
                    "auc_se": se,
                    "direction": "higher_score_more_nonresponse" if auc >= 0.5 else "reversed",
                    "i2": np.nan,
                    "tau2": np.nan,
                }
            )
        score_rows = pd.DataFrame([r for r in rows if r["analysis"] == "cohort" and r["score"] == score])
        pooled = random_effects(score_rows)
        rows.append(
            {
                "analysis": "pooled_random_effects",
                "score": score,
                "cohort": "POOLED",
                "n_patients": int(score_rows["n_patients"].sum()) if len(score_rows) else 0,
                "n_R": int(score_rows["n_R"].sum()) if len(score_rows) else 0,
                "n_NR": int(score_rows["n_NR"].sum()) if len(score_rows) else 0,
                "therapy": "mixed",
                "auc": pooled["auc"],
                "ci_low": pooled["ci_low"],
                "ci_high": pooled["ci_high"],
                "auc_se": np.nan,
                "direction": "higher_score_more_nonresponse" if pooled["auc"] >= 0.5 else "reversed",
                "i2": pooled["i2"],
                "tau2": pooled["tau2"],
            }
        )
        for drop in score_rows["cohort"]:
            loco = random_effects(score_rows[score_rows["cohort"] != drop])
            rows.append(
                {
                    "analysis": "leave_one_cohort_out",
                    "score": score,
                    "cohort": f"without_{drop}",
                    "n_patients": int(score_rows.loc[score_rows["cohort"] != drop, "n_patients"].sum()),
                    "n_R": int(score_rows.loc[score_rows["cohort"] != drop, "n_R"].sum()),
                    "n_NR": int(score_rows.loc[score_rows["cohort"] != drop, "n_NR"].sum()),
                    "therapy": "mixed",
                    "auc": loco["auc"],
                    "ci_low": loco["ci_low"],
                    "ci_high": loco["ci_high"],
                    "auc_se": np.nan,
                    "direction": "higher_score_more_nonresponse" if loco["auc"] >= 0.5 else "reversed",
                    "i2": loco["i2"],
                    "tau2": loco["tau2"],
                }
            )

    out = pd.DataFrame(rows)
    root_file = P.outputs / "multicohort_auc.tsv"
    script_file = OUT / "multicohort_auc.tsv"
    out.to_csv(root_file, sep="\t", index=False)
    out.to_csv(script_file, sep="\t", index=False)

    sen = out[(out["score"] == "senmayo") & (out["analysis"] == "cohort")].copy()
    pooled = out[(out["score"] == "senmayo") & (out["analysis"] == "pooled_random_effects")].iloc[0]
    sen = pd.concat([sen, pooled.to_frame().T], ignore_index=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    y = np.arange(len(sen))
    colors = ["#2f6f8f" if a == "cohort" else "#b03a2e" for a in sen["analysis"]]
    ax.errorbar(
        sen["auc"].astype(float),
        y,
        xerr=[sen["auc"].astype(float) - sen["ci_low"].astype(float), sen["ci_high"].astype(float) - sen["auc"].astype(float)],
        fmt="none",
        ecolor="#444444",
        elinewidth=1,
        capsize=3,
    )
    ax.scatter(sen["auc"].astype(float), y, s=46, color=colors, zorder=3)
    labels = [f"{r.cohort} (n={int(r.n_patients)})" for r in sen.itertuples()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.axvline(0.5, color="#888888", linestyle="--", linewidth=1)
    ax.set_xlim(0.25, 1.0)
    ax.set_xlabel("AUC for baseline SenMayo predicting non-response")
    ax.set_title("Task 3 multi-cohort bulk validation")
    fig.tight_layout()
    fig_file = P.outputs / "Fig_task3_multicohort.png"
    fig.savefig(fig_file, dpi=200)
    fig.savefig(OUT / "Fig_task3_multicohort.png", dpi=200)
    print(out[out["score"] == "senmayo"].to_string(index=False))
    print(f"\nwrote {root_file}")
    print(f"wrote {fig_file}")


if __name__ == "__main__":
    main()
