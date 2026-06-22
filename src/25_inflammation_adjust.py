"""Task 3 inflammation-adjusted prediction models."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut

from paths import P


OUT = P.out("25_inflammation_adjust")
COHORTS = ["GSE16879", "GSE73661", "GSE12251", "GSE23597", "GSE92415"]
THERAPY_KEEP = {
    "GSE16879": {"IFX"},
    "GSE73661": {"IFX", "VDZ"},
    "GSE12251": {"IFX"},
    "GSE23597": {"IFX"},
    "GSE92415": {"GLM"},
}
ADJUSTMENTS = {
    "crp_signature": ["inflam_crp"],
    "neutrophil_proxy": ["neutrophil_proxy"],
    "gene_inflammation": ["inflammation_score"],
    "mayo": ["mayo"],
    "gene_inflammation_plus_mayo": ["inflammation_score", "mayo"],
}


def patient_baseline(cohort: str) -> pd.DataFrame:
    sc = pd.read_csv(P.interim / f"{cohort}_scored.tsv", sep="\t")
    b = sc[(sc["timepoint"] == "baseline") & (sc["response"].isin(["R", "NR"]))].copy()
    keep = THERAPY_KEEP.get(cohort)
    if keep is not None:
        b = b[b["therapy"].isin(keep)]
    wanted = [
        "cohort",
        "patient",
        "response",
        "therapy",
        "senmayo",
        "inflam_crp",
        "neutrophil_proxy",
        "inflammation_score",
        "mayo",
    ]
    for col in wanted:
        if col not in b:
            b[col] = np.nan
    b = b[wanted]
    agg = {
        "response": ("response", "first"),
        "therapy": ("therapy", lambda x: ",".join(sorted(set(map(str, x))))),
        "senmayo": ("senmayo", "mean"),
        "inflam_crp": ("inflam_crp", "mean"),
        "neutrophil_proxy": ("neutrophil_proxy", "mean"),
        "inflammation_score": ("inflammation_score", "mean"),
        "mayo": ("mayo", "mean"),
    }
    out = b.groupby("patient", as_index=False).agg(**agg)
    out["cohort"] = cohort
    return out


def complete(df: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    cols = ["cohort", "patient", "response", "therapy", "senmayo"] + features
    return df[cols].dropna(subset=["senmayo"] + features).copy()


def z_train_test(train: pd.DataFrame, test: pd.DataFrame, features: list[str], cohort_fe: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    xtr = train[features].astype(float).copy()
    xte = test[features].astype(float).copy()
    mu = xtr.mean()
    sd = xtr.std().replace(0, 1).fillna(1)
    xtr = (xtr - mu) / sd
    xte = (xte - mu) / sd
    if cohort_fe:
        dtr = pd.get_dummies(train["cohort"], prefix="cohort", drop_first=True, dtype=float)
        dte = pd.get_dummies(test["cohort"], prefix="cohort", drop_first=True, dtype=float).reindex(columns=dtr.columns, fill_value=0)
        xtr = pd.concat([xtr.reset_index(drop=True), dtr.reset_index(drop=True)], axis=1)
        xte = pd.concat([xte.reset_index(drop=True), dte.reset_index(drop=True)], axis=1)
    return xtr, xte


def loo_auc(df: pd.DataFrame, features: list[str], cohort_fe: bool = False) -> float:
    if len(df) < 6:
        return np.nan
    y = (df["response"] == "NR").astype(int).to_numpy()
    if len(np.unique(y)) < 2:
        return np.nan
    preds = np.full(len(df), np.nan)
    for tr, te in LeaveOneOut().split(df):
        if len(np.unique(y[tr])) < 2:
            continue
        train, test = df.iloc[tr], df.iloc[te]
        xtr, xte = z_train_test(train, test, features, cohort_fe)
        try:
            lr = LogisticRegression(max_iter=1000, solver="liblinear")
            lr.fit(xtr, y[tr])
            preds[te] = lr.predict_proba(xte)[:, 1]
        except Exception:
            continue
    ok = ~np.isnan(preds)
    if ok.sum() < 3 or len(np.unique(y[ok])) < 2:
        return np.nan
    return float(roc_auc_score(y[ok], preds[ok]))


def logit_coef(df: pd.DataFrame, features: list[str], cohort_fe: bool = False) -> tuple[float, float, str]:
    y = (df["response"] == "NR").astype(int)
    if len(np.unique(y)) < 2:
        return np.nan, np.nan, "one_class"
    x = df[features].astype(float).copy()
    x = (x - x.mean()) / x.std().replace(0, 1).fillna(1)
    if cohort_fe:
        x = pd.concat([x.reset_index(drop=True), pd.get_dummies(df["cohort"], prefix="cohort", drop_first=True, dtype=float).reset_index(drop=True)], axis=1)
    x = sm.add_constant(x, has_constant="add")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = sm.Logit(y.to_numpy(), x.to_numpy()).fit(disp=0, maxiter=200)
        names = list(x.columns)
        idx = names.index("senmayo")
        return float(fit.params[idx]), float(fit.pvalues[idx]), "ok"
    except Exception as exc:
        return np.nan, np.nan, type(exc).__name__


def summarize(df: pd.DataFrame, label: str, adj_name: str, covars: list[str], cohort_fe: bool) -> dict[str, object] | None:
    dat = complete(df, covars)
    if len(dat) < 6 or dat["response"].nunique() < 2:
        return None
    y = (dat["response"] == "NR").astype(int)
    auc_sen = loo_auc(dat, ["senmayo"], cohort_fe=cohort_fe)
    auc_adj = loo_auc(dat, covars, cohort_fe=cohort_fe)
    auc_both = loo_auc(dat, covars + ["senmayo"], cohort_fe=cohort_fe)
    coef, pval, status = logit_coef(dat, covars + ["senmayo"], cohort_fe=cohort_fe)
    return {
        "analysis": label,
        "adjustment": adj_name,
        "covariates": ",".join(covars),
        "n_patients": len(dat),
        "n_R": int((y == 0).sum()),
        "n_NR": int((y == 1).sum()),
        "therapy": ",".join(sorted(set(",".join(dat["therapy"]).split(",")))),
        "auc_senmayo": auc_sen,
        "auc_adjustment": auc_adj,
        "auc_adjusted": auc_both,
        "delta_auc_vs_adjustment": auc_both - auc_adj if np.isfinite(auc_both) and np.isfinite(auc_adj) else np.nan,
        "senmayo_coef_std": coef,
        "senmayo_p": pval,
        "fit_status": status,
    }


def main() -> None:
    base = pd.concat([patient_baseline(c) for c in COHORTS], ignore_index=True)
    rows: list[dict[str, object]] = []
    for cohort, sub in base.groupby("cohort"):
        for adj_name, covars in ADJUSTMENTS.items():
            row = summarize(sub, cohort, adj_name, covars, cohort_fe=False)
            if row:
                rows.append(row)
    for adj_name, covars in ADJUSTMENTS.items():
        row = summarize(base, "pooled_with_cohort_FE", adj_name, covars, cohort_fe=True)
        if row:
            rows.append(row)

    out = pd.DataFrame(rows)
    root_file = P.outputs / "inflammation_adjusted.tsv"
    script_file = OUT / "inflammation_adjusted.tsv"
    out.to_csv(root_file, sep="\t", index=False)
    out.to_csv(script_file, sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"\nwrote {root_file}")


if __name__ == "__main__":
    main()
