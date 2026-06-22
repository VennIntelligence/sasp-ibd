"""Single-cell signature NNLS deconvolution of bulk response cohorts."""
from __future__ import annotations

import math
import os
import time
from pathlib import Path

_THREADS = "30"
for _key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ.setdefault(_key, _THREADS)

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
import statsmodels.api as sm
from joblib import Parallel, delayed
from pydantic import BaseModel, ConfigDict, Field
from scipy import stats
from scipy.optimize import nnls
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paths import P


OUT = P.out("deconv")
SC01 = P.out("sc_01")
COHORTS = ("GSE12251", "GSE16879", "GSE23597", "GSE73661", "GSE92415")
NEUTROPHIL_GENES = ("S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR1", "CXCR2", "MMP8", "MMP9", "ELANE", "MPO", "CEACAM8", "LCN2")
REFRACTORY_MODULE = (
    "OSM", "OSMR", "TREM1", "IL13RA2", "CXCR2", "CCL8", "IL11", "IL24",
    "CXCL1", "CXCL2", "CXCL3", "CXCL8", "MMP1", "MMP3", "MMP8", "MMP9",
    "MMP10", "MMP12", "S100A8", "S100A9", "S100A12", "FCGR2A", "FCGR3B",
    "CSF3R", "C5AR1", "FPR1", "TREM1", "BCL2A1", "AQP9", "PROK2",
    "COL1A1", "COL1A2", "COL3A1", "FAP", "POSTN", "CXCL10", "IL1B", "TNF",
)
FOCUS_GENES = tuple(dict.fromkeys(NEUTROPHIL_GENES + REFRACTORY_MODULE + ("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "LYZ", "LST1", "CD68")))


class Config(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    n_jobs: int = Field(default=30, ge=1)
    chunk_cells: int = Field(default=10_000, ge=1000)
    min_cells: int = Field(default=80, ge=10)
    markers_per_celltype: int = Field(default=180, ge=20)
    max_markers_total: int = Field(default=1600, ge=100)
    signature_paths: tuple[Path, ...] = (
        SC01 / "smillie_uc_qc.h5ad",
        SC01 / "martin_cd_qc.h5ad",
    )


CFG = Config()


class FractionRow(BaseModel):
    model_config = ConfigDict(extra="forbid")

    cohort: str
    gsm: str
    patient: str
    response: str
    therapy: str
    cell_type: str
    fraction: float
    nnls_residual_rel: float
    n_marker_genes: int


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def broad_cell_type(label: str, dataset: str) -> str:
    s = str(label)
    low = s.lower()
    if dataset == "Martin_CD":
        return {
            "T_NK": "T_NK",
            "B_Plasma": "B_Plasma",
            "Epithelial": "Epithelial",
            "Myeloid": "Myeloid",
            "Fibroblast": "Fibroblast",
            "Endothelial": "Endothelial",
            "Mast": "Mast",
        }.get(s, "Unknown")
    if any(x in low for x in ["macrophage", "monocyte", "dc1", "dc2"]):
        return "Myeloid"
    if any(x in low for x in ["fibro", "wnt2b", "wnt5b", "myofibro", "rspo3"]):
        return "Fibroblast"
    if any(x in low for x in ["enterocyte", "goblet", "tuft", "stem", "ta", "m cell", "enteroendocrine"]):
        return "Epithelial"
    if any(x in low for x in ["cd4", "cd8", "treg", "nks", "ilcs", "cycling t"]):
        return "T_NK"
    if any(x in low for x in ["plasma", "follicular", "cycling b", "gc"]):
        return "B_Plasma"
    if any(x in low for x in ["endothelial", "venule", "microvascular"]):
        return "Endothelial"
    if "mast" in low:
        return "Mast"
    if any(x in low for x in ["pericyte", "glia"]):
        return "Stromal_other"
    return "Unknown"


def symbol_index(adata: ad.AnnData) -> pd.Index:
    if "gene_symbol" in adata.var:
        vals = adata.var["gene_symbol"].astype(str)
    else:
        vals = pd.Index(adata.var_names.astype(str))
    vals = vals.str.upper().str.strip()
    return pd.Index(vals)


def group_mean_backed(adata: ad.AnnData, idx: np.ndarray, chunk_cells: int) -> np.ndarray:
    total = np.zeros(adata.n_vars, dtype=np.float64)
    n = 0
    idx = np.asarray(idx, dtype=np.int64)
    for start in range(0, len(idx), chunk_cells):
        x = adata.X[idx[start : start + chunk_cells], :]
        vals = np.asarray(x.sum(axis=0)).ravel() if sp.issparse(x) else np.asarray(x).sum(axis=0)
        total += vals
        n += x.shape[0]
    return (total / max(n, 1)).astype(np.float32)


def build_signature(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    log(f"building signature from {path.name}")
    adata = ad.read_h5ad(path, backed="r")
    dataset = str(adata.obs["dataset"].iloc[0]) if "dataset" in adata.obs else path.stem
    obs = adata.obs[["cell_type"]].copy()
    obs["dataset"] = dataset
    obs["broad_cell_type"] = [broad_cell_type(x, dataset) for x in obs["cell_type"]]
    counts = (
        obs.groupby(["dataset", "broad_cell_type", "cell_type"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["dataset", "broad_cell_type", "n_cells"], ascending=[True, True, False])
    )
    groups = obs.groupby("broad_cell_type", observed=True).indices
    symbols = symbol_index(adata)
    means = {}
    for ct, idx in groups.items():
        if ct == "Unknown" or len(idx) < CFG.min_cells:
            continue
        means[ct] = group_mean_backed(adata, np.asarray(idx), CFG.chunk_cells)
    adata.file.close()
    sig = pd.DataFrame(means, index=symbols)
    sig = sig[(sig.index != "") & ~sig.index.str.startswith("NAN")]
    sig = sig.groupby(level=0, sort=False).mean()
    sig.columns.name = "cell_type"
    sig.index.name = "gene"
    return sig, counts


def build_all_signatures() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sigs, counts = {}, []
    for path in CFG.signature_paths:
        sig, ct = build_signature(path)
        sigs[path.stem] = sig
        counts.append(ct)
        sig.to_csv(OUT / f"{path.stem}_broad_signature.tsv", sep="\t")
    long = []
    for source, sig in sigs.items():
        x = sig.reset_index().melt(id_vars="gene", var_name="cell_type", value_name="mean_logexpr")
        x["signature_source"] = source
        long.append(x)
    pd.concat(long, ignore_index=True).to_csv(OUT / "broad_signature_long.tsv", sep="\t", index=False)
    combo = pd.concat(sigs, names=["signature_source", "gene"])
    combo = combo.groupby("gene", sort=False).mean()
    preferred = ["Epithelial", "Fibroblast", "Myeloid", "T_NK", "B_Plasma", "Endothelial", "Mast", "Stromal_other"]
    combo = combo[[c for c in preferred if c in combo.columns]].fillna(0.0)
    combo.to_csv(OUT / "ibd_combined_broad_signature.tsv", sep="\t")
    cell_counts = pd.concat(counts, ignore_index=True)
    cell_counts.to_csv(OUT / "signature_celltype_counts.tsv", sep="\t", index=False)
    return combo, cell_counts, pd.concat(long, ignore_index=True)


def marker_genes(sig: pd.DataFrame) -> pd.DataFrame:
    rows = []
    vals = sig.fillna(0.0)
    for ct in vals.columns:
        other = vals.drop(columns=ct)
        spec = vals[ct] - other.max(axis=1)
        x = pd.DataFrame({"gene": vals.index, "cell_type": ct, "mean_logexpr": vals[ct].values, "specificity": spec.values})
        x = x[(x["specificity"] > 0) & (x["mean_logexpr"] > 0.02)]
        rows.append(x.sort_values(["specificity", "mean_logexpr"], ascending=False).head(CFG.markers_per_celltype))
    markers = pd.concat(rows, ignore_index=True)
    if FOCUS_GENES:
        focus = [g for g in FOCUS_GENES if g in vals.index]
        for gene in focus:
            ct = vals.loc[gene].idxmax()
            markers = pd.concat(
                [
                    markers,
                    pd.DataFrame(
                        [{"gene": gene, "cell_type": ct, "mean_logexpr": float(vals.loc[gene, ct]), "specificity": float(vals.loc[gene, ct] - vals.loc[gene].drop(ct).max())}]
                    ),
                ],
                ignore_index=True,
            )
    markers = markers.sort_values(["cell_type", "specificity", "mean_logexpr"], ascending=False).drop_duplicates(["gene", "cell_type"])
    keep = markers.drop_duplicates("gene").sort_values(["specificity", "mean_logexpr"], ascending=False)
    if len(keep) > CFG.max_markers_total:
        keep = keep.head(CFG.max_markers_total)
    markers = markers[markers["gene"].isin(set(keep["gene"]))].copy()
    markers.to_csv(OUT / "signature_marker_genes.tsv", sep="\t", index=False)
    return markers


def collapse_expr(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr.copy()
    expr.index = expr.index.astype(str).str.upper().str.strip()
    expr = expr[(expr.index != "") & ~expr.index.str.startswith("NAN")]
    if np.nanmax(expr.to_numpy(dtype=float)) > 100:
        expr = np.log2(expr.clip(lower=0) + 1)
    return expr.groupby(level=0, sort=False).mean()


def rank01(v: np.ndarray) -> np.ndarray:
    r = stats.rankdata(v, method="average")
    return ((r - 1) / max(len(r) - 1, 1)).astype(np.float64)


def scaled_score(expr: pd.DataFrame, genes: tuple[str, ...], samples: list[str]) -> pd.Series:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=samples)
    sub = expr.loc[present, samples].astype(float)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def deconvolve_cohort(cohort: str, sig: pd.DataFrame, markers: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = collapse_expr(pd.read_parquet(P.interim / f"{cohort}_expr.parquet"))
    scored = pd.read_csv(P.interim / f"{cohort}_scored.tsv", sep="\t")
    base = scored[
        scored["timepoint"].eq("baseline")
        & scored["response"].isin(["R", "NR"])
        & ~scored["therapy"].astype(str).str.lower().isin(["placebo", "none"])
        & scored["gsm"].isin(expr.columns)
    ].copy()
    genes = [g for g in markers["gene"].drop_duplicates() if g in expr.index and g in sig.index]
    if len(base) == 0 or len(genes) < 20:
        return pd.DataFrame(), pd.DataFrame()
    raw_a = sig.loc[genes].fillna(0.0)
    span = raw_a.max(axis=1) - raw_a.min(axis=1)
    keep = span > 1e-6
    raw_a = raw_a.loc[keep]
    genes = raw_a.index.tolist()
    a = raw_a.sub(raw_a.min(axis=1), axis=0).div(span.loc[genes], axis=0).to_numpy(np.float64)
    norms = np.linalg.norm(a, axis=0)
    norms[norms == 0] = 1.0
    a_fit = a / norms
    rows = []
    for sample in base["gsm"]:
        b = rank01(expr.loc[genes, sample].to_numpy(float))
        coef, resid = nnls(a_fit, b, maxiter=20_000)
        coef = coef / norms
        frac = coef / coef.sum() if coef.sum() > 0 else np.repeat(1 / len(coef), len(coef))
        rel_resid = float(resid / (np.linalg.norm(b) + 1e-12))
        meta = base.loc[base["gsm"].eq(sample)].iloc[0]
        for ct, value in zip(raw_a.columns, frac, strict=True):
            rows.append(
                FractionRow(
                    cohort=cohort,
                    gsm=str(sample),
                    patient=str(meta.get("patient", "")),
                    response=str(meta["response"]),
                    therapy=str(meta.get("therapy", "")),
                    cell_type=str(ct),
                    fraction=float(value),
                    nnls_residual_rel=rel_resid,
                    n_marker_genes=len(genes),
                ).model_dump()
            )
    scores = base[["gsm", "cohort", "patient", "response", "therapy", "senmayo", "inflam_crp", "neutrophil_proxy", "inflammation_score"]].copy()
    scores["cxcr2_expr_z"] = scaled_score(expr, ("CXCR2",), scores["gsm"].tolist()).reindex(scores["gsm"]).to_numpy()
    scores["refractory_module_score"] = scaled_score(expr, REFRACTORY_MODULE, scores["gsm"].tolist()).reindex(scores["gsm"]).to_numpy()
    scores["neutrophil_marker_score"] = scaled_score(expr, NEUTROPHIL_GENES, scores["gsm"].tolist()).reindex(scores["gsm"]).to_numpy()
    scores["n_refractory_module_genes"] = sum(g in expr.index for g in REFRACTORY_MODULE)
    scores["n_neutrophil_marker_genes"] = sum(g in expr.index for g in NEUTROPHIL_GENES)
    return pd.DataFrame(rows), scores


def logistic_beta(y: np.ndarray, x: np.ndarray) -> tuple[float, float, float, str]:
    if len(np.unique(y)) < 2 or np.nanstd(x) == 0:
        return np.nan, np.nan, np.nan, "not_estimable"
    z = (x - np.nanmean(x)) / (np.nanstd(x, ddof=1) or 1.0)
    try:
        fit = sm.Logit(y, sm.add_constant(z, has_constant="add")).fit(disp=False, maxiter=300)
        return float(fit.params[1]), float(fit.bse[1]), float(fit.pvalues[1]), "ok"
    except Exception as exc:
        return np.nan, np.nan, np.nan, f"logit_failed: {type(exc).__name__}"


def feature_tests(feature_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (feature_class, feature), d in feature_df.groupby(["feature_class", "feature"], observed=True):
        for cohort, sub in d.groupby("cohort", observed=True):
            sub = sub.dropna(subset=["value", "nr"])
            if sub["nr"].nunique() < 2 or len(sub) < 6:
                continue
            y, x = sub["nr"].to_numpy(int), sub["value"].to_numpy(float)
            beta, se, p_logit, status = logistic_beta(y, x)
            try:
                auc = float(roc_auc_score(y, x))
            except ValueError:
                auc = np.nan
            r = sub.loc[sub["response"].eq("R"), "value"]
            nr = sub.loc[sub["response"].eq("NR"), "value"]
            mw = stats.mannwhitneyu(nr, r, alternative="two-sided") if len(r) and len(nr) else None
            rows.append(
                {
                    "feature_class": feature_class,
                    "feature": feature,
                    "cohort": cohort,
                    "n": len(sub),
                    "n_NR": int(y.sum()),
                    "n_R": int((1 - y).sum()),
                    "auc_NR_high": auc,
                    "median_NR": float(nr.median()) if len(nr) else np.nan,
                    "median_R": float(r.median()) if len(r) else np.nan,
                    "delta_median_NR_minus_R": float(nr.median() - r.median()) if len(r) and len(nr) else np.nan,
                    "logit_beta_per_sd": beta,
                    "logit_se": se,
                    "logit_OR_per_sd": math.exp(beta) if pd.notna(beta) else np.nan,
                    "logit_p": p_logit,
                    "mannwhitney_p": float(mw.pvalue) if mw else np.nan,
                    "status": status,
                }
            )
    out = pd.DataFrame(rows)
    meta_rows = []
    ok = out[out["status"].eq("ok") & out["logit_se"].notna() & (out["logit_se"] > 0)]
    for (feature_class, feature), d in ok.groupby(["feature_class", "feature"], observed=True):
        beta = d["logit_beta_per_sd"].to_numpy(float)
        var = np.square(d["logit_se"].to_numpy(float))
        w = 1 / var
        fixed = float(np.sum(w * beta) / np.sum(w))
        q = float(np.sum(w * np.square(beta - fixed)))
        df = len(beta) - 1
        c = float(np.sum(w) - np.sum(np.square(w)) / np.sum(w))
        tau2 = max(0.0, (q - df) / c) if c > 0 and df > 0 else 0.0
        wr = 1 / (var + tau2)
        re = float(np.sum(wr * beta) / np.sum(wr))
        se = float(np.sqrt(1 / np.sum(wr)))
        p = float(2 * stats.norm.sf(abs(re / se))) if se > 0 else np.nan
        meta_rows.append(
            {
                "feature_class": feature_class,
                "feature": feature,
                "cohort": "random_effects",
                "n": int(d["n"].sum()),
                "n_NR": int(d["n_NR"].sum()),
                "n_R": int(d["n_R"].sum()),
                "auc_NR_high": np.nan,
                "median_NR": np.nan,
                "median_R": np.nan,
                "delta_median_NR_minus_R": np.nan,
                "logit_beta_per_sd": re,
                "logit_se": se,
                "logit_OR_per_sd": math.exp(re),
                "logit_p": p,
                "mannwhitney_p": np.nan,
                "meta_tau2": tau2,
                "meta_Q": q,
                "meta_I2": max(0.0, (q - df) / q) if q > 0 and df > 0 else 0.0,
                "status": "ok_random_effects",
            }
        )
    if meta_rows:
        out = pd.concat([out, pd.DataFrame(meta_rows)], ignore_index=True)
    out["logit_fdr"] = np.nan
    ix = out["logit_p"].notna()
    if ix.any():
        out.loc[ix, "logit_fdr"] = multipletests(out.loc[ix, "logit_p"], method="fdr_bh")[1]
    return out.sort_values(["feature_class", "feature", "cohort"])


def pooled_increment(feature_df: pd.DataFrame) -> pd.DataFrame:
    base_scores = feature_df[feature_df["feature"].eq("senmayo")][["cohort", "gsm", "value"]].rename(columns={"value": "senmayo_value"})
    rows = []
    for feature in sorted(feature_df.loc[feature_df["feature_class"].eq("cell_fraction"), "feature"].unique()):
        d = feature_df[(feature_df["feature_class"].eq("cell_fraction")) & (feature_df["feature"].eq(feature))]
        m = d.merge(base_scores, on=["cohort", "gsm"], how="inner").dropna(subset=["value", "senmayo_value", "nr"])
        if len(m) < 12 or m["nr"].nunique() < 2:
            continue
        y = m["nr"].to_numpy(int)
        cohort_dummies = pd.get_dummies(m["cohort"], drop_first=True, dtype=float)
        z_s = (m["senmayo_value"] - m["senmayo_value"].mean()) / (m["senmayo_value"].std() or 1.0)
        z_f = (m["value"] - m["value"].mean()) / (m["value"].std() or 1.0)
        x0 = sm.add_constant(pd.concat([z_s.rename("senmayo"), cohort_dummies], axis=1), has_constant="add")
        x1 = sm.add_constant(pd.concat([z_s.rename("senmayo"), z_f.rename(feature), cohort_dummies], axis=1), has_constant="add")
        try:
            f0 = sm.Logit(y, x0).fit(disp=False, maxiter=300)
            f1 = sm.Logit(y, x1).fit(disp=False, maxiter=300)
            lrt = 2 * (f1.llf - f0.llf)
            p_lrt = float(stats.chi2.sf(lrt, 1))
            auc0 = float(roc_auc_score(y, f0.predict(x0)))
            auc1 = float(roc_auc_score(y, f1.predict(x1)))
            rows.append(
                {
                    "feature": feature,
                    "n": len(m),
                    "n_NR": int(y.sum()),
                    "senmayo_auc": auc0,
                    "senmayo_plus_feature_auc": auc1,
                    "delta_auc": auc1 - auc0,
                    "feature_beta_adjusted_for_senmayo": float(f1.params[feature]),
                    "feature_p_adjusted_for_senmayo": float(f1.pvalues[feature]),
                    "lrt_p_increment": p_lrt,
                    "status": "ok",
                }
            )
        except Exception as exc:
            rows.append({"feature": feature, "n": len(m), "status": f"failed: {type(exc).__name__}"})
    return pd.DataFrame(rows)


def plot_results(stats_df: pd.DataFrame, inc: pd.DataFrame) -> Path:
    keep = stats_df[
        stats_df["cohort"].isin(COHORTS)
        & (
            stats_df["feature"].isin(["Myeloid", "Fibroblast", "neutrophil_marker_score", "refractory_module_score", "cxcr2_expr_z", "senmayo"])
        )
    ].copy()
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
    mat = keep.pivot_table(index="feature", columns="cohort", values="auc_NR_high", aggfunc="mean")
    sns.heatmap(mat, vmin=0.25, vmax=0.75, center=0.5, cmap="vlag", annot=True, fmt=".2f", ax=axes[0], cbar_kws={"shrink": 0.7})
    axes[0].set_title("Baseline non-response AUC")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("")
    pooled = stats_df[stats_df["cohort"].eq("random_effects") & stats_df["feature"].isin(mat.index)].copy()
    pooled["err"] = 1.96 * pooled["logit_se"]
    order = pooled.sort_values("logit_beta_per_sd")["feature"]
    axes[1].barh(pooled.set_index("feature").loc[order].index, pooled.set_index("feature").loc[order, "logit_beta_per_sd"], color="#4c78a8")
    axes[1].errorbar(
        pooled.set_index("feature").loc[order, "logit_beta_per_sd"],
        np.arange(len(order)),
        xerr=pooled.set_index("feature").loc[order, "err"],
        fmt="none",
        ecolor="#333333",
        lw=1,
    )
    axes[1].axvline(0, color="#777777", lw=0.8)
    axes[1].set_title("Random-effects log OR per SD")
    axes[1].set_xlabel("higher feature -> NR")
    fig.tight_layout()
    p = OUT / "Fig_deconv_response.png"
    fig.savefig(p, dpi=220)
    plt.close(fig)
    return p


def write_summary(stats_df: pd.DataFrame, inc: pd.DataFrame, cell_counts: pd.DataFrame) -> None:
    pooled = stats_df[stats_df["cohort"].eq("random_effects")].copy()
    top_cell = pooled[pooled["feature_class"].eq("cell_fraction")].sort_values("logit_p").head(5)
    focus = pooled[pooled["feature"].isin(["Myeloid", "Fibroblast", "neutrophil_marker_score", "refractory_module_score", "cxcr2_expr_z", "senmayo"])]
    neut_present = cell_counts["cell_type"].astype(str).str.contains("Neut", case=False, na=False).any()
    lines = [
        "# Single-cell deconvolution response summary",
        "",
        "## Direct answer",
    ]
    if not top_cell.empty:
        best = top_cell.iloc[0]
        direction = "higher in non-response" if best["logit_beta_per_sd"] > 0 else "lower in non-response"
        lines.append(
            f"- Strongest NNLS cell-fraction signal: **{best['feature']}** ({direction}; random-effects OR/SD={best['logit_OR_per_sd']:.2f}, p={best['logit_p']:.3g})."
        )
    for name in ["Myeloid", "Fibroblast", "neutrophil_marker_score", "refractory_module_score", "cxcr2_expr_z", "senmayo"]:
        r = focus[focus["feature"].eq(name)]
        if len(r):
            x = r.iloc[0]
            lines.append(f"- {name}: OR/SD={x['logit_OR_per_sd']:.2f}, p={x['logit_p']:.3g}, I2={x.get('meta_I2', np.nan):.2f}.")
    lines += [
        "",
        "## Triangulation",
        f"- Explicit neutrophil cluster in the local scRNA references: **{'yes' if neut_present else 'no'}**. Therefore neutrophils are not interpreted as a learned NNLS fraction; they are tested as a bulk CXCR2/neutrophil marker score.",
        "- CXCR2 remains genetically protective in neutrophil/blood contexts, so a high neutrophil/CXCR2 mucosal marker signal is interpreted as disease-state/refractory biology, not proof that blocking CXCR2 is beneficial.",
        "- The refractory module score is a targeted bulk triangulation of the neutrophil + OSM-fibroblast + myeloid secretion module rather than a cell fraction.",
        "",
        "## Increment over bulk SASP",
    ]
    if len(inc):
        show = inc.sort_values("delta_auc", ascending=False).head(8)
        lines += ["```tsv", show.to_csv(sep="\t", index=False).strip(), "```"]
    else:
        lines.append("- Incremental models were not estimable.")
    lines += [
        "",
        "## Method caveats",
        "- NNLS uses broad single-cell signatures and marker-restricted rank targets to reduce cross-platform scale mismatch. Fractions should be read as relative cell-type enrichment estimates, not absolute histologic percentages.",
        "- Smillie UC has detailed author cell labels; Martin CD labels are broad marker-derived labels from local 10x matrices.",
        "- Missing neutrophils in the scRNA reference limits direct neutrophil deconvolution; neutrophil evidence comes from bulk marker/CXCR2 scoring.",
    ]
    (OUT / "SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_status(stats_df: pd.DataFrame) -> None:
    pooled = stats_df[stats_df["cohort"].eq("random_effects")]
    focus = pooled[pooled["feature"].isin(["Myeloid", "Fibroblast", "neutrophil_marker_score", "refractory_module_score", "senmayo"])]
    bullets = []
    for r in focus.itertuples(index=False):
        bullets.append(f"  - {r.feature}: OR/SD={r.logit_OR_per_sd:.2f}, p={r.logit_p:.3g}")
    block = (
        "\n## 2026-06-23 JST - single-cell deconvolution response\n\n"
        "- Ran `src/29_deconv_response.py` CPU-only with marker-restricted NNLS from `outputs/sc_01/*.h5ad` into five baseline bulk response cohorts.\n"
        "- Wrote `outputs/deconv/{deconv_proportions.tsv,celltype_fraction_vs_response.tsv,targeted_scores.tsv,Fig_deconv_response.png,SUMMARY.md}` and promoted final table/figure copies into `results/`.\n"
        "- Key random-effects response associations:\n"
        + "\n".join(bullets)
        + "\n- Caveat: local scRNA references do not contain an explicit neutrophil cluster, so neutrophil evidence is marker/CXCR2 based rather than an NNLS neutrophil fraction.\n"
    )
    for path in [P.journal / "status" / "task3_status.md", P.journal / "status" / "overnight_autorun_log.md"]:
        with path.open("a", encoding="utf-8") as fh:
            fh.write(block)


def main() -> None:
    t0 = time.perf_counter()
    log("start deconvolution")
    sig, cell_counts, _ = build_all_signatures()
    markers = marker_genes(sig)
    results = Parallel(n_jobs=min(CFG.n_jobs, len(COHORTS)), prefer="processes")(
        delayed(deconvolve_cohort)(cohort, sig, markers) for cohort in COHORTS
    )
    prop = pd.concat([x[0] for x in results if len(x[0])], ignore_index=True)
    scores = pd.concat([x[1] for x in results if len(x[1])], ignore_index=True)
    prop.to_csv(OUT / "deconv_proportions.tsv", sep="\t", index=False)
    scores.to_csv(OUT / "targeted_scores.tsv", sep="\t", index=False)

    frac_features = prop.rename(columns={"cell_type": "feature", "fraction": "value"})
    frac_features["feature_class"] = "cell_fraction"
    frac_features = frac_features[["cohort", "gsm", "patient", "response", "therapy", "feature_class", "feature", "value"]]
    score_features = scores.melt(
        id_vars=["cohort", "gsm", "patient", "response", "therapy"],
        value_vars=["senmayo", "inflam_crp", "neutrophil_proxy", "inflammation_score", "cxcr2_expr_z", "refractory_module_score", "neutrophil_marker_score"],
        var_name="feature",
        value_name="value",
    )
    score_features["feature_class"] = np.where(score_features["feature"].isin(["senmayo", "inflam_crp", "neutrophil_proxy", "inflammation_score"]), "bulk_score", "bulk_marker_score")
    features = pd.concat([frac_features, score_features], ignore_index=True)
    features["nr"] = features["response"].eq("NR").astype(int)
    features.to_csv(OUT / "deconv_features_long.tsv", sep="\t", index=False)
    stats_df = feature_tests(features)
    stats_df.to_csv(OUT / "celltype_fraction_vs_response.tsv", sep="\t", index=False)
    inc = pooled_increment(features)
    inc.to_csv(OUT / "deconv_incremental_vs_senmayo.tsv", sep="\t", index=False)
    fig = plot_results(stats_df, inc)
    write_summary(stats_df, inc, cell_counts)

    for f in [
        OUT / "celltype_fraction_vs_response.tsv",
        OUT / "deconv_proportions.tsv",
        OUT / "targeted_scores.tsv",
        OUT / "deconv_incremental_vs_senmayo.tsv",
    ]:
        P.promote_table(f)
    P.promote_figure(fig)
    append_status(stats_df)
    log(f"done in {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
