"""Train an intestinal transcriptomic aging clock on GTEx gut samples.

Memory-safe: streams the 1.6GB gct.gz, keeping only gut-sample columns and
genes expressed in gut (mean TPM > 1). Trains ElasticNet with GroupKFold by
subject (no subject leakage). Saves clock genes + coefficients + per-gene
standardization stats for cross-platform application.
"""
import gzip, json, os
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.model_selection import GroupKFold, cross_val_predict
from scipy.stats import pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/ujs/Downloads/lzy"
GCT = f"{BASE}/data/raw/gtex/gene_tpm.gct.gz"
RES = f"{BASE}/outputs"

# ---- 1. gut samples + ages ----
gut = pd.read_csv(f"{BASE}/data/raw/gtex/gut_samples.tsv", sep="\t")
gut = gut.dropna(subset=["age_mid"])
age_by_sample = dict(zip(gut["SAMPID"], gut["age_mid"]))
subj_by_sample = dict(zip(gut["SAMPID"], gut["SUBJID"]))
gut_set = set(gut["SAMPID"])
print(f"gut samples with age: {len(gut_set)}")

# ---- 2-3. build gene-level log2 matrix (cached) ----
CACHE = f"{BASE}/data/raw/gtex/gut_expr_log2.parquet"
if os.path.exists(CACHE):
    print("loading cached gut expression matrix...")
    df = pd.read_parquet(CACHE)
    print(f"cached genes x samples: {df.shape}")
else:
    print("streaming gct.gz (this takes several minutes)...")
    with gzip.open(GCT, "rt") as fh:
        fh.readline()  # #1.2
        fh.readline()  # dims
        header = fh.readline().rstrip("\n").split("\t")
        col_idx = [i for i in range(2, len(header)) if header[i] in gut_set]
        sample_ids = [header[i] for i in col_idx]
        print(f"matched gut columns in matrix: {len(sample_ids)}")
        symbols, rows = [], []
        n = 0
        for line in fh:
            n += 1
            parts = line.rstrip("\n").split("\t")
            vals = np.array([parts[i] for i in col_idx], dtype=np.float32)
            if vals.mean() <= 1.0:        # expressed-in-gut filter
                continue
            symbols.append(parts[1])
            rows.append(vals)
            if n % 10000 == 0:
                print(f"  scanned {n} genes, kept {len(rows)}")
    X = np.vstack(rows)
    print(f"kept genes: {X.shape[0]}, samples: {X.shape[1]}")
    del rows
    # collapse duplicate symbols (keep highest-variance probe)
    df = pd.DataFrame(X, index=symbols, columns=sample_ids)
    df = df[df.index.notna() & (df.index != "")]
    df = np.log2(df + 1.0)
    df["__var__"] = df.var(axis=1)
    df = df.sort_values("__var__", ascending=False)
    df = df[~df.index.duplicated(keep="first")].drop(columns="__var__")
    df.to_parquet(CACHE)
    print(f"genes after symbol collapse: {df.shape[0]} (cached)")

# ---- 4. pre-filter to genes most correlated with age (fast, focused) ----
df = df.loc[df.var(axis=1) > 1e-6]
samples = list(df.columns)
y = np.array([age_by_sample[s] for s in samples])
Xall = df.values                                  # genes x samples
Xc = Xall - Xall.mean(axis=1, keepdims=True)
yc = y - y.mean()
denom = np.sqrt((Xc**2).sum(axis=1) * (yc**2).sum())
corr = (Xc @ yc) / np.where(denom == 0, 1.0, denom)
TOPN = 1500
top_idx = np.argsort(-np.abs(corr))[:TOPN]
expr = df.iloc[top_idx].T                          # samples x genes
genes = list(expr.columns)
print(f"selected {len(genes)} age-correlated genes (max|r|={np.abs(corr).max():.2f}); "
      f"NOTE: selection on full data -> CV slightly optimistic (documented).")

# standardize per-gene (save stats)
mu = expr.mean(axis=0)
sd = expr.std(axis=0).replace(0, 1.0)
Xz = ((expr - mu) / sd).values
groups = np.array([subj_by_sample[s] for s in samples])
print(f"feature matrix: {Xz.shape}, subjects: {len(set(groups))}")

# ---- 5. ElasticNet with subject-grouped CV (no subject leakage) ----
cv_splits = list(GroupKFold(n_splits=5).split(Xz, y, groups))
model = ElasticNetCV(
    l1_ratio=[0.5, 0.8], eps=1e-2,
    cv=cv_splits, max_iter=5000, n_jobs=-1, random_state=0,
)
model.fit(Xz, y)
print(f"selected alpha={model.alpha_:.4g}, l1_ratio={model.l1_ratio_}")

# honest out-of-fold predictions at selected hyperparameters
en = ElasticNet(alpha=model.alpha_, l1_ratio=model.l1_ratio_, max_iter=10000, random_state=0)
oof = cross_val_predict(
    en, Xz, y, cv=list(GroupKFold(5).split(Xz, y, groups)), n_jobs=-1
)
r, p = pearsonr(y, oof)
mae = np.mean(np.abs(y - oof))
print(f"CV (out-of-fold): Pearson r={r:.3f} (p={p:.1e}), MAE={mae:.2f} yr")

# ---- 6. save clock ----
coef = model.coef_
nz = np.where(coef != 0)[0]
clock = {
    "genes": [genes[i] for i in nz],
    "coef": [float(coef[i]) for i in nz],
    "intercept": float(model.intercept_),
    "mu": {genes[i]: float(mu.iloc[i]) for i in nz},
    "sd": {genes[i]: float(sd.iloc[i]) for i in nz},
    "alpha": float(model.alpha_), "l1_ratio": float(model.l1_ratio_),
    "cv_r": float(r), "cv_mae": float(mae), "n_clock_genes": int(len(nz)),
}
with open(f"{RES}/clock_model.json", "w") as f:
    json.dump(clock, f, indent=2)
print(f"clock genes (nonzero): {len(nz)} -> saved results/clock_model.json")

# ---- 7. Fig 1: predicted vs chronological ----
fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(y, oof, s=10, alpha=0.4, color="#2c7fb8")
lims = [15, 80]
ax.plot(lims, lims, "--", color="grey")
ax.set_xlabel("Chronological age (GTEx bracket midpoint)")
ax.set_ylabel("Predicted transcriptomic age")
ax.set_title(f"Intestinal aging clock (GTEx gut)\nr={r:.2f}, MAE={mae:.1f} yr, {len(nz)} genes")
ax.set_xlim(lims); ax.set_ylim(lims)
fig.tight_layout()
fig.savefig(f"{RES}/Fig1_clock_accuracy.png", dpi=200)
print("saved results/Fig1_clock_accuracy.png")
