"""Link the aging clock to cellular senescence / SASP (SenMayo).

- SenMayo single-sample score (mean of within-cohort z-scored member genes).
- Show clock age-acceleration correlates with SenMayo  (clock captures senescence).
- SenMayo elevated in IBD baseline and reverses in responders (parallels clock).
- In GTEx normal gut: does SenMayo rise with age? (senescence-aging link).
- Leading SASP genes that fall with treatment response.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/ujs/Downloads/lzy"
PROC = f"{BASE}/data/interim"
RES = f"{BASE}/outputs"

sets = json.load(open(f"{BASE}/data/external/genesets/senescence_sets.json"))
SENMAYO = sets["SenMayo"]


def ss_score(expr, genes):
    g = [x for x in genes if x in expr.index]
    sub = expr.loc[g]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0)
    return z.mean(axis=0), len(g)

out = {}
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))

# ---- GTEx normal gut: SenMayo vs age ----
gtex = pd.read_parquet(f"{BASE}/data/raw/gtex/gut_expr_log2.parquet")
gut = pd.read_csv(f"{BASE}/data/raw/gtex/gut_samples.tsv", sep="\t").dropna(subset=["age_mid"])
gtex = gtex[[c for c in gut["SAMPID"] if c in gtex.columns]]
sm_gtex, n_sm = ss_score(gtex, SENMAYO)
age_map = dict(zip(gut["SAMPID"], gut["age_mid"]))
ages = np.array([age_map[c] for c in gtex.columns])
rho, pg = stats.spearmanr(ages, sm_gtex.values)
out["GTEx_SenMayo_vs_age_rho"] = float(rho); out["GTEx_SenMayo_vs_age_p"] = float(pg)
print(f"GTEx normal gut: SenMayo vs age Spearman rho={rho:.3f}, p={pg:.2e} ({n_sm} genes)")
ax = axes[0]
ax.scatter(ages, sm_gtex.values, s=10, alpha=0.4, color="#756bb1")
ax.set_xlabel("Chronological age"); ax.set_ylabel("SenMayo score")
ax.set_title(f"GTEx normal gut\nSenMayo rises with age (rho={rho:.2f})")

# ---- per cohort: clock vs SenMayo, and reversal ----
for j, gse in enumerate(["GSE16879", "GSE73661"], start=1):
    expr = pd.read_parquet(f"{PROC}/{gse}_expr.parquet")
    sc = pd.read_csv(f"{PROC}/{gse}_scored.tsv", sep="\t", index_col=0)
    sm, ng = ss_score(expr, SENMAYO)
    sc["senmayo"] = sm.reindex(sc.index)
    sc.to_csv(f"{PROC}/{gse}_scored.tsv", sep="\t")  # add senmayo column back
    r, p = stats.spearmanr(sc["age_accel"], sc["senmayo"], nan_policy="omit")
    out[f"{gse}_clock_vs_senmayo_rho"] = float(r)
    print(f"{gse}: clock age-accel vs SenMayo Spearman rho={r:.3f}, p={p:.2e}")
    # baseline group difference
    base = sc[sc["timepoint"].isin(["baseline", "control"])]
    bg = base.groupby("group")["senmayo"].median()
    print(f"  SenMayo median by group: {bg.to_dict()}")

    if j == 1:  # plot clock-senmayo concordance for primary cohort
        ax = axes[1]
        for grp, col in [("Control", "#1f77b4"), ("UC", "#ff7f0e"), ("CD", "#2ca02c")]:
            s = sc[sc["group"] == grp]
            ax.scatter(s["senmayo"], s["age_accel"], s=14, alpha=0.6, label=grp, color=col)
        ax.set_xlabel("SenMayo score"); ax.set_ylabel("Clock age acceleration")
        ax.legend(fontsize=8)
        ax.set_title(f"GSE16879: clock vs senescence\nrho={r:.2f}")

# ---- leading SASP genes reversing with response (GSE16879) ----
expr = pd.read_parquet(f"{PROC}/GSE16879_expr.parquet")
sc = pd.read_csv(f"{PROC}/GSE16879_scored.tsv", sep="\t", index_col=0)
def grp_samples(cond):
    return sc[cond].index.intersection(expr.columns)
sm_genes = [g for g in SENMAYO if g in expr.index]
deltas = {}
for g in sm_genes:
    rb = expr.loc[g, grp_samples((sc.response == "R") & (sc.timepoint == "baseline"))].mean()
    rp = expr.loc[g, grp_samples((sc.response == "R") & (sc.timepoint == "post"))].mean()
    deltas[g] = rp - rb
ds = pd.Series(deltas).sort_values()
out["top_down_SASP_in_responders"] = ds.head(12).round(3).to_dict()
print("\nTop SASP genes DOWN in responders post-treatment:")
print(ds.head(12).to_string())
ax = axes[2]
top = ds.head(12)
ax.barh(range(len(top)), top.values, color="#2ca02c")
ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontsize=8)
ax.invert_yaxis(); ax.axvline(0, color="grey", lw=0.8)
ax.set_xlabel("Δ expression (post − baseline, responders)")
ax.set_title("SASP genes resolved by\nresponse (GSE16879)")

fig.suptitle("Cellular senescence / SASP underlies the reversible aging signal",
             fontweight="bold")
fig.tight_layout()
fig.savefig(f"{RES}/Fig4_senescence.png", dpi=200)
json.dump(out, open(f"{RES}/senescence_stats.json", "w"), indent=2)
print("\nsaved results/Fig4_senescence.png and senescence_stats.json")
