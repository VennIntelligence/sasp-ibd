"""Apply the GTEx-trained intestinal aging clock to IBD cohorts.

Cross-platform strategy: within each cohort, z-score every clock gene across
samples, form the linear predictor with the clock coefficients, then express
each sample's transcriptomic age acceleration RELATIVE TO HEALTHY CONTROLS
(control-anchored z). All key comparisons are within-cohort, so absolute
cross-platform calibration is not required.
"""
import json
import numpy as np
import pandas as pd

BASE = "/Users/ujs/Downloads/lzy"
PROC = f"{BASE}/data/interim"

clock = json.load(open(f"{BASE}/outputs/clock_model.json"))
cg = clock["genes"]
coef = pd.Series(clock["coef"], index=cg)
print(f"clock genes: {len(cg)}")


def normalize_meta_16879(meta):
    m = pd.DataFrame(index=meta.index)
    m["cohort"] = "GSE16879"
    m["patient"] = (meta["title"]
                    .str.replace("_beforeT", "", regex=False)
                    .str.replace("_afterT", "", regex=False))
    m["group"] = meta["disease"].map({"Control": "Control", "UC": "UC", "CD": "CD"})
    m["tissue"] = meta["tissue"]
    tp = meta["before or after first infliximab treatment"]
    m["timepoint"] = np.where(tp.str.startswith("Before"), "baseline",
                       np.where(tp.str.startswith("After"), "post", "control"))
    r = meta["response to infliximab"]
    m["response"] = np.where(r == "Yes", "R", np.where(r == "No", "NR", "NA"))
    m["therapy"] = np.where(m["group"] == "Control", "none", "IFX")
    return m


def normalize_meta_73661(meta):
    m = pd.DataFrame(index=meta.index)
    m["cohort"] = "GSE73661"
    m["patient"] = "p" + meta["study individual number"].astype(str)
    wk = meta["week (w)"].str.upper()
    is_ctrl = (wk == "CO") | (meta["mayo endoscopic subscore"] == "CO")
    m["group"] = np.where(is_ctrl, "Control", "UC")
    m["tissue"] = "Colon"
    m["timepoint"] = np.where(is_ctrl, "control",
                       np.where(wk == "W0", "baseline", "post"))
    m["week"] = wk
    th = meta["induction therapy_maintenance therapy"]
    m["therapy"] = np.where(is_ctrl, "none",
                    np.where(th == "IFX", "IFX",
                     np.where(th.str.startswith("vdz"), "VDZ",
                      np.where(th.str.startswith("plac"), "placebo", "other"))))
    mayo = pd.to_numeric(meta["mayo endoscopic subscore"], errors="coerce")
    m["mayo"] = mayo
    # response per patient = mucosal healing (mayo<=1) at any post timepoint
    healed = {}
    for pid, sub in m.assign(mayo=mayo).groupby("patient"):
        post = sub[sub["timepoint"] == "post"]
        if len(post) and post["mayo"].notna().any():
            healed[pid] = "R" if post["mayo"].min() <= 1 else "NR"
    m["response"] = m["patient"].map(healed).fillna("NA")
    return m


NORM = {"GSE16879": normalize_meta_16879, "GSE73661": normalize_meta_73661}

for gse in ["GSE16879", "GSE73661"]:
    expr = pd.read_parquet(f"{PROC}/{gse}_expr.parquet")
    meta = pd.read_csv(f"{PROC}/{gse}_meta.tsv", sep="\t", index_col=0)
    m = NORM[gse](meta)

    present = [g for g in cg if g in expr.index]
    cov = len(present) / len(cg)
    print(f"\n{gse}: clock-gene coverage {len(present)}/{len(cg)} = {cov:.1%}")

    sub = expr.loc[present]                       # genes x samples
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0)
    lp = z.mul(coef[present], axis=0).sum(axis=0)  # linear predictor per sample

    out = m.copy()
    out["lp"] = lp.reindex(out.index)
    ctrl = out.loc[out["group"] == "Control", "lp"]
    mu_c, sd_c = ctrl.mean(), (ctrl.std() or 1.0)
    out["age_accel"] = (out["lp"] - mu_c) / sd_c   # control-anchored z
    out.to_csv(f"{PROC}/{gse}_scored.tsv", sep="\t")
    print(f"  saved {gse}_scored.tsv  (n={len(out)})")
    # quick sanity: baseline IBD vs control
    base = out[out["timepoint"].isin(["baseline", "control"])]
    print(base.groupby("group")["age_accel"].agg(["mean", "count"]).to_string())
