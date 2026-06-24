"""Phenome-wide safety map for the CXCL8 / IL-8-receptor (CXCR1/CXCR2) axis.

The manuscript argues that genetics cautions against blocking the IL-8-receptor
axis in IBD because higher receptor expression is *protective* (blood eQTLGen MR
OR 0.753, coloc 0.951; FinnGen IBD/UC replication). A reviewer will ask: what
else does this locus do? If the axis is a broad neutrophil/infection/immune
homeostasis hub, blunt antagonism is more likely to carry off-target risk than
to help IBD -- which is exactly the safety message.

This script builds that evidence cheaply and reproducibly: it pulls the FinnGen
R12 *pheweb* phenome-wide scan for the axis lead variant(s) -- a single ~1 MB
JSON per variant, ~2,470 disease endpoints each, NO large summary-stat downloads.

  blood-axis lead : rs62183956  hg38 chr2:218,181,399  (eQTLGen CXCR1/CXCR2 lead; the IBD coloc variant)
  neutrophil lead : rs6737563   hg38 chr2:218,080,951  (BLUEPRINT neutrophil eQTL lead; robustness)

Everything is oriented to the IBD-PROTECTIVE allele (the direction a drug that
*raised* axis activity would mimic), using FinnGen's own IBD endpoint as the
anchor, so a negative oriented beta = "the protective direction also lowers this
trait" and positive = "the protective direction raises it".

Honesty: this is a single-variant, Finnish-only PheWAS. It maps pleiotropy /
shared genetic architecture, NOT proven causality for every endpoint. The IBD
anchor must reproduce the known protective direction or orientation is suspect.

Run:  .venv/bin/python src/49_phewas_axis.py   (from repo root)
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict
from statsmodels.stats.multitest import multipletests

from paths import P

FINNGEN_API = "https://r12.finngen.fi/api/variant"
FINNGEN_N_ENDPOINTS = 2470  # full R12 catalog; API returns only p<0.05 hits but denominator must be total

# Axis instruments. `fg_ids` are FinnGen hg38 chr-pos-ref-alt candidates tried in
# order (the locus is multiallelic; the first that returns results is used).
VARIANTS = [
    {
        "label": "blood_lead_rs62183956",
        "rsid": "rs62183956",
        "fg_ids": ["2-218181399-C-T", "2-218181399-C-A"],
        "note": "eQTLGen CXCR1/CXCR2 blood lead; the IBD-colocalising variant (PP4 0.951).",
    },
    {
        "label": "neutrophil_lead_rs6737563",
        "rsid": "rs6737563",
        "fg_ids": ["2-218080951-T-C"],
        "note": "BLUEPRINT neutrophil eQTL lead (cell-type-matched); robustness.",
    },
]

# FinnGen IBD endpoint phenocodes used as the orientation anchor (first present wins).
IBD_ANCHORS = ["K11_IBD_STRICT", "K11_IBD", "K11_CD_STRICT2", "K11_UC_STRICT2", "CD_STRICT", "UC_STRICT"]

# Mechanistic category buckets that, if enriched among significant hits, support
# the "neutrophil / infection / immune homeostasis axis" reading.
MECH_BUCKETS = {
    "infection": ["infectious", "AB1_", "infection", "pneumonia", "sepsis", "abscess"],
    "respiratory": ["respiratory", "J10_", "asthma", "copd", "bronch", "lung", "pulmonary"],
    "autoimmune_inflammatory": ["autoimmun", "rheumat", "psoriasis", "spondyl", "inflammatory",
                                "L12_", "M13_", "K11_", "sarcoid", "connective"],
    "blood_immune": ["blood", "haematolog", "hematolog", "D3_", "neutrophil", "leukocyte",
                     "immune mechanism", "immunodef"],
    "skin": ["skin", "L12_", "dermat", "hidradenitis", "acne"],
}


class Inputs(BaseModel):
    model_config = ConfigDict(frozen=True)
    out_dir: Path = P.out("49_phewas_axis")


def fetch_phewas(fg_ids: list[str], cache: Path) -> tuple[list[dict], str]:
    if cache.exists() and cache.stat().st_size > 0:
        d = json.loads(cache.read_text())
        return d.get("results", []), d.get("_fg_id_used", "")
    last_err = None
    for fg_id in fg_ids:
        url = f"{FINNGEN_API}/{fg_id}"
        for attempt in range(4):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "research-phewas/1.0"})
                with urllib.request.urlopen(req, timeout=90) as r:
                    d = json.load(r)
                if d.get("results"):
                    d["_fg_id_used"] = fg_id
                    cache.write_text(json.dumps(d))
                    return d["results"], fg_id
                break  # valid response but no results -> try next candidate id
            except urllib.error.HTTPError as exc:
                last_err = exc
                if exc.code == 404:
                    break  # try next candidate id
                time.sleep(5 * (attempt + 1))
            except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
                last_err = exc
                time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"no FinnGen PheWAS for {fg_ids}: {last_err}")


def to_frame(results: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(results)
    for col in ["beta", "sebeta", "pval", "mlogp", "n_case", "n_control", "maf"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def orient(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Flip all betas so the reference direction is the IBD-protective allele."""
    anchor_row, anchor_code = None, None
    for code in IBD_ANCHORS:
        hit = df[df["phenocode"].astype(str) == code]
        if len(hit) and pd.notna(hit.iloc[0].get("beta")):
            anchor_row, anchor_code = hit.iloc[0], code
            break
    if anchor_row is None:
        # fall back to any phenostring containing "inflammatory bowel"
        hit = df[df["phenostring"].astype(str).str.contains("inflammatory bowel", case=False, na=False)]
        if len(hit):
            anchor_row, anchor_code = hit.sort_values("pval").iloc[0], str(hit.sort_values("pval").iloc[0]["phenocode"])
    if anchor_row is None:
        flip = 1.0
        anchor_info = {"anchor": "none", "anchor_beta": np.nan, "anchor_p": np.nan, "flip": flip}
    else:
        # protective = lowers IBD risk -> we want anchor oriented beta < 0
        flip = -1.0 if anchor_row["beta"] > 0 else 1.0
        anchor_info = {
            "anchor": anchor_code,
            "anchor_beta_reported": float(anchor_row["beta"]),
            "anchor_p": float(anchor_row["pval"]),
            "flip": flip,
            "anchor_oriented_beta": float(anchor_row["beta"] * flip),
        }
    out = df.copy()
    out["oriented_beta"] = out["beta"] * flip
    return out, anchor_info


def bucket(row: pd.Series) -> str:
    text = (str(row.get("category", "")) + " " + str(row.get("phenostring", "")) + " "
            + str(row.get("phenocode", ""))).lower()
    for name, keys in MECH_BUCKETS.items():
        if any(k.lower() in text for k in keys):
            return name
    return "other"


def annotate(df: pd.DataFrame) -> pd.DataFrame:
    out = df.dropna(subset=["pval"]).copy()
    out["mech_bucket"] = out.apply(bucket, axis=1)
    out["bonf_sig"] = out["pval"] < (0.05 / FINNGEN_N_ENDPOINTS)
    out["fdr"] = multipletests(out["pval"], method="fdr_bh")[1]
    out["fdr_sig"] = out["fdr"] < 0.05
    return out.sort_values("pval")


def per_variant(v: dict, cfg: Inputs) -> tuple[pd.DataFrame, dict]:
    cache = cfg.out_dir / f"phewas_{v['rsid']}.json"
    results, fg_id = fetch_phewas(v["fg_ids"], cache)
    df = to_frame(results)
    df, anchor = orient(df)
    df = annotate(df)
    df.insert(0, "variant", v["label"])
    df.insert(1, "fg_id", fg_id)
    anchor["variant"] = v["label"]
    anchor["fg_id"] = fg_id
    anchor["n_tested"] = int(len(df))
    anchor["n_bonf_sig"] = int(df["bonf_sig"].sum())
    anchor["n_fdr_sig"] = int(df["fdr_sig"].sum())
    return df, anchor


def plot_phewas(df: pd.DataFrame, anchor: dict, path: Path) -> None:
    d = df.dropna(subset=["mlogp"]).copy()
    cats = (d.groupby("mech_bucket")["mlogp"].max().sort_values(ascending=False).index.tolist())
    order = [c for c in ["infection", "respiratory", "autoimmune_inflammatory", "blood_immune", "skin", "other"]
             if c in cats]
    palette = {"infection": "#E64B35", "respiratory": "#4DBBD5", "autoimmune_inflammatory": "#00A087",
               "blood_immune": "#3C5488", "skin": "#F39B7F", "other": "#B0B0B0"}
    d["_cat"] = pd.Categorical(d["mech_bucket"], categories=order, ordered=True)
    d = d.sort_values(["_cat", "pval"])
    d["x"] = np.arange(len(d))
    bonf = -np.log10(0.05 / FINNGEN_N_ENDPOINTS)

    fig, ax = plt.subplots(figsize=(11, 4.6))
    for cat in order:
        sub = d[d["mech_bucket"] == cat]
        ax.scatter(sub["x"], sub["mlogp"], s=14, c=palette[cat], label=cat.replace("_", "/"),
                   edgecolor="none", alpha=0.85)
    ax.axhline(bonf, ls="--", lw=0.8, color="#444",
               label=f"Bonferroni 0.05/{FINNGEN_N_ENDPOINTS}")
    # annotate top significant hits + the IBD anchor
    top = d[d["bonf_sig"]].sort_values("pval").head(14)
    for _, r in top.iterrows():
        ax.annotate(str(r["phenostring"])[:34], (r["x"], r["mlogp"]),
                    fontsize=6.2, rotation=32, ha="left", va="bottom",
                    xytext=(1, 2), textcoords="offset points")
    ax.set_ylabel(r"$-\log_{10}P$  (FinnGen R12)")
    ax.set_xlabel("phenotype endpoints, grouped by mechanistic category")
    ax.set_xticks([])
    ttl = (f"CXCR1/CXCR2 axis PheWAS  ({anchor.get('variant','')}, {anchor.get('fg_id','')})  "
           f"IBD anchor {anchor.get('anchor','?')} p={anchor.get('anchor_p', float('nan')):.1e}")
    ax.set_title(ttl, fontsize=9)
    ax.legend(fontsize=7, frameon=False, ncol=3, loc="upper right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_summary(frames: list[pd.DataFrame], anchors: list[dict], out_dir: Path) -> None:
    lines = ["# CXCL8 / IL-8-receptor (CXCR1/CXCR2) axis -- phenome-wide safety map", ""]
    lines.append("Source: FinnGen R12 pheweb per-variant PheWAS API (single JSON per variant; no "
                 "large summary-stat downloads). Betas oriented to the IBD-protective allele via the "
                 "FinnGen IBD endpoint, so a negative oriented beta means the IBD-protective direction "
                 "also lowers that trait.\n")
    for df, anchor in zip(frames, anchors):
        lines.append(f"## {anchor['variant']}  (FinnGen {anchor['fg_id']})")
        lines.append(f"- IBD anchor = {anchor.get('anchor')}, reported beta="
                     f"{anchor.get('anchor_beta_reported', float('nan')):.3f}, p="
                     f"{anchor.get('anchor_p', float('nan')):.2e} -> orientation flip={anchor.get('flip')}. "
                     f"(Oriented IBD beta should be <0 = protective.)")
        lines.append(f"- Endpoints tested: {anchor['n_tested']}; Bonferroni-significant: "
                     f"{anchor['n_bonf_sig']}; FDR<0.05: {anchor['n_fdr_sig']}.")
        sig = df[df["bonf_sig"]]
        if len(sig):
            by_cat = sig["mech_bucket"].value_counts()
            lines.append(f"- Significant-hit categories: "
                         + ", ".join(f"{k}={v}" for k, v in by_cat.items()) + ".")
            lines.append("- Top Bonferroni hits (oriented to IBD-protective allele):")
            for _, r in sig.sort_values("pval").head(12).iterrows():
                arrow = "down" if r["oriented_beta"] < 0 else "up"
                lines.append(f"    - {r['phenostring']} ({r['phenocode']}): p={r['pval']:.2e}, "
                             f"oriented beta={r['oriented_beta']:+.3f} ({arrow}); cat={r['mech_bucket']}")
        else:
            lines.append("- No Bonferroni-significant endpoints (locus pleiotropy is narrow here).")
        lines.append("")
    lines.append("## Interpretation (honest)")
    lines.append("- If the significant hits concentrate in infection / respiratory / autoimmune / "
                 "blood-immune categories, the axis behaves as a broad neutrophil/immune homeostasis "
                 "hub -- consistent with the paper's warning that blunt IL-8-receptor antagonism is "
                 "genetically cautioned rather than a clean IBD opportunity.")
    lines.append("- Caveats: single-variant, Finnish-only PheWAS; maps shared genetic architecture, "
                 "not proven causality for every endpoint; orientation depends on the IBD anchor "
                 "reproducing the protective direction (checked above).")
    (out_dir / "SUMMARY.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    cfg = Inputs()
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    frames, anchors = [], []
    for v in VARIANTS:
        df, anchor = per_variant(v, cfg)
        df.to_csv(cfg.out_dir / f"phewas_{v['label']}.tsv", sep="\t", index=False)
        plot_phewas(df, anchor, cfg.out_dir / f"Fig_phewas_{v['label']}.png")
        frames.append(df)
        anchors.append(anchor)
        print(f"{anchor['variant']}: tested={anchor['n_tested']} bonf_sig={anchor['n_bonf_sig']} "
              f"fdr_sig={anchor['n_fdr_sig']} anchor={anchor.get('anchor')} "
              f"anchor_p={anchor.get('anchor_p', float('nan')):.2e}")

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(cfg.out_dir / "phewas_axis_all.tsv", sep="\t", index=False)
    sig = combined[combined["bonf_sig"]].sort_values(["variant", "pval"])
    sig.to_csv(cfg.out_dir / "phewas_axis_significant.tsv", sep="\t", index=False)
    pd.DataFrame(anchors).to_csv(cfg.out_dir / "phewas_anchors.tsv", sep="\t", index=False)
    write_summary(frames, anchors, cfg.out_dir)

    # promote the primary (blood-lead) figure + the significant-hit table
    P.promote_figure(cfg.out_dir / f"Fig_phewas_{VARIANTS[0]['label']}.png")
    P.promote_table(cfg.out_dir / "phewas_axis_significant.tsv")
    print(f"\nwrote outputs to {cfg.out_dir}")


if __name__ == "__main__":
    main()
