"""Task 3 parse GEO matrices and score senescence/inflammation signatures."""
from __future__ import annotations

import gzip
import json
import os
import re
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from pydantic import BaseModel, ConfigDict

from paths import P


RAW_GEO = P.raw / "geo"
OUT = P.out("23_score_all")


class CohortConfig(BaseModel):
    model_config = ConfigDict(frozen=True)

    accession: str
    platform: str
    therapy_filter: tuple[str, ...] = ()


COHORTS = [
    CohortConfig(accession="GSE16879", platform="GPL570"),
    CohortConfig(accession="GSE73661", platform="GPL6244"),
    CohortConfig(accession="GSE12251", platform="GPL570"),
    CohortConfig(accession="GSE23597", platform="GPL570"),
    CohortConfig(accession="GSE92415", platform="GPL13158"),
]

SEN_UP = ["CDKN1A", "CDKN2A", "CDKN2B", "GLB1", "SERPINE1"]
SEN_DOWN = ["MKI67", "LMNB1"]
INFLAM_CRP_LIKE = [
    "IL1B",
    "IL6",
    "TNF",
    "CXCL1",
    "CXCL2",
    "CXCL3",
    "CXCL8",
    "CCL2",
    "PTGS2",
    "NFKBIA",
    "ICAM1",
    "SELE",
    "VCAM1",
    "S100A8",
    "S100A9",
    "LCN2",
    "REG1A",
    "REG3A",
    "DUOX2",
]
NEUTROPHIL_PROXY = [
    "S100A8",
    "S100A9",
    "FCGR3B",
    "CSF3R",
    "CXCR1",
    "CXCR2",
    "MMP8",
    "MMP9",
    "ELANE",
    "MPO",
    "CEACAM8",
    "LCN2",
]


def first_symbol(raw: str) -> str:
    for part in re.split(r"///|//|;|,", str(raw)):
        sym = part.strip()
        if sym and sym not in {"---", "NA", "nan"}:
            return sym
    return ""


def load_probe2gene(gpl: str) -> dict[str, str]:
    path = RAW_GEO / f"{gpl}.annot.gz"
    if not path.exists():
        raise FileNotFoundError(path)
    rows: dict[str, str] = {}
    with gzip.open(path, "rt", errors="ignore") as fh:
        header = None
        for line in fh:
            if not line.strip() or line[0] in "^#!":
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                lower = [h.lower() for h in header]
                i_id = lower.index("id")
                try:
                    i_sym = lower.index("gene symbol")
                except ValueError:
                    i_sym = lower.index("gene_symbol")
                continue
            if len(parts) <= max(i_id, i_sym):
                continue
            sym = first_symbol(parts[i_sym])
            if sym:
                rows[parts[i_id]] = sym
    return rows


def parse_series_matrix(gse: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = RAW_GEO / f"{gse}_series_matrix.txt.gz"
    if not path.exists():
        raise FileNotFoundError(path)
    gsm: list[str] = []
    titles: list[str] = []
    source: list[str] = []
    char_rows: list[tuple[str, list[str]]] = []
    table: list[list[str]] = []
    in_table = False
    with gzip.open(path, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                gsm = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                source = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                field = vals[0].split(":", 1)[0].strip() if vals and ":" in vals[0] else "char"
                cleaned = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
                char_rows.append((field, cleaned))
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif line.startswith("!series_matrix_table_end"):
                in_table = False
            elif in_table:
                table.append(line.rstrip("\n").split("\t"))

    if not table:
        raise ValueError(f"{gse}: no expression table found")
    thead = [x.strip('"') for x in table[0]]
    probes = [r[0].strip('"') for r in table[1:]]
    data = np.array(
        [[np.nan if v in {"", '""'} else float(v.strip('"')) for v in r[1:]] for r in table[1:]],
        dtype=np.float32,
    )
    expr = pd.DataFrame(data, index=probes, columns=thead[1:])
    meta = pd.DataFrame({"gsm": gsm, "title": titles})
    if source:
        meta["source"] = source
    seen: dict[str, int] = {}
    for field, vals in char_rows:
        col = field
        if col in seen:
            seen[col] += 1
            col = f"{field}_{seen[col]}"
        else:
            seen[col] = 0
        meta[col] = vals
    meta = meta.set_index("gsm")
    expr = expr[[c for c in meta.index if c in expr.columns]]
    return expr, meta


def collapse_to_genes(expr: pd.DataFrame, p2g: dict[str, str]) -> pd.DataFrame:
    if np.nanmax(expr.values) > 100:
        expr = np.log2(expr.clip(lower=0) + 1)
    expr = expr[expr.index.isin(p2g)].copy()
    expr.index = [p2g[p] for p in expr.index]
    expr["__mean__"] = expr.mean(axis=1)
    expr = expr.sort_values("__mean__", ascending=False)
    return expr[~expr.index.duplicated(keep="first")].drop(columns="__mean__")


def map_response(s: pd.Series) -> pd.Series:
    vals = s.astype(str).str.strip().str.lower()
    return np.where(vals.isin(["yes", "y", "responder", "response", "r"]), "R", np.where(vals.isin(["no", "n", "non-responder", "nonresponder", "nr"]), "NR", "NA"))


def normalize_meta(gse: str, meta: pd.DataFrame) -> pd.DataFrame:
    m = pd.DataFrame(index=meta.index)
    m["cohort"] = gse
    if gse == "GSE16879":
        m["patient"] = meta["title"].str.replace("_beforeT", "", regex=False).str.replace("_afterT", "", regex=False)
        m["group"] = meta["disease"].map({"Control": "Control", "UC": "UC", "CD": "CD"})
        m["tissue"] = meta["tissue"]
        tp = meta["before or after first infliximab treatment"]
        m["timepoint"] = np.where(tp.str.startswith("Before"), "baseline", np.where(tp.str.startswith("After"), "post", "control"))
        m["response"] = map_response(meta["response to infliximab"])
        m["therapy"] = np.where(m["group"] == "Control", "none", "IFX")
    elif gse == "GSE73661":
        m["patient"] = "p" + meta["study individual number"].astype(str)
        wk = meta["week (w)"].astype(str).str.upper()
        is_ctrl = (wk == "CO") | (meta["mayo endoscopic subscore"] == "CO")
        m["group"] = np.where(is_ctrl, "Control", "UC")
        m["tissue"] = "Colon"
        m["timepoint"] = np.where(is_ctrl, "control", np.where(wk == "W0", "baseline", "post"))
        m["week"] = wk
        th = meta["induction therapy_maintenance therapy"].astype(str)
        m["therapy"] = np.where(is_ctrl, "none", np.where(th == "IFX", "IFX", np.where(th.str.startswith("vdz"), "VDZ", np.where(th.str.startswith("plac"), "placebo", "other"))))
        m["mayo"] = pd.to_numeric(meta["mayo endoscopic subscore"], errors="coerce")
        healed: dict[str, str] = {}
        for pid, sub in m.groupby("patient"):
            post = sub[sub["timepoint"] == "post"]
            if len(post) and post["mayo"].notna().any():
                healed[pid] = "R" if post["mayo"].min() <= 1 else "NR"
        m["response"] = m["patient"].map(healed).fillna("NA")
    elif gse == "GSE12251":
        m["patient"] = meta["title"].str.extract(r"^(P\d+)", expand=False).fillna(meta["title"])
        m["group"] = "UC"
        m["tissue"] = "Colon"
        m["timepoint"] = "baseline"
        m["response"] = map_response(meta["WK8RSPHM"])
        m["therapy"] = "IFX"
        m["dose"] = meta["title"].str.extract(r"/([^/]+mg/kg)/", expand=False)
    elif gse == "GSE23597":
        m["patient"] = meta["title"].str.extract(r"^(P\d+)", expand=False).fillna(meta["title"])
        m["group"] = "UC"
        m["tissue"] = "Colon"
        m["week"] = meta["time"].astype(str).str.upper()
        m["timepoint"] = np.where(m["week"] == "W0", "baseline", "post")
        m["dose"] = meta["dose"]
        m["therapy"] = np.where(meta["dose"].astype(str).str.lower() == "placebo", "placebo", "IFX")
        m["response"] = map_response(meta["wk8 response"])
        m["wk30_response"] = map_response(meta["wk30 response"])
    elif gse == "GSE92415":
        disease = meta["disease"].astype(str)
        m["patient"] = meta["title"]
        m["group"] = np.where(disease.str.contains("Healthy", case=False, na=False), "Control", "UC")
        m["tissue"] = "Colon"
        visit = meta["visit"].astype(str)
        m["timepoint"] = np.where(visit.eq("Week 0"), "baseline", np.where(visit.eq("Week 6"), "post", "control"))
        treatment = meta["treatment"].astype(str)
        m["therapy"] = np.where(treatment.str.contains("golimumab", case=False, na=False), "GLM", np.where(treatment.str.contains("placebo", case=False, na=False), "placebo", "none"))
        m["response"] = map_response(meta["wk6response"].fillna(""))
        m["mayo"] = pd.to_numeric(meta["mayo score"], errors="coerce")
    else:
        raise KeyError(gse)
    return m


def zscore_rows(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.DataFrame, list[str]]:
    present = [g for g in genes if g in expr.index]
    sub = expr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0)
    return z, present


def mean_z_score(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, int]:
    z, present = zscore_rows(expr, genes)
    return z.mean(axis=0), len(present)


def scaled(s: pd.Series) -> pd.Series:
    sd = s.std()
    return (s - s.mean()) / (sd if sd and not np.isnan(sd) else 1.0)


def score_cohort(gse: str, expr: pd.DataFrame, meta: pd.DataFrame, sets: dict[str, list[str] | str]) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    m = normalize_meta(gse, meta)
    m = m.loc[[c for c in expr.columns if c in m.index]].copy()
    expr = expr[m.index]

    old_path = P.interim / f"{gse}_scored.tsv"
    if old_path.exists():
        old = pd.read_csv(old_path, sep="\t", index_col=0)
        for col in ["lp", "age_accel"]:
            if col in old and col not in m:
                m[col] = old[col].reindex(m.index)

    coverage: list[dict[str, object]] = []
    signatures = {
        "senmayo": list(sets["SenMayo"]),
        "inflam_crp": INFLAM_CRP_LIKE,
        "neutrophil_proxy": NEUTROPHIL_PROXY,
    }
    for name, genes in signatures.items():
        score, n = mean_z_score(expr, genes)
        m[name] = score.reindex(m.index)
        coverage.append({"cohort": gse, "score": name, "n_present": n, "n_total": len(genes)})

    z_up, g_up = zscore_rows(expr, SEN_UP)
    z_down, g_down = zscore_rows(expr, SEN_DOWN)
    m["core_sen"] = (z_up.mean(axis=0) - z_down.mean(axis=0)).reindex(m.index)
    coverage.append({"cohort": gse, "score": "core_sen_up", "n_present": len(g_up), "n_total": len(SEN_UP)})
    coverage.append({"cohort": gse, "score": "core_sen_down", "n_present": len(g_down), "n_total": len(SEN_DOWN)})

    for marker in ["CDKN1A", "CDKN2A", "CDKN2B", "MKI67"]:
        if marker in expr.index:
            m[marker] = scaled(expr.loc[marker]).reindex(m.index)

    m["inflammation_score"] = pd.concat([scaled(m["inflam_crp"]), scaled(m["neutrophil_proxy"])], axis=1).mean(axis=1)
    return m, coverage


def process_cohort(cfg: CohortConfig) -> tuple[str, list[dict[str, object]], str, float]:
    t0 = time.perf_counter()
    p2g = load_probe2gene(cfg.platform)
    expr_raw, meta = parse_series_matrix(cfg.accession)
    expr = collapse_to_genes(expr_raw, p2g)
    expr.to_parquet(P.interim / f"{cfg.accession}_expr.parquet")
    meta.to_csv(P.interim / f"{cfg.accession}_meta.tsv", sep="\t", index_label="gsm")
    sets = json.load((P.external / "genesets" / "senescence_sets.json").open())
    scored, cov = score_cohort(cfg.accession, expr, meta, sets)
    scored.to_csv(P.interim / f"{cfg.accession}_scored.tsv", sep="\t", index_label="gsm")
    base = scored[(scored["timepoint"] == "baseline") & (scored["response"].isin(["R", "NR"]))]
    counts = base.groupby(["therapy", "response"]).size().to_string() if len(base) else "no labeled baseline rows"
    msg = f"{cfg.accession}: gene-level expr {expr.shape}; scored {scored.shape}\n{counts}"
    return cfg.accession, cov, msg, time.perf_counter() - t0


def main() -> None:
    sets = json.load((P.external / "genesets" / "senescence_sets.json").open())
    coverage_rows: list[dict[str, object]] = []
    n_jobs = min(int(os.environ.get("TASK3_N_JOBS", "4")), len(COHORTS), os.cpu_count() or 1)
    t0 = time.perf_counter()
    results = Parallel(n_jobs=n_jobs, prefer="processes")(delayed(process_cohort)(cfg) for cfg in COHORTS)
    for accession, cov, msg, elapsed in results:
        coverage_rows.extend(cov)
        print(f"\n===== {accession} =====")
        print(msg)
        print(f"elapsed_seconds={elapsed:.1f}")

    cov = pd.DataFrame(coverage_rows)
    cov["coverage"] = cov["n_present"] / cov["n_total"]
    cov_file = OUT / "task3_score_coverage.tsv"
    root_file = P.outputs / "task3_score_coverage.tsv"
    cov.to_csv(cov_file, sep="\t", index=False)
    cov.to_csv(root_file, sep="\t", index=False)
    print(f"\nwrote {cov_file}")
    print(f"wrote {root_file}")
    print(f"total_elapsed_seconds={time.perf_counter() - t0:.1f}; n_jobs={n_jobs}; senescence_sets_loaded={len(sets)}")


if __name__ == "__main__":
    main()
