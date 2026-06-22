"""Parse GEO series matrices -> gene-level expression + sample metadata.

Maps Affy probes to gene symbols via GEO platform .annot.gz, collapses probes
to genes by max mean expression. Saves per-cohort expr (genes x samples) and
metadata. Also prints all characteristic field names for inspection.
"""
import gzip, sys, re
import numpy as np
import pandas as pd

BASE = "/Users/ujs/Downloads/lzy/data/raw/geo"
OUT = "/Users/ujs/Downloads/lzy/data/interim"
import os; os.makedirs(OUT, exist_ok=True)

COHORTS = {
    "GSE16879": "GPL570",
    "GSE73661": "GPL6244",
}


def load_probe2gene(gpl):
    """GEO .annot.gz: skip comment/caret lines, columns ID, Gene symbol."""
    rows = {}
    with gzip.open(f"{BASE}/{gpl}.annot.gz", "rt", errors="ignore") as fh:
        header = None
        for line in fh:
            if line[0] in "^#!":
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                i_id = header.index("ID")
                i_sym = header.index("Gene symbol")
                continue
            if len(parts) <= max(i_id, i_sym):
                continue
            sym = parts[i_sym].split("///")[0].strip()
            if sym:
                rows[parts[i_id]] = sym
    return rows


def parse_series_matrix(gse):
    path = f"{BASE}/{gse}_series_matrix.txt.gz"
    gsm, titles, source = [], [], []
    char_rows = []  # list of (field, [values])
    table = []
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
                field = vals[0].split(":")[0].strip() if ":" in vals[0] else "char"
                cleaned = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
                char_rows.append((field, cleaned))
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            elif line.startswith("!series_matrix_table_end"):
                in_table = False
            elif in_table:
                table.append(line.rstrip("\n").split("\t"))

    # expression table
    thead = [x.strip('"') for x in table[0]]
    probes = [r[0].strip('"') for r in table[1:]]
    data = np.array([[np.nan if v in ("", '""') else float(v.strip('"'))
                      for v in r[1:]] for r in table[1:]], dtype=np.float32)
    expr = pd.DataFrame(data, index=probes, columns=thead[1:])

    # metadata
    meta = pd.DataFrame({"gsm": gsm, "title": titles})
    if source:
        meta["source"] = source
    seen = {}
    for field, vals in char_rows:
        col = field
        if col in seen:
            seen[col] += 1; col = f"{field}_{seen[col]}"
        else:
            seen[col] = 0
        meta[col] = vals
    meta = meta.set_index("gsm")
    # align expr columns (GSM) with meta
    expr = expr[[c for c in meta.index if c in expr.columns]]
    return expr, meta


for gse, gpl in COHORTS.items():
    print(f"\n===== {gse} ({gpl}) =====")
    p2g = load_probe2gene(gpl)
    print(f"probe->gene map: {len(p2g)}")
    expr, meta = parse_series_matrix(gse)
    print(f"expr probes x samples: {expr.shape}; meta fields: {list(meta.columns)}")
    for c in meta.columns:
        if c in ("title",):
            continue
        u = meta[c].unique()
        if len(u) <= 12:
            print(f"  [{c}] -> {list(u)}")

    # log2 if looks linear
    if np.nanmax(expr.values) > 100:
        expr = np.log2(expr.clip(lower=0) + 1)
        print("  applied log2")

    # map to genes, collapse by max mean
    expr = expr[expr.index.isin(p2g)]
    expr.index = [p2g[p] for p in expr.index]
    expr["__m__"] = expr.mean(axis=1)
    expr = expr.sort_values("__m__", ascending=False)
    expr = expr[~expr.index.duplicated(keep="first")].drop(columns="__m__")
    print(f"  gene-level expr: {expr.shape}")

    expr.to_parquet(f"{OUT}/{gse}_expr.parquet")
    meta.to_csv(f"{OUT}/{gse}_meta.tsv", sep="\t")
    print(f"  saved {OUT}/{gse}_expr.parquet + meta")
