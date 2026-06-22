"""Shared helpers for Task 2 single-cell scripts."""
from __future__ import annotations

import json
import pickle
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from paths import P


RNG = np.random.default_rng(20260622)

HEALTH_ORDER = ["Healthy", "Non-inflamed", "Inflamed"]


def read_json(path: Path):
    with path.open() as fh:
        return json.load(fh)


def write_tsv(df: pd.DataFrame, path: Path, index: bool = False) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=index)
    return path


def load_senmayo_symbols() -> list[str]:
    sets = read_json(P.external / "genesets" / "senescence_sets.json")
    return [g.upper() for g in sets["SenMayo"]]


def normalize_copy(adata: ad.AnnData, target_sum: float = 1e4) -> ad.AnnData:
    x = adata.copy()
    if "counts" in x.layers:
        x.X = x.layers["counts"].copy()
    sc.pp.normalize_total(x, target_sum=target_sum)
    sc.pp.log1p(x)
    return x


def symbol_index(adata: ad.AnnData) -> pd.Index:
    if "gene_symbol" in adata.var:
        vals = adata.var["gene_symbol"].astype(str).str.upper()
    else:
        vals = pd.Index(adata.var_names.astype(str)).str.upper()
    return pd.Index(vals)


def present_genes(adata: ad.AnnData, genes: Iterable[str]) -> list[str]:
    present = set(symbol_index(adata))
    return [g for g in (x.upper() for x in genes) if g in present]


def add_gene_scores(
    adata: ad.AnnData,
    genes: list[str],
    score_name: str,
    ctrl_size: int = 50,
) -> int:
    genes = present_genes(adata, genes)
    if not genes:
        adata.obs[score_name] = np.nan
        return 0
    old = adata.var_names.copy()
    symbols = pd.Series(symbol_index(adata).astype(str))
    dup_n = symbols.groupby(symbols, sort=False).cumcount()
    tmp_names = symbols.where(dup_n.eq(0), symbols + "__dup" + dup_n.astype(str))
    adata.var_names = pd.Index(tmp_names)
    genes = [g for g in genes if g in adata.var_names]
    sc.tl.score_genes(
        adata,
        gene_list=genes,
        score_name=score_name,
        ctrl_size=min(ctrl_size, max(1, adata.n_vars - len(genes))),
        random_state=20260622,
        use_raw=False,
    )
    adata.var_names = old
    return len(genes)


def sparse_col_frame(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    genes = [g.upper() for g in genes]
    sym = symbol_index(adata)
    cols = {}
    for gene in genes:
        loc = np.flatnonzero(sym == gene)
        if len(loc) == 0:
            cols[gene] = np.zeros(adata.n_obs, dtype=np.float32)
            continue
        x = adata.X[:, loc[0]]
        if sp.issparse(x):
            x = x.toarray().ravel()
        else:
            x = np.asarray(x).ravel()
        cols[gene] = x.astype(np.float32, copy=False)
    return pd.DataFrame(cols, index=adata.obs_names)


def group_expression(
    adata: ad.AnnData,
    genes: list[str],
    by: list[str],
    min_expr: float = 0.0,
) -> pd.DataFrame:
    expr = sparse_col_frame(adata, genes)
    meta = adata.obs[by].reset_index(drop=True).copy()
    expr = expr.reset_index(drop=True)
    rows = []
    for keys, idx in meta.groupby(by, observed=True).indices.items():
        if not isinstance(keys, tuple):
            keys = (keys,)
        sub = expr.iloc[idx]
        means = sub.mean(axis=0)
        fracs = (sub > min_expr).mean(axis=0)
        for gene in genes:
            row = dict(zip(by, keys, strict=True))
            row.update(
                gene=gene,
                mean_expr=float(means[gene]),
                frac_expr=float(fracs[gene]),
                n_cells=int(len(idx)),
            )
            rows.append(row)
    return pd.DataFrame(rows)


def stratified_sample_obs(
    obs: pd.DataFrame,
    n: int | None,
    strata: list[str],
    seed: int = 20260622,
) -> np.ndarray:
    if not n or n <= 0 or n >= len(obs):
        return np.arange(len(obs))
    rng = np.random.default_rng(seed)
    keep: list[int] = []
    grouped = obs.reset_index(drop=True).groupby(strata, dropna=False, observed=True).indices
    per = max(1, int(np.ceil(n / max(1, len(grouped)))))
    for idx in grouped.values():
        idx = np.asarray(idx)
        k = min(len(idx), per)
        keep.extend(rng.choice(idx, size=k, replace=False).tolist())
    if len(keep) > n:
        keep = rng.choice(np.asarray(keep), size=n, replace=False).tolist()
    return np.asarray(sorted(keep), dtype=int)


def geneformer_paths() -> dict[str, Path]:
    base = P.root / "models" / "Geneformer"
    dictionaries = base / "geneformer" / "gene_dictionaries_30m"
    return {
        "model": base / "Geneformer-V1-10M",
        "token_dict": dictionaries / "token_dictionary_gc30M.pkl",
        "median_dict": dictionaries / "gene_median_dictionary_gc30M.pkl",
        "name_to_id": dictionaries / "gene_name_id_dict_gc30M.pkl",
    }


def load_pickle(path: Path):
    with path.open("rb") as fh:
        return pickle.load(fh)


def append_status(text: str) -> None:
    path = P.journal / "status" / "task2_status.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(text.rstrip() + "\n\n")
