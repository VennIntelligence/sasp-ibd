"""Task 2 GPU hardening: null-calibrated Geneformer perturbation statistics.

Run one independent model x dataset job:
    python src/sc_05_perturb_stats.py --job --model-key gf_v1_10m --dataset smillie --device cuda:0

Run the planned two-GPU schedule and combine outputs:
    python src/sc_05_perturb_stats.py --schedule
"""
from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from pydantic import BaseModel, ConfigDict, Field
from sklearn.metrics import accuracy_score, average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from torch import nn
from torch.utils.data import DataLoader, TensorDataset
from tqdm.auto import tqdm
from transformers import AutoModel

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paths import P
from sc_task2_utils import append_status, load_pickle, stratified_sample_obs, symbol_index, write_tsv


OUT = P.out("sc_05_perturb_stats")
JOB_DIR = OUT / "jobs"
LOG_DIR = OUT / "logs"
SC01 = P.out("sc_01")

TARGET_GENES = ["CCL8", "CXCR2"]
CONTROL_GENES = ["CXCL8", "MMP3", "MMP9", "CXCL10"]
PANEL_GENES = TARGET_GENES + CONTROL_GENES
PERTURBATIONS = ["delete", "overexpress"]
EXPECTED_DIRECTION = {
    ("CCL8", "delete"): -1,
    ("CCL8", "overexpress"): 1,
    ("CXCR2", "delete"): 1,
    ("CXCR2", "overexpress"): -1,
}


class RunConfig(BaseModel):
    model_config = ConfigDict(frozen=True)

    model_key: str
    dataset: str
    max_cells_per_dataset: int = Field(default=1000, ge=100)
    batch_size: int = Field(default=128, ge=1)
    null_genes: int = Field(default=100, ge=10)
    n_bootstrap: int = Field(default=200, ge=0)
    seed: int = 20260622
    device: str = "cuda"


@dataclass(frozen=True)
class ModelSpec:
    key: str
    label: str
    model_path: Path
    token_dict_path: Path
    median_dict_path: Path
    name_to_id_path: Path
    max_len: int
    special_token: bool
    emb_mode: str


class LinearHead(nn.Module):
    def __init__(self, n_features: int) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_features, 64),
            nn.ReLU(),
            nn.Dropout(0.15),
            nn.Linear(64, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(1)


@dataclass
class ClassifierBundle:
    model: LinearHead
    mean: np.ndarray
    std: np.ndarray
    device: str
    auc: float
    auprc: float
    accuracy: float
    train_n: int
    val_n: int


def log(msg: str) -> None:
    print(f"[sc_05] {msg}", flush=True)


def stable_int(text: str) -> int:
    return int(hashlib.md5(text.encode("utf-8")).hexdigest()[:8], 16)


def model_specs() -> dict[str, ModelSpec]:
    base = P.root / "models" / "Geneformer"
    gf = base / "geneformer"
    return {
        "gf_v1_10m": ModelSpec(
            key="gf_v1_10m",
            label="Geneformer V1-10M",
            model_path=base / "Geneformer-V1-10M",
            token_dict_path=gf / "gene_dictionaries_30m" / "token_dictionary_gc30M.pkl",
            median_dict_path=gf / "gene_dictionaries_30m" / "gene_median_dictionary_gc30M.pkl",
            name_to_id_path=gf / "gene_dictionaries_30m" / "gene_name_id_dict_gc30M.pkl",
            max_len=2048,
            special_token=False,
            emb_mode="cell",
        ),
        "gf_v2_104m": ModelSpec(
            key="gf_v2_104m",
            label="Geneformer V2-104M",
            model_path=base / "Geneformer-V2-104M",
            token_dict_path=gf / "token_dictionary_gc104M.pkl",
            median_dict_path=gf / "gene_median_dictionary_gc104M.pkl",
            name_to_id_path=gf / "gene_name_id_dict_gc104M.pkl",
            max_len=4096,
            special_token=True,
            emb_mode="cls",
        ),
    }


def dataset_path(key: str) -> Path:
    paths = {
        "smillie": SC01 / "smillie_uc_qc.h5ad",
        "martin": SC01 / "martin_cd_qc.h5ad",
    }
    if key not in paths:
        raise KeyError(f"Unknown dataset {key}. Valid: {sorted(paths)}")
    if not paths[key].exists():
        raise FileNotFoundError(paths[key])
    return paths[key]


def load_assets(spec: ModelSpec) -> dict:
    for path in [spec.model_path, spec.token_dict_path, spec.median_dict_path, spec.name_to_id_path]:
        if not path.exists():
            raise FileNotFoundError(path)
    token_dict = load_pickle(spec.token_dict_path)
    return {
        "token_dict": token_dict,
        "median_dict": load_pickle(spec.median_dict_path),
        "name_to_id": {str(k).upper(): str(v) for k, v in load_pickle(spec.name_to_id_path).items()},
        "pad_id": int(token_dict.get("<pad>", 0)),
        "cls_id": token_dict.get("<cls>"),
        "eos_id": token_dict.get("<eos>"),
    }


def gene_ids_for_adata(adata: ad.AnnData, name_to_id: dict[str, str]) -> np.ndarray:
    if "ensembl_id" in adata.var:
        raw = adata.var["ensembl_id"].astype(str).str.upper().to_numpy()
        symbols = symbol_index(adata).astype(str).str.upper().to_numpy()
        return np.array([x if x.startswith("ENSG") else name_to_id.get(s, "") for x, s in zip(raw, symbols, strict=True)])
    return np.array([name_to_id.get(s, "") for s in symbol_index(adata).astype(str).str.upper()])


def counts_layer(adata: ad.AnnData) -> sp.csr_matrix:
    x = adata.layers["counts"] if "counts" in adata.layers else adata.X
    return x.tocsr() if sp.issparse(x) else sp.csr_matrix(x)


def variable_table(
    adata: ad.AnnData,
    rows: np.ndarray,
    token_dict: dict[str, int],
    median_dict: dict[str, float],
    name_to_id: dict[str, str],
) -> pd.DataFrame:
    x = counts_layer(adata)[rows]
    csc = x.tocsc()
    mean_counts = np.asarray(csc.mean(axis=0)).ravel().astype(np.float64)
    detected = np.diff(csc.indptr).astype(np.float64)
    ensg = gene_ids_for_adata(adata, name_to_id)
    symbols = symbol_index(adata).astype(str).str.upper().to_numpy()
    tokens = np.array([int(token_dict.get(g, 0)) if g else 0 for g in ensg], dtype=np.int32)
    med = np.array([float(median_dict.get(g, np.nan)) if g else np.nan for g in ensg], dtype=np.float64)
    out = pd.DataFrame(
        {
            "var_index": np.arange(adata.n_vars, dtype=np.int32),
            "gene": symbols,
            "ensembl_id": ensg,
            "token": tokens,
            "median": med,
            "mean_counts": mean_counts,
            "detected_fraction": detected / max(1, x.shape[0]),
            "detected_cells": detected.astype(int),
        }
    )
    out = out[(out["token"] > 0) & np.isfinite(out["median"]) & (out["median"] > 0)]
    out = out.sort_values(["gene", "detected_cells", "mean_counts"], ascending=[True, False, False])
    return out.drop_duplicates("gene", keep="first").reset_index(drop=True)


def augment_rows_for_target_detection(
    adata: ad.AnnData,
    rows: np.ndarray,
    genes: list[str],
    seed: int,
    max_per_gene: int = 500,
) -> np.ndarray:
    x = counts_layer(adata)
    sym = symbol_index(adata).astype(str).str.upper().to_numpy()
    keep = set(int(i) for i in rows)
    rng = np.random.default_rng(seed)
    for gene in genes:
        loc = np.flatnonzero(sym == gene)
        if loc.size == 0:
            continue
        detected = x[:, int(loc[0])].nonzero()[0]
        if detected.size == 0:
            continue
        take_n = min(max_per_gene, detected.size)
        take = rng.choice(detected, size=take_n, replace=False)
        keep.update(int(i) for i in take)
    return np.asarray(sorted(keep), dtype=int)


def build_tokens(
    adata: ad.AnnData,
    rows: np.ndarray,
    cfg: RunConfig,
    spec: ModelSpec,
    assets: dict,
) -> tuple[list[list[int]], pd.DataFrame, pd.DataFrame]:
    x = counts_layer(adata)
    ensg = gene_ids_for_adata(adata, assets["name_to_id"])
    tokens_by_var = np.array([int(assets["token_dict"].get(g, 0)) if g else 0 for g in ensg], dtype=np.int32)
    medians = np.array([float(assets["median_dict"].get(g, np.nan)) if g else np.nan for g in ensg], dtype=np.float32)
    ok = (tokens_by_var > 0) & np.isfinite(medians) & (medians > 0)

    tokenized: list[list[int]] = []
    kept_rows: list[int] = []
    for row in tqdm(rows, desc=f"{cfg.dataset} tokenize {spec.key}", unit="cell"):
        r = x.getrow(int(row))
        total = float(r.sum())
        if total <= 0:
            continue
        mask = ok[r.indices]
        idx = r.indices[mask]
        if idx.size == 0:
            continue
        scaled = (r.data[mask].astype(np.float32) / total * 10000.0) / medians[idx]
        order = np.argsort(-scaled, kind="stable")
        seq = tokens_by_var[idx][order].astype(int).tolist()
        if spec.special_token:
            seq = [int(assets["cls_id"])] + seq[: spec.max_len - 2] + [int(assets["eos_id"])]
        else:
            seq = seq[: spec.max_len]
        tokenized.append(seq)
        kept_rows.append(int(row))

    meta = adata.obs.iloc[kept_rows].copy().reset_index(names="cell")
    meta["row_index"] = kept_rows
    meta["tokenized_length"] = [len(x) for x in tokenized]
    genes = variable_table(adata, np.asarray(kept_rows, dtype=int), assets["token_dict"], assets["median_dict"], assets["name_to_id"])
    return tokenized, meta, genes


def select_null_panel(genes: pd.DataFrame, cfg: RunConfig) -> pd.DataFrame:
    candidates = genes[
        (genes["detected_cells"] >= 3)
        & (~genes["gene"].isin(PANEL_GENES))
        & (~genes["gene"].str.startswith("MT-"))
    ].copy()
    if candidates.empty:
        raise ValueError("No candidate genes available for null matching")

    feature = candidates[["mean_counts", "detected_fraction"]].copy()
    feature["log_mean"] = np.log1p(feature["mean_counts"])
    scale = feature[["log_mean", "detected_fraction"]].std().replace(0, 1).to_numpy()
    rows = []
    for gene in PANEL_GENES:
        hit = genes[genes["gene"].eq(gene)]
        if hit.empty:
            log(f"skip {gene}: absent from tokenized gene table")
            continue
        row = hit.iloc[0].to_dict()
        row["gene_role"] = "target" if gene in TARGET_GENES else "control"
        row["null_parent"] = gene if gene in TARGET_GENES else ""
        row["null_rank"] = 0
        rows.append(row)

    for parent in TARGET_GENES:
        target = genes[genes["gene"].eq(parent)]
        if target.empty:
            continue
        t = target.iloc[0]
        center = np.array([np.log1p(float(t["mean_counts"])), float(t["detected_fraction"])])
        dist = np.sqrt((((feature[["log_mean", "detected_fraction"]].to_numpy() - center) / scale) ** 2).sum(axis=1))
        pool_n = min(len(candidates), max(cfg.null_genes, cfg.null_genes * 5))
        nearest = candidates.assign(match_distance=dist).nsmallest(pool_n, "match_distance").copy()
        rng = np.random.default_rng(cfg.seed + stable_int(f"{cfg.dataset}:{cfg.model_key}:{parent}"))
        if len(nearest) > cfg.null_genes:
            weights = 1.0 / (nearest["match_distance"].to_numpy() + 1e-3)
            weights = weights / weights.sum()
            take = rng.choice(np.arange(len(nearest)), size=cfg.null_genes, replace=False, p=weights)
            nearest = nearest.iloc[np.sort(take)].copy()
        nearest["gene_role"] = "null"
        nearest["null_parent"] = parent
        nearest["null_rank"] = np.arange(1, len(nearest) + 1)
        rows.extend(nearest.to_dict("records"))

    panel = pd.DataFrame(rows).drop_duplicates(["gene", "gene_role", "null_parent"])
    return panel.sort_values(["gene_role", "null_parent", "null_rank", "gene"]).reset_index(drop=True)


def pad_batch(batch: list[list[int]], pad_id: int, device: str) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    max_len = max(len(x) for x in batch)
    arr = np.full((len(batch), max_len), pad_id, dtype=np.int64)
    mask = np.zeros((len(batch), max_len), dtype=np.int64)
    lengths = np.zeros(len(batch), dtype=np.int64)
    for i, seq in enumerate(batch):
        arr[i, : len(seq)] = seq
        mask[i, : len(seq)] = 1
        lengths[i] = len(seq)
    return torch.from_numpy(arr).to(device), torch.from_numpy(mask).to(device), torch.from_numpy(lengths).to(device)


@torch.inference_mode()
def embed_tokens(
    model: AutoModel,
    tokens: list[list[int]],
    cfg: RunConfig,
    spec: ModelSpec,
    pad_id: int,
    desc: str,
) -> np.ndarray:
    out = []
    use_amp = cfg.device.startswith("cuda")
    for start in tqdm(range(0, len(tokens), cfg.batch_size), desc=desc, unit="batch"):
        batch = tokens[start : start + cfg.batch_size]
        ids, mask, lengths = pad_batch(batch, pad_id, cfg.device)
        with torch.autocast(device_type="cuda", dtype=torch.float16, enabled=use_amp):
            pred = model(input_ids=ids, attention_mask=mask)
        hidden = pred.last_hidden_state.float()
        if spec.emb_mode == "cls":
            pooled = hidden[:, 0, :]
        else:
            if spec.special_token:
                use_hidden = hidden[:, 1:, :]
                use_lengths = torch.clamp(lengths - 2, min=1)
                valid = torch.arange(use_hidden.shape[1], device=cfg.device).view(1, -1) < use_lengths.view(-1, 1)
                pooled = (use_hidden * valid.unsqueeze(-1)).sum(1) / use_lengths.view(-1, 1)
            else:
                denom = mask.sum(1).clamp(min=1).view(-1, 1)
                pooled = (hidden * mask.unsqueeze(-1)).sum(1) / denom
        out.append(pooled.cpu().numpy().astype(np.float32))
    return np.vstack(out)


def perturb_one(seq: list[int], token: int, mode: str, spec: ModelSpec, assets: dict) -> list[int]:
    if mode == "delete":
        if spec.special_token:
            body = [x for x in seq[1:-1] if x != token]
            return [seq[0]] + body[: spec.max_len - 2] + [seq[-1]]
        return [x for x in seq if x != token][: spec.max_len]

    if mode != "overexpress":
        raise ValueError(mode)

    if spec.special_token:
        body = [x for x in seq[1:-1] if x != token]
        out = [int(assets["cls_id"]), token] + body
        if len(out) > spec.max_len - 1:
            out = out[: spec.max_len - 1]
        return out + [int(assets["eos_id"])]
    body = [x for x in seq if x != token]
    return ([token] + body)[: spec.max_len]


def projection(delta: np.ndarray, direction: np.ndarray) -> np.ndarray:
    denom = float(np.linalg.norm(direction))
    if denom == 0:
        return np.zeros(delta.shape[0], dtype=np.float32)
    return delta @ (direction / denom)


def cosine_to(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    den = np.linalg.norm(a, axis=1) * max(float(np.linalg.norm(b)), 1e-8)
    return (a @ b) / np.maximum(den, 1e-8)


def bootstrap_ci(values: np.ndarray, cfg: RunConfig) -> tuple[float, float]:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size <= 1 or cfg.n_bootstrap <= 0:
        val = float(values.mean()) if values.size else np.nan
        return val, val
    rng = np.random.default_rng(cfg.seed + values.size)
    idx = rng.integers(0, values.size, size=(cfg.n_bootstrap, values.size))
    means = values[idx].mean(axis=1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return float(lo), float(hi)


def summarize_effect(
    values: np.ndarray,
    probs: np.ndarray | None,
    compute_ci: bool,
    cfg: RunConfig,
) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    lo, hi = bootstrap_ci(values, cfg) if compute_ci else (np.nan, np.nan)
    out = {
        "n_cells": int(values.size),
        "projection_toward_inflamed_mean": float(np.mean(values)) if values.size else np.nan,
        "projection_toward_inflamed_sem": float(np.std(values, ddof=1) / np.sqrt(values.size)) if values.size > 1 else 0.0,
        "projection_toward_inflamed_ci_low": lo,
        "projection_toward_inflamed_ci_high": hi,
    }
    if probs is None:
        out.update(
            delta_p_inflamed_mean=np.nan,
            delta_p_inflamed_sem=np.nan,
            delta_p_inflamed_ci_low=np.nan,
            delta_p_inflamed_ci_high=np.nan,
        )
    else:
        probs = np.asarray(probs, dtype=np.float64)
        p_lo, p_hi = bootstrap_ci(probs, cfg) if compute_ci else (np.nan, np.nan)
        out.update(
            delta_p_inflamed_mean=float(np.mean(probs)) if probs.size else np.nan,
            delta_p_inflamed_sem=float(np.std(probs, ddof=1) / np.sqrt(probs.size)) if probs.size > 1 else 0.0,
            delta_p_inflamed_ci_low=p_lo,
            delta_p_inflamed_ci_high=p_hi,
        )
    return out


def completed_effect_keys(path: Path) -> set[tuple[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return set()
    df = pd.read_csv(path, sep="\t", usecols=["gene", "perturbation"])
    return set(zip(df["gene"].astype(str), df["perturbation"].astype(str), strict=True))


def upsert_tsv(path: Path, rows: pd.DataFrame, key_cols: list[str]) -> None:
    if rows.empty:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = rows.copy()
    if path.exists() and path.stat().st_size > 0:
        old = pd.read_csv(path, sep="\t")
        key_new = pd.MultiIndex.from_frame(rows[key_cols].astype(str))
        key_old = pd.MultiIndex.from_frame(old[key_cols].astype(str))
        old = old.loc[~key_old.isin(key_new)]
        rows = pd.concat([old, rows], ignore_index=True)
    write_tsv(rows, path)


def split_by_subject(meta: pd.DataFrame, y: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray]:
    groups = meta["Subject"].astype(str).to_numpy() if "Subject" in meta.columns else meta["cell"].astype(str).to_numpy()
    n_splits = min(5, len(np.unique(groups)))
    if n_splits < 2:
        idx = np.arange(len(y))
        rng = np.random.default_rng(seed)
        rng.shuffle(idx)
        cut = max(1, int(len(idx) * 0.8))
        return idx[:cut], idx[cut:]
    splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for train, val in splitter.split(np.zeros(len(y)), y, groups):
        if len(np.unique(y[train])) == 2 and len(np.unique(y[val])) == 2:
            return train, val
    raise ValueError("No valid subject-level split with both classes in train and validation")


def train_classifier(base: np.ndarray, meta: pd.DataFrame, cfg: RunConfig) -> ClassifierBundle | None:
    y = meta["Health"].eq("Inflamed").astype(np.float32).to_numpy()
    if len(np.unique(y)) < 2:
        log("classifier skipped: only one class")
        return None
    train, val = split_by_subject(meta, y.astype(int), cfg.seed)
    x_train = base[train].astype(np.float32)
    x_val = base[val].astype(np.float32)
    mean = x_train.mean(axis=0, keepdims=True)
    std = x_train.std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    x_train = (x_train - mean) / std
    x_val = (x_val - mean) / std
    y_train = y[train]
    y_val = y[val]

    device = cfg.device if cfg.device.startswith("cuda") else "cpu"
    model = LinearHead(x_train.shape[1]).to(device)
    pos_weight = torch.tensor([(len(y_train) - y_train.sum()) / max(float(y_train.sum()), 1.0)], device=device)
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=1e-4)
    ds = TensorDataset(torch.from_numpy(x_train), torch.from_numpy(y_train))
    loader = DataLoader(ds, batch_size=min(512, len(ds)), shuffle=True, drop_last=False)

    best_auc = -np.inf
    best_state = None
    no_improve = 0
    for epoch in tqdm(range(120), desc="classifier head", unit="epoch"):
        model.train()
        for xb, yb in loader:
            xb = xb.to(device)
            yb = yb.to(device)
            opt.zero_grad(set_to_none=True)
            loss = loss_fn(model(xb), yb)
            loss.backward()
            opt.step()
        if epoch % 5:
            continue
        model.eval()
        with torch.no_grad():
            pred = torch.sigmoid(model(torch.from_numpy(x_val).to(device))).cpu().numpy()
        try:
            auc = roc_auc_score(y_val, pred)
        except ValueError:
            auc = np.nan
        score = -np.inf if not np.isfinite(auc) else auc
        if score > best_auc:
            best_auc = score
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            no_improve = 0
        else:
            no_improve += 1
        if no_improve >= 12:
            break
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    pred = predict_prob(model, x_val * std + mean, mean, std, device)
    auc = float(roc_auc_score(y_val, pred))
    auprc = float(average_precision_score(y_val, pred))
    acc = float(accuracy_score(y_val, pred >= 0.5))
    return ClassifierBundle(
        model=model,
        mean=mean.astype(np.float32),
        std=std.astype(np.float32),
        device=device,
        auc=auc,
        auprc=auprc,
        accuracy=acc,
        train_n=int(len(train)),
        val_n=int(len(val)),
    )


@torch.inference_mode()
def predict_prob(model: LinearHead, emb: np.ndarray, mean: np.ndarray, std: np.ndarray, device: str) -> np.ndarray:
    x = ((emb.astype(np.float32) - mean) / std).astype(np.float32)
    out = []
    for start in range(0, len(x), 4096):
        xb = torch.from_numpy(x[start : start + 4096]).to(device)
        out.append(torch.sigmoid(model(xb)).cpu().numpy())
    return np.concatenate(out).reshape(-1)


def run_job(cfg: RunConfig) -> None:
    spec = model_specs()[cfg.model_key]
    assets = load_assets(spec)
    path = dataset_path(cfg.dataset)
    log(f"job start model={cfg.model_key} dataset={cfg.dataset} device={cfg.device}")
    adata = ad.read_h5ad(path)
    rows = stratified_sample_obs(
        adata.obs,
        cfg.max_cells_per_dataset,
        ["Health", "cell_type"] if "cell_type" in adata.obs else ["Health"],
        seed=cfg.seed + stable_int(f"{cfg.model_key}:{cfg.dataset}"),
    )
    rows = augment_rows_for_target_detection(
        adata,
        rows,
        TARGET_GENES,
        seed=cfg.seed + stable_int(f"target-detection:{cfg.model_key}:{cfg.dataset}"),
    )
    tokens, meta, genes = build_tokens(adata, rows, cfg, spec, assets)
    if not tokens:
        raise ValueError(f"{path}: no cells tokenized")
    panel = select_null_panel(genes, cfg)
    write_tsv(panel, JOB_DIR / f"{cfg.model_key}_{cfg.dataset}_null_panel.tsv")

    log(f"tokenized {len(tokens):,} cells; panel genes={len(panel)}")
    model = AutoModel.from_pretrained(spec.model_path, output_hidden_states=False).to(cfg.device).eval()
    if cfg.device.startswith("cuda"):
        model = model.half()
    base = embed_tokens(model, tokens, cfg, spec, assets["pad_id"], desc=f"{cfg.model_key} {cfg.dataset} base")
    meta["embedding_norm"] = np.linalg.norm(base, axis=1)

    inflamed = meta["Health"].eq("Inflamed").to_numpy()
    quiet = meta["Health"].isin(["Healthy", "Non-inflamed"]).to_numpy()
    if inflamed.sum() == 0 or quiet.sum() == 0:
        raise ValueError("Need both inflamed and healthy/non-inflamed cells")
    direction = base[inflamed].mean(axis=0) - base[quiet].mean(axis=0)
    inflamed_centroid = base[inflamed].mean(axis=0)
    quiet_centroid = base[quiet].mean(axis=0)

    clf = train_classifier(base, meta, cfg)
    clf_row = {
        "model": cfg.model_key,
        "model_label": spec.label,
        "dataset": str(meta["dataset"].iloc[0]) if "dataset" in meta.columns else cfg.dataset,
        "dataset_key": cfg.dataset,
        "train_n": 0,
        "val_n": 0,
        "auc": np.nan,
        "auprc": np.nan,
        "accuracy": np.nan,
    }
    if clf is not None:
        clf_row.update(train_n=clf.train_n, val_n=clf.val_n, auc=clf.auc, auprc=clf.auprc, accuracy=clf.accuracy)
        log(f"classifier AUC={clf.auc:.3f} AUPRC={clf.auprc:.3f}")

    effect_partial = JOB_DIR / f"{cfg.model_key}_{cfg.dataset}_effects.partial.tsv"
    cell_partial = JOB_DIR / f"{cfg.model_key}_{cfg.dataset}_cell_targets.partial.tsv"
    completed = completed_effect_keys(effect_partial)
    if completed:
        log(f"resume: skipping {len(completed)} completed gene perturbations from {effect_partial.relative_to(P.root)}")

    for item in tqdm(panel.to_dict("records"), desc=f"{cfg.model_key} {cfg.dataset} perturb genes", unit="gene"):
        gene = str(item["gene"])
        token = int(item["token"])
        role = str(item["gene_role"])
        null_parent = str(item.get("null_parent", ""))
        has_token = np.fromiter((token in seq for seq in tokens), dtype=bool, count=len(tokens))
        for mode in PERTURBATIONS:
            if (gene, mode) in completed:
                continue
            use = has_token if mode == "delete" else np.ones(len(tokens), dtype=bool)
            if use.sum() == 0:
                continue
            idx = np.flatnonzero(use)
            perturbed = [perturb_one(tokens[i], token, mode, spec, assets) for i in idx]
            pert = embed_tokens(
                model,
                perturbed,
                cfg,
                spec,
                assets["pad_id"],
                desc=f"{cfg.model_key} {cfg.dataset} {gene} {mode}",
            )
            delta = pert - base[idx]
            proj = projection(delta, direction)
            d_inflamed = cosine_to(pert, inflamed_centroid) - cosine_to(base[idx], inflamed_centroid)
            d_quiet = cosine_to(pert, quiet_centroid) - cosine_to(base[idx], quiet_centroid)
            delta_prob = None
            if clf is not None:
                delta_prob = predict_prob(clf.model, pert, clf.mean, clf.std, clf.device) - predict_prob(
                    clf.model, base[idx], clf.mean, clf.std, clf.device
                )
            row = {
                "model": cfg.model_key,
                "model_label": spec.label,
                "dataset": str(meta["dataset"].iloc[0]) if "dataset" in meta.columns else cfg.dataset,
                "dataset_key": cfg.dataset,
                "gene": gene,
                "perturbation": mode,
                "gene_role": role,
                "null_parent": null_parent,
                "null_rank": int(item.get("null_rank", 0)),
                "target_detected_fraction": float(has_token.mean()),
                "mean_counts": float(item["mean_counts"]),
                "detected_fraction": float(item["detected_fraction"]),
                "detected_cells": int(item["detected_cells"]),
                "delta_cosine_inflamed_mean": float(np.mean(d_inflamed)),
                "delta_cosine_quiescent_mean": float(np.mean(d_quiet)),
            }
            row.update(summarize_effect(proj, delta_prob, compute_ci=(role != "null"), cfg=cfg))
            if role != "null":
                sub = meta.iloc[idx][["cell", "dataset", "Health", "cell_type", "Subject"]].copy()
                sub["model"] = cfg.model_key
                sub["gene"] = gene
                sub["perturbation"] = mode
                sub["projection_toward_inflamed"] = proj
                sub["delta_cosine_inflamed"] = d_inflamed
                sub["delta_cosine_quiescent"] = d_quiet
                sub["delta_p_inflamed"] = delta_prob if delta_prob is not None else np.nan
                upsert_tsv(cell_partial, sub, ["model", "gene", "perturbation", "cell"])
            upsert_tsv(effect_partial, pd.DataFrame([row]), ["model", "dataset_key", "gene", "perturbation"])
            completed.add((gene, mode))

    if not effect_partial.exists():
        raise FileNotFoundError(f"No checkpointed effect rows were written for {cfg.model_key} {cfg.dataset}")
    effects = pd.read_csv(effect_partial, sep="\t")
    cells = pd.read_csv(cell_partial, sep="\t") if cell_partial.exists() else pd.DataFrame()
    write_tsv(effects, JOB_DIR / f"{cfg.model_key}_{cfg.dataset}_effects.tsv")
    write_tsv(cells, JOB_DIR / f"{cfg.model_key}_{cfg.dataset}_cell_targets.tsv")
    write_tsv(pd.DataFrame([clf_row]), JOB_DIR / f"{cfg.model_key}_{cfg.dataset}_classifier.tsv")
    log(f"job done {cfg.model_key} {cfg.dataset}: {effects.shape[0]} effect rows")


def add_null_statistics(effects: pd.DataFrame) -> pd.DataFrame:
    out = effects.copy()
    stat_cols = [
        "projection_null_n",
        "projection_null_mean",
        "projection_null_sd",
        "projection_null_z",
        "projection_empirical_p_abs0",
        "delta_p_null_n",
        "delta_p_null_mean",
        "delta_p_null_sd",
        "delta_p_null_z",
        "delta_p_empirical_p_abs0",
    ]
    for col in stat_cols:
        out[col] = np.nan
    target_mask = out["gene_role"].eq("target")
    for idx, row in out[target_mask].iterrows():
        nulls = out[
            out["gene_role"].eq("null")
            & out["model"].eq(row["model"])
            & out["dataset_key"].eq(row["dataset_key"])
            & out["perturbation"].eq(row["perturbation"])
            & out["null_parent"].eq(row["gene"])
        ]
        if nulls.empty:
            continue
        for prefix, value_col in [
            ("projection", "projection_toward_inflamed_mean"),
            ("delta_p", "delta_p_inflamed_mean"),
        ]:
            vals = nulls[value_col].dropna().to_numpy(dtype=float)
            if vals.size == 0 or not np.isfinite(row[value_col]):
                continue
            mean = float(vals.mean())
            sd = float(vals.std(ddof=1)) if vals.size > 1 else np.nan
            z = float((row[value_col] - mean) / sd) if sd and np.isfinite(sd) and sd > 0 else np.nan
            p_abs = float((np.sum(np.abs(vals) >= abs(float(row[value_col]))) + 1) / (vals.size + 1))
            out.loc[idx, f"{prefix}_null_n"] = vals.size
            out.loc[idx, f"{prefix}_null_mean"] = mean
            out.loc[idx, f"{prefix}_null_sd"] = sd
            out.loc[idx, f"{prefix}_null_z"] = z
            out.loc[idx, f"{prefix}_empirical_p_abs0"] = p_abs
    return out


def make_consensus(stats: pd.DataFrame) -> pd.DataFrame:
    rows = []
    target = stats[stats["gene_role"].eq("target")].copy()
    for (gene, mode), sub in target.groupby(["gene", "perturbation"], observed=True):
        expected = EXPECTED_DIRECTION.get((gene, mode), np.nan)
        row = {
            "gene": gene,
            "perturbation": mode,
            "expected_direction": expected,
            "n_model_dataset": int(len(sub)),
        }
        for prefix, value_col, p_col in [
            ("projection", "projection_toward_inflamed_mean", "projection_empirical_p_abs0"),
            ("delta_p", "delta_p_inflamed_mean", "delta_p_empirical_p_abs0"),
        ]:
            vals = sub[value_col].to_numpy(dtype=float)
            signs = np.sign(vals)
            match = signs == expected
            sig = sub[p_col].to_numpy(dtype=float) < 0.05
            row[f"{prefix}_mean_effect"] = float(np.nanmean(vals))
            row[f"{prefix}_direction_match_n"] = int(np.nansum(match))
            row[f"{prefix}_significant_n"] = int(np.nansum(sig))
            row[f"{prefix}_significant_match_n"] = int(np.nansum(match & sig))
            row[f"{prefix}_models_with_sig_match"] = int(sub.loc[match & sig, "model"].nunique())
            row[f"{prefix}_consistent_significant"] = bool(row[f"{prefix}_models_with_sig_match"] >= 2)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["gene", "perturbation"]).reset_index(drop=True)


def plot_stats(stats: pd.DataFrame) -> Path:
    target = stats[stats["gene_role"].eq("target")].copy()
    target["label"] = (
        target["dataset_key"].str.title()
        + " | "
        + target["model"].str.replace("gf_", "GF-", regex=False)
        + " | "
        + target["gene"]
        + " "
        + target["perturbation"]
    )
    target = target.sort_values(["gene", "perturbation", "model", "dataset_key"]).reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, max(5.5, 0.33 * len(target))), sharey=True)
    for ax, metric, pcol, title in [
        (axes[0], "projection_toward_inflamed_mean", "projection_empirical_p_abs0", "Embedding projection"),
        (axes[1], "delta_p_inflamed_mean", "delta_p_empirical_p_abs0", "Classifier probability"),
    ]:
        vals = target[metric].to_numpy(dtype=float)
        colors = np.where(vals >= 0, "#b4432c", "#1f6f8b")
        y = np.arange(len(target))
        ax.barh(y, vals, color=colors, alpha=0.9)
        ax.axvline(0, color="#333333", lw=0.8)
        ax.set_title(title)
        ax.set_xlabel("Move toward inflamed state" if metric.startswith("projection") else "Delta P(inflamed)")
        for yi, val, p in zip(y, vals, target[pcol].to_numpy(dtype=float), strict=True):
            if np.isfinite(p) and p < 0.05:
                ax.text(val, yi, " *", va="center", ha="left" if val >= 0 else "right", fontsize=9)
    axes[0].set_yticks(np.arange(len(target)))
    axes[0].set_yticklabels(target["label"], fontsize=8)
    fig.suptitle("Null-calibrated Geneformer in-silico perturbation", fontweight="bold")
    fig.tight_layout()
    out = OUT / "Fig_task2_perturb_stats.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def write_summary(stats: pd.DataFrame, consensus: pd.DataFrame, classifiers: pd.DataFrame) -> Path:
    target = stats[stats["gene_role"].eq("target")].copy()

    def md_tsv(df: pd.DataFrame) -> str:
        return "```tsv\n" + df.to_csv(sep="\t", index=False).strip() + "\n```"

    lines = [
        "# Task2 perturbation hardening summary",
        "",
        "Positive projection or Delta P means movement toward the inflamed state. Null statistics use expression/detection-matched genes sampled separately for CCL8 and CXCR2.",
        "",
        "## Classifier validation",
        md_tsv(classifiers[["model", "dataset_key", "train_n", "val_n", "auc", "auprc", "accuracy"]]),
        "",
        "## Target perturbations",
        md_tsv(
            target[
                [
                    "model",
                    "dataset_key",
                    "gene",
                    "perturbation",
                    "projection_toward_inflamed_mean",
                    "projection_null_z",
                    "projection_empirical_p_abs0",
                    "delta_p_inflamed_mean",
                    "delta_p_null_z",
                    "delta_p_empirical_p_abs0",
                ]
            ]
        ),
        "",
        "## Consensus",
        md_tsv(consensus),
        "",
        "## Interpretation",
    ]
    for _, row in consensus.iterrows():
        gene = row["gene"]
        mode = row["perturbation"]
        proj_ok = bool(row["projection_consistent_significant"])
        prob_ok = bool(row["delta_p_consistent_significant"])
        lines.append(
            f"- {gene} {mode}: projection_sig_match_models={row['projection_models_with_sig_match']}, "
            f"deltaP_sig_match_models={row['delta_p_models_with_sig_match']}; "
            f"multi-model significant support={proj_ok or prob_ok}."
        )
    out = OUT / "task2_perturb_stats_summary.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def combine_outputs() -> None:
    effect_files = sorted(JOB_DIR.glob("*_effects.tsv"))
    classifier_files = sorted(JOB_DIR.glob("*_classifier.tsv"))
    if not effect_files:
        raise FileNotFoundError(f"No job effects in {JOB_DIR}")
    effects = pd.concat((pd.read_csv(p, sep="\t") for p in effect_files), ignore_index=True)
    classifiers = pd.concat((pd.read_csv(p, sep="\t") for p in classifier_files), ignore_index=True)
    stats = add_null_statistics(effects)
    consensus = make_consensus(stats)
    fig = plot_stats(stats)
    summary = write_summary(stats, consensus, classifiers)
    write_tsv(stats, OUT / "perturb_stats.tsv")
    write_tsv(consensus, OUT / "perturb_consensus.tsv")
    write_tsv(classifiers, OUT / "classifier_metrics.tsv")
    P.promote_table(OUT / "perturb_stats.tsv")
    P.promote_table(OUT / "perturb_consensus.tsv")
    P.promote_table(OUT / "classifier_metrics.tsv")
    P.promote_figure(fig)
    append_status(
        "## sc_05 perturbation hardening\n"
        f"- combined {len(effect_files)} model x dataset jobs from `outputs/sc_05_perturb_stats/jobs/`.\n"
        "- outputs: `outputs/sc_05_perturb_stats/perturb_stats.tsv`, `perturb_consensus.tsv`, "
        "`classifier_metrics.tsv`, `Fig_task2_perturb_stats.png`.\n"
        f"- promoted final table/figure copies to `results/`; summary: `{summary.relative_to(P.root)}`."
    )
    log(f"combined outputs; figure={fig.relative_to(P.root)}")


def schedule(args: argparse.Namespace) -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    JOB_DIR.mkdir(parents=True, exist_ok=True)
    for old in JOB_DIR.glob("*.tsv"):
        old.unlink()
    waves = [
        [
            ("gf_v2_104m", "smillie", 0, args.batch_size_v2),
            ("gf_v2_104m", "martin", 1, args.batch_size_v2),
        ],
        [
            ("gf_v1_10m", "smillie", 0, args.batch_size),
            ("gf_v1_10m", "martin", 1, args.batch_size),
        ],
    ]
    started = time.time()
    failures = []
    for wave_i, combos in enumerate(waves, start=1):
        procs: list[tuple[str, subprocess.Popen, object]] = []
        for model_key, dataset, gpu, batch_size in combos:
            name = f"{model_key}_{dataset}"
            cmd = [
                sys.executable,
                str(P.src / "sc_05_perturb_stats.py"),
                "--job",
                "--model-key",
                model_key,
                "--dataset",
                dataset,
                "--device",
                "cuda:0",
                "--max-cells-per-dataset",
                str(args.max_cells_per_dataset),
                "--batch-size",
                str(batch_size),
                "--null-genes",
                str(args.null_genes),
                "--n-bootstrap",
                str(args.n_bootstrap),
            ]
            env = os.environ.copy()
            env["CUDA_VISIBLE_DEVICES"] = str(gpu)
            env["PYTHONPATH"] = str(P.src)
            env.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
            log_path = LOG_DIR / f"{name}.log"
            fh = log_path.open("w", encoding="utf-8")
            log(
                f"launch wave{wave_i} {name} on physical GPU{gpu} batch={batch_size}; "
                f"log={log_path.relative_to(P.root)}"
            )
            procs.append((name, subprocess.Popen(cmd, cwd=P.root, env=env, stdout=fh, stderr=subprocess.STDOUT), fh))

        while procs:
            still = []
            for name, proc, fh in procs:
                rc = proc.poll()
                if rc is None:
                    still.append((name, proc, fh))
                else:
                    fh.close()
                    if rc:
                        failures.append((name, rc))
                        log(f"FAILED {name} rc={rc}")
                    else:
                        log(f"finished {name}")
            procs = still
            if procs:
                time.sleep(30)
        if failures:
            break
    elapsed = (time.time() - started) / 60.0
    if failures:
        raise RuntimeError(f"Job failures: {failures}")
    append_status(
        "## sc_05 two-GPU schedule\n"
        f"- launched Geneformer V2-104M as one job per GPU after measured 24 GB memory pressure, "
        f"then Geneformer V1-10M as one job per GPU; all jobs used `CUDA_VISIBLE_DEVICES`.\n"
        f"- max_cells_per_dataset={args.max_cells_per_dataset}, null_genes={args.null_genes}, "
        f"v1_batch_size={args.batch_size}, v2_batch_size={args.batch_size_v2}; elapsed={elapsed:.1f} minutes."
    )
    combine_outputs()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--job", action="store_true")
    mode.add_argument("--combine", action="store_true")
    mode.add_argument("--schedule", action="store_true")
    p.add_argument("--model-key", choices=sorted(model_specs()), default="gf_v1_10m")
    p.add_argument("--dataset", choices=["smillie", "martin"], default="smillie")
    p.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    p.add_argument("--max-cells-per-dataset", type=int, default=1000)
    p.add_argument("--batch-size", type=int, default=128)
    p.add_argument("--batch-size-v2", type=int, default=48)
    p.add_argument("--null-genes", type=int, default=100)
    p.add_argument("--n-bootstrap", type=int, default=200)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.job:
        run_job(
            RunConfig(
                model_key=args.model_key,
                dataset=args.dataset,
                max_cells_per_dataset=args.max_cells_per_dataset,
                batch_size=args.batch_size,
                null_genes=args.null_genes,
                n_bootstrap=args.n_bootstrap,
                device=args.device,
            )
        )
    elif args.combine:
        combine_outputs()
    elif args.schedule:
        schedule(args)


if __name__ == "__main__":
    main()
