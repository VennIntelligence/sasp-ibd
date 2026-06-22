"""Task 2 Geneformer zero-shot embeddings and in-silico perturbation."""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from pydantic import BaseModel, ConfigDict, Field
from transformers import AutoModel
from tqdm.auto import tqdm

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paths import P
from sc_task2_utils import (
    append_status,
    geneformer_paths,
    load_pickle,
    stratified_sample_obs,
    symbol_index,
    write_tsv,
)


OUT = P.out("sc_03")
SC01 = P.out("sc_01")
TARGET_GENES = ["CCL8", "CXCR2"]
CONTROL_GENES = ["CXCL8", "MMP3", "MMP9", "CXCL10"]


class PerturbConfig(BaseModel):
    model_config = ConfigDict(frozen=True)

    max_cells_per_dataset: int = Field(default=12000, ge=0)
    batch_size: int = Field(default=48, ge=1)
    max_len: int = Field(default=2048, ge=128)
    device: str = "cuda"


def log(msg: str) -> None:
    print(f"[sc_03] {msg}", flush=True)


def input_paths(use_subsample: bool) -> list[Path]:
    names = (
        ("smillie_uc_subsample.h5ad", "martin_cd_subsample.h5ad")
        if use_subsample
        else ("smillie_uc_qc.h5ad", "martin_cd_qc.h5ad")
    )
    return [SC01 / x for x in names if (SC01 / x).exists()]


def gene_id_map(adata: ad.AnnData, name_to_id: dict[str, str]) -> list[str | None]:
    if "ensembl_id" in adata.var:
        raw = adata.var["ensembl_id"].astype(str).tolist()
        return [x if x.startswith("ENSG") else name_to_id.get(x.upper()) for x in raw]
    ids = []
    for sym in symbol_index(adata):
        ids.append(name_to_id.get(str(sym).upper()))
    return ids


def build_tokens(
    adata: ad.AnnData,
    rows: np.ndarray,
    cfg: PerturbConfig,
    token_dict: dict[str, int],
    median_dict: dict[str, float],
    name_to_id: dict[str, str],
) -> tuple[list[list[int]], pd.DataFrame, dict[str, int]]:
    x = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if not sp.issparse(x):
        x = sp.csr_matrix(x)
    else:
        x = x.tocsr()
    ensg = gene_id_map(adata, name_to_id)
    tokens = np.array([token_dict.get(g, 0) if g else 0 for g in ensg], dtype=np.int32)
    medians = np.array([median_dict.get(g, np.nan) if g else np.nan for g in ensg], dtype=np.float32)
    ok = (tokens > 0) & np.isfinite(medians) & (medians > 0)
    target_ensg = {g: name_to_id.get(g, g) for g in TARGET_GENES + CONTROL_GENES}
    target_tokens = {g: token_dict[e] for g, e in target_ensg.items() if e in token_dict}
    out: list[list[int]] = []
    kept_obs = []
    for row in rows:
        r = x.getrow(row)
        n_counts = float(r.sum())
        if n_counts <= 0:
            continue
        mask = ok[r.indices]
        idx = r.indices[mask]
        if len(idx) == 0:
            continue
        vals = (r.data[mask].astype(np.float32) / n_counts * 10000.0) / medians[idx]
        order = np.argsort(-vals, kind="stable")
        out.append(tokens[idx][order][: cfg.max_len].astype(int).tolist())
        kept_obs.append(row)
    meta = adata.obs.iloc[kept_obs].copy()
    meta["tokenized_length"] = [len(t) for t in out]
    return out, meta, target_tokens


def pad_batch(batch: list[list[int]], pad_id: int, device: str) -> tuple[torch.Tensor, torch.Tensor]:
    max_len = max(len(x) for x in batch)
    arr = np.full((len(batch), max_len), pad_id, dtype=np.int64)
    mask = np.zeros((len(batch), max_len), dtype=np.int64)
    for i, toks in enumerate(batch):
        arr[i, : len(toks)] = toks
        mask[i, : len(toks)] = 1
    return torch.from_numpy(arr).to(device), torch.from_numpy(mask).to(device)


@torch.inference_mode()
def embed_tokens(
    model,
    tokens: list[list[int]],
    cfg: PerturbConfig,
    pad_id: int,
    desc: str = "embedding",
) -> np.ndarray:
    embs = []
    use_amp = cfg.device.startswith("cuda")
    for start in tqdm(range(0, len(tokens), cfg.batch_size), desc=desc, unit="batch"):
        batch = tokens[start : start + cfg.batch_size]
        ids, mask = pad_batch(batch, pad_id, cfg.device)
        with torch.autocast(device_type="cuda", dtype=torch.float16, enabled=use_amp):
            out = model(input_ids=ids, attention_mask=mask)
        hidden = out.last_hidden_state.float()
        denom = mask.sum(1).clamp(min=1).view(-1, 1)
        pooled = (hidden * mask.unsqueeze(-1)).sum(1) / denom
        embs.append(pooled.cpu().numpy().astype(np.float32))
    return np.vstack(embs)


def perturb_one(tokens: list[int], token: int, mode: str, max_len: int) -> list[int]:
    if mode == "delete":
        return [t for t in tokens if t != token]
    if mode == "overexpress":
        rest = [t for t in tokens if t != token]
        return ([token] + rest)[:max_len]
    raise ValueError(mode)


def projection(delta: np.ndarray, direction: np.ndarray) -> np.ndarray:
    denom = np.linalg.norm(direction)
    if denom == 0:
        return np.zeros(delta.shape[0], dtype=np.float32)
    return delta @ (direction / denom)


def cosine_to(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    den = np.linalg.norm(a, axis=1) * max(np.linalg.norm(b), 1e-8)
    return (a @ b) / np.maximum(den, 1e-8)


def run_dataset(path: Path, cfg: PerturbConfig, model, dictionaries: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    log(f"load {path.name}")
    adata = ad.read_h5ad(path)
    rows = stratified_sample_obs(
        adata.obs,
        cfg.max_cells_per_dataset,
        ["Health", "cell_type"] if "cell_type" in adata.obs else ["Health"],
    )
    tokens, meta, target_tokens = build_tokens(
        adata,
        rows,
        cfg,
        dictionaries["token_dict"],
        dictionaries["median_dict"],
        dictionaries["name_to_id"],
    )
    if not tokens:
        raise ValueError(f"{path}: no cells tokenized")
    log(f"{path.stem}: tokenized {len(tokens):,} cells; target tokens={target_tokens}")
    base = embed_tokens(model, tokens, cfg, pad_id=dictionaries["pad_id"], desc=f"{path.stem} base")
    meta = meta.reset_index(names="cell")
    meta["emb_norm"] = np.linalg.norm(base, axis=1)
    inflamed = meta["Health"].eq("Inflamed").values
    quiescent = meta["Health"].isin(["Healthy", "Non-inflamed"]).values
    if inflamed.sum() == 0 or quiescent.sum() == 0:
        raise ValueError(f"{path}: need both inflamed and non-inflamed/healthy cells")
    inflamed_centroid = base[inflamed].mean(axis=0)
    quiescent_centroid = base[quiescent].mean(axis=0)
    direction = inflamed_centroid - quiescent_centroid
    rows_out = []
    cell_rows = []
    for gene, token in tqdm(target_tokens.items(), desc=f"{path.stem} genes", unit="gene"):
        for mode in tqdm(["delete", "overexpress"], desc=f"{gene} perturb", leave=False, unit="mode"):
            has_token = np.array([token in t for t in tokens])
            use = has_token if mode == "delete" else np.ones(len(tokens), dtype=bool)
            if use.sum() == 0:
                continue
            idx = np.flatnonzero(use)
            perturbed = [perturb_one(tokens[i], token, mode, cfg.max_len) for i in idx]
            pert = embed_tokens(
                model,
                perturbed,
                cfg,
                pad_id=dictionaries["pad_id"],
                desc=f"{path.stem} {gene} {mode}",
            )
            delta = pert - base[idx]
            proj = projection(delta, direction)
            d_inflamed = cosine_to(pert, inflamed_centroid) - cosine_to(base[idx], inflamed_centroid)
            d_quiet = cosine_to(pert, quiescent_centroid) - cosine_to(base[idx], quiescent_centroid)
            sub = meta.iloc[idx].copy()
            sub["gene"] = gene
            sub["perturbation"] = mode
            sub["projection_toward_inflamed"] = proj
            sub["delta_cosine_inflamed"] = d_inflamed
            sub["delta_cosine_quiescent"] = d_quiet
            sub["target_detected"] = has_token[idx]
            cell_rows.append(sub)
            summary = (
                sub.groupby(["dataset", "disease", "Health", "cell_type", "gene", "perturbation"], observed=True)
                .agg(
                    n_cells=("projection_toward_inflamed", "size"),
                    projection_toward_inflamed_mean=("projection_toward_inflamed", "mean"),
                    projection_toward_inflamed_sem=(
                        "projection_toward_inflamed",
                        lambda s: float(s.std(ddof=1) / math.sqrt(len(s))) if len(s) > 1 else 0.0,
                    ),
                    delta_cosine_inflamed_mean=("delta_cosine_inflamed", "mean"),
                    delta_cosine_quiescent_mean=("delta_cosine_quiescent", "mean"),
                    target_detected_fraction=("target_detected", "mean"),
                )
                .reset_index()
            )
            rows_out.append(summary)
    return pd.concat(rows_out, ignore_index=True), pd.concat(cell_rows, ignore_index=True)


def plot_perturb(summary: pd.DataFrame) -> Path:
    compact = (
        summary.groupby(["dataset", "gene", "perturbation"], observed=True)
        .agg(projection=("projection_toward_inflamed_mean", "mean"))
        .reset_index()
    )
    compact["label"] = compact["gene"] + " " + compact["perturbation"]
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    labels = compact["dataset"] + " | " + compact["label"]
    colors = np.where(compact["projection"] >= 0, "#d55e00", "#0072b2")
    ax.barh(np.arange(len(compact)), compact["projection"], color=colors)
    ax.axvline(0, color="#333333", lw=0.8)
    ax.set_yticks(np.arange(len(compact)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Geneformer embedding projection toward inflamed centroid")
    ax.set_title("In-silico perturbation direction")
    fig.tight_layout()
    out = OUT / "insilico_perturbation.png"
    fig.savefig(out, dpi=190)
    plt.close(fig)
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--subsample", action="store_true")
    p.add_argument("--max-cells-per-dataset", type=int, default=12000)
    p.add_argument("--batch-size", type=int, default=48)
    p.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = PerturbConfig(
        max_cells_per_dataset=args.max_cells_per_dataset,
        batch_size=args.batch_size,
        device=args.device,
    )
    paths = input_paths(args.subsample)
    if not paths:
        raise FileNotFoundError("No sc_01 h5ad files found. Run src/sc_01_qc.py first.")
    gp = geneformer_paths()
    dictionaries = {
        "token_dict": load_pickle(gp["token_dict"]),
        "median_dict": load_pickle(gp["median_dict"]),
        "name_to_id": {str(k).upper(): v for k, v in load_pickle(gp["name_to_id"]).items()},
    }
    dictionaries["pad_id"] = dictionaries["token_dict"].get("<pad>", 0)
    model = AutoModel.from_pretrained(gp["model"]).to(cfg.device).eval()
    if cfg.device.startswith("cuda"):
        model = model.half()
    summaries = []
    cells = []
    for path in paths:
        s, c = run_dataset(path, cfg, model, dictionaries)
        summaries.append(s)
        cells.append(c)
    summary = pd.concat(summaries, ignore_index=True)
    cell_level = pd.concat(cells, ignore_index=True)
    write_tsv(summary, OUT / "insilico_perturbation.tsv")
    write_tsv(cell_level, OUT / "insilico_perturbation_celllevel.tsv")
    fig = plot_perturb(summary)
    P.promote_table(OUT / "insilico_perturbation.tsv")
    append_status(
        "## sc_03 Geneformer perturbation\n"
        f"- inputs: {', '.join(str(p.relative_to(P.root)) for p in paths)}\n"
        f"- model: {gp['model'].relative_to(P.root)}, max_cells_per_dataset={cfg.max_cells_per_dataset}, "
        f"batch_size={cfg.batch_size}, device={cfg.device}\n"
        f"- outputs: outputs/sc_03/insilico_perturbation.tsv, {fig.relative_to(P.root)}\n"
        "- metric: positive projection_toward_inflamed means the perturbation moved the Geneformer "
        "cell embedding toward the inflamed centroid; negative means away from inflammation."
    )
    print("wrote", OUT / "insilico_perturbation.tsv")


if __name__ == "__main__":
    main()
