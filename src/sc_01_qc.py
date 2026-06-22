"""Task 2 single-cell QC/preprocessing for Smillie UC and Martin CD."""
from __future__ import annotations

import argparse
import gzip
import os
import re
import tarfile
import time
from pathlib import Path

_DEFAULT_THREADS = str(max(1, min(24, os.cpu_count() or 1)))
for _key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ.setdefault(_key, _DEFAULT_THREADS)

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
from joblib import Parallel, delayed
from pydantic import BaseModel, ConfigDict, Field
from tqdm.auto import tqdm

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paths import P
from sc_task2_utils import append_status, write_tsv


OUT = P.out("sc_01")
SCRNA = P.data / "scrna"
SCRIPT_T0 = time.perf_counter()


class QCConfig(BaseModel):
    model_config = ConfigDict(frozen=True)

    min_genes: int = Field(default=200, ge=1)
    min_counts: int = Field(default=500, ge=1)
    max_mito_pct: float = Field(default=30.0, gt=0)
    n_top_hvg: int = Field(default=3000, ge=100)
    n_pcs: int = Field(default=50, ge=2)
    scrublet_max_cells: int = Field(default=25000, ge=0)
    n_jobs: int = Field(default=1, ge=1)


SMILLIE_PARTS = {
    "Epi": "5cdc540d328cee7a2efc2348",
    "Fib": "5cdc540d328cee7a2efc2349",
    "Imm": "5cdc540d328cee7a2efc234a",
}

MARTIN_SAMPLES = {
    "GSM3972009": ("69", "Inflamed"),
    "GSM3972010": ("68", "Non-inflamed"),
    "GSM3972011": ("122", "Inflamed"),
    "GSM3972012": ("123", "Non-inflamed"),
    "GSM3972013": ("128", "Inflamed"),
    "GSM3972014": ("129", "Non-inflamed"),
    "GSM3972015": ("135", "Non-inflamed"),
    "GSM3972016": ("138", "Inflamed"),
    "GSM3972017": ("158", "Inflamed"),
    "GSM3972018": ("159", "Non-inflamed"),
    "GSM3972019": ("180", "Non-inflamed"),
    "GSM3972020": ("181", "Inflamed"),
    "GSM3972021": ("186", "Non-inflamed"),
    "GSM3972022": ("187", "Inflamed"),
    "GSM3972023": ("189", "Non-inflamed"),
    "GSM3972024": ("190", "Inflamed"),
    "GSM3972025": ("192", "Non-inflamed"),
    "GSM3972026": ("193", "Inflamed"),
    "GSM3972027": ("195", "Non-inflamed"),
    "GSM3972028": ("196", "Inflamed"),
    "GSM3972029": ("208", "Non-inflamed"),
    "GSM3972030": ("209", "Inflamed"),
}

BROAD_MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "VIL1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"],
    "Myeloid": ["LYZ", "LST1", "CD68", "FCGR3A", "MS4A7"],
    "T_NK": ["CD3D", "CD3E", "TRAC", "NKG7", "KLRD1"],
    "B_Plasma": ["MS4A1", "CD79A", "MZB1", "JCHAIN", "XBP1"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "RAMP2"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
}


def read_lines(path: Path) -> list[str]:
    return pd.read_csv(path, sep="\t", header=None, dtype=str)[0].tolist()


def log(msg: str) -> None:
    elapsed = time.perf_counter() - SCRIPT_T0
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


class stage:
    def __init__(self, msg: str):
        self.msg = msg
        self.t0 = 0.0

    def __enter__(self):
        self.t0 = time.perf_counter()
        log(f"START {self.msg}")
        return self

    def __exit__(self, exc_type, exc, tb):
        dt = time.perf_counter() - self.t0
        if exc_type is None:
            log(f"DONE  {self.msg} ({dt:.1f}s)")
        else:
            log(f"FAIL  {self.msg} after {dt:.1f}s: {exc}")
        return False


def subsample_indices(n: int, keep: int | None, seed: int) -> np.ndarray | slice:
    if not keep or keep <= 0 or keep >= n:
        return slice(None)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, size=keep, replace=False))


def stable_seed(label: str) -> int:
    return int.from_bytes(label.encode("utf-8"), "little", signed=False) % 2**32


def configure_threads(n_jobs: int) -> None:
    n_jobs = max(1, min(n_jobs, os.cpu_count() or 1))
    for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
        os.environ.setdefault(key, str(n_jobs))
    sc.settings.n_jobs = n_jobs
    try:
        import numba

        numba.set_num_threads(n_jobs)
    except Exception:
        pass
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(limits=n_jobs)
    except Exception:
        pass


def read_mtx(path: Path) -> sp.csr_matrix:
    with stage(f"read Matrix Market {path.name}"):
        x = scipy.io.mmread(path).T.tocsr()
    x.data = x.data.astype(np.float32, copy=False)
    return x


def load_smillie_part(part: str, cfg: QCConfig, subsample: int | None) -> ad.AnnData:
    folder = SCRNA / "SCP259" / "SCP259" / "expression" / SMILLIE_PARTS[part]
    with stage(f"load Smillie {part} gene/barcode tables"):
        genes = read_lines(folder / f"{part}.genes.tsv")
        barcodes = read_lines(folder / f"{part}.barcodes2.tsv")
    x = read_mtx(folder / f"gene_sorted-{part}.matrix.mtx")
    idx = subsample_indices(x.shape[0], subsample, seed=stable_seed(part))
    x = x[idx, :]
    obs_names = np.asarray(barcodes, dtype=object)[idx]
    var = pd.DataFrame(index=pd.Index(genes, name="gene_symbol"))
    var["gene_symbol"] = var.index.astype(str)
    obs = pd.DataFrame(index=pd.Index(obs_names, name="NAME"))
    obs["dataset"] = "Smillie_UC"
    obs["disease"] = "UC"
    obs["data_compartment"] = part
    a = ad.AnnData(X=x, obs=obs, var=var)
    a.var_names_make_unique()
    with stage(f"basic QC Smillie {part}"):
        sc.pp.filter_cells(a, min_genes=cfg.min_genes)
        sc.pp.filter_cells(a, min_counts=cfg.min_counts)
    log(f"Smillie {part}: {a.n_obs:,} cells x {a.n_vars:,} genes after basic QC")
    return a


def attach_smillie_metadata(adata: ad.AnnData) -> ad.AnnData:
    with stage("attach Smillie metadata"):
        meta = pd.read_csv(
            SCRNA / "SCP259" / "SCP259" / "metadata" / "all.meta2.txt",
            sep="\t",
            skiprows=[1],
        ).set_index("NAME")
    missing = adata.obs_names.difference(meta.index)
    if len(missing):
        raise ValueError(f"Smillie metadata missing {len(missing)} barcodes")
    cols = ["Cluster", "nGene", "nUMI", "Subject", "Health", "Location", "Sample"]
    adata.obs = adata.obs.join(meta[cols], how="left")
    adata.obs["cell_type"] = adata.obs["Cluster"].astype(str)
    return adata


def tar_member(prefix: str, suffix: str) -> str:
    return f"{prefix}_{suffix}"


def extract_martin_tar() -> Path:
    raw_dir = OUT / "martin_raw"
    done = raw_dir / ".extract_done"
    if done.exists():
        return raw_dir
    raw_dir.mkdir(parents=True, exist_ok=True)
    with stage("extract Martin ileal 10x tar members"):
        tar = tarfile.open(SCRNA / "GSE134809_RAW.tar")
        wanted = []
        for gsm, (num, _) in MARTIN_SAMPLES.items():
            prefix = f"{gsm}_{num}"
            wanted.extend(
                [
                    tar_member(prefix, "barcodes.tsv.gz"),
                    tar_member(prefix, "genes.tsv.gz"),
                    tar_member(prefix, "matrix.mtx.gz"),
                ]
            )
        members = [m for m in tar.getmembers() if m.name in set(wanted)]
        tar.extractall(raw_dir, members=members)
        tar.close()
    done.write_text("ok\n")
    return raw_dir


def read_gzip_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n") for line in fh]


def load_martin_sample(gsm: str, num: str, health: str, raw_dir: Path, cfg: QCConfig) -> ad.AnnData:
    prefix = raw_dir / f"{gsm}_{num}"
    genes_df = pd.read_csv(
        f"{prefix}_genes.tsv.gz",
        sep="\t",
        header=None,
        names=["ensembl_id", "gene_symbol"],
        dtype=str,
    )
    barcodes = read_gzip_lines(Path(f"{prefix}_barcodes.tsv.gz"))
    x = scipy.io.mmread(f"{prefix}_matrix.mtx.gz").T.tocsr()
    x.data = x.data.astype(np.float32, copy=False)
    n_genes = np.asarray((x > 0).sum(axis=1)).ravel()
    n_counts = np.asarray(x.sum(axis=1)).ravel()
    keep = (n_genes >= cfg.min_genes) & (n_counts >= cfg.min_counts)
    x = x[keep, :]
    kept_barcodes = np.asarray(barcodes, dtype=object)[keep]
    obs = pd.DataFrame(index=pd.Index([f"{gsm}_{b}" for b in kept_barcodes], name="cell"))
    obs["dataset"] = "Martin_CD"
    obs["disease"] = "CD"
    obs["Subject"] = num
    obs["Sample"] = f"{gsm}_{num}"
    obs["Health"] = health
    obs["Location"] = "Ileum"
    obs["data_compartment"] = "Ileum"
    var = pd.DataFrame(index=genes_df["ensembl_id"].astype(str).values)
    var["gene_symbol"] = genes_df["gene_symbol"].astype(str).values
    var["ensembl_id"] = genes_df["ensembl_id"].astype(str).values
    a = ad.AnnData(X=x, obs=obs, var=var)
    a.var_names_make_unique()
    sc.pp.filter_cells(a, min_genes=cfg.min_genes)
    return a


def add_mito_and_filter(adata: ad.AnnData, cfg: QCConfig) -> ad.AnnData:
    with stage("mitochondrial QC and raw count layer"):
        if "gene_symbol" not in adata.var:
            adata.var["gene_symbol"] = adata.var_names.astype(str)
        adata.var["mt"] = adata.var["gene_symbol"].astype(str).str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
        before = adata.n_obs
        adata = adata[adata.obs["pct_counts_mt"] <= cfg.max_mito_pct].copy()
        adata.layers["counts"] = adata.X.copy()
        adata.obs["n_counts"] = adata.obs["total_counts"].astype(np.float32)
    log(f"mito filter retained {adata.n_obs:,}/{before:,} cells")
    return adata


def maybe_scrublet(adata: ad.AnnData, cfg: QCConfig) -> ad.AnnData:
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False
    if cfg.scrublet_max_cells == 0:
        return adata
    try:
        import scrublet as scr
    except Exception:
        return adata
    for label, idx in tqdm(
        adata.obs.groupby("dataset", observed=True).indices.items(),
        desc="scrublet groups",
        unit="group",
    ):
        idx = np.asarray(idx)
        run_idx = idx
        if len(idx) > cfg.scrublet_max_cells:
            run_idx = np.random.default_rng(20260622).choice(
                idx, size=cfg.scrublet_max_cells, replace=False
            )
        counts = adata[run_idx].layers["counts"]
        scrub = scr.Scrublet(counts, expected_doublet_rate=0.06)
        scores, calls = scrub.scrub_doublets(verbose=False)
        adata.obs.iloc[run_idx, adata.obs.columns.get_loc("doublet_score")] = scores
        adata.obs.iloc[run_idx, adata.obs.columns.get_loc("predicted_doublet")] = calls
        print(f"scrublet {label}: ran {len(run_idx):,} / {len(idx):,} cells")
    return adata[~adata.obs["predicted_doublet"].fillna(False)].copy()


def preprocess_embedding(adata: ad.AnnData, cfg: QCConfig) -> ad.AnnData:
    with stage("normalize_total + log1p"):
        adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    with stage("highly variable genes"):
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=min(cfg.n_top_hvg, adata.n_vars - 1),
            flavor="seurat",
            subset=False,
        )
    with stage("scale sparse matrix"):
        sc.pp.scale(adata, zero_center=False, max_value=10)
    with stage("PCA"):
        sc.tl.pca(adata, n_comps=min(cfg.n_pcs, adata.n_vars - 1), use_highly_variable=True)
    with stage("neighbors graph"):
        sc.pp.neighbors(adata, n_pcs=min(cfg.n_pcs, adata.obsm["X_pca"].shape[1]))
    with stage("UMAP"):
        sc.tl.umap(adata)
    with stage("Leiden clustering"):
        sc.tl.leiden(
            adata,
            key_added="leiden",
            resolution=1.0,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    return adata


def mean_marker_scores(adata: ad.AnnData, markers: dict[str, list[str]]) -> pd.DataFrame:
    symbols = pd.Index(adata.var["gene_symbol"].astype(str).str.upper())
    frames = {}
    for label, genes in markers.items():
        loc = [np.flatnonzero(symbols == g)[0] for g in genes if len(np.flatnonzero(symbols == g))]
        if not loc:
            frames[label] = np.zeros(adata.n_obs, dtype=np.float32)
            continue
        sub = adata[:, loc].X
        vals = np.asarray(sub.mean(axis=1)).ravel()
        frames[label] = vals.astype(np.float32, copy=False)
    return pd.DataFrame(frames, index=adata.obs_names)


def annotate_martin_celltypes(adata: ad.AnnData) -> ad.AnnData:
    work = adata.copy()
    if "counts" in work.layers:
        work.X = work.layers["counts"].copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    scores = mean_marker_scores(work, BROAD_MARKERS)
    labels = scores.idxmax(axis=1)
    labels[scores.max(axis=1) <= 0] = "Unknown"
    adata.obs["cell_type"] = labels.reindex(adata.obs_names).values
    adata.obs["Cluster"] = adata.obs["cell_type"]
    return adata


def qc_plot(adata: ad.AnnData, stem: str) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    axes[0].hist(adata.obs["total_counts"], bins=80, color="#4477aa")
    axes[0].set_title("UMI counts")
    axes[0].set_xlabel("total_counts")
    axes[1].hist(adata.obs["n_genes_by_counts"], bins=80, color="#66aa55")
    axes[1].set_title("Detected genes")
    axes[1].set_xlabel("n_genes")
    if "X_umap" in adata.obsm:
        color = adata.obs["Health"].astype(str).map(
            {"Healthy": "#4c78a8", "Non-inflamed": "#59a14f", "Inflamed": "#e15759"}
        ).fillna("#999999")
        axes[2].scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1], s=1, c=color)
        axes[2].set_title("UMAP by Health")
        axes[2].set_xticks([])
        axes[2].set_yticks([])
    fig.tight_layout()
    out = OUT / f"{stem}_qc.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def build_smillie(cfg: QCConfig, subsample_per_compartment: int | None) -> Path:
    t0 = time.time()
    log("building Smillie UC")
    parts = [
        load_smillie_part(k, cfg, subsample_per_compartment)
        for k in tqdm(SMILLIE_PARTS, desc="Smillie compartments", unit="part")
    ]
    with stage("concat Smillie compartments"):
        adata = ad.concat(parts, axis=0, join="outer", merge="same", fill_value=0)
        log(f"Smillie concat: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    adata = attach_smillie_metadata(adata)
    adata = add_mito_and_filter(adata, cfg)
    adata = maybe_scrublet(adata, cfg)
    adata = preprocess_embedding(adata, cfg)
    out = OUT / ("smillie_uc_subsample.h5ad" if subsample_per_compartment else "smillie_uc_qc.h5ad")
    with stage(f"write {out.name}"):
        adata.write_h5ad(out, compression="gzip")
    with stage("write Smillie QC plot and obs table"):
        qc_plot(adata, "smillie_uc_subsample" if subsample_per_compartment else "smillie_uc")
        write_tsv(adata.obs.reset_index(names="cell"), OUT / out.with_suffix(".obs.tsv").name)
    print(f"wrote {out} cells={adata.n_obs:,} genes={adata.n_vars:,} in {(time.time()-t0)/60:.1f} min")
    return out


def build_martin(cfg: QCConfig, subsample_per_sample: int | None) -> Path:
    t0 = time.time()
    log("building Martin CD")
    raw_dir = extract_martin_tar()
    samples = []
    def one(item: tuple[str, tuple[str, str]]) -> ad.AnnData:
        gsm, (num, health) = item
        a = load_martin_sample(gsm, num, health, raw_dir, cfg)
        if subsample_per_sample and subsample_per_sample < a.n_obs:
            idx = subsample_indices(a.n_obs, subsample_per_sample, seed=int(num))
            a = a[idx, :].copy()
        print(f"Martin {gsm}_{num} {health}: {a.n_obs:,} cells")
        return a

    items = list(MARTIN_SAMPLES.items())
    if cfg.n_jobs > 1 and len(items) > 1:
        samples = Parallel(n_jobs=min(cfg.n_jobs, len(items)), prefer="processes")(
            delayed(one)(item) for item in tqdm(items, desc="Martin samples", unit="sample")
        )
    else:
        for item in tqdm(items, desc="Martin samples", unit="sample"):
            samples.append(one(item))
    with stage("concat Martin samples"):
        adata = ad.concat(samples, axis=0, join="inner", merge="same")
        log(f"Martin concat: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    adata = add_mito_and_filter(adata, cfg)
    adata = maybe_scrublet(adata, cfg)
    adata = annotate_martin_celltypes(adata)
    adata = preprocess_embedding(adata, cfg)
    out = OUT / ("martin_cd_subsample.h5ad" if subsample_per_sample else "martin_cd_qc.h5ad")
    with stage(f"write {out.name}"):
        adata.write_h5ad(out, compression="gzip")
    with stage("write Martin QC plot and obs table"):
        qc_plot(adata, "martin_cd_subsample" if subsample_per_sample else "martin_cd")
        write_tsv(adata.obs.reset_index(names="cell"), OUT / out.with_suffix(".obs.tsv").name)
    print(f"wrote {out} cells={adata.n_obs:,} genes={adata.n_vars:,} in {(time.time()-t0)/60:.1f} min")
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", choices=["smillie", "martin", "all"], default="all")
    p.add_argument("--subsample-per-compartment", type=int, default=0)
    p.add_argument("--subsample-per-martin-sample", type=int, default=0)
    p.add_argument("--min-genes", type=int, default=200)
    p.add_argument("--min-counts", type=int, default=500)
    p.add_argument("--max-mito-pct", type=float, default=30.0)
    p.add_argument("--scrublet-max-cells", type=int, default=25000)
    p.add_argument("--n-jobs", type=int, default=max(1, min(24, os.cpu_count() or 1)))
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = QCConfig(
        min_genes=args.min_genes,
        min_counts=args.min_counts,
        max_mito_pct=args.max_mito_pct,
        scrublet_max_cells=args.scrublet_max_cells,
        n_jobs=args.n_jobs,
    )
    configure_threads(cfg.n_jobs)
    log(f"configured n_jobs={cfg.n_jobs}, scanpy_n_jobs={sc.settings.n_jobs}")
    outputs = []
    if args.dataset in {"smillie", "all"}:
        outputs.append(build_smillie(cfg, args.subsample_per_compartment or None))
    if args.dataset in {"martin", "all"}:
        outputs.append(build_martin(cfg, args.subsample_per_martin_sample or None))
    append_status(
        "## sc_01 QC\n"
        f"- outputs: {', '.join(str(x.relative_to(P.root)) for x in outputs)}\n"
        "- Martin Health mapping follows GEO GSE134809 sample titles: Ileal Involved=Inflamed, "
        "Ileal Uninvolved=Non-inflamed; PBMC samples were excluded from gut replication."
    )


if __name__ == "__main__":
    main()
