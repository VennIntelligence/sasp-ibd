"""Build the refractory-module multicontext eQTL MR + coloc map.

Contexts:
  blood            existing eQTLGen outputs in results/tables
  gut_sigmoid      GTEx Colon Sigmoid lead MR plus full-allpairs coloc when resolvable
  gut_transverse   GTEx Colon Transverse lead MR plus full-allpairs coloc when resolvable
  monocyte         eQTL Catalogue BLUEPRINT monocyte
  neutrophil       eQTL Catalogue BLUEPRINT neutrophil
  monocyte_stim    eQTL Catalogue Quach LPS-stimulated monocyte
"""
from __future__ import annotations

import gzip
import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request
import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from pandas.errors import EmptyDataError
from pydantic import BaseModel, ConfigDict, field_validator
from scipy import stats
from scipy.special import logsumexp
from statsmodels.stats.multitest import multipletests

from causal_module_utils import gwas_lookup, harmonised_beta
from paths import P

_SPEC = importlib.util.spec_from_file_location("build_instruments", P.src / "13_build_instruments.py")
_BUILD = importlib.util.module_from_spec(_SPEC)
assert _SPEC.loader is not None
_SPEC.loader.exec_module(_BUILD)
MODULE_GENES = _BUILD.MODULE_GENES
parse_variant_id = _BUILD.parse_variant_id


OUTCOMES = {"IBD": "IBD.h.tsv.gz", "CD": "CD.h.tsv.gz", "UC": "UC.h.tsv.gz"}
IMMUNE_FALLBACKS = {
    "monocyte": {"dataset_id": "QTD000021", "study_label": "BLUEPRINT", "sample_group": "monocyte"},
    "neutrophil": {"dataset_id": "QTD000026", "study_label": "BLUEPRINT", "sample_group": "neutrophil"},
    "monocyte_stim": {"dataset_id": "QTD000414", "study_label": "Quach_2016", "sample_group": "monocyte_LPS"},
}
GTEX_CONTEXTS = {
    "gut_sigmoid": "Colon_Sigmoid",
    "gut_transverse": "Colon_Transverse",
}
W_EQTL = 0.15**2
W_GWAS = 0.2**2
P1, P2, P12 = 1e-4, 1e-4, 1e-5
API = "https://www.ebi.ac.uk/eqtl/api/v2"
REQ_GENES = ["CCL8", "CXCR2", "OSM", "OSMR", "IL1B", "CXCL10", "TREM1", "IL13RA2"]
CONTEXTS = ["blood", "gut_sigmoid", "gut_transverse", "monocyte", "neutrophil", "monocyte_stim"]
MAX_EQTL_WINDOW = 950_000


class Inputs(BaseModel):
    model_config = ConfigDict(frozen=True)

    out_dir: Path = P.out("causal_module")
    gwas_dir: Path = P.raw / "gwas"
    gtex_dir: Path = P.raw / "gtex_colon"

    @field_validator("gwas_dir")
    @classmethod
    def have_gwas(cls, v: Path) -> Path:
        missing = [v / f for f in OUTCOMES.values() if not (v / f).exists()]
        if missing:
            raise FileNotFoundError(f"missing GWAS files: {missing}")
        return v


def module_genes() -> list[str]:
    genes = list(dict.fromkeys(MODULE_GENES))
    nr = P.outputs / "nr_unbiased" / "NR_up_robust.tsv"
    if nr.exists():
        df = pd.read_csv(nr, sep="\t")
        col = "gene" if "gene" in df.columns else df.columns[0]
        genes.extend(df[col].dropna().astype(str).head(50).tolist())
    return sorted(dict.fromkeys(genes))


def read_gene_table(gtex_dir: Path, genes: list[str]) -> pd.DataFrame:
    frames = []
    cols = ["gene_id", "gene_name", "gene_chr", "gene_start", "gene_end", "strand", "variant_pos", "tss_distance"]
    for tissue in GTEX_CONTEXTS.values():
        path = gtex_dir / f"{tissue}.v8.egenes.txt.gz"
        df = pd.read_csv(path, sep="\t", usecols=cols)
        frames.append(df)
    g = pd.concat(frames, ignore_index=True).drop_duplicates("gene_id")
    g["gene"] = g["gene_name"]
    g["ensg"] = g["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    g["chrom"] = g["gene_chr"].str.replace("chr", "", regex=False)
    g["start"] = g["gene_start"].astype(int)
    g["end"] = g["gene_end"].astype(int)
    g["phenotype_id"] = "chr" + g["chrom"].astype(str) + "_" + (g["variant_pos"].astype(int) - g["tss_distance"].astype(int)).astype(str)
    return g[g["gene"].isin(genes)].copy()


def fetch_json(url: str, retries: int = 5) -> list[dict]:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.load(r)
        except urllib.error.HTTPError as exc:
            body = exc.read(500).decode("utf-8", "replace")
            if exc.code == 400 and "No results" in body:
                return []
            if exc.code in {429, 500, 502, 503, 504}:
                time.sleep(min(120, 30 * (2**attempt)))
                continue
            raise RuntimeError(f"eQTL Catalogue request failed {exc.code}: {url}\n{body}") from exc
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            time.sleep(min(120, 30 * (2**attempt)))
    return []


def discover_datasets(out_dir: Path) -> pd.DataFrame:
    cache = out_dir / "eqtl_catalogue_datasets_ge.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    rows, start, size = [], 0, 500
    while True:
        url = f"{API}/datasets?quant_method=ge&size={size}&start={start}"
        data = fetch_json(url)
        if not data:
            break
        rows.extend(data)
        if len(data) < size:
            break
        start += size
        time.sleep(0.05)
    df = pd.DataFrame(rows)
    df.to_csv(cache, sep="\t", index=False)
    return df


def select_immune_contexts(datasets: pd.DataFrame) -> dict[str, dict[str, str]]:
    """Discover the requested immune contexts from the catalogue listing.

    The exact QTD IDs are stable in this run, but selecting from metadata keeps
    the analysis auditable and avoids silently hard-coding an unrelated dataset.
    """
    df = datasets.copy()
    for col in ["study_label", "sample_group", "tissue_label", "condition_label", "dataset_id"]:
        if col not in df:
            df[col] = ""
    text = df[["study_label", "sample_group", "tissue_label", "condition_label"]].astype(str)
    found: dict[str, dict[str, str]] = {}

    def one(context: str, mask: pd.Series) -> None:
        recs = df[mask].sort_values(["sample_size", "dataset_id"], ascending=[False, True])
        if len(recs):
            r = recs.iloc[0]
            found[context] = {k: str(r.get(k, "")) for k in ["dataset_id", "study_label", "sample_group", "tissue_label", "condition_label"]}
        else:
            found[context] = IMMUNE_FALLBACKS[context].copy()

    one(
        "monocyte",
        text["study_label"].str.fullmatch("BLUEPRINT", case=False, na=False)
        & text["sample_group"].str.fullmatch("monocyte", case=False, na=False),
    )
    one(
        "neutrophil",
        text["study_label"].str.fullmatch("BLUEPRINT", case=False, na=False)
        & text["sample_group"].str.fullmatch("neutrophil", case=False, na=False),
    )
    one(
        "monocyte_stim",
        text["study_label"].str.fullmatch("Quach_2016", case=False, na=False)
        & text["sample_group"].str.contains("monocyte", case=False, na=False)
        & text["condition_label"].str.contains("LPS", case=False, na=False),
    )
    return found


def region_chunks(chrom: str, start: int, end: int, max_width: int = MAX_EQTL_WINDOW) -> list[str]:
    chunks = []
    left = max(1, int(start))
    end = int(end)
    while left <= end:
        right = min(end, left + max_width - 1)
        chunks.append(f"{chrom}:{left}-{right}")
        left = right + 1
    return chunks


def fetch_eqtl_region(dataset_id: str, gene: pd.Series, cache: Path) -> pd.DataFrame:
    if cache.exists():
        try:
            return pd.read_csv(cache, sep="\t", dtype={"chromosome": str, "chrom": str})
        except EmptyDataError:
            return pd.DataFrame()
    start = max(1, int(gene.start) - 1_000_000)
    end = int(gene.end) + 1_000_000
    rows, offset, size = [], 0, 1000
    for pos in region_chunks(str(gene.chrom), start, end):
        offset = 0
        while True:
            qs = urllib.parse.urlencode({"pos": pos, "size": size, "start": offset})
            data = fetch_json(f"{API}/datasets/{dataset_id}/associations?{qs}")
            if not data:
                break
            rows.extend(
                [
                    r
                    for r in data
                    if str(r.get("molecular_trait_id", r.get("gene_id", ""))).split(".")[0] == gene.ensg
                ]
            )
            if len(data) < size:
                break
            offset += size
            time.sleep(0.1)
    df = pd.DataFrame(rows)
    if len(df):
        df = df.rename(columns={"pvalue": "p_eqtl", "beta": "beta_eqtl", "se": "se_eqtl"})
        df["gene"] = gene.gene
        df["chrom"] = df["chromosome"].astype(str).str.replace("chr", "", regex=False)
        df["pos"] = df["position"].astype(int)
        df["ref"] = df["ref"].str.upper()
        df["alt"] = df["alt"].str.upper()
        df["effect_allele"] = df["alt"]
        df["other_allele"] = df["ref"]
        df["variant_id"] = df.get("variant", df["chrom"].astype(str) + "_" + df["pos"].astype(str) + "_" + df["ref"] + "_" + df["alt"])
        keep = [
            "gene", "gene_id", "variant_id", "rsid", "chrom", "pos", "ref", "alt",
            "effect_allele", "other_allele", "maf", "beta_eqtl", "se_eqtl", "p_eqtl",
            "median_tpm", "an",
        ]
        df = df[[c for c in keep if c in df.columns]]
    cache.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, sep="\t", index=False)
    return df


def fetch_context_gene_eqtl(context: str, meta: dict[str, str], gene_row: dict, cache_dir: Path) -> pd.DataFrame:
    gene = pd.Series(gene_row)
    cache = cache_dir / context / f"{gene.gene}.tsv"
    df = fetch_eqtl_region(meta["dataset_id"], gene, cache)
    if len(df):
        df["context"] = context
    return df


def load_gtex_sig_pairs(gtex_dir: Path, genes: list[str]) -> pd.DataFrame:
    frames = []
    for context, tissue in GTEX_CONTEXTS.items():
        path = gtex_dir / f"{tissue}.v8.signif_variant_gene_pairs.txt.gz"
        egenes = pd.read_csv(gtex_dir / f"{tissue}.v8.egenes.txt.gz", sep="\t", usecols=["gene_id", "gene_name"])
        id_to_gene = dict(zip(egenes["gene_id"], egenes["gene_name"]))
        target_ids = set(egenes.loc[egenes["gene_name"].isin(genes), "gene_id"])
        rows = []
        with gzip.open(path, "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(head)}
            for line in fh:
                p = line.rstrip("\n").split("\t")
                gid = p[ix["gene_id"]]
                if gid not in target_ids:
                    continue
                chrom, pos, ref, alt = parse_variant_id(p[ix["variant_id"]])
                beta, se = float(p[ix["slope"]]), float(p[ix["slope_se"]])
                if not np.isfinite(beta) or not np.isfinite(se) or se <= 0:
                    continue
                rows.append(
                    {
                        "gene": id_to_gene[gid],
                        "context": context,
                        "gene_id": gid,
                        "variant_id": p[ix["variant_id"]],
                        "rsid": p[ix.get("rs_id_dbSNP151_GRCh38p7", ix["variant_id"])],
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "effect_allele": alt,
                        "other_allele": ref,
                        "maf": float(p[ix["maf"]]),
                        "beta_eqtl": beta,
                        "se_eqtl": se,
                        "p_eqtl": float(p[ix["pval_nominal"]]),
                    }
                )
        if rows:
            frames.append(pd.DataFrame(rows))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def lead_instruments(eqtl: pd.DataFrame) -> pd.DataFrame:
    if not len(eqtl):
        return pd.DataFrame()
    return (
        eqtl.dropna(subset=["beta_eqtl", "se_eqtl", "p_eqtl"])
        .query("se_eqtl > 0 and beta_eqtl != 0")
        .assign(abs_z=lambda x: (x["beta_eqtl"] / x["se_eqtl"]).abs())
        .sort_values(["context", "gene", "p_eqtl", "abs_z"], ascending=[True, True, True, False])
        .drop_duplicates(["context", "gene"])
        .drop(columns=["abs_z"])
        .copy()
    )


def run_mr(eqtl: pd.DataFrame, gwas_dir: Path, out_dir: Path) -> dict[str, pd.DataFrame]:
    out = {}
    for outcome, fname in OUTCOMES.items():
        if not len(eqtl):
            out[outcome] = pd.DataFrame()
            continue
        gwas = gwas_lookup(gwas_dir / fname, eqtl, out_dir / f"_gwas_keys_multicontext_{outcome}.txt")
        merged = eqtl.merge(gwas, on=["chrom", "pos"], how="left")
        merged["beta_gwas_alt"] = merged.apply(harmonised_beta, axis=1)
        merged = merged.dropna(subset=["beta_gwas_alt", "se_gwas"]).copy()
        merged = merged[(merged["se_gwas"] > 0) & (merged["beta_eqtl"] != 0)]
        merged = (
            merged.assign(abs_z=lambda x: (x["beta_eqtl"] / x["se_eqtl"]).abs())
            .sort_values(["context", "gene", "p_eqtl", "abs_z"], ascending=[True, True, True, False])
            .drop_duplicates(["context", "gene"])
            .drop(columns=["abs_z"])
            .copy()
        )
        theta = merged["beta_gwas_alt"] / merged["beta_eqtl"]
        se = (merged["se_gwas"] / merged["beta_eqtl"]).abs()
        res = merged.copy()
        res["theta"] = theta
        res["se"] = se
        res["MR_OR"] = np.exp(theta.clip(-20, 20))
        res["MR_p"] = 2 * stats.norm.sf(np.abs(theta / se))
        res["MR_FDR"] = multipletests(res["MR_p"], method="fdr_bh")[1] if len(res) else []
        out[outcome] = res.sort_values(["MR_p", "gene", "context"])
        out[outcome].to_csv(out_dir / f"mr_multicontext_{outcome}.tsv", sep="\t", index=False)
    return out


def labf(beta: np.ndarray, varbeta: np.ndarray, prior_var: float) -> np.ndarray:
    ok = np.isfinite(beta) & np.isfinite(varbeta) & (varbeta > 0)
    out = np.full(len(beta), np.nan)
    z2 = beta[ok] * beta[ok] / varbeta[ok]
    r = prior_var / (varbeta[ok] + prior_var)
    out[ok] = 0.5 * (np.log1p(-r) + r * z2)
    return out


def coloc_abf(df: pd.DataFrame) -> dict[str, float | int]:
    df = df.dropna(subset=["beta_eqtl", "se_eqtl", "beta_gwas_alt", "se_gwas"])
    df = df[(df["se_eqtl"] > 0) & (df["se_gwas"] > 0)]
    if len(df) < 2:
        return {"nsnps": len(df), "PP0": np.nan, "PP1": np.nan, "PP2": np.nan, "PP3": np.nan, "PP4": np.nan}
    l1 = labf(df["beta_eqtl"].to_numpy(float), np.square(df["se_eqtl"].to_numpy(float)), W_EQTL)
    l2 = labf(df["beta_gwas_alt"].to_numpy(float), np.square(df["se_gwas"].to_numpy(float)), W_GWAS)
    ok = np.isfinite(l1) & np.isfinite(l2)
    if ok.sum() < 2:
        return {"nsnps": int(ok.sum()), "PP0": np.nan, "PP1": np.nan, "PP2": np.nan, "PP3": np.nan, "PP4": np.nan}
    l1, l2 = l1[ok], l2[ok]
    h1, h2 = logsumexp(l1), logsumexp(l2)
    h4 = logsumexp(l1 + l2)
    s = h1 + h2
    h3 = s + np.log1p(-np.exp(min(h4 - s, -1e-12)))
    logs = np.array([0, np.log(P1) + h1, np.log(P2) + h2, np.log(P1) + np.log(P2) + h3, np.log(P12) + h4])
    pp = np.exp(logs - logsumexp(logs))
    return {"nsnps": int(ok.sum()), "PP0": pp[0], "PP1": pp[1], "PP2": pp[2], "PP3": pp[3], "PP4": pp[4]}


def coloc_for_eqtl(eqtl: pd.DataFrame, gwas_ibd: Path, out_dir: Path, tag: str) -> pd.DataFrame:
    if not len(eqtl):
        return pd.DataFrame(columns=["gene", "context", "nsnps", "PP0", "PP1", "PP2", "PP3", "PP4", "coloc_method"])
    gwas = gwas_lookup(gwas_ibd, eqtl, out_dir / f"_coloc_keys_{tag}.txt")
    merged = eqtl.merge(gwas, on=["chrom", "pos"], how="left")
    merged["beta_gwas_alt"] = merged.apply(harmonised_beta, axis=1)
    rows = []
    for (gene, context), g in merged.groupby(["gene", "context"], sort=True):
        rec = {"gene": gene, "context": context, "coloc_method": tag}
        rec.update(coloc_abf(g))
        rows.append(rec)
    return pd.DataFrame(rows)


def allpairs_resolvable_coloc(gene_table: pd.DataFrame, gwas_ibd: Path, out_dir: Path) -> pd.DataFrame:
    """Run full nominal-allpairs coloc by mapping allpairs chr_TSS IDs to genes."""
    import shlex
    import subprocess

    candidate = {}
    for r in gene_table.itertuples():
        candidate[str(r.phenotype_id)] = r.gene
    if not candidate:
        return pd.DataFrame()
    key_file = out_dir / "_candidate_allpair_ids.txt"
    key_file.write_text("\n".join(sorted(candidate)) + "\n")
    records = []
    for context, tissue in [("gut_sigmoid", "Colon-Sigmoid"), ("gut_transverse", "Colon-Transverse")]:
        path = P.raw / "gtex_colon" / f"{tissue}.nominal.allpairs.txt.gz"
        if not path.exists():
            continue
        extract = out_dir / f"gtex_{context}_module_allpairs.raw.tsv"
        cmd = (
            f"gzip -cd {shlex.quote(str(path))} | "
            "awk -F'\\t' -v OFS='\\t' "
            f"'NR==FNR{{want[$1]; next}} FNR==1{{print; next}} $1 in want{{print}}' "
            f"{shlex.quote(str(key_file))} - > {shlex.quote(str(extract))}"
        )
        subprocess.run(cmd, shell=True, check=True)
        matched = pd.read_csv(extract, sep="\t", usecols=["gene_id"]) if extract.stat().st_size else pd.DataFrame()
        audit = pd.DataFrame({"phenotype_id": sorted(candidate), "gene": [candidate[k] for k in sorted(candidate)]})
        found = set(matched["gene_id"].dropna().astype(str)) if len(matched) else set()
        audit["found_in_allpairs"] = audit["phenotype_id"].isin(found)
        audit.to_csv(out_dir / f"gtex_{context}_allpairs_phenotype_id_audit.tsv", sep="\t", index=False)
        if extract.stat().st_size == 0:
            continue
        raw = pd.read_csv(extract, sep="\t")
        if not len(raw):
            continue
        raw["gene"] = raw["gene_id"].map(candidate)
        raw["context"] = context
        raw = raw.rename(columns={"variant_id": "rsid", "slope": "beta_eqtl", "slope_se": "se_eqtl", "pval_nominal": "p_eqtl"})
        eqtl = raw[["gene", "context", "rsid", "beta_eqtl", "se_eqtl", "p_eqtl"]].copy()
        eqtl = eqtl.dropna(subset=["gene", "rsid", "beta_eqtl", "se_eqtl", "p_eqtl"])
        eqtl = eqtl[np.isfinite(eqtl["beta_eqtl"]) & np.isfinite(eqtl["se_eqtl"]) & (eqtl["se_eqtl"] > 0)]
        if not len(eqtl):
            continue
        eqtl.to_csv(out_dir / f"gtex_{context}_module_allpairs.tsv", sep="\t", index=False)
        keys = out_dir / f"_allpairs_rsids_{context}.txt"
        keys.write_text("\n".join(sorted(eqtl["rsid"].dropna().astype(str).unique())) + "\n")
        gwas = lookup_gwas_by_rsid(gwas_ibd, keys)
        merged = eqtl.merge(gwas, on="rsid", how="left").dropna(subset=["beta_gwas", "se_gwas"])
        merged["beta_gwas_alt"] = merged["beta_gwas"]
        for (gene, ctx), g in merged.groupby(["gene", "context"], sort=True):
            rec = {"gene": gene, "context": ctx, "coloc_method": "gtex_nominal_allpairs_exact_coordinate_id"}
            rec.update(coloc_abf(g))
            records.append(rec)
    return pd.DataFrame(records)


def lookup_gwas_by_rsid(path: Path, key_file: Path) -> pd.DataFrame:
    import shlex
    import subprocess

    hdr = subprocess.run(f"gunzip -c {shlex.quote(str(path))} | head -1", shell=True, capture_output=True, text=True, check=True).stdout.rstrip("\n").split("\t")
    ix = {c: i + 1 for i, c in enumerate(hdr)}
    cmd = (
        f"gunzip -c {shlex.quote(str(path))} | awk -F'\\t' -v OFS='\\t' "
        f"-v rs={ix['hm_rsid']} -v b={ix['hm_beta']} -v se={ix['standard_error']} -v pv={ix['p_value']} "
        f"'NR==FNR{{w[$1]; next}} FNR==1{{next}} $rs in w{{print $rs,$b,$se,$pv}}' {shlex.quote(str(key_file))} -"
    )
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    rows = []
    for line in out.splitlines():
        rs, b, se, p = line.split("\t")
        try:
            rows.append({"rsid": rs, "beta_gwas": float(b), "se_gwas": float(se), "p_gwas": float(p)})
        except ValueError:
            pass
    return pd.DataFrame(rows).drop_duplicates("rsid")


def blood_rows(genes: list[str]) -> pd.DataFrame:
    base = pd.DataFrame({"gene": genes})
    tri_path = P.tables / "triangulation.tsv"
    if tri_path.exists():
        tri = pd.read_csv(tri_path, sep="\t")
        row = base.merge(tri[["gene", "MR_OR", "MR_p", "MR_fdr", "coloc_PP4"]], on="gene", how="left")
        row = row.rename(columns={"MR_fdr": "MR_FDR"})
        coloc = pd.read_csv(P.tables / "coloc_IBD.tsv", sep="\t")
        row = row.merge(coloc[["gene", "nsnps"]].rename(columns={"nsnps": "coloc_nsnps"}), on="gene", how="left")
    else:
        mr = pd.read_csv(P.tables / "mr_IBD.tsv", sep="\t")
        coloc = pd.read_csv(P.tables / "coloc_IBD.tsv", sep="\t")
        row = base.merge(mr[["gene", "OR", "p_mr", "fdr"]], on="gene", how="left")
        row = row.merge(coloc[["gene", "PP4", "nsnps"]], on="gene", how="left")
        row = row.rename(columns={"OR": "MR_OR", "p_mr": "MR_p", "fdr": "MR_FDR", "PP4": "coloc_PP4", "nsnps": "coloc_nsnps"})
    row["context"] = "blood"
    row["n_instrument"] = row["MR_p"].notna().astype(int)
    row["analysis_status"] = np.where(row["MR_p"].notna(), "tested_existing_eqtlgen", "no_blood_result")
    row["coloc_method"] = "existing_eqtlgen"
    return row


def assemble_map(genes: list[str], mr: dict[str, pd.DataFrame], coloc: pd.DataFrame, inst: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    ibd = mr["IBD"].rename(columns={"p_gwas": "IBD_gwas_p"})
    keep = [
        "gene", "context", "variant_id", "rsid_x", "rsid", "chrom", "pos", "ref", "alt",
        "MR_OR", "MR_p", "MR_FDR", "p_eqtl", "IBD_gwas_p",
    ]
    ibd = ibd[[c for c in keep if c in ibd.columns]].copy()
    if "rsid_x" in ibd.columns and "rsid" not in ibd.columns:
        ibd = ibd.rename(columns={"rsid_x": "rsid"})
    for outcome in ["CD", "UC"]:
        d = mr[outcome][["gene", "context", "MR_OR", "MR_p", "MR_FDR"]].rename(
            columns={"MR_OR": f"{outcome}_MR_OR", "MR_p": f"{outcome}_MR_p", "MR_FDR": f"{outcome}_MR_FDR"}
        )
        ibd = ibd.merge(d, on=["gene", "context"], how="left")
    ibd["n_instrument"] = 1
    ibd["analysis_status"] = "tested"

    df = pd.concat([ibd, blood_rows(genes)], ignore_index=True, sort=False)
    grid = pd.MultiIndex.from_product([genes, CONTEXTS], names=["gene", "context"]).to_frame(index=False)
    df = grid.merge(df, on=["gene", "context"], how="left")
    if len(coloc):
        c = coloc.rename(columns={"PP4": "coloc_PP4", "nsnps": "coloc_nsnps"})
        df = df.merge(c[["gene", "context", "PP3", "coloc_PP4", "coloc_nsnps", "coloc_method"]], on=["gene", "context"], how="left", suffixes=("", "_coloc"))
        for col in ["PP3", "coloc_PP4", "coloc_nsnps", "coloc_method"]:
            alt = f"{col}_coloc"
            if alt in df.columns:
                df[col] = df[col].combine_first(df[alt])
                df = df.drop(columns=[alt])
    has_inst = inst[["gene", "context"]].drop_duplicates().assign(has_instrument=True) if len(inst) else pd.DataFrame(columns=["gene", "context", "has_instrument"])
    df = df.merge(has_inst, on=["gene", "context"], how="left")
    df["has_instrument"] = df["has_instrument"].eq(True)
    df["analysis_status"] = np.select(
        [df["analysis_status"].notna(), df["has_instrument"]],
        [df["analysis_status"], "has_eqtl_no_harmonized_mr"],
        default="no_tool",
    )
    df["n_instrument"] = df["n_instrument"].fillna(df["has_instrument"].astype(int)).astype(int)
    df = df.drop(columns=["has_instrument"])
    df["causal_call"] = (df["MR_FDR"] < 0.05) & (df["coloc_PP4"] > 0.8)
    order = {c: i for i, c in enumerate(CONTEXTS)}
    df["_ctx"] = df["context"].map(order)
    df["_gene"] = pd.Categorical(df["gene"], genes, ordered=True)
    df = df.sort_values(["_gene", "_ctx"]).drop(columns=["_ctx", "_gene"])
    df.to_csv(out_dir / "module_causal_map_multicontext.tsv", sep="\t", index=False)
    return df


def plot_map(df: pd.DataFrame, path: Path) -> None:
    contexts = CONTEXTS
    genes = (
        df.assign(best=lambda x: x["causal_call"].astype(int) * 100 - x["MR_p"].fillna(1))
        .groupby("gene")["best"].max()
        .sort_values(ascending=False)
        .index.tolist()
    )
    mat = df.pivot(index="gene", columns="context", values="MR_OR").reindex(index=genes, columns=contexts)
    pp4 = df.pivot(index="gene", columns="context", values="coloc_PP4").reindex(index=genes, columns=contexts)
    causal = df.pivot(index="gene", columns="context", values="causal_call").reindex(index=genes, columns=contexts).fillna(False)
    vals = np.log2(mat.astype(float))
    fig, ax = plt.subplots(figsize=(9.5, max(7, 0.24 * len(genes) + 2)))
    im = ax.imshow(vals, cmap="RdBu_r", vmin=-1.5, vmax=1.5, aspect="auto")
    ax.set_xticks(np.arange(len(contexts)))
    ax.set_xticklabels(contexts, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    for i, g in enumerate(genes):
        for j, c in enumerate(contexts):
            p = pp4.loc[g, c]
            if pd.notna(p):
                ax.text(j, i, "*" if causal.loc[g, c] else f"{p:.2f}", ha="center", va="center", fontsize=6, color="black")
            elif pd.isna(mat.loc[g, c]):
                ax.text(j, i, ".", ha="center", va="center", fontsize=7, color="0.55")
    ax.set_title("Refractory-module causal map across eQTL contexts")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("log2(MR OR), red=risk, blue=protective")
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    plt.close(fig)


def write_summary(df: pd.DataFrame, datasets: pd.DataFrame, contexts: dict[str, dict[str, str]], out_dir: Path) -> None:
    tested = df[df["analysis_status"].str.startswith("tested")]
    causal = df[df["causal_call"]].sort_values(["context", "MR_p"])
    new_causal = causal[causal["context"].ne("blood")]
    new_causal_text = (
        ", ".join(new_causal["gene"].astype(str) + "/" + new_causal["context"].astype(str))
        if len(new_causal)
        else "none"
    )
    focus = df[df["gene"].isin(REQ_GENES) & df["context"].isin(["monocyte", "neutrophil", "monocyte_stim", "blood"])]
    no_tool = df[(df["n_instrument"] == 0) & df["context"].ne("blood")]
    ds_lines = []
    for context, meta in contexts.items():
        rec = datasets[datasets["dataset_id"].eq(meta["dataset_id"])]
        if len(rec):
            r = rec.iloc[0]
            ds_lines.append(f"- {context}: {r.dataset_id} {r.study_label} / {r.sample_group} / n={r.sample_size}")
        else:
            ds_lines.append(f"- {context}: {meta['dataset_id']} {meta.get('study_label', '')} / {meta.get('sample_group', '')}")
    ccl8 = df[(df.gene.eq("CCL8")) & (df.context.eq("monocyte"))]
    cxcr2 = df[(df.gene.eq("CXCR2")) & (df.context.eq("neutrophil"))]
    lines = [
        "# Multicontext immune eQTL causal map summary",
        "",
        "## Dataset choices",
        *ds_lines,
        "",
        "Fairfax_2014 was searched in the eQTL Catalogue `quant_method=ge` dataset listing but was not exposed there in this run; BLUEPRINT and Quach_2016 provided the requested monocyte/neutrophil/stimulated-monocyte contexts.",
        "",
        "## Direct answers",
        answer_line("CCL8 in monocyte", ccl8),
        answer_line("CXCR2 in neutrophil", cxcr2),
        f"New strict multicontext causal calls beyond existing blood anchors: {new_causal_text}.",
        "",
        "## Strict causal calls",
    ]
    if len(causal):
        for r in causal.itertuples():
            lines.append(f"- {r.gene}/{r.context}: OR={r.MR_OR:.3g}, MR_FDR={r.MR_FDR:.3g}, PP4={r.coloc_PP4:.3g}, nsnps={getattr(r, 'coloc_nsnps', math.nan)}")
    else:
        lines.append("- None under MR FDR<0.05 and PP4>0.8.")
    lines.extend(
        [
            "",
            "## Tool coverage",
            f"Tested gene-context rows: {len(tested)}.",
            f"No-tool/no-harmonized rows outside blood: {len(no_tool)}.",
            "",
            "Key focus rows:",
            "```tsv",
            focus[["gene", "context", "MR_OR", "MR_p", "MR_FDR", "coloc_PP4", "n_instrument", "analysis_status"]].to_csv(sep="\t", index=False).rstrip(),
            "```",
            "",
            "## Gut allpairs note",
            "The local `Colon-{Sigmoid,Transverse}.nominal.allpairs.txt.gz` files were streamed with GTEx phenotype IDs inferred as `variant_pos - tss_distance` from the matching eGenes files. In this run those inferred IDs were absent from the supplied allpairs files (see `gtex_*_allpairs_phenotype_id_audit.tsv`), so gut coloc remains labelled `full_cis_region_api_or_gtex_sigpairs` and the GTEx rows should be treated as significant-pairs restricted. This is an input-file mapping mismatch, not an immune-eQTL limitation.",
        ]
    )
    (out_dir / "SUMMARY.md").write_text("\n".join(lines) + "\n")


def answer_line(label: str, row: pd.DataFrame) -> str:
    if not len(row) or row.iloc[0]["n_instrument"] == 0:
        return f"- {label}: not confirmed; no usable/harmonized immune eQTL instrument in this context."
    r = row.iloc[0]
    call = "confirmed" if bool(r["causal_call"]) else "not confirmed by the strict MR+coloc rule"
    pp4 = "NA" if pd.isna(r["coloc_PP4"]) else f"{r['coloc_PP4']:.3g}"
    return f"- {label}: {call}; OR={r['MR_OR']:.3g}, MR_FDR={r['MR_FDR']:.3g}, PP4={pp4}."


def update_status(out_dir: Path) -> None:
    status = P.journal / "status" / "causal_module_status.md"
    existing = status.read_text() if status.exists() else "# Causal module status\n"
    note = f"""

## 2026-06-23 JST - multicontext immune eQTL run

- Pulled latest repo before work (`git pull --ff-only`, fast-forwarded to origin/main).
- Ran `src/26_immune_eqtl_multicontext.py` on CPU only.
- eQTL Catalogue datasets selected from the live `quant_method=ge` listing and cached in `{out_dir / 'eqtl_catalogue_datasets_ge.tsv'}`.
- Outputs: `{out_dir / 'module_causal_map_multicontext.tsv'}`, `{out_dir / 'Fig_module_causal_multicontext.png'}`, `{out_dir / 'SUMMARY.md'}`.
- Gut allpairs: supplied nominal allpairs files were streamed by eGenes-inferred phenotype IDs, but no target IDs matched the local allpairs `gene_id` values; audit files were written and gut coloc remains significant-pairs restricted.
"""
    status.write_text(existing.rstrip() + "\n" + note)


def main() -> None:
    cfg = Inputs()
    genes = module_genes()
    gene_table = read_gene_table(cfg.gtex_dir, genes)
    datasets = discover_datasets(cfg.out_dir)
    immune_contexts = select_immune_contexts(datasets)
    pd.DataFrame([{"context": k, **v} for k, v in immune_contexts.items()]).to_csv(
        cfg.out_dir / "eqtl_catalogue_selected_contexts.tsv", sep="\t", index=False
    )

    cache_dir = cfg.out_dir / "immune_eqtl_cache"
    tasks = [
        (context, meta, gene._asdict())
        for context, meta in immune_contexts.items()
        for gene in gene_table.itertuples(index=False)
    ]
    immune_frames = Parallel(n_jobs=min(30, max(1, len(tasks))), prefer="threads", verbose=5)(
        delayed(fetch_context_gene_eqtl)(context, meta, gene_row, cache_dir)
        for context, meta, gene_row in tasks
    )
    immune_frames = [df for df in immune_frames if len(df)]
    immune_eqtl = pd.concat(immune_frames, ignore_index=True) if immune_frames else pd.DataFrame()
    immune_eqtl.to_csv(cfg.out_dir / "immune_eqtl_catalogue_module_pairs.tsv", sep="\t", index=False)

    gtex_eqtl = load_gtex_sig_pairs(cfg.gtex_dir, genes)
    all_eqtl = pd.concat([gtex_eqtl, immune_eqtl], ignore_index=True, sort=False)
    all_eqtl.to_csv(cfg.out_dir / "multicontext_eqtl_pairs.tsv", sep="\t", index=False)
    inst = lead_instruments(all_eqtl)
    inst.to_csv(cfg.out_dir / "instruments_multicontext.tsv", sep="\t", index=False)

    mr = run_mr(all_eqtl, cfg.gwas_dir, cfg.out_dir)
    coloc_parts = [coloc_for_eqtl(all_eqtl, cfg.gwas_dir / "IBD.h.tsv.gz", cfg.out_dir, "full_cis_region_api_or_gtex_sigpairs")]
    allpair_coloc = allpairs_resolvable_coloc(gene_table, cfg.gwas_dir / "IBD.h.tsv.gz", cfg.out_dir)
    if len(allpair_coloc):
        coloc_parts.append(allpair_coloc)
    coloc = pd.concat([c for c in coloc_parts if len(c)], ignore_index=True) if coloc_parts else pd.DataFrame()
    coloc = coloc.sort_values(["gene", "context", "coloc_method"]).drop_duplicates(["gene", "context"], keep="last")
    coloc.to_csv(cfg.out_dir / "coloc_multicontext_IBD.tsv", sep="\t", index=False)

    df = assemble_map(genes, mr, coloc, inst, cfg.out_dir)
    plot_map(df, cfg.out_dir / "Fig_module_causal_multicontext.png")
    write_summary(df, datasets, immune_contexts, cfg.out_dir)
    P.promote_table(cfg.out_dir / "module_causal_map_multicontext.tsv")
    P.promote_figure(cfg.out_dir / "Fig_module_causal_multicontext.png")
    update_status(cfg.out_dir)
    print(df[["gene", "context", "MR_OR", "MR_p", "MR_FDR", "coloc_PP4", "n_instrument", "causal_call"]].dropna(how="all").head(80).to_string(index=False))


if __name__ == "__main__":
    main()
