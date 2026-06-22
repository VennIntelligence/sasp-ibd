"""FinnGen R12 replication for CCL8/CXCR2 cis-MR and coloc.

The script uses the public R12 manifest to resolve exact endpoint URLs. It does
not auto-download endpoint GWAS files by default; missing files are recorded in
the output with their verified URL so downloads remain deliberate.
"""
from __future__ import annotations

import argparse
import gzip
import shlex
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict
from scipy import stats
from scipy.special import logsumexp

from paths import P


class Config(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    out_dir: Path = P.out("mr")
    manifest: Path = P.raw / "gwas" / "finngen_R12_manifest.tsv"
    manifest_url: str = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv"
    cis_candidates: Path = P.raw / "eqtlgen" / "cis_full_candidates.tsv"
    eqtl_af: Path = P.raw / "eqtlgen" / "snp_af.txt.gz"
    instruments: Path = P.out("mr") / "instruments.tsv"
    endpoints: dict[str, str] = {
        "IBD": "K11_IBD_STRICT",
        "CD": "K11_CD_STRICT2",
        "UC": "K11_UC_STRICT2",
    }
    genes: tuple[str, ...] = ("CCL8", "CXCR2")
    w_eqtl: float = 0.15**2
    w_gwas: float = 0.2**2
    p1: float = 1e-4
    p2: float = 1e-4
    p12: float = 1e-5


CFG = Config()


def run_text(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def download(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url} -> {path}", flush=True)
    cmd = ["curl", "-L", "--fail", "--continue-at", "-", "--progress-bar", "--output", str(path), url]
    subprocess.run(cmd, check=True)


def endpoint_path(phenocode: str) -> Path:
    return P.raw / "gwas" / f"finngen_R12_{phenocode}.gz"


def manifest_endpoints(download_missing: bool = False) -> pd.DataFrame:
    if not CFG.manifest.exists():
        if download_missing:
            download(CFG.manifest_url, CFG.manifest)
        else:
            raise FileNotFoundError(f"missing FinnGen manifest: {CFG.manifest_url}")
    manifest = pd.read_csv(CFG.manifest, sep="\t")
    rows = []
    for label, phenocode in CFG.endpoints.items():
        hit = manifest[manifest["phenocode"] == phenocode]
        if hit.empty:
            rows.append({"outcome": label, "phenocode": phenocode, "status": "phenocode not found in R12 manifest"})
            continue
        r = hit.iloc[0].to_dict()
        rows.append(
            {
                "outcome": label,
                "phenocode": phenocode,
                "phenotype": r["phenotype"],
                "num_cases": r["num_cases"],
                "num_controls": r["num_controls"],
                "path_https": r["path_https"],
                "local_path": endpoint_path(phenocode),
                "status": "ok",
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(CFG.out_dir / "finngen_R12_task1_endpoints.tsv", sep="\t", index=False)
    return out


def header(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").lstrip("#").split("\t")


def pick(cols: list[str], names: list[str]) -> str:
    lower = {c.lower(): c for c in cols}
    for name in names:
        if name.lower() in lower:
            return lower[name.lower()]
    raise ValueError(f"none of {names} found in columns: {cols[:20]}")


def finngen_lookup(path: Path, rsids: set[str]) -> pd.DataFrame:
    cols = header(path)
    idx = {c: i + 1 for i, c in enumerate(cols)}
    rs_col = pick(cols, ["rsids", "rsid"])
    beta_col = pick(cols, ["beta"])
    se_col = pick(cols, ["sebeta", "se", "standard_error"])
    p_col = pick(cols, ["pval", "p_value", "p"])
    ea_col = pick(cols, ["alt", "effect_allele"])
    oa_col = pick(cols, ["ref", "other_allele"])
    rsfile = CFG.out_dir / "_finngen_lookup_rsids.txt"
    rsfile.write_text("\n".join(sorted(rsids)) + "\n")
    cmd = (
        f"gunzip -c {shlex.quote(str(path))} | "
        "awk -F'\\t' -v OFS='\\t' "
        f"-v rs={idx[rs_col]} -v b={idx[beta_col]} -v se={idx[se_col]} -v pv={idx[p_col]} "
        f"-v ea={idx[ea_col]} -v oa={idx[oa_col]} "
        "'NR==FNR{want[$1]=1; next} FNR==1{next} "
        "{n=split($rs,a,\",\"); for(i=1;i<=n;i++){if(a[i] in want){print a[i],$b,$se,$pv,$ea,$oa; break}}}' "
        f"{shlex.quote(str(rsfile))} -"
    )
    rows = []
    for line in run_text(cmd).splitlines():
        rsid, beta, se, p, ea, oa = line.split("\t")
        try:
            rows.append({"rsid": rsid, "beta": float(beta), "se": float(se), "p": float(p), "ea": ea.upper(), "oa": oa.upper()})
        except ValueError:
            continue
    return pd.DataFrame(rows)


def align_beta(beta: float, ea: str, oa: str, assessed: str, other: str) -> float | None:
    if ea == assessed and oa == other:
        return beta
    if ea == other and oa == assessed:
        return -beta
    return None


def wald_rows(outcome: str, phenocode: str, path: Path, url: str) -> list[dict]:
    inst = pd.read_csv(CFG.instruments, sep="\t")
    inst = inst[inst["gene"].isin(CFG.genes)].copy()
    fg = finngen_lookup(path, set(inst["rsid"]))
    rows = []
    for r in inst.merge(fg, on="rsid", how="left").itertuples(index=False):
        if pd.isna(r.beta):
            rows.append(base_row(outcome, phenocode, r.gene, "wald", "lead rsid absent from FinnGen endpoint", url))
            continue
        bg = align_beta(r.beta, r.ea, r.oa, r.assessed, r.other)
        if bg is None:
            rows.append(base_row(outcome, phenocode, r.gene, "wald", "allele mismatch", url))
            continue
        theta = bg / r.beta_eqtl
        se = abs(r.se / r.beta_eqtl)
        p = float(2 * stats.norm.sf(abs(theta / se)))
        rows.append(
            {
                "outcome": outcome,
                "phenocode": phenocode,
                "gene": r.gene,
                "method": "wald",
                "rsid": r.rsid,
                "theta": theta,
                "se": se,
                "OR": float(np.exp(theta)),
                "p": p,
                "PP4": np.nan,
                "nsnps": 1,
                "endpoint_url": url,
                "status": "ok",
            }
        )
    return rows


def base_row(outcome: str, phenocode: str, gene: str, method: str, status: str, url: str) -> dict:
    return {
        "outcome": outcome,
        "phenocode": phenocode,
        "gene": gene,
        "method": method,
        "rsid": "",
        "theta": np.nan,
        "se": np.nan,
        "OR": np.nan,
        "p": np.nan,
        "PP4": np.nan,
        "nsnps": 0,
        "endpoint_url": url,
        "status": status,
    }


def lookup_af(rsids: set[str]) -> dict[str, tuple[str, float]]:
    rsfile = CFG.out_dir / "_finngen_coloc_rsids.txt"
    rsfile.write_text("\n".join(sorted(rsids)) + "\n")
    cmd = (
        f"gunzip -c {shlex.quote(str(CFG.eqtl_af))} | "
        "awk -F'\\t' -v OFS='\\t' "
        "'NR==FNR{want[$1]; next} FNR==1{next} ($1 in want){print $1,$5,$9}' "
        f"{shlex.quote(str(rsfile))} -"
    )
    af = {}
    for line in run_text(cmd).splitlines():
        rsid, allele_b, freq_b = line.split("\t")
        try:
            af[rsid] = (allele_b.upper(), float(freq_b))
        except ValueError:
            continue
    return af


def labf(z: float, v: float, w: float) -> float:
    r = w / (v + w)
    return 0.5 * (np.log(1 - r) + r * z * z)


def coloc_rows(outcome: str, phenocode: str, path: Path, url: str) -> list[dict]:
    cols = ["Pvalue", "SNP", "Zscore", "AssessedAllele", "OtherAllele", "GeneSymbol", "NrSamples"]
    cis = pd.read_csv(CFG.cis_candidates, sep="\t", usecols=cols)
    cis = cis[cis["GeneSymbol"].isin(CFG.genes)].rename(
        columns={
            "Pvalue": "p_eqtl",
            "SNP": "rsid",
            "Zscore": "z",
            "AssessedAllele": "assessed",
            "OtherAllele": "other",
            "GeneSymbol": "gene",
            "NrSamples": "n",
        }
    )
    for c in ["p_eqtl", "z", "n"]:
        cis[c] = pd.to_numeric(cis[c], errors="coerce")
    cis = cis.dropna(subset=["p_eqtl", "z", "n"])
    af = lookup_af(set(cis["rsid"]))
    fg = finngen_lookup(path, set(cis["rsid"]))
    rows = []
    for gene, g in cis.groupby("gene"):
        d = g.merge(fg, on="rsid", how="inner")
        l1, l2 = [], []
        for r in d.itertuples(index=False):
            if r.rsid not in af or r.se <= 0:
                continue
            allele_b, freq_b = af[r.rsid]
            assessed, other = r.assessed.upper(), r.other.upper()
            if allele_b == assessed:
                f = freq_b
            elif allele_b == other:
                f = 1 - freq_b
            else:
                continue
            if f <= 0 or f >= 1:
                continue
            v_e = 1 / (2 * r.n * f * (1 - f))
            z_g = r.beta / r.se
            l1.append(labf(r.z, v_e, CFG.w_eqtl))
            l2.append(labf(z_g, r.se**2, CFG.w_gwas))
        if len(l1) < 5:
            rows.append(base_row(outcome, phenocode, gene, "coloc", "fewer than 5 overlapping coloc SNPs", url))
            continue
        l1, l2 = np.array(l1), np.array(l2)
        h1, h2, h4 = logsumexp(l1), logsumexp(l2), logsumexp(l1 + l2)
        s = h1 + h2
        h3 = s + np.log1p(-np.exp(min(h4 - s, -1e-12)))
        logs = np.array([0, np.log(CFG.p1) + h1, np.log(CFG.p2) + h2, np.log(CFG.p1) + np.log(CFG.p2) + h3, np.log(CFG.p12) + h4])
        pp = np.exp(logs - logsumexp(logs))
        rows.append(
            {
                "outcome": outcome,
                "phenocode": phenocode,
                "gene": gene,
                "method": "coloc",
                "rsid": "",
                "theta": np.nan,
                "se": np.nan,
                "OR": np.nan,
                "p": np.nan,
                "PP4": float(pp[4]),
                "nsnps": len(l1),
                "endpoint_url": url,
                "status": "ok",
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true", help="download missing FinnGen endpoint files before analysis")
    args = parser.parse_args()

    CFG.out_dir.mkdir(parents=True, exist_ok=True)
    try:
        endpoints = manifest_endpoints(download_missing=args.download)
    except FileNotFoundError as exc:
        rows = [
            base_row(outcome, phenocode, gene, method, str(exc), CFG.manifest_url)
            for outcome, phenocode in CFG.endpoints.items()
            for gene in CFG.genes
            for method in ["wald", "coloc"]
        ]
        pd.DataFrame(rows).to_csv(CFG.out_dir / "replication_finngen.tsv", sep="\t", index=False)
        print(exc)
        return

    rows = []
    for r in endpoints.itertuples(index=False):
        if r.status != "ok":
            rows.extend(base_row(r.outcome, r.phenocode, gene, method, r.status, "") for gene in CFG.genes for method in ["wald", "coloc"])
            continue
        path = Path(r.local_path)
        if not path.exists() and args.download:
            download(r.path_https, path)
        if not path.exists():
            status = f"missing local FinnGen endpoint; download {r.path_https} to {path}"
            rows.extend(base_row(r.outcome, r.phenocode, gene, method, status, r.path_https) for gene in CFG.genes for method in ["wald", "coloc"])
            continue
        rows.extend(wald_rows(r.outcome, r.phenocode, path, r.path_https))
        rows.extend(coloc_rows(r.outcome, r.phenocode, path, r.path_https))

    res = pd.DataFrame(rows)
    res.to_csv(CFG.out_dir / "replication_finngen.tsv", sep="\t", index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
