"""Task 3 GEO cohort harvest.

Lightweight GEO verification first, then download only cohorts with baseline
mucosal bulk expression, response labels, and usable platform annotation.
"""
from __future__ import annotations

import gzip
import shutil
import tempfile
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd
from pydantic import BaseModel, ConfigDict

from paths import P


RAW_GEO = P.raw / "geo"
OUT = P.out("22_harvest_cohorts")


class CohortPlan(BaseModel):
    model_config = ConfigDict(frozen=True)

    accession: str
    platform: str | None
    drug: str
    decision: str
    reason: str
    label_field: str
    baseline_field: str
    download: bool = False
    duplicate_of: str = ""


COHORTS = [
    CohortPlan(
        accession="GSE16879",
        platform="GPL570",
        drug="IFX",
        decision="include_existing",
        reason="existing mucosal IBD IFX pre/post cohort with response labels",
        label_field="response to infliximab",
        baseline_field="before or after first infliximab treatment",
    ),
    CohortPlan(
        accession="GSE73661",
        platform="GPL6244",
        drug="VDZ/IFX",
        decision="include_existing",
        reason="existing mucosal UC biologic cohort with W0 baseline, follow-up Mayo, and response derivable from mucosal healing",
        label_field="mayo endoscopic subscore at post-treatment visits",
        baseline_field="week (w)=W0",
    ),
    CohortPlan(
        accession="GSE12251",
        platform="GPL570",
        drug="IFX",
        decision="include_download",
        reason="baseline colonic biopsy; week-8 endoscopic/histologic response label",
        label_field="WK8RSPHM",
        baseline_field="title contains W0",
        download=True,
    ),
    CohortPlan(
        accession="GSE23597",
        platform="GPL570",
        drug="IFX",
        decision="include_download",
        reason="baseline colonic biopsy; infliximab dose, time, and wk8/wk30 response fields available",
        label_field="wk8 response",
        baseline_field="time=W0",
        download=True,
    ),
    CohortPlan(
        accession="GSE92415",
        platform="GPL13158",
        drug="GLM",
        decision="include_download",
        reason="PURSUIT-SC colon mucosa; Week 0 baseline, golimumab/placebo treatment field, wk6 response, Mayo score",
        label_field="wk6response",
        baseline_field="visit=Week 0",
        download=True,
    ),
    CohortPlan(
        accession="GSE14580",
        platform="GPL570",
        drug="IFX",
        decision="exclude_duplicate",
        reason="same GSM accession block as the UC baseline subset already contained in GSE16879",
        label_field="response to infliximab",
        baseline_field="before first infliximab treatment",
        duplicate_of="GSE16879",
    ),
    CohortPlan(
        accession="GSE52746",
        platform="GPL17996",
        drug="anti-TNF",
        decision="exclude_no_baseline_response",
        reason="colon bulk data exist, but metadata distinguish active/inactive under anti-TNF rather than treatment-naive baseline response",
        label_field="biopsy",
        baseline_field="not available",
    ),
    CohortPlan(
        accession="GSE111761",
        platform="GPL13497",
        drug="anti-TNF",
        decision="exclude_not_mucosal_baseline_bulk",
        reason="isolated LPMC expression from ongoing anti-TNF responders/non-responders; no baseline mucosal biopsy prediction setting",
        label_field="patient",
        baseline_field="not available",
    ),
]


@dataclass
class GeoHeader:
    accession: str
    title: str
    platforms: str
    n_samples: int
    matrix_bytes: int | None
    char_fields: str
    matrix_url: str
    geo_url: str


def geo_bucket(accession: str) -> str:
    return f"{accession[:-3]}nnn"


def matrix_url(accession: str) -> str:
    bucket = geo_bucket(accession)
    return (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{bucket}/{accession}/matrix/{accession}_series_matrix.txt.gz"
    )


def platform_url(gpl: str) -> str:
    bucket = geo_bucket(gpl)
    return (
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
        f"{bucket}/{gpl}/annot/{gpl}.annot.gz"
    )


def content_length(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            val = resp.headers.get("Content-Length")
            return int(val) if val else None
    except Exception:
        return None


def read_geo_header(accession: str) -> GeoHeader:
    url = matrix_url(accession)
    gsm: list[str] = []
    title = ""
    platforms: list[str] = []
    char_fields: list[str] = []
    with urllib.request.urlopen(url, timeout=120) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            for raw in gz:
                line = raw.decode("utf-8", "ignore").rstrip("\n")
                if line.startswith("!series_matrix_table_begin"):
                    break
                parts = [x.strip().strip('"') for x in line.split("\t")]
                tag, vals = parts[0], parts[1:]
                if tag == "!Series_title" and vals:
                    title = vals[0]
                elif tag == "!Series_platform_id":
                    platforms.extend(vals)
                elif tag == "!Sample_geo_accession":
                    gsm = vals
                elif tag == "!Sample_characteristics_ch1" and vals:
                    field = vals[0].split(":", 1)[0].strip() if ":" in vals[0] else "char"
                    char_fields.append(field)
    return GeoHeader(
        accession=accession,
        title=title,
        platforms=";".join(sorted(set(platforms))),
        n_samples=len(gsm),
        matrix_bytes=content_length(url),
        char_fields=";".join(dict.fromkeys(char_fields)),
        matrix_url=url,
        geo_url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
    )


def download(url: str, dest: Path) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return False
    with tempfile.NamedTemporaryFile(dir=dest.parent, delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        with urllib.request.urlopen(url, timeout=300) as resp, tmp_path.open("wb") as fh:
            shutil.copyfileobj(resp, fh, length=1024 * 1024)
        tmp_path.replace(dest)
        return True
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def main() -> None:
    RAW_GEO.mkdir(parents=True, exist_ok=True)
    rows = []
    for plan in COHORTS:
        header = read_geo_header(plan.accession)
        matrix_path = RAW_GEO / f"{plan.accession}_series_matrix.txt.gz"
        platform_path = RAW_GEO / f"{plan.platform}.annot.gz" if plan.platform else None
        matrix_downloaded = False
        platform_downloaded = False
        if plan.download:
            matrix_downloaded = download(header.matrix_url, matrix_path)
            if plan.platform and platform_path is not None:
                platform_downloaded = download(platform_url(plan.platform), platform_path)
        rows.append(
            {
                **plan.model_dump(),
                **asdict(header),
                "matrix_path": str(matrix_path) if matrix_path.exists() else "",
                "platform_path": str(platform_path) if platform_path and platform_path.exists() else "",
                "matrix_downloaded": matrix_downloaded,
                "platform_downloaded": platform_downloaded,
            }
        )

    out = pd.DataFrame(rows)
    out_file = OUT / "task3_cohort_harvest.tsv"
    root_file = P.outputs / "task3_cohort_harvest.tsv"
    out.to_csv(out_file, sep="\t", index=False)
    out.to_csv(root_file, sep="\t", index=False)
    print(out[["accession", "decision", "n_samples", "platforms", "matrix_downloaded", "platform_downloaded"]].to_string(index=False))
    print(f"\nwrote {out_file}")
    print(f"wrote {root_file}")


if __name__ == "__main__":
    main()
