# Task 3 Status

Status: local CPU deliverables complete

Owner: task3 worker

Scope:
- Multi-cohort bulk validation and inflammation-adjusted prediction.
- Write new scripts as `src/22_*`, `src/23_*`, `src/24_*`, `src/25_*` or similarly clear task3 names.
- Write downloaded GEO raw files under `data/raw/geo/`, derived files under `data/interim/`, and task3 outputs under `outputs/`.

Progress log:
- 2026-06-22: status file initialized by parent agent.
- 2026-06-22: task3 worker read shared context, task3 plan, paths, and existing GEO/senescence scripts. No task1/task2 files touched.
- 2026-06-22: Lightweight GEO metadata check before download:
  - Include/download: GSE12251 (IFX baseline colonic biopsy, WK8RSPHM response), GSE23597 (IFX baseline W0 biopsy, wk8 response), GSE92415 (PURSUIT golimumab Week 0 colon mucosa, wk6response, Mayo).
  - Include existing: GSE16879, GSE73661.
  - Exclude: GSE14580 duplicate GSM block already in GSE16879; GSE52746 lacks baseline response design; GSE111761 is isolated LPMC/on-treatment responder status, not baseline mucosal bulk prediction.
- 2026-06-22: Added task3-only scripts `src/22_harvest_cohorts.py`, `src/23_score_all.py`, `src/24_loco_meta.py`, `src/25_inflammation_adjust.py`; no shared-script edits.
- 2026-06-22: Ran `src/22_harvest_cohorts.py`; downloaded GSE12251, GSE23597, GSE92415 series matrices and GPL13158 annotation into `data/raw/geo/`. Existing GPL570/GPL6244 and existing matrices reused. Runtime about 3.2 min; slow point was network download, kept sequential/low-concurrency to avoid GEO failures.
- 2026-06-22: Optimized task3 parsing/scoring in `src/23_score_all.py` with limited cohort-level `joblib.Parallel` (`TASK3_N_JOBS`, default max 4). Each worker writes only its own cohort parquet/meta/scored files, then the parent process writes coverage. Runtime with 5 cohorts: total 2.5 s; per cohort 0.5-1.7 s. This avoids serial reparse/score overhead as more GEO cohorts are added.
- 2026-06-22: Ran `src/23_score_all.py`; regenerated existing GSE16879/GSE73661 and created `data/interim/GSE12251_{expr.parquet,meta.tsv,scored.tsv}`, `GSE23597_*`, `GSE92415_*`. Baseline labeled biologic samples after excluding placebo downstream: GSE12251 22, GSE23597 31, GSE92415 59.
- 2026-06-22: Ran `src/24_loco_meta.py`; wrote `outputs/multicohort_auc.tsv` and `outputs/Fig_task3_multicohort.png`. Bootstrap CI/meta was the slowest CPU step (~30 s) but acceptable for current 5 cohorts; no large-table joins used.
  - SenMayo AUCs: GSE16879 0.852, GSE73661 0.741 (reproduces existing scale), GSE12251 0.933, GSE23597 0.762, GSE92415 0.671.
  - Random-effects pooled SenMayo AUC 0.782 (95% CI 0.677-0.860), I2 46.2%; LOCO pooled AUCs remain 0.750-0.814.
- 2026-06-22: Ran `src/25_inflammation_adjust.py`; wrote `outputs/inflammation_adjusted.tsv`. Runtime ~14 s. Pooled cohort-FE SenMayo remains positive after CRP-like adjustment (coef 0.74, p=0.022), neutrophil proxy (coef 0.87, p=0.0026), and Mayo where available (coef 0.75, p=0.0009). Stronger combined gene-inflammation+Mayo adjustment attenuates SenMayo (coef 0.51, p=0.23), so interpretation should be "partly independent, partly overlapping with inflammation."
- 2026-06-22: Promoted current task3 tables because these are the final local CPU validation outputs for this task: `results/tables/multicohort_auc.tsv`, `results/tables/inflammation_adjusted.tsv`, `results/tables/task3_cohort_harvest.tsv`, `results/tables/task3_score_coverage.tsv`.

## 2026-06-23 JST - single-cell deconvolution response

- Ran `src/29_deconv_response.py` CPU-only with marker-restricted NNLS from `outputs/sc_01/*.h5ad` into five baseline bulk response cohorts.
- Wrote `outputs/deconv/{deconv_proportions.tsv,celltype_fraction_vs_response.tsv,targeted_scores.tsv,Fig_deconv_response.png,SUMMARY.md}` and promoted final table/figure copies into `results/`.
- Key random-effects response associations:
  - neutrophil_marker_score: OR/SD=3.16, p=0.000192
  - refractory_module_score: OR/SD=3.31, p=0.000137
  - senmayo: OR/SD=3.22, p=4.2e-05
  - Fibroblast: OR/SD=1.58, p=0.115
  - Myeloid: OR/SD=2.65, p=4.83e-07
- Caveat: local scRNA references do not contain an explicit neutrophil cluster, so neutrophil evidence is marker/CXCR2 based rather than an NNLS neutrophil fraction.
