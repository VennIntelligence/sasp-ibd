# Task2 single-cell foundation perturbation status

Date: 2026-06-22

## Runs completed

- `src/sc_01_qc.py --dataset smillie --scrublet-max-cells 5000 --n-jobs 24`
  - Output: `outputs/sc_01/smillie_uc_qc.h5ad`
  - Retained 320,199 Smillie UC cells and 21,784 genes after basic QC and mitochondrial filtering.
  - Runtime was about 4.5 minutes after adding stage logs/tqdm. The main cost was UMAP; kNN and Leiden were not bottlenecks.
- `src/sc_01_qc.py --dataset martin --scrublet-max-cells 5000 --n-jobs 24`
  - Output: `outputs/sc_01/martin_cd_qc.h5ad`
  - Retained 103,849 Martin CD ileal cells and 33,694 genes.
  - Martin Health mapping follows GEO GSE134809 sample titles: Ileal Involved=Inflamed, Ileal Uninvolved=Non-inflamed. PBMC samples were excluded.
- `src/sc_02_senescence.py`
  - Outputs: `outputs/sc_02/senescence_per_celltype.tsv`, `outputs/sc_02/senescence_by_health_celltype.tsv`, per-cell score TSVs and score figures.
  - Gene coverage: Smillie SenMayo 118 genes; Martin SenMayo 123 genes; arrest/proliferation marker coverage complete.
  - Strict candidate definition: SASP z>=0.5, arrest z>=0.5, positive absolute SASP/arrest scores, non-positive proliferation score, and no MKI67 expression.
- `src/sc_03_foundation_perturb.py --max-cells-per-dataset 12000 --batch-size 128`
  - Outputs: `outputs/sc_03/insilico_perturbation.tsv`, `outputs/sc_03/insilico_perturbation_celllevel.tsv`, `outputs/sc_03/insilico_perturbation.png`.
  - Model: `models/Geneformer/Geneformer-V1-10M`, fp16 on GPU0. GPU utilization was ~98%, memory ~10 GB.
  - Metric: positive `projection_toward_inflamed` means movement toward the inflamed centroid.
- `src/sc_04_cellcomm.py`
  - Outputs: `outputs/sc_04/cellcomm_SASP_axis.tsv`, `outputs/sc_04/cellcomm_SASP_axis.png`, `outputs/sc_04/Fig_task2_singlecell.png`, `outputs/sc_04/task2_honest_summary.md`.
  - Promoted final artifacts to `results/figures/Fig_task2_singlecell.png` and result tables.

## Main interpretation

- Strict bona fide senescence candidates are present but not dominant. The strongest cell-type-level fractions are about 8-11%, not a wholesale senescent conversion of mucosa.
- Smillie UC points most strongly to non-proliferative epithelial/stromal-like compartments: M cells, Goblet cells, immature enterocytes, Best4+ enterocytes and WNT2B/WNT5B fibroblast-like groups. Inflamed epithelial groups show higher candidate fractions than their global averages.
- Martin CD broadly replicates the direction in inflamed fibroblast/endothelial/myeloid/epithelial compartments, but Martin cell types are broad marker-derived labels because the local tar has raw 10x matrices and no author per-cell annotation.
- Geneformer perturbation is strongest for the CXCR2 protective interpretation: CXCR2 deletion moves toward inflammation and CXCR2 overexpression moves away. This matches MR/coloc direction and supports caution about naive CXCR2 blockade.
- CCL8 overexpression moves toward inflammation, supporting risk biology, but CCL8 deletion is near-null overall and differs by dataset. Treat this as partial support, not a clean knockout rescue.
- Focused ligand-receptor scoring is dominated by IL1B/TNF inflammatory SASP axes. CCL8->CCR and CXCL chemokine->CXCR2 axes are present and auditable but are not uniformly the top interactions.

## Caveats

- The single-cell datasets do not include biologic response labels, so response triangulation uses existing bulk-response results (`results/tables/multicohort_auc.tsv`) rather than single-cell response modeling.
- `sc_04` intentionally uses a focused curated SASP ligand-receptor score instead of installing CellPhoneDB into the main environment. This preserves the required environment isolation principle and keeps the CCL8/CXCR2 axis auditable, but it is not a full CellPhoneDB database-wide permutation analysis.
