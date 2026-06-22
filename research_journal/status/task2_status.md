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
## task_gpu_fix stage0
- `git pull --ff-only`: already up to date.
- Read task plan, AGENTS, `src/sc_03_foundation_perturb.py`, and `outputs/sc_04/task2_honest_summary.md`.
- Downloaded `Geneformer-V2-104M` weights from the official `ctheodoris/Geneformer` Hugging Face repo into `models/Geneformer/Geneformer-V2-104M`.
- D reconnaissance found GSE282122 as a public anti-TNF longitudinal IBD single-cell atlas with remission outcomes; processed filtered archive is 2.8 GB, so it is a viable follow-on dataset after A/B/C hardening.

## sc_05 benchmark note
- Geneformer V2-104M with batch_size=96 on 512 Martin cells OOMed on a 24 GB RTX 3090 (process used about 17 GB and needed another 6 GB).
- Formal schedule will use V2 one job per GPU with smaller batch, and V1 jobs separately; this is a measured memory adaptation rather than the optimistic 4-slot plan.

## scGPT environment attempt
- Created isolated `.venv-scgpt` and installed `scgpt==0.2.4` without touching the main `.venv`.
- Added missing `IPython`, but `import scgpt` fails because `torchtext/libtorchtext.so` has an undefined symbol against the installed torch ABI.
- Per task fallback, scGPT is not used in the main matrix; cross-model evidence uses Geneformer V1-10M and V2-104M.

## sc_05 GPU0 failure
- During the V2-104M wave, physical GPU0 disappeared from `nvidia-smi` with `Unable to determine the device handle ... Unknown Error`.
- The `gf_v2_104m_smillie` job failed at 6/205 perturb genes with CUDA `unspecified launch failure`; no completed effect table was written.
- Physical GPU1 remained healthy and continued `gf_v2_104m_martin` at ~99% utilization.
- Recovery plan: do not reboot/reset from the script; allow GPU1 job to finish, then continue remaining jobs on GPU1 with conservative batch size.

## sc_05 CUDA unavailable after GPU failure
- After both V2 jobs exited, `nvidia-smi` still reported GPU0 as `Unable to determine the device handle ... Unknown Error`; GPU1 appeared idle.
- A fresh PyTorch process with `CUDA_VISIBLE_DEVICES=1` could not initialize CUDA (`cuda=False`, `CUDA unknown error`), so no remaining GPU work can proceed safely without driver/device reset or reboot.
- No completed `*_effects.tsv` tables were written for the V2 wave; only null panel files exist. The next executable step after hardware recovery is to rerun remaining jobs conservatively on a healthy GPU.

## sc_05 perturbation hardening
- combined 4 model x dataset jobs from `outputs/sc_05_perturb_stats/jobs/`.
- outputs: `outputs/sc_05_perturb_stats/perturb_stats.tsv`, `perturb_consensus.tsv`, `classifier_metrics.tsv`, `Fig_task2_perturb_stats.png`.
- promoted final table/figure copies to `results/`; summary: `outputs/sc_05_perturb_stats/task2_perturb_stats_summary.md`.

