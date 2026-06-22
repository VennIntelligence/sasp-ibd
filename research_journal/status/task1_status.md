# Task 1 Status

Status: FinnGen replication complete; LTL external large-file blocker remains

Owner: task1 worker

Scope:
- Genetic causal separation: MVMR/reverse MR/LTL MR/FinnGen replication where feasible.
- Write new scripts as `src/18_*`, `src/19_*`, `src/20_*`, `src/21_*` or similarly clear task1 names.
- Write task1 outputs under `outputs/mr/` first; promote curated tables to `results/tables/` only when final.

Progress log:
- 2026-06-22: status file initialized by parent agent.
- 2026-06-22: Re-ran `src/14_mr.py` with the optimized awk hash-join path. Reproduced current MR counts: IBD 77 genes tested / 9 FDR<0.05, CD 77 / 7, UC 77 / 5. CCL8 and CXCR2 remain the two high-PP4 coloc candidates; TNFRSF1A remains IBD risk-direction positive control.
- 2026-06-22: External task1 data sources checked before download:
  - LTL exposure: Codd 2021 leukocyte telomere length, GWAS Catalog `GCST90002398`, harmonised file `32888494-GCST90002398-EFO_0004833.h.tsv.gz`, expected size ~1.4G from GWAS Catalog FTP. Download status: deferred after CRP took 21m46s at ~212KB/s; expected LTL download at this throughput is >100 min. Script/output blocker is in place.
  - CRP exposure: Said 2022 C-reactive protein, GWAS Catalog `GCST90029070`, harmonised file `35459240-GCST90029070-EFO_0004458.h.tsv.gz`, expected size ~272M from GWAS Catalog FTP. Download status: complete as `data/raw/gwas/CRP_GCST90029070.h.tsv.gz`; gzip test passed.
  - FinnGen replication: task plan requested R12; official R12 manifest is `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv`. Download status: manifest complete; endpoint summary stats not downloaded after URL/size check (~2.4G total).
- 2026-06-22: Performance audit for task1 scripts: `rg "grep -f|grep"` over `src/13_build_instruments.py`, `src/14_mr.py`, `src/15_coloc.py`, `src/16_triangulate.py`, and new `src/18-21_*` found no remaining grep-based large-file lookups. Existing `src/14_mr.py` already uses awk hash join; new task1 scripts use awk hash joins for rsid lookups and chunked pandas only for p-threshold scans.
- 2026-06-22: Added task1 scripts:
  - `src/18_mvmr.py`: local two-exposure MVMR for target-gene cis-eQTL effects adjusted for CRP; writes explicit blocker/status rows if CRP is absent, partial, or local instruments are rank-insufficient.
  - `src/19_reverse_mr.py`: reverse MR using disease genome-wide significant distance-clumped instruments and eQTLGen cis outcomes; eQTL outcome lookup is a single awk hash scan over target rsids/genes.
  - `src/20_ltl_mr.py`: Codd 2021 LTL -> IBD/CD/UC with IVW, MR-Egger, weighted median, and leave-one-out; writes blocker rows if LTL data is missing/partial.
  - `src/21_replicate_finngen.py`: FinnGen R12 manifest-driven CCL8/CXCR2 replication script. R12 endpoint names confirmed as `K11_IBD_STRICT`, `K11_CD_STRICT2`, `K11_UC_STRICT2`; full endpoint GWAS downloads are not automatic unless `--download` is passed.
- 2026-06-22: Ran `src/21_replicate_finngen.py` without endpoint downloads. It wrote `outputs/mr/finngen_R12_task1_endpoints.tsv` and `outputs/mr/replication_finngen.tsv` with explicit missing-file statuses. HEAD checks show R12 endpoint sizes are ~814M (`K11_IBD_STRICT`), ~808M (`K11_CD_STRICT2`), and ~813M (`K11_UC_STRICT2`), ~2.4G total; not downloaded yet because current network throughput is slow and CRP/LTL are higher-priority task1 inputs.
- 2026-06-22: Ran `src/20_ltl_mr.py` before LTL download; it wrote `outputs/mr/mr_LTL.tsv` with explicit blocker rows pointing to the verified GCST90002398 harmonised URL and expected local path.
- 2026-06-22: Ran `src/18_mvmr.py` after CRP download. Initial 250kb distance clump left CXCR2 with only 2 instruments; changed only the new task1 script default to 100kb, giving CXCR2 5 strong cis-eQTL instruments while CCL8 remains non-estimable with 1 instrument. Output: `outputs/mr/mvmr_results.tsv` and `outputs/mr/mvmr_instruments.tsv`. CXCR2 expression remains protective after CRP adjustment for IBD/CD/UC (IBD OR~0.749, p~1.1e-10), but model status flags weak CRP conditional F (~2.09), so this is supportive/exploratory rather than strong MVMR validation.
- 2026-06-22: Ran `src/19_reverse_mr.py` with awk eQTLGen lookup. Runtime ~327s. Outputs: `outputs/mr/reverse_mr.tsv`, `outputs/mr/reverse_mr_disease_instruments.tsv`, `outputs/mr/reverse_mr_eqtl_outcomes.tsv`, `outputs/mr/reverse_mr_instruments.tsv`. Disease instruments after distance clump: IBD 147, CD 103, UC 69. Only 13 eQTL outcome rows overlap after AF/allele alignment, so most reverse estimates are single-SNP. CCL8 has reverse overlap for IBD/CD; CXCR2 has strong reverse overlap for IBD only. Interpret as possible bidirectional/local-locus signal, not definitive reverse causality.
- 2026-06-22: Final syntax check passed for `src/18_mvmr.py`, `src/19_reverse_mr.py`, `src/20_ltl_mr.py`, and `src/21_replicate_finngen.py`.
- 2026-06-22: FinnGen R12 task1 replication was checked in `tmux lzy:4` and local files. The three endpoints completed download and passed `gzip -t`: `data/raw/gwas/finngen_R12_K11_IBD_STRICT.gz` (776M), `data/raw/gwas/finngen_R12_K11_CD_STRICT2.gz` (771M), and `data/raw/gwas/finngen_R12_K11_UC_STRICT2.gz` (775M). The existing `outputs/mr/replication_finngen.tsv` has all CCL8/CXCR2 Wald and coloc rows with `status=ok`, so no re-download/re-run was needed.
- 2026-06-22: Wrote `outputs/mr/finngen_replication_summary.tsv` comparing FinnGen R12 against de Lange discovery MR (`outputs/mr/mr_{IBD,CD,UC}.tsv`), IBD coloc (`outputs/mr/coloc_IBD.tsv`), and `results/tables/triangulation.tsv`. CCL8 strongly replicates for IBD (de Lange OR 2.35/FDR 0.030/PP4 0.955; FinnGen OR 4.37, p=2.9e-6, PP4 0.950). FinnGen UC also shows a strong CCL8 risk/shared-causal signal (OR 6.28, p=1.9e-6, PP4 0.952), but de Lange UC MR was not FDR-significant, so this is FinnGen-positive extension rather than formal UC replication. CD CCL8 does not replicate (OR 0.87, p=0.84, PP4 0.016).
- 2026-06-22: CXCR2 Wald direction replicates for IBD and UC (FinnGen IBD OR 0.80, p=1.3e-5; UC OR 0.74, p=1.6e-6), but FinnGen coloc does not support shared causal variation at PP4>0.75 (IBD PP4 0.619; UC PP4 0.438). CD CXCR2 is near-null (OR 0.97, p=0.77, PP4 0.010). Lead SNPs were present and allele-aligned for all Wald rows; non-replication is not a missing-SNP or allele-mismatch artifact.
