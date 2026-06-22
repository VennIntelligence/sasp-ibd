# TNP Shared Context

Updated: 2026-06-22

Purpose: shared read-only context for task-specific workers. Do not use this file for frequent progress updates; each task writes its own status file to avoid edit conflicts.

Project state:
- Current repo root: `/Users/ujs/Downloads/lzy`.
- Code lives in `src/`; the directory has been reorganized from the older `scripts/` layout.
- Raw data lives in `data/raw/`; derived matrices and scored tables live in `data/interim/`.
- Machine outputs live in `outputs/`; curated paper-facing tables and figures live in `results/`.
- Path registry exists at `src/paths.py`, but not all scripts have been converted to it yet.

Existing high-signal results:
- Bulk SASP/SenMayo is elevated in active IBD and predicts non-response in GSE16879 and GSE73661.
- Baseline SenMayo AUC: about 0.85 in GSE16879 and 0.74 in GSE73661.
- GSE73661 baseline Mayo severity is weakly correlated with SenMayo, and SenMayo adds predictive value over severity.
- MR plus coloc leaves CCL8 as risk and CXCR2 as protective shared-causal candidates; many other MR hits look like LD artifacts.

Current data footprint:
- `data/raw/eqtlgen/cis_full.txt.gz`: about 4.3 GB.
- `data/raw/eqtlgen/cis_sig.txt.gz`: about 308 MB.
- `data/raw/eqtlgen/snp_af.txt.gz`: about 229 MB.
- `data/raw/gwas/{IBD,CD,UC}.h.tsv.gz`: about 283-288 MB each.
- Existing GEO bulk files are small: GSE16879 about 13 MB, GSE73661 about 31 MB.

Efficiency rules:
- For large GWAS/eQTLGen lookups, use awk hash joins. Do not use `grep -f` on large compressed files.
- Use rsid joins between eQTLGen and GWAS because builds differ.
- Keep task1 and task3 writes separate where possible.
- `src/14_mr.py` was updated on 2026-06-22 to replace `grep -f` with an awk hash join. Validation run completed in about 61.5 seconds for IBD/CD/UC.

Task status files:
- Task 1 writes progress to `research_journal/status/task1_status.md`.
- Task 3 writes progress to `research_journal/status/task3_status.md`.
