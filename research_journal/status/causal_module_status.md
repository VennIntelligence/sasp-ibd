2026-06-23T00:28:17+09:00
# Causal module gut eQTL MR+coloc status

- Started after git pull on 2026-06-23 JST.
- Confirmed task requires CPU only; no GPU calls used.
- Local data/raw was absent; downloading de Lange GWAS and GTEx colon eQTL from public sources.
- GTEx individual tissue URLs from brief returned 404; verified public object is adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar. all-associations direct anonymous access is unavailable so coloc may need a significant-pairs approximation unless an allpairs source becomes accessible.
- Downloaded and gzip-validated de Lange IBD/CD/UC GWAS.
- Extracted GTEx v8 Colon_Sigmoid/Colon_Transverse eGenes and significant variant-gene pairs from adult-gtex GTEx_Analysis_v8_eQTL.tar.
- Pitfall: a GTEx editing-qtl/all_associations download appeared under data/raw/gtex_colon; its columns were not eQTL gene associations, so it was stopped, removed, and not used.
- Full GTEx eQTL allpairs remained unavailable via anonymous direct eQTL paths; coloc is labelled significant-pairs restricted.
- Built GTEx colon instruments: 11 gene-tissue instruments covering 10 genes; 28/38 module genes had no significant GTEx colon cis-eQTL instrument.
- Ran MR for IBD/CD/UC. IBD: SELE passed MR FDR<0.05 (OR 0.862, p=0.00186, FDR=0.0186) but failed coloc PP4 (0.033), so it is not a strict causal call.
- TNFRSF1A had a GTEx colon instrument but no harmonised de Lange GWAS allele row at the lead variant.
- CCL8 and CXCR2 remain blood-only causal anchors in the old eQTLGen results because GTEx colon significant-pairs yielded no gut instrument for them.
- Final strict gut causal calls under MR FDR<0.05 + coloc PP4>0.8: none.


## 2026-06-23 JST - multicontext immune eQTL run

- Pulled latest repo before work (`git pull --ff-only`, fast-forwarded to origin/main).
- Ran `src/26_immune_eqtl_multicontext.py` on CPU only.
- eQTL Catalogue datasets selected from the live `quant_method=ge` listing and cached in `/home/ujs/mycode/sasp-ibd/outputs/causal_module/eqtl_catalogue_datasets_ge.tsv`.
- Outputs: `/home/ujs/mycode/sasp-ibd/outputs/causal_module/module_causal_map_multicontext.tsv`, `/home/ujs/mycode/sasp-ibd/outputs/causal_module/Fig_module_causal_multicontext.png`, `/home/ujs/mycode/sasp-ibd/outputs/causal_module/SUMMARY.md`.
- Gut allpairs: supplied nominal allpairs files were streamed by eGenes-inferred phenotype IDs, but no target IDs matched the local allpairs `gene_id` values; audit files were written and gut coloc remains significant-pairs restricted.


## 2026-06-23 JST - causal hardening

- Ran `src/27_causal_hardening.py` on CPU only with `n_jobs=30` for parallel-safe steps.
- Downloaded/resumed CRP, FinnGen R12 endpoints, and eQTLGen full cis/allele-frequency inputs; eQTLGen HTTPS certificate is expired, so `src/27` scopes `curl --insecure` only to that host.
- Wrote hardening outputs under `/home/ujs/mycode/sasp-ibd/outputs/causal_hardening` and promoted final table/figure copies into `results/`.
- Reverse-MR is not a clean pass: significant one-SNP reverse rows appear for CD->CCL8, IBD->CCL8, and IBD->CXCR2, so reverse causality/disease-linked expression remains a caveat.
- Steiger supports exposure->outcome for the blood CCL8/CXCR2 lead instruments. CXCR2 FinnGen R12 replicates for IBD and UC but not CD. CXCR2 MVMR-CRP keeps the protective direction, but CRP conditional F is weak; CCL8 MVMR is not estimable with one instrument.
- deCODE CCL8 pQTL coloc is best-effort blocked because the public deCODE proteomics folder requires an interactive/token workflow and no local CCL8 pQTL file is present.
- OSM/OSMR/TREM1/IL13RA2 are formalised as bystanders in `bystander_triage.tsv`.


## 2026-06-23 JST - causal hardening 2

- Ran `src/28_causal_hardening2.py` CPU-only with `n_jobs=30`.
- Wrote strict bidirectional reverse-MR with genome-wide disease instruments, 1Mb distance clumping, and target-gene cis exclusion to `outputs/causal_hardening2/reverse_mr_proper.tsv`.
- Built relaxed CCL8 cis instruments from eQTLGen (`p<0.005`, 100kb distance clump), then ran CCL8 MVMR-CRP and IVW/Egger/weighted-median/leave-one-out sensitivity.
- Audited CCL8/MCP-2 plasma pQTL sources including SCALLOP CVD1 Zenodo probes, SCALLOP-INF GWAS Catalog accession GCST90274822, and the deCODE summary-data landing page; all source attempts are recorded in `pqtl_ccl8_v2.tsv`.
- Promoted the hardening2 final figure and requested tables into `results/`.
