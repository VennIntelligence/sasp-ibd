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
