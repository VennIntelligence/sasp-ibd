# Single-cell deconvolution response summary

## Direct answer
- Strongest NNLS cell-fraction signal: **Myeloid** (higher in non-response; random-effects OR/SD=2.65, p=4.83e-07).
- Myeloid: OR/SD=2.65, p=4.83e-07, I2=0.12.
- Fibroblast: OR/SD=1.58, p=0.115, I2=0.69.
- neutrophil_marker_score: OR/SD=3.16, p=0.000192, I2=0.50.
- refractory_module_score: OR/SD=3.31, p=0.000137, I2=0.50.
- cxcr2_expr_z: OR/SD=2.93, p=0.00103, I2=0.56.
- senmayo: OR/SD=3.22, p=4.2e-05, I2=0.46.

## Triangulation
- Explicit neutrophil cluster in the local scRNA references: **no**. Therefore neutrophils are not interpreted as a learned NNLS fraction; they are tested as a bulk CXCR2/neutrophil marker score.
- CXCR2 remains genetically protective in neutrophil/blood contexts, so a high neutrophil/CXCR2 mucosal marker signal is interpreted as disease-state/refractory biology, not proof that blocking CXCR2 is beneficial.
- The refractory module score is a targeted bulk triangulation of the neutrophil + OSM-fibroblast + myeloid secretion module rather than a cell fraction.

## Increment over bulk SASP
```tsv
feature	n	n_NR	senmayo_auc	senmayo_plus_feature_auc	delta_auc	feature_beta_adjusted_for_senmayo	feature_p_adjusted_for_senmayo	lrt_p_increment	status
B_Plasma	239	117	0.7927000140114894	0.80369903320723	0.010999019195740645	-0.27297817932622387	0.10776559366313002	0.10230903297346118	ok
Myeloid	239	117	0.7927000140114894	0.8004063331932184	0.007706319181729038	0.5036492355228369	0.07278489192138388	0.06874032558860035	ok
Fibroblast	239	117	0.7927000140114894	0.7945215076362618	0.0018214936247724633	-0.19054216771741805	0.4302117293948019	0.43058843929137236	ok
Stromal_other	239	117	0.7927000140114894	0.7943813927420484	0.0016813787305590688	0.12738721253782598	0.4775599032920679	0.4768208974165249	ok
Endothelial	239	117	0.7927000140114894	0.7940311055065153	0.0013310914950259711	0.21728543154748114	0.2877859156034356	0.28742456271278427	ok
T_NK	239	117	0.7927000140114894	0.7931904161412358	0.0004904021297463812	-0.014203196524027163	0.9454175089475213	0.9454175034067063	ok
Epithelial	239	117	0.7927000140114894	0.7918593246462099	-0.0008406893652794789	-0.030295736489532467	0.9126517934856933	0.9125952481041852	ok
Mast	239	117	0.7927000140114894	0.7913689225164635	-0.00133109149502586	-0.1710119398763415	0.5465298105740235	0.5459498545755177	ok
```

## Method caveats
- NNLS uses broad single-cell signatures and marker-restricted rank targets to reduce cross-platform scale mismatch. Fractions should be read as relative cell-type enrichment estimates, not absolute histologic percentages.
- Smillie UC has detailed author cell labels; Martin CD labels are broad marker-derived labels from local 10x matrices.
- Missing neutrophils in the scRNA reference limits direct neutrophil deconvolution; neutrophil evidence comes from bulk marker/CXCR2 scoring.
