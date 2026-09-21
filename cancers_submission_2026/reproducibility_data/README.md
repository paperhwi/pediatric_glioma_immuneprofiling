# Public reproducibility data

This directory contains non-identifiable derived data needed to inspect or reproduce the principal
analysis outputs. Sample-level files use pseudonymous Kids First biospecimen identifiers only.

| Path | Rows x data columns | Description |
|---|---:|---|
| `ecotype_assignments.tsv` | 349 x 2 | Locked k = 3 cluster number and ecotype label |
| `feature_matrices/main_cohort_LM22_ssGSEA_z.tsv` | 349 x 46 | Standardized 22-feature LM22 plus 24-feature ssGSEA matrix used for the main clustering |
| `feature_matrices/TPM_harmonization_quanTIseq_ssGSEA_z.tsv` | 349 x 29 | Standardized 10-feature quanTIseq plus 19-feature retained ssGSEA TPM sensitivity matrix |
| `embeddings/pca_coordinates.tsv` | 349 x 3 | PCA coordinates and locked ecotype |
| `embeddings/umap_coordinates.tsv` | 349 x 3 | UMAP coordinates and locked ecotype |
| `clustering/` | aggregate | k = 2-10 selection criteria and 1,000-resample consensus summaries |
| `differential_expression/` | 55,758 rows | Molecular-group-adjusted gene-level statistics for the three ecotype contrasts |
| `contrast_structure/` | aggregate/gene-level | Contrast overlap and ordered shared immune-presence-axis summaries |
| `single_cell_spatial/` | aggregate | Signature AUROC, cell-type source, spatial assignment, autocorrelation, and methods summaries |
| `survival_robustness/` | aggregate | OS evaluability, PERMANOVA/Cox-interaction summary, and IPW sensitivity summary |
| `signatures/` | gene sets | Brain immune-program GMT used by the revision workflows |

`T5_analysis_set_reconciliation.tsv` documents the analysis-set counts and
`T5_software_environment.tsv` records the original software environment. `SHA256SUMS.tsv` covers
every deposited file in this directory except itself.

The TPM harmonization matrix contains 29 retained features (10 quanTIseq and 19 ssGSEA), matching
the deposited file exactly. It is distinct from the 46-feature main LM22 + ssGSEA matrix.

Excluded by design: raw TPM/count matrices, individual-level clinical variables, survival time or
status, histology records, and serialized RDS/H5AD objects.
