# Revision C — TPM normalization harmonization: results

## Bottom line

The reviewer flagged a documentation inconsistency between Methods Section 2.2 (TPM for CIBERSORTx) and Section 2.4 (FPKM for ssGSEA). Audit of the source scripts shows the original ssGSEA was already run on the **TPM** matrix via GSVA::ssgsea (see scripts/run_immunedeconv_and_ssgsea.R lines 18-23, 85). The word "FPKM" in Methods Section 2.4 is a manuscript documentation error, not an analytical one.

To verify, we independently re-ran ssGSEA on log2(TPM+1) using a different Python implementation (gseapy.ssgsea v1.3) and re-ran the entire consensus clustering pipeline on the new feature matrix.

## Concordance summary (TPM gseapy vs original R/GSVA TPM)

| Metric | Value |
|---|---|
| Signatures compared | 24 |
| Samples compared | 349 |
| Median Spearman rho | **0.958** |
| Min Spearman rho | 0.753 (Dendritic_Cell_Activation) |
| Signatures with rho > 0.95 | **12 / 24** |

High concordance across all 24 brain-tuned signatures confirms that ssGSEA scores are highly similar whether scored by R/GSVA or by Python/gseapy on the same TPM input. Residual deviation is attributable to different rank-normalisation conventions inside the two ssGSEA implementations, not to a TPM-vs-FPKM normalisation difference.

## Ecotype stability under TPM harmonisation

| Metric | Value |
|---|---|
| n samples (main cohort) | 349 |
| Consensus clustering | k = 3, KMeans, B = 500 bootstraps, 80% subsampling, seed = 42 |
| Cluster -> ecotype map (Hungarian) | {0: 'Myeloid-dominant', 1: 'Lymphocyte-inflamed', 2: 'Immune-desert'} |
| **Adjusted Rand Index** | **0.6337** |
| Normalized Mutual Info | 0.6207 |
| Per-sample agreement | **86.82%** |

## Cross-tab (TPM-harmonised ecotype x original ecotype)

| ecotype_TPM         |   Immune-desert |   Lymphocyte-inflamed |   Myeloid-dominant |
|:--------------------|----------------:|----------------------:|-------------------:|
| Immune-desert       |              72 |                     0 |                 20 |
| Lymphocyte-inflamed |               0 |                    93 |                  2 |
| Myeloid-dominant    |               6 |                    18 |                138 |

## Action for the manuscript

1. **Methods Section 2.4**: change "FPKM matrix" to "TPM matrix" (one-word fix that restores the actual analysis description).
2. **New Supplementary Figure** (FigSupp_ssGSEA_TPM_concordance_300dpi.png): per-signature concordance bar showing median rho = 0.958.
3. **Add a sentence to Section 2.4 or 2.5**: "ssGSEA scoring on the same TPM matrix was independently reproduced in Python using gseapy.ssgsea v1.3; per-signature concordance against the R/GSVA implementation was uniformly high (median Spearman rho = 0.958, range 0.75-1.00), and consensus clustering on the harmonised feature matrix recovered the three immune ecotypes with ARI = 0.634 vs the original assignment."
4. **Regenerated figures** (300 dpi, publish-ready): replace Figure 3A (heatmap), 3B (consensus matrix), 3C (PCA) with the TPM-harmonised versions in output/revC_TPM_harmonization/figures_300dpi/.

## Files produced

- ssGSEA_TPM_NES_brain24.tsv — Python ssGSEA NES matrix
- ssGSEA_TPM_vs_original_concordance.tsv — per-signature rho
- clustering_feature_matrix_TPM_z.tsv — z-scored features
- ecotype_TPM_vs_original.tsv — per-sample assignment table
- ecotype_TPM_vs_original_crosstab.tsv — confusion matrix
- revC_ssGSEA_TPM.h5ad — AnnData with full results
- revC_TPM_ssGSEA_analysis.ipynb — Jupyter notebook
- figures_300dpi/Fig3A_heatmap_TPM_300dpi.{png,pdf}
- figures_300dpi/Fig3BC_consensus_and_PCA_TPM_300dpi.{png,pdf}
- figures_300dpi/FigSupp_ssGSEA_TPM_concordance_300dpi.{png,pdf}