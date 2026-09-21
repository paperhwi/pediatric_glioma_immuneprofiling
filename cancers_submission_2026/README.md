# Cancers final-submission reproducibility bundle

Prepared 22 September 2026 for the submission of

> Integrated transcriptomic immune profiling identifies survival-associated immune ecotypes
> in pediatric diffuse high-grade glioma

## Contents

| Path | Contents |
|---|---|
| `analysis_code/core_pipeline/` | cohort construction, immune profiling, clustering, pathway, differential-expression, survival, and sensitivity scripts |
| `analysis_code/revision_scripts/` | final cluster-number, contrast, selection-bias, single-cell, and spatial workflows |
| `analysis_code/notebooks/` | executed revision/finalization notebooks and retained executable sensitivity notebooks |
| `reproducibility_data/` | non-identifiable feature matrices, ecotype assignments, summary tables, signatures, software versions, parameters, seeds, and checksums |
| `main_figures/` | Figures 1-7 as vector PDF and 600 dpi PNG |
| `figure_separate_files/` | each main-figure panel as PDF and PNG, with a manifest |
| `supplementary_figures/` | Figures S1-S15 as vector PDF and 600 dpi PNG |
| `legends/` | final figure and table legends |
| `supplementary_material/` | combined supplementary document and the 25-table workbook |
| `scripts/` | final display-item assembly and verification code |

## Reproducibility data

The directly reusable public inputs include:

- locked n = 349 ecotype assignments (111 Lymphocyte-inflamed, 160 Myeloid-dominant, and
  78 Immune-desert);
- the 46-feature main LM22 plus ssGSEA standardized matrix;
- the 34-feature TPM sensitivity matrix (10 quanTIseq plus 24 ssGSEA features);
- PCA/UMAP coordinates and k = 2-10 selection/stability summaries;
- molecular-group-adjusted differential-expression and contrast summaries;
- single-cell, spatial, survival-robustness, software-environment, parameter, and seed summaries.

Detailed file descriptions are in [`reproducibility_data/README.md`](reproducibility_data/README.md).

## Analysis order

The original workflow is filename-driven. The core order is:

```text
cohort construction and TPM extraction
  -> immune deconvolution and ssGSEA
  -> feature integration
  -> locked k = 3 consensus clustering
  -> cluster-number and resampling-stability evaluation
  -> ecotype annotation and clinical comparisons
  -> pathway and differential-expression analyses
  -> survival and sensitivity analyses
  -> single-cell/spatial support analyses
  -> final figure, table, legend, and supplementary-document assembly
```

Scripts preserve their analysis-time paths for provenance. Before rerunning, point their input and
output constants at local copies of the required OpenPedCan/OpenPBTA resources and the derived
files documented here. The final locked settings include random seed 20260722, 500 bootstrap
replicates for the locked k = 3 assignment, 1,000 replicates across k = 2-10 for cluster-number
evaluation/resampling stability, 4,999 PERMANOVA permutations, and 9,999 contingency Monte Carlo
replicates.

## Validation

Run from this directory:

```bash
python analysis_code/validate_public_bundle.py
```

## Public-data boundary

Raw TPM/count matrices, restricted clinical tables, individual survival times or vital status,
histology records, and RDS/H5AD objects are excluded. Sample identifiers in the public derived
matrices and ecotype assignment are pseudonymous Kids First biospecimen identifiers.
