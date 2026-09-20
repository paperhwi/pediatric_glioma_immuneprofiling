# Cancers submission — display items and deposited analyses

Prepared 21 September 2026 for the submission of

> Integrated transcriptomic immune profiling identifies survival-associated immune ecotypes
> in pediatric diffuse high-grade glioma

to *Cancers*. The manuscript text is not in this repository. Everything here is either a
display item of that submission or an analysis that was removed from it and is deposited
so that the record stays complete.

## Layout

| Path | Contents |
|---|---|
| `main_figures/` | Figure 1-7, vector PDF and 600 dpi PNG |
| `supplementary_figures/` | Figure S1-S16, vector PDF |
| `legends/` | all figure and table legends |
| `supplementary_material/` | the combined supplementary document and the 25-table workbook |
| `scripts/` | every script that generates or assembles the above, plus the verification pass |
| `archived_analyses/` | the analyses removed from the submission, with their source tables |

## Main figures

| # | Panels |
|---|---|
| 1 | A study workflow; B derivation of the analysis sets |
| 2 | A LM22 quality control; B cross-method concordance; C cross-method estimates by ecotype |
| 3 | A PAC; B silhouette; C relative dAUC; D gap statistic; E minimum cluster size; F PCA and UMAP; G defining-feature heatmap |
| 4 | A anatomical location; B integrated molecular group; C within-group comparison; D sequential PERMANOVA |
| 5 | A single-cell AUROC; B cell type expressing the ecotype-axis genes; C spot-to-centroid assignment in three sections; D spatial coherence; E Moran's I |
| 6 | A-C molecular-group-adjusted volcanoes; D contrast overlap; E gradient; F paired effect sizes |
| 7 | A overall Kaplan-Meier; B multivariable Cox |

Supplementary Figures S1-S16 and Supplementary Tables S1-S25 are renumbered from scratch for
this submission. The `Crosswalk` sheet of `supplementary_material/Supplementary_Tables_Cancers.xlsx`
maps every item of the previous version onto this one, including the items that were moved
here instead of being submitted.

## Reproducing the figures

`scripts/` is self-contained apart from its inputs, which are the locked pipeline outputs in
`output/` of this repository and the revision tables listed at the top of each script. Order:

```
fig1.py            Figure 1
fig3_ksel.py       Figure 3A-E
fig5.py            Figure 5A, 5B, 5D, 5E
fig5_spatial.py    Figure 5C
fig6.py            Figure 6A-F
build_main.py      assembles Figure 1-7
build_suppl.py     assembles Figure S1-S16
build_tables.py    builds the supplementary table workbook
build_docs.py      builds the legend and supplementary documents
verify.py          checks files, panel letters and cross-references
append_audit.py    writes the citation checklist into the workbook
```

`compose.py` places existing single-panel PDFs into one multi-panel vector PDF; `patch_text.py`
rewrites text inside a reused panel without rasterising it; `legends.py` holds the legend text.

## Deposited analyses

`archived_analyses/` holds the pooled-cohort sensitivity analyses, the BRAF/RTK projection,
the checkpoint and CAR-T transcript panels, the differential-expression audit copy, the very
small survival strata and the hematopoietic progenitor analysis. See its `README.txt` for what
each one shows and why it is not in the submission.
