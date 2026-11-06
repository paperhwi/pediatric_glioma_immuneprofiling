# Revision bundle for the final cohort n = 343

This directory contains the current November 2026 revision. The original n = 349 September submission in the parent directory remains a historical record.

## Current definitions

Six specimens outside the pediatric-type diffuse high-grade glioma cohort are excluded from the final analyses, including clustering. The retained 46-feature standardized input is restricted to the 343 eligible specimens; this is not a claim that all original features were newly standardized after exclusion.

Immune-state labels were re-derived with the locked consensus procedure and assigned by their mode over seeds 42–66: **Immune-high 122, Immune-intermediate 144, Immune-low 77**. Historical label strings in source-loading or comparison code are provenance identifiers, not alternative current state names.

The primary survival specification remains immune state plus molecular group, age, and sex. One origin-review rule is applied throughout the cohort, with four records re-timed. Aggregate sensitivity results are retained regardless of significance. Individual survival durations are not included in this repository.

## Contents

- `analysis_code/`: the delivered analysis sources, including shared definitions and numeric pipeline stages. Original absolute paths require adaptation and some stages depend on intermediate outputs not bundled here.
- `reproducibility_data/`: final assignments, aggregate results, differential-expression tables, clustering diagnostics, seed records, and the software/parameter workbook.
- `figures/`: delivered main and supplementary PDFs with their original analysis filenames.
- `validate_revision_bundle.py`: checks consistency of the stored public outputs.

## Publication numbering crosswalk

The cleaned manuscript package uses these delivery labels; repository analysis filenames remain stable for provenance.

| Submission label | Repository analysis figure |
|---|---|
| Figure S1 | `FigureS_survival.pdf` |
| Figure S2 | `FigureS_state_structure.pdf` |
| Figure S3 | `FigureS_immune_features.pdf` |
| Figure S4 | `figures/FigureS4_ScPCA_axis_linkage.pdf`; final public SVG rendering source in `final_figures_20261105/` |
| Supplementary Table S1 | Detailed n343 analysis workbook supplied with the manuscript |
| Supplementary Table S2 | Software/parameter workbook retained here under its original `Supplementary_Table_S25_software_parameters_seeds_n343.xlsx` filename |
| Supplementary Table S3 | [`additional_analyses_20261103/supplementary/Supplementary_Table_S3_gene_sets_and_sources.xlsx`](additional_analyses_20261103/supplementary/Supplementary_Table_S3_gene_sets_and_sources.xlsx) |
| Supplementary Table S4 | Single-cell per-patient summaries supplied with the manuscript |
| Figure S5 | [`figures/FigureS5_standardization_consistency.pdf`](figures/FigureS5_standardization_consistency.pdf) |
| Figure S6 | [`figures/FigureS6_axis_library_survival.pdf`](figures/FigureS6_axis_library_survival.pdf) |
| Supplementary Table S6 | [`additional_analyses_20261103/supplementary/Supplementary_Table_S6_standardization_library_survival.xlsx`](additional_analyses_20261103/supplementary/Supplementary_Table_S6_standardization_library_survival.xlsx) |
| Supplementary Table S5 | Long-observation and survival-origin-correction records supplied with the manuscript; individual survival records are not redistributed here |

## Validation and analytic execution

```bash
python cancers_submission_2026/revision_2026_11_n343/validate_revision_bundle.py
```

This validator checks stored primary and follow-up results and does not fit models. The two additional sensitivity analyses were executed separately; the original primary results and gene membership remain unchanged.

Full analytic execution requires source OpenPedCan v15 expression and clinical files, GSE268577, ScPCA SCPCP000001, and intermediate files identified by the scripts. Do not assume that this folder is a standalone end-to-end rerun environment. Before corrected-OS survival stages, the uniform audit in `s02_os_audit.py` must derive the local, git-ignored `restricted_os_fix.tsv`; subsequent survival scripts must run in a fresh process so `common.py` reads the corrections. Without the required inputs and corrections, stored-output validation is the supported check.

The additional re-standardization and continuous-axis library-control results are supplied with the revised manuscript and available in [`additional_analyses_20261103/`](additional_analyses_20261103/). Re-standardization preserves 334/343 modal labels (97.4%; ARI 0.915), with 116/149/78 sensitivity counts. The library-adjusted axis estimate attenuates to HR 0.861 (95% CI 0.726–1.022), P 0.088. Primary labels remain 122/144/77.

Corrected-OS loading now fails explicitly if the four locally derived corrections are absent. Run the uniform audit and restart the analysis process; released durations must not silently stand in for corrected outcomes.

## Terminology consistency

Use [TERMINOLOGY.md](TERMINOLOGY.md) and the current [figure legends](FIGURE_LEGENDS.md). Canonical state names are Immune-high, Immune-intermediate, and Immune-low throughout current figures and Table S6. Updated primary SVG/PNG sources accompany the PDFs; Figure S4 is [figures/FigureS4_ScPCA_axis_linkage.pdf](figures/FigureS4_ScPCA_axis_linkage.pdf). Historical n349 material retains its original wording.


## Final reviewer follow-up and version 2.0

See [reviewer_followup_20261103](reviewer_followup_20261103/) for the re-standardized-partition primary Cox model, simultaneous autopsy/library categorical models, restored IPW analysis, and updated submitted Tables 2–5. The original feature-scaling reference is the initial 702-specimen pool (mean and sample SD, ddof=1). Current gene definitions and the exact-match limitations are documented in the [signature provenance audit](reviewer_followup_20261103/SIGNATURE_PROVENANCE.md).

All seven main and seven supplementary figures are provided in SVG, PDF and PNG; the existing analysis filenames for S1–S3 are retained. Figure 7 and the subtype KM panels display five years while using full follow-up. Figure S7 restores the outcome-availability/IPW sensitivity. Additional supplements S7 and S8 are in `reviewer_followup_20261103/supplementary/`. Gene-source S3 and standardization S6 are also updated in the prior addition.

## Final November 6 version v2.0.4

The final Data Availability link is https://github.com/paperhwi/pediatric_glioma_immuneprofiling/releases/tag/v2.0.4 . Use [final_figures_20261105](final_figures_20261105/) for the current SVG/PDF/PNG figures, Arial 12-pt panel-letter specifications, revised legends, reproducible rendering code and unchanged-geometry checks. The active figures retain the stable analysis filenames listed above.

Figure 5A marks Cytotoxic Lymphocyte Curated with an asterisk: originally myeloid, displayed after reinterpretation against annotated CD4/CD8 T cells (pooled AUROC 0.925; specimen range 0.927–0.976). **Original-lineage concordance remains 15/16**, with original-target median AUROC 0.868. The November 4 retargeted summaries are historical provenance; their 16/16 count is not the prespecified-lineage result. Pooling adds between-specimen cell comparisons and can yield an AUROC below each within-specimen estimate.

S5/S6 overall titles and S6 footer explanations are now in [FIGURE_LEGENDS.md](FIGURE_LEGENDS.md). All 14 figure sets use Arial and consistent panel letters. Primary model results, clustering assignments, full-follow-up inference and source geometry are unchanged. Earlier tags and archived analysis outputs are preserved.


## November 6 layout corrections

The final release is v2.0.4. S2C now has a full-width panel, separated percentage/ARI annotations and readable 0/25/50/75/100 ticks. The adjacent reference-deconvolution panel has its own row. Figure 1B centers the sensitivity heading and places delayed-entry n = 187 and non-autopsy n = 189 on separate lines inside the box. Figure 4 uses “Monte Carlo P” and consistent “n = value” spacing.

This release also incorporates the subsequent November 5 fixes: violin distributions follow the data-axis transform without marker-size compensation, and the at-risk labels and year-zero counts are independently aligned in Figure 7A and S1B/C. The shared rendering rule applies consistently to filled distributions. All source data paths, estimates, assignments and underlying analyses are preserved. Regression checks now verify the actual violin centers, table columns, S2C tick/annotation gaps and Figure 1B text containment.

Figure 1B spacing correction in v2.0.4: the red missing-OS box is separated from the blue survival-analysis box by about 10 pt at the supplied figure width. Its original count of 98 and wording are retained across three lines, and the connector is repositioned between the boxes.
