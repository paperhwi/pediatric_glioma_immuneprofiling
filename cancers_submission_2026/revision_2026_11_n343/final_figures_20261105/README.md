# Final figure typography and interpretation, November 5, 2026

This is the current figure supplement for **v2.0.3**. It supersedes the November 4 display interpretation and typography. Primary n = 343 labels (122/144/77), the 116/149/78 sensitivity partition, all fitted models, effect estimates, Kaplan–Meier curves and full-follow-up inference are unchanged.

## Figure presentation

All 54 panel letters use black (#000000) Arial Bold at 12 pt, outside the upper-left plot-label boundary. Panel titles use Arial 10 pt and ordinary plot text Arial 8 pt, with smaller mathematical superscripts/subscripts. Numeric heatmap labels and small embedded subheaders use 9–10 pt. Main and supplementary PDFs have a physical width of 178 mm, except the single-panel S5 heatmap at 89 mm. These font sizes apply when placed at the supplied width; resizing the whole figure changes them. Dense figures are taller to preserve legibility. SVG text stays editable; PDF text uses true vector outlines of the same Arial glyphs. PNGs are rendered at 400 dpi.

S5/S6 overall titles and S6 explanatory footer text have moved to FIGURE_LEGENDS.md. Plot keys, axis captions and data annotations use Arial. Axis-label spacing and row layouts were adjusted without changing data paths. S1 subtype models occupy their own row. S7A uses molecular-group abbreviations defined in its legend and the manuscript (DHG, DMG, IHG, pHGG).

## Figure 5A: keep the prespecified result separate from reinterpretation

**15 of 16 directional signatures were concordant with their originally specified lineages** (original-target median AUROC 0.868). Cytotoxic Lymphocyte Curated was the original myeloid-target exception. Its displayed bar is marked with an asterisk and reinterpreted against annotated CD4/CD8 T cells: pooled AUROC 0.925, specimen-specific range 0.927–0.976. NK cells were not independently annotated. The reinterpretation does not change original-lineage concordance to 16/16. Historical November 4 retargeted summaries remain in the earlier folder and tag solely as provenance.

`tables/Figure5A_final_display_source.tsv` distinguishes `display_target_verdict`, `prespecified_expected`, `prespecified_verdict`, and `prespecified_lineage_concordant`. Whiskers show specimen minimum–maximum, not confidence intervals. Values retain the source three-decimal precision; no unavailable per-cell score is reconstructed.

Pooled AUROC includes between-specimen as well as within-specimen positive–negative cell pairs, so specimen-specific score shifts can place the pooled value below all within-specimen estimates. There are 65,410,080 total pairs, of which 16,047,522 are within specimen and 49,362,558 are between specimens. The reported summaries are mathematically compatible with a between-specimen AUROC of 0.9084–0.9243. This is a compatibility bound, not a newly estimated study AUROC. The statistical basis is the pairwise ranking interpretation of AUROC ([Hanley and McNeil, 1982](https://pubmed.ncbi.nlm.nih.gov/7063747/)).

## Reproduce the presentation without restricted inputs

The 14 SVGs in `inputs/baseline_v2.0.1/` are the previously public figure geometry, not raw expression or patient-level clinical data. The deterministic script preserves their scientific paths, point coordinates, image payloads and markers while rendering text and geometry independently. It does not refit models. Arial regular/bold must be installed, and Cairo must be available for CairoSVG.

From this repository root:

```bash
python cancers_submission_2026/revision_2026_11_n343/final_figures_20261105/code/render_final_figures.py --baseline cancers_submission_2026/revision_2026_11_n343/final_figures_20261105/inputs/baseline_v2.0.1 --output /path/to/output/01_figures
python cancers_submission_2026/revision_2026_11_n343/validate_revision_bundle.py
```

The renderer emits publication-numbered files in main/supplementary subfolders. The repository's active `figures/` keeps the existing S1–S3 analysis filenames (`FigureS_survival`, `FigureS_state_structure`, `FigureS_immune_features`). Stored-output validation compares all scientific SVG geometry with its baseline, confirms the panel-letter style, the 15/16 original-lineage result and AUROC summary, and verifies checksums. This validates presentation and stored outputs; it does not rerun models on restricted clinical inputs.

## Public contents

This addition contains figure geometry, presentation code, aggregate AUROC checks and documentation only. No raw matrix, restricted clinical record, individual survival duration or individual correction value is added. Existing v1.0, v2.0, v2.0.1 and v2.0.2 records remain preserved. The repository-wide README lists all earlier n = 343 code, public assignments, exclusion lists, aggregate correction records and reviewer-requested sensitivities.


## November 6 layout corrections

The final release is v2.0.3. S2C now has a full-width panel, separated percentage/ARI annotations and readable 0/25/50/75/100 ticks. The adjacent reference-deconvolution panel has its own row. Figure 1B centers the sensitivity heading and places delayed-entry n = 187 and non-autopsy n = 189 on separate lines inside the box. Figure 4 uses “Monte Carlo P” and consistent “n = value” spacing.

This release also incorporates the subsequent November 5 fixes: violin distributions follow the data-axis transform without marker-size compensation, and the at-risk labels and year-zero counts are independently aligned in Figure 7A and S1B/C. The shared rendering rule applies consistently to filled distributions. All source data paths, estimates, assignments and underlying analyses are preserved. Regression checks now verify the actual violin centers, table columns, S2C tick/annotation gaps and Figure 1B text containment.
