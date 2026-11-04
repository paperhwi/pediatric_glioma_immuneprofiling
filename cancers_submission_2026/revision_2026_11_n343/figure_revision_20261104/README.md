# November 4 figure and target corrections

Version **v2.0.1** corrects the November 3 **v2.0** n343 revision. v2.0 and the submitted **v1.0** snapshot remain immutable. Primary 122/144/77 assignments, 116/149/78 scaling sensitivity, survival inputs, Cox/IPW fits, and signature gene memberships are unchanged.

## Figure 5A target

`tables/Figure5A_AUROC_current_targets.tsv` is the active display source. The curated cytotoxic signature is evaluated against annotated CD4/CD8 T cells, giving pooled AUROC **0.925**, four-specimen range **0.927–0.976**. NK cells were not separately annotated. The original notebook's 18,619 eligible cells partition into 13,920 Myeloid and 4,699 CD4/CD8 cells; 93 unclear cells are excluded. The unchanged score therefore permits the exact complementary-label identity AUROC(complement) = 1 − AUROC(original), including ties, at the original table's three-decimal precision. Range endpoints reverse. No unavailable cell-level score precision is invented. See `TARGET_CORRECTION.md` and `tables/Figure5A_target_audit.json`. The other 15 target AUROCs are unchanged; current lineage concordance is 16/16 and median AUROC is 0.879.

The earlier `reviewer_followup_20261103/tables/Figure5A_source_AUROC.tsv` is a **historical original-target audit**, not the current Figure 5A source. The old target is mentioned in one sentence in the revised Figure 5 legend.

## Figure presentation

All 14 main/supplementary figures have PNG, vector PDF, and editable SVG files in `../figures/`. Figure 2B uses horizontal bars and Figure 2C a 2×2 layout; Figure 5B uses horizontal cell-class bars. Panel letters are independent, large and bold. American spelling, the S1 prior-label phrase, S3C capitalization, S4 header/legend separation and S7 clinical factor labels are corrected. Displayed P values use three decimals for P ≥ 0.001 and one mantissa decimal with a superscript power of ten below 0.001. Source tables retain their existing numerical precision. Five-year KM displays and full-follow-up modeling are unchanged.

## Reproduction

`code/render_updated_figures.py` regenerates Figure 2, Figure 5 and Figure S7 from source data and updates the other SVGs. Supply the immutable **v2.0 figures directory** as `--baseline-figures`, an authorized local JSON matching `local_inputs.example.json`, and a separate output directory. Do not use v2.0.1 SVGs as the baseline for a second layout pass.

```text
python code/render_updated_figures.py --input-config /authorized-inputs/local_inputs.json --baseline-figures /v2.0-snapshot/cancers_submission_2026/revision_2026_11_n343/figures --output-dir /local-output/figure-update
```

Input matrices and cell/spot-level files remain outside this repository. The renderer exports figures and text-edit logs, not clinical observations. It uses the final n343 Visium H5AD for current centroids, axis scores, coordinates and deconvolution shares. The historical Figure 5 source layout remains in `analysis_code/s23_fig5.py`, now pointed at the corrected current-target AUROC source; use the new renderer for the actual publication layout. Existing survival plotting code also uses the corrected P formatter and prior-label wording.

Requested clinical supplements, full manuscript DOCX, and the authors' delivery ZIP are not uploaded. Individual requested tables and references [1]–[54], unchanged in numbering, are delivered locally to the authors. `../DATA_AVAILABILITY.md` provides the current and preceding release links.
