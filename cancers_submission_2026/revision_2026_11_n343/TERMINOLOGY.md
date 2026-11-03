# Immune-state terminology

Use `immune state` for the categorical study variable and `immune-state` as an attributive adjective. The three display names are `Immune-high`, `Immune-intermediate`, and `Immune-low`. Use the full names in figures, risk tables, captions, workbook headings, and prose. Line breaks after `Immune-` are permitted in compact figures.

Keep the continuous immune axis distinct from the categorical states. Thorsson immune subtypes, Bagaev classes, external published immune clusters, diagnostic high-grade/low-grade labels, and biological cell states retain their own meanings. Machine-readable columns such as `n_High` remain stable provenance identifiers.

The November 3 terminology pass preserves all numeric figure text, plotted paths/markers/raster panels, counts, assignments, model estimates, and formula values. It updates Figures 2–5 and 7, Supplementary Figures S1/S2/S4/S5, Table S6 headings, and the current figure legends. Edited SVGs retain selectable text; the PDF exports outline text using the source fonts and retain vector scientific graphics. Final PNGs preserve the delivered convention of omitting panel letters.

`analysis_code/standardize_immune_state_labels.py` can normalize legacy SVG labels and regenerate vector PDF/PNG outputs from a package root containing `figures/`. This is a graphics-only operation; it does not load clinical records or fit models. The additional-analysis plotting script now emits full state names directly. The older primary renderer can use the graphics normalization step after generation; its stable `SHORT` dictionary is also used in stored data-column names and must not be globally renamed.
