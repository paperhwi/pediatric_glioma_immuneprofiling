# Analysis code

This directory contains the analysis code retained for the September 2026 *Cancers* submission.

- `revision_scripts/` contains the Python scripts used for cluster-number evaluation,
  molecular-group-adjusted differential-expression contrasts, selection-bias checks,
  hematopoietic-progenitor analyses, and single-cell/spatial support analyses.
- `notebooks/` contains the corresponding revision notebooks plus the PERMANOVA interaction,
  TPM harmonization, and inverse-probability-weighting sensitivity notebooks.
- `validate_public_bundle.py` checks the public data boundary and the locked analysis invariants.

All notebook outputs and execution counters were removed before deposit. This prevents cached raw
or clinical-level records from being embedded while preserving executable code and Markdown
documentation. Scripts retain their original provenance and may refer to the project layout used
during analysis. Controlled OpenPBTA expression/clinical inputs and large public-accession objects
must be supplied separately; they are not redistributed by this repository.

The main inputs that can be shared safely are in `../reproducibility_data/`. Figure assembly code
is kept separately in `../scripts/`.
