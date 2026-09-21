# Analysis code

This directory contains the code retained for the September 2026 *Cancers* submission.

- `core_pipeline/` contains the main cohort, deconvolution, integration, clustering, pathway,
  differential-expression, survival, and sensitivity scripts.
- `revision_scripts/` contains the final cluster-number, contrast, selection-bias, single-cell,
  and spatial support workflows.
- `notebooks/` contains executed revision and finalization notebooks. A small number of retained
  sensitivity notebooks are executable source notebooks without cached output because their
  controlled inputs cannot be redistributed.
- `validate_public_bundle.py` checks the public data boundary and locked analysis invariants.

The public inputs that can be shared safely are in `../reproducibility_data/`. Controlled
OpenPedCan/OpenPBTA expression and clinical inputs and large public-accession objects must be
obtained separately. Scripts preserve their analysis-time paths for provenance and may require
local path configuration before rerunning.

Figure and document assembly code is kept separately in `../scripts/`.
