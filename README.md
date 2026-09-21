# Pediatric glioma immune profiling

Reproducibility materials for the *Cancers* submission:

> Integrated transcriptomic immune profiling identifies survival-associated immune ecotypes
> in pediatric diffuse high-grade glioma

The publication link will be added after acceptance.

## Final Cancers submission bundle

[`cancers_submission_2026/`](cancers_submission_2026/) contains the analysis code, executed
notebooks, non-identifiable derived feature matrices, locked ecotype assignments, software and
package versions, analysis parameters and random seeds, final figures, legends, and supplementary
files used for the September 2026 submission.

The repository is intentionally limited to analyses reported in that submission. Raw expression
data, restricted clinical data, individual survival times or status, histology records, and large
serialized analysis objects are not redistributed. Sample-level derived files use pseudonymous
Kids First biospecimen identifiers only.

## Quick validation

From the repository root:

```bash
python cancers_submission_2026/analysis_code/validate_public_bundle.py
```

This checks the locked cohort size (n = 349), ecotype counts (111/160/78), feature-matrix
dimensions, identifier agreement, checksums, and the public-data boundary.

See the [submission bundle README](cancers_submission_2026/README.md) for the directory map,
analysis order, and input requirements.
