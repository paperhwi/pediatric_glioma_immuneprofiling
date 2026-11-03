"""Connect non-public outcome data without storing it in the Git checkout."""
from pathlib import Path
import os
import pandas as pd
import numpy as np
from analysis_paths import OUT

REPO = next((p for p in OUT.parents if (p / ".git").exists()), None)

def outside_checkout(value):
    p = Path(value).expanduser().resolve()
    if REPO is not None and p.is_relative_to(REPO):
        raise ValueError("Restricted survival inputs and generated clinical data must be outside this Git checkout.")
    return p

def private_output_dir():
    source = os.environ.get("PGI_RESTRICTED_SURVIVAL_INPUT")
    if not source:
        raise RuntimeError("Set PGI_RESTRICTED_SURVIVAL_INPUT to the non-public, corrected 343-row survival/covariate TSV outside this checkout.")
    default = outside_checkout(source).parent / ".n343_private_outputs"
    target = outside_checkout(os.environ.get("PGI_RESTRICTED_OUTPUT_DIR", str(default)))
    target.mkdir(parents=True, exist_ok=True)
    return target

def load_survival_observations():
    source = os.environ.get("PGI_RESTRICTED_SURVIVAL_INPUT")
    if not source:
        raise RuntimeError("Set PGI_RESTRICTED_SURVIVAL_INPUT to the non-public corrected 343-row survival/covariate TSV. Individual durations and status are not redistributed.")
    a = pd.read_csv(outside_checkout(source), sep="\t", index_col=0)
    public = pd.read_csv(OUT / "inputs/public_annotations_n343.tsv", sep="\t", index_col=0)
    required = {"OS_days_corrected", "OS_event", "OS_corrected_record", "in_survival_analysis", "age_at_diagnosis_years", "sex", "autopsy_specimen", "molecular_group", "RNA_library", "immune_axis"}
    missing = required.difference(a.columns)
    if missing:
        raise ValueError("Restricted input is missing required columns: " + ", ".join(sorted(missing)))
    if not a.index.is_unique or set(a.index) != set(public.index):
        raise ValueError("Restricted input must contain exactly the 343 public biospecimen IDs, once each.")
    a = a.loc[public.index]
    for field in ["molecular_group", "RNA_library"]:
        if not a[field].equals(public[field]):
            raise ValueError("Restricted/public annotation mismatch for " + field)
    if not np.allclose(a.immune_axis, public.immune_axis, rtol=0, atol=1e-12):
        raise ValueError("Restricted input must use the fixed published n343 immune axis.")
    for field in ["OS_corrected_record", "in_survival_analysis", "autopsy_specimen"]:
        if not a[field].dropna().isin([True, False, 0, 1]).all() or a[field].isna().any():
            raise ValueError("Invalid boolean field " + field)
        a[field] = a[field].astype(bool)
    return a
