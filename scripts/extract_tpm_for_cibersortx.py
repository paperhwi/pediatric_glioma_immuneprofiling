"""
Extract TPM expression subset for CIBERSORTx upload.

Reads gene-expression-rsem-tpm-collapsed.rds (gene x sample),
subsets to the 702 cohort biospecimen IDs (Main 349 + BRAF_main 31 + BRAF_LGG 322),
and writes a CIBERSORTx-ready tab-delimited file:

    GeneSymbol  Sample1  Sample2  ...
    A1BG        2.34     0.0      ...
    A2M         12.5     8.7      ...
"""
import pyreadr
import pandas as pd
import numpy as np
from pathlib import Path
import sys

ROOT = Path(r"C:\Users\scott\Downloads\Open PBTA")
DATA = ROOT / "data"
OUT  = ROOT / "output"

RDS_PATH = DATA / "gene-expression-rsem-tpm-collapsed.rds"
COHORT_PATH = OUT / "cohort_for_deconvolution.tsv"

print(f"[1/5] Reading RDS: {RDS_PATH}")
print(f"      File size: {RDS_PATH.stat().st_size / 1e6:.1f} MB")

res = pyreadr.read_r(str(RDS_PATH))
print(f"      Keys in RDS: {list(res.keys())}")

# The collapsed RDS contains a single matrix/dataframe
key = list(res.keys())[0]
df = res[key]
print(f"      Shape (raw): {df.shape}")
print(f"      Column types: {df.dtypes.iloc[:5].to_dict()}")
print(f"      First 3 columns: {list(df.columns[:3])}")
print(f"      First 3 row indices: {list(df.index[:3])}")

# OpenPedCan collapsed RDS: rows = genes (HGNC symbol as rownames), columns = biospecimens
# Confirm by checking column names look like Kids_First_Biospecimen_IDs (BS_XXXXXXXX)
sample_cols = [c for c in df.columns if str(c).startswith("BS_")]
print(f"      BS_*-style columns: {len(sample_cols)} / {df.shape[1]}")

# Gene symbols are in the index (e.g. 'RN7SL1', 'RN7SL2', ...).
# Move them to a "GeneSymbol" column.
df = df.reset_index()
old_first = df.columns[0]
df = df.rename(columns={old_first: "GeneSymbol"})

print(f"      Shape after gene-col promotion: {df.shape}")
print(f"      First 3 gene symbols: {list(df['GeneSymbol'].head(3))}")

# -------- Load cohort biospecimen IDs --------
print(f"\n[2/5] Loading cohort biospecimen IDs from {COHORT_PATH}")
cohort = pd.read_csv(COHORT_PATH, sep="\t", dtype=str)
bs_ids = cohort["Kids_First_Biospecimen_ID"].unique().tolist()
print(f"      Cohort biospecimen IDs: {len(bs_ids)}")

# -------- Subset --------
available = [b for b in bs_ids if b in df.columns]
missing = [b for b in bs_ids if b not in df.columns]
print(f"\n[3/5] Matching against expression matrix columns")
print(f"      Available in expression matrix: {len(available)} / {len(bs_ids)}")
print(f"      Missing (no expression data): {len(missing)}")
if missing[:5]:
    print(f"      First missing IDs: {missing[:5]}")

# Final subset
sub = df[["GeneSymbol"] + available].copy()

# -------- Sanity QC --------
print(f"\n[4/5] QC checks")
print(f"      Shape: {sub.shape}")
print(f"      Total NA values: {sub.iloc[:, 1:].isna().sum().sum()}")
print(f"      Duplicate gene symbols: {sub['GeneSymbol'].duplicated().sum()}")
neg = (sub.iloc[:, 1:] < 0).sum().sum()
print(f"      Negative values (should be 0 for TPM): {neg}")
print(f"      Median per-sample non-zero TPM: {sub.iloc[:, 1:].replace(0, np.nan).median().median():.2f}")

# Drop duplicate gene symbols if any (keep first)
if sub["GeneSymbol"].duplicated().any():
    sub = sub.drop_duplicates(subset="GeneSymbol", keep="first")
    print(f"      Deduplicated: {sub.shape}")

# Drop rows with all-zero expression (no info)
nonzero = (sub.iloc[:, 1:].sum(axis=1) > 0)
n_drop = (~nonzero).sum()
sub = sub[nonzero].reset_index(drop=True)
print(f"      Dropped {n_drop} all-zero genes -> {sub.shape}")

# -------- Save --------
out_tpm = OUT / "tpm_for_cibersortx.tsv"
sub.to_csv(out_tpm, sep="\t", index=False)
print(f"\n[5/5] Wrote: {out_tpm}")
print(f"      Size: {out_tpm.stat().st_size / 1e6:.1f} MB")

# Also save a manifest of which biospecimens are included with their cohort_group
manifest = cohort[cohort["Kids_First_Biospecimen_ID"].isin(available)][
    ["Kids_First_Biospecimen_ID", "Kids_First_Participant_ID", "cohort_group",
     "harmonized_diagnosis", "molecular_subtype", "location_class",
     "age_years", "age_dev_group", "reported_gender", "has_survival"]
].copy()
manifest.to_csv(OUT / "tpm_manifest.tsv", sep="\t", index=False)
print(f"      Manifest: {OUT / 'tpm_manifest.tsv'} ({len(manifest)} samples)")

# Per-group sample count in final TPM
print(f"\n=== Final samples in TPM matrix by cohort_group ===")
g_counts = manifest["cohort_group"].value_counts()
print(g_counts)
print(f"Total: {len(manifest)}")
