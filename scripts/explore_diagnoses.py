"""Quick exploration of OpenPedCan v15 to scope the cohort."""
import pandas as pd
from pathlib import Path

DATA = Path(r"C:\Users\scott\Downloads\Open PBTA\data")
H = pd.read_csv(DATA / "histologies.tsv", sep="\t", low_memory=False, dtype=str)

print(f"Total rows: {len(H):,}")
print(f"Unique participants: {H['Kids_First_Participant_ID'].nunique():,}")
print(f"Unique biospecimens:  {H['Kids_First_Biospecimen_ID'].nunique():,}")

print("\n--- cohort ---")
print(H["cohort"].value_counts(dropna=False).head(20))

print("\n--- experimental_strategy ---")
print(H["experimental_strategy"].value_counts(dropna=False))

print("\n--- broad_histology (top 20) ---")
print(H["broad_histology"].value_counts(dropna=False).head(20))

# Focus on RNA-seq, PBTA cohort
rna = H[(H["experimental_strategy"] == "RNA-Seq") & (H["cohort"] == "PBTA")]
print(f"\n=== PBTA RNA-Seq rows: {len(rna):,} ===")

print("\n--- pathology_diagnosis (RNA-seq, PBTA) ---")
print(rna["pathology_diagnosis"].value_counts(dropna=False).head(30))

print("\n--- harmonized_diagnosis (RNA-seq, PBTA, top 30) ---")
print(rna["harmonized_diagnosis"].value_counts(dropna=False).head(30))

print("\n--- molecular_subtype (RNA-seq, PBTA, top 40) ---")
print(rna["molecular_subtype"].value_counts(dropna=False).head(40))
