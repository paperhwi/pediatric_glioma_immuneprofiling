"""
Build pediatric glioma cohort from OpenPedCan v15.

Main cohort (per guideline):
  1. DMG_K27   - Diffuse midline glioma, H3 K27/K28-altered
  2. DHG_G34   - Diffuse hemispheric glioma, H3 G34/G35-mutant
  3. pHGG_WT   - Pediatric-type HGG, H3/IDH-wildtype
  4. IHG       - Infant-type hemispheric glioma (age<18mo + hemispheric + RTK/BRAF fusion)

Secondary:
  5. BRAF_ALT  - BRAF V600E, KIAA1549-BRAF fusion, other MAPK LGG/GNG (supplementary)

Excluded from main: LGG (other), ependymoma, medulloblastoma, ATRT, embryonal, chordoma, etc.
"""
import pandas as pd
import numpy as np
from pathlib import Path

DATA = Path(r"C:\Users\scott\Downloads\Open PBTA\data")
OUT  = Path(r"C:\Users\scott\Downloads\Open PBTA\output")
OUT.mkdir(exist_ok=True)

H   = pd.read_csv(DATA / "histologies.tsv", sep="\t", low_memory=False, dtype=str)
FUS = pd.read_csv(DATA / "fusion-putative-oncogenic.tsv", sep="\t", low_memory=False, dtype=str)
IND = pd.read_csv(DATA / "independent-specimens.rnaseqpanel.primary-plus.tsv", sep="\t", low_memory=False, dtype=str)

# -------- 1. Base filter: PBTA + RNA-Seq + Tumor --------
rna = H[
    (H["cohort"] == "PBTA") &
    (H["experimental_strategy"] == "RNA-Seq") &
    (H["sample_type"] == "Tumor")
].copy()

# Numeric conversions
rna["age_at_diagnosis_days"] = pd.to_numeric(rna["age_at_diagnosis_days"], errors="coerce")
rna["age_years"]   = rna["age_at_diagnosis_days"] / 365.25
rna["age_months"]  = rna["age_at_diagnosis_days"] / 30.4375
rna["OS_days_num"] = pd.to_numeric(rna["OS_days"], errors="coerce")

print(f"[Step 1] PBTA RNA-Seq tumor specimens: {len(rna):,}")

# -------- 2. Diagnosis grouping --------
hd = rna["harmonized_diagnosis"].fillna("")
ms = rna["molecular_subtype"].fillna("")
pd_dx = rna["pathology_diagnosis"].fillna("")

# DMG H3 K27 (K28 in IUPAC) – includes K28-altered, with or without TP53
is_dmg_k27 = (
    hd.str.contains("Diffuse midline glioma, H3 K28-altered", case=False, na=False) |
    ms.str.startswith("DMG, H3 K28", na=False)
)

# DHG H3 G34 (G35 in IUPAC)
is_dhg_g34 = (
    hd.str.contains("Diffuse hemispheric glioma, H3 G35-mutant", case=False, na=False) |
    ms.str.startswith("DHG, H3 G35", na=False)
)

# Pediatric HGG, H3/IDH-WT (exclude IDH-mutant adult-type)
is_phgg_wt = (
    hd.str.contains("High-grade glioma, IDH-wildtype and H3-wildtype", case=False, na=False) |
    ms.str.startswith("HGG, H3 wildtype", na=False)
)

# Infant-type hemispheric glioma: operational definition
#   age at dx < 18 months + hemispheric + RTK/BRAF fusion (NTRK/ALK/ROS1/MET/BRAF)
hemispheric_terms = ["Frontal", "Parietal", "Temporal", "Occipital", "Hemispheric"]
is_hemi = rna["CNS_region"].fillna("").isin(["Hemispheric"]) | \
          rna["primary_site"].fillna("").str.contains("|".join(hemispheric_terms), case=False, na=False)

# Build per-participant fusion calls (for IHG and BRAF_ALT classification)
FUS["FusionName"] = FUS["FusionName"].fillna("")
ihg_fusion_pattern = r"(NTRK[123]|ALK|ROS1|MET|^BRAF--|--BRAF)"
fus_ihg = FUS[FUS["FusionName"].str.contains(ihg_fusion_pattern, regex=True, na=False)]
ihg_fusion_samples = set(fus_ihg["Sample"].dropna().unique())

braf_fusion = FUS[FUS["FusionName"].str.contains(r"BRAF", regex=True, na=False)]
braf_fusion_samples = set(braf_fusion["Sample"].dropna().unique())

# Match by Kids_First_Biospecimen_ID OR sample_id (fusion uses biospecimen IDs typically)
rna["has_ihg_fusion"]  = rna["Kids_First_Biospecimen_ID"].isin(ihg_fusion_samples) | \
                        rna["sample_id"].isin(ihg_fusion_samples)
rna["has_braf_fusion"] = rna["Kids_First_Biospecimen_ID"].isin(braf_fusion_samples) | \
                        rna["sample_id"].isin(braf_fusion_samples)

# IHG definition: age <18mo + hemispheric + (RTK/BRAF fusion OR molecular_subtype hints)
is_ihg = (
    (rna["age_months"] < 18) &
    is_hemi &
    (rna["has_ihg_fusion"] | ms.str.contains("RTK", case=False, na=False))
)

# Secondary: BRAF-altered / MAPK pediatric gliomas (LGG + GNG + PXA, BRAF/MAPK)
is_braf_alt = (
    ms.str.contains("BRAF V600E", case=False, na=False) |
    ms.str.contains("KIAA1549-BRAF", case=False, na=False) |
    ms.str.contains("other MAPK", case=False, na=False) |
    ms.str.startswith("LGG, RTK", na=False) |
    ms.str.startswith("LGG, FGFR", na=False) |
    ms.str.startswith("HGG, PXA", na=False) |
    rna["has_braf_fusion"]
) & ~(is_dmg_k27 | is_dhg_g34 | is_phgg_wt | is_ihg)

# Assign cohort_group with priority: IHG > DMG_K27 > DHG_G34 > pHGG_WT > BRAF_ALT
rna["cohort_group"] = "Other"
rna.loc[is_braf_alt, "cohort_group"] = "BRAF_ALT"
rna.loc[is_phgg_wt,  "cohort_group"] = "pHGG_WT"
rna.loc[is_dhg_g34,  "cohort_group"] = "DHG_G34"
rna.loc[is_dmg_k27,  "cohort_group"] = "DMG_K27"
rna.loc[is_ihg,      "cohort_group"] = "IHG"

# -------- 3. Midline vs Hemispheric classification --------
midline_regions = ["Mixed", "Midline", "Suprasellar", "Optic pathway", "Pineal", "Ventricles"]
midline_sites   = ["Brain Stem", "Pons", "Thalamus", "Spinal Cord", "Spine", "Hypothalamus", "Brainstem", "Diencephalon", "Midbrain", "Medulla"]

def classify_location(row):
    cns = row["CNS_region"] if isinstance(row["CNS_region"], str) else ""
    ps  = row["primary_site"] if isinstance(row["primary_site"], str) else ""
    if cns == "Hemispheric": return "Hemispheric"
    if any(t in cns for t in midline_regions): return "Midline"
    if any(t in ps for t in midline_sites): return "Midline"
    if any(t in ps for t in hemispheric_terms): return "Hemispheric"
    if "Cerebell" in ps or "Posterior Fossa" in cns: return "Posterior_fossa"
    return "Other/NA"

rna["location_class"] = rna.apply(classify_location, axis=1)

# -------- 4. Age developmental group --------
def age_group(months):
    if pd.isna(months): return "NA"
    if months < 18: return "infant"   # <18mo
    if months < 60: return "young_child"  # 18mo–5y
    if months < 144: return "child"   # 5–12y
    return "adolescent"               # 12–21y
rna["age_dev_group"] = rna["age_months"].apply(age_group)

# -------- 5. Survival availability --------
rna["has_survival"] = rna["OS_days_num"].notna() & rna["OS_status"].notna()

# -------- 6. Independent-specimen flag --------
indep_ids = set(IND["Kids_First_Biospecimen_ID"].unique())
rna["is_independent_primary_plus"] = rna["Kids_First_Biospecimen_ID"].isin(indep_ids)

# -------- 7. Compile annotation table --------
keep_cols = [
    "Kids_First_Biospecimen_ID", "Kids_First_Participant_ID", "sample_id", "aliquot_id",
    "cohort_group", "pathology_diagnosis", "harmonized_diagnosis", "molecular_subtype",
    "broad_histology", "primary_site", "CNS_region", "location_class",
    "tumor_descriptor", "RNA_library",
    "reported_gender", "age_at_diagnosis_days", "age_years", "age_months", "age_dev_group",
    "OS_days", "OS_status", "EFS_days", "EFS_event_type", "has_survival",
    "has_braf_fusion", "has_ihg_fusion",
    "is_independent_primary_plus",
]
ann = rna[keep_cols].copy()

# -------- 8. Save --------
ann.to_csv(OUT / "cohort_annotation_full.tsv", sep="\t", index=False)

MAIN = ["DMG_K27", "DHG_G34", "pHGG_WT", "IHG"]
main_cohort = ann[ann["cohort_group"].isin(MAIN)].copy()
main_cohort.to_csv(OUT / "cohort_main.tsv", sep="\t", index=False)

secondary = ann[ann["cohort_group"] == "BRAF_ALT"].copy()
secondary.to_csv(OUT / "cohort_secondary_BRAF.tsv", sep="\t", index=False)

# -------- 9. Summary --------
print("\n=== Cohort group breakdown (all PBTA RNA-Seq tumors) ===")
print(ann["cohort_group"].value_counts())

print("\n=== Main cohort (n biospecimens) ===")
print(main_cohort["cohort_group"].value_counts())
print(f"Main cohort unique patients: {main_cohort['Kids_First_Participant_ID'].nunique()}")

print("\n=== Main cohort × location_class ===")
print(pd.crosstab(main_cohort["cohort_group"], main_cohort["location_class"]))

print("\n=== Main cohort × age_dev_group ===")
print(pd.crosstab(main_cohort["cohort_group"], main_cohort["age_dev_group"]))

print("\n=== Survival data availability (main cohort) ===")
print(main_cohort.groupby("cohort_group")["has_survival"].agg(["sum", "count"]))

print("\n=== Independent primary-plus flag (main cohort) ===")
print(main_cohort.groupby("cohort_group")["is_independent_primary_plus"].agg(["sum", "count"]))

print(f"\n=== Sex distribution (main cohort) ===")
print(pd.crosstab(main_cohort["cohort_group"], main_cohort["reported_gender"]))

# Median age per group
print("\n=== Age (years) - median [IQR] per group ===")
for g in MAIN:
    sub = main_cohort[main_cohort["cohort_group"] == g]["age_years"].dropna()
    if len(sub):
        print(f"  {g}: n={len(sub)}, median={sub.median():.1f}, IQR=[{sub.quantile(0.25):.1f}, {sub.quantile(0.75):.1f}], range=[{sub.min():.1f}, {sub.max():.1f}]")

print(f"\nFiles written to: {OUT}")
print("  - cohort_annotation_full.tsv (every PBTA RNA-Seq tumor with cohort_group tag)")
print("  - cohort_main.tsv (DMG_K27 + DHG_G34 + pHGG_WT + IHG only)")
print("  - cohort_secondary_BRAF.tsv (BRAF-altered/MAPK secondary cohort)")
