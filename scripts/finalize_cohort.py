"""
Finalize main cohort + secondary BRAF cohort and generate Table 1.

Decisions (from PI confirmation):
  - Independent-primary-plus filter applied to main analysis
  - pHGG_WT: age < 21 yr cap
  - Pineal / Suprasellar / Optic pathway / Ventricles -> 'Ambiguous'
  - BRAF_ALT split: HGG/PXA/infant-fusion vs LGG
"""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(r"C:\Users\scott\Downloads\Open PBTA")
OUT  = ROOT / "output"
OUT.mkdir(exist_ok=True)

ann = pd.read_csv(OUT / "cohort_annotation_full.tsv", sep="\t", low_memory=False)

# Numeric
ann["age_years"]  = pd.to_numeric(ann["age_years"], errors="coerce")
ann["age_months"] = pd.to_numeric(ann["age_months"], errors="coerce")
ann["is_independent_primary_plus"] = ann["is_independent_primary_plus"].astype(str).str.lower() == "true"
ann["has_braf_fusion"] = ann["has_braf_fusion"].astype(str).str.lower() == "true"
ann["has_ihg_fusion"]  = ann["has_ihg_fusion"].astype(str).str.lower() == "true"
ann["has_survival"]    = ann["has_survival"].astype(str).str.lower() == "true"

MAIN = ["DMG_K27", "DHG_G34", "pHGG_WT", "IHG"]

# -------- 1. Ambiguous location reclassification --------
ambiguous_sites = ["Pineal", "Suprasellar", "Optic pathway", "Ventricles",
                   "Sellar", "Hypophyseal", "Pituitary", "Septum"]
ambiguous_regions = ["Mixed", "Optic pathway", "Suprasellar", "Pineal", "Ventricles"]

def reclassify_loc(row):
    cns = row["CNS_region"] if isinstance(row["CNS_region"], str) else ""
    ps  = row["primary_site"] if isinstance(row["primary_site"], str) else ""
    if any(t in cns for t in ambiguous_regions): return "Ambiguous"
    if any(t in ps for t in ambiguous_sites): return "Ambiguous"
    return row["location_class"]

ann["location_class"] = ann.apply(reclassify_loc, axis=1)

# -------- 2. Main cohort filters --------
main = ann[ann["cohort_group"].isin(MAIN)].copy()
main_n_pre = len(main)
print(f"[Pre-filter] Main cohort biospecimens: {main_n_pre}")

# Apply independent filter
main = main[main["is_independent_primary_plus"]].copy()
print(f"[Step A] After independent-primary-plus filter: {len(main)}")

# Apply pHGG_WT age cap (<21yr)
before = len(main)
mask_drop = (main["cohort_group"] == "pHGG_WT") & (main["age_years"] >= 21)
main = main[~mask_drop].copy()
print(f"[Step B] After pHGG_WT age<21 cap: {len(main)}  (dropped {before - len(main)} adult-age pHGG_WT)")

# -------- 3. Secondary BRAF_ALT split --------
braf = ann[ann["cohort_group"] == "BRAF_ALT"].copy().reset_index(drop=True)

ms = braf["molecular_subtype"].fillna("")
hd = braf["harmonized_diagnosis"].fillna("")
pd_dx = braf["pathology_diagnosis"].fillna("")

# BRAF_ALT_main: HGG / PXA / infant-type fusion-driven
is_hgg = (
    ms.str.startswith("HGG") |
    hd.str.contains("high-grade glioma", case=False, na=False) |
    pd_dx.str.contains("high-grade", case=False, na=False)
)
is_pxa = ms.str.contains("PXA", case=False, na=False) | \
         hd.str.contains("pleomorphic xanthoastrocytoma", case=False, na=False) | \
         pd_dx.str.contains("pleomorphic xanthoastrocytoma", case=False, na=False)
is_infant_fusion = (braf["age_months"].astype(float) < 18) & (braf["has_braf_fusion"] | braf["has_ihg_fusion"])

is_braf_main = is_hgg | is_pxa | is_infant_fusion
braf["braf_subgroup"] = np.where(is_braf_main, "BRAF_ALT_main", "BRAF_ALT_LGG")

braf_main = braf[braf["braf_subgroup"] == "BRAF_ALT_main"].copy()
braf_lgg  = braf[braf["braf_subgroup"] == "BRAF_ALT_LGG"].copy()

# Apply independent filter to BRAF as well
braf_main_indep = braf_main[braf_main["is_independent_primary_plus"]].copy()
braf_lgg_indep  = braf_lgg[braf_lgg["is_independent_primary_plus"]].copy()

print(f"\n[BRAF split]")
print(f"  BRAF_ALT_main (HGG/PXA/infant-fusion): {len(braf_main)}  (independent: {len(braf_main_indep)})")
print(f"  BRAF_ALT_LGG (LGG-only):              {len(braf_lgg)}   (independent: {len(braf_lgg_indep)})")

# -------- 4. Save finalized cohorts --------
main.to_csv(OUT / "cohort_main_final.tsv", sep="\t", index=False)
braf_main_indep.to_csv(OUT / "cohort_BRAF_ALT_main.tsv", sep="\t", index=False)
braf_lgg_indep.to_csv(OUT / "cohort_BRAF_ALT_LGG.tsv", sep="\t", index=False)

# All samples that will be sent to CIBERSORTx (main + BRAF_ALT_main + BRAF_ALT_LGG)
deconv_input = pd.concat([main, braf_main_indep, braf_lgg_indep], ignore_index=True)
deconv_input.to_csv(OUT / "cohort_for_deconvolution.tsv", sep="\t", index=False)
print(f"\n[Deconvolution input] Total unique biospecimens: {deconv_input['Kids_First_Biospecimen_ID'].nunique()}")

# -------- 5. Summary breakdown --------
print(f"\n=== Final main cohort breakdown ===")
print(main["cohort_group"].value_counts())

print(f"\n=== Main × location_class ===")
print(pd.crosstab(main["cohort_group"], main["location_class"]))

print(f"\n=== Main × age_dev_group ===")
print(pd.crosstab(main["cohort_group"], main["age_dev_group"]))

# -------- 6. Table 1 — Clinical/molecular characteristics --------
def fmt_age(s):
    s = s.dropna()
    if len(s) == 0: return "—"
    return f"{s.median():.1f} ({s.quantile(0.25):.1f}-{s.quantile(0.75):.1f})"

def fmt_pct(n, total):
    return f"{n} ({100*n/total:.1f}%)" if total else "—"

rows = []
groups = [("DMG_K27", "DMG, H3 K27/K28-altered"),
          ("DHG_G34", "DHG, H3 G34/G35-mutant"),
          ("pHGG_WT", "pHGG, H3/IDH-wildtype"),
          ("IHG",     "Infant-type hemispheric")]

for g, label in groups:
    sub = main[main["cohort_group"] == g]
    n = len(sub)
    n_pt = sub["Kids_First_Participant_ID"].nunique()
    row = {
        "Cohort": label,
        "n biospecimen": n,
        "n unique patient": n_pt,
        "Age yr median (IQR)": fmt_age(sub["age_years"]),
        "Sex Male": fmt_pct((sub["reported_gender"] == "Male").sum(), n),
        "Sex Female": fmt_pct((sub["reported_gender"] == "Female").sum(), n),
        "Midline": fmt_pct((sub["location_class"] == "Midline").sum(), n),
        "Hemispheric": fmt_pct((sub["location_class"] == "Hemispheric").sum(), n),
        "Posterior fossa": fmt_pct((sub["location_class"] == "Posterior_fossa").sum(), n),
        "Ambiguous/Other": fmt_pct(sub["location_class"].isin(["Ambiguous", "Other/NA"]).sum(), n),
        "BRAF/RTK fusion+": fmt_pct((sub["has_braf_fusion"] | sub["has_ihg_fusion"]).sum(), n),
        "OS data available": fmt_pct(sub["has_survival"].sum(), n),
    }
    rows.append(row)

# Total row
n_total = len(main)
total_row = {
    "Cohort": "Total (Main)",
    "n biospecimen": n_total,
    "n unique patient": main["Kids_First_Participant_ID"].nunique(),
    "Age yr median (IQR)": fmt_age(main["age_years"]),
    "Sex Male": fmt_pct((main["reported_gender"] == "Male").sum(), n_total),
    "Sex Female": fmt_pct((main["reported_gender"] == "Female").sum(), n_total),
    "Midline": fmt_pct((main["location_class"] == "Midline").sum(), n_total),
    "Hemispheric": fmt_pct((main["location_class"] == "Hemispheric").sum(), n_total),
    "Posterior fossa": fmt_pct((main["location_class"] == "Posterior_fossa").sum(), n_total),
    "Ambiguous/Other": fmt_pct(main["location_class"].isin(["Ambiguous", "Other/NA"]).sum(), n_total),
    "BRAF/RTK fusion+": fmt_pct((main["has_braf_fusion"] | main["has_ihg_fusion"]).sum(), n_total),
    "OS data available": fmt_pct(main["has_survival"].sum(), n_total),
}
rows.append(total_row)

table1 = pd.DataFrame(rows)
table1.to_csv(OUT / "Table1_cohort_characteristics.tsv", sep="\t", index=False)

print(f"\n=== Table 1 - Cohort characteristics ===")
with pd.option_context("display.max_columns", None, "display.width", 200):
    print(table1.T.to_string())

# Also save BRAF-secondary summary
braf_rows = []
for g, df in [("BRAF_ALT_main (HGG/PXA/infant-fusion)", braf_main_indep),
              ("BRAF_ALT_LGG (LGG-only)", braf_lgg_indep)]:
    n = len(df)
    braf_rows.append({
        "Subgroup": g,
        "n biospecimen": n,
        "n unique patient": df["Kids_First_Participant_ID"].nunique(),
        "Age yr median (IQR)": fmt_age(df["age_years"]),
        "BRAF/RTK fusion+": fmt_pct((df["has_braf_fusion"] | df["has_ihg_fusion"]).sum(), n),
        "OS data available": fmt_pct(df["has_survival"].sum(), n),
    })
pd.DataFrame(braf_rows).to_csv(OUT / "Table1_supp_BRAF_secondary.tsv", sep="\t", index=False)

print(f"\nFiles written to: {OUT}")
print("  - cohort_main_final.tsv")
print("  - cohort_BRAF_ALT_main.tsv")
print("  - cohort_BRAF_ALT_LGG.tsv")
print("  - cohort_for_deconvolution.tsv  (CIBERSORTx 입력 sample 목록)")
print("  - Table1_cohort_characteristics.tsv")
print("  - Table1_supp_BRAF_secondary.tsv")
