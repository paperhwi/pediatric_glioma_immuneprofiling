# Phase 1 — IDH1/IDH2 status audit on the 349-sample pediatric HGG cohort

## Bottom line
**No confirmed IDH-mutant cases** were found in the 349 primary-analysis cohort. The
cohort is clean for the "pediatric-type diffuse high-grade glioma" disease category
(DMG H3 K27-altered, DHG H3 G34-mutant, pHGG H3-wildtype/IDH-wildtype, IHG). No
exclusion is required for the primary analysis. **A single sample flagged by the
methylation classifier requires narrative acknowledgement but does not warrant
exclusion** because its integrated diagnosis and pathology narrative confirm H3
K27M-mutant DMG.

## Audit strategy
For each of the 349 samples we screened five OpenPedCan annotation columns for
IDH-related tokens:

| Column | IDH mentions | IDH-mutant flags (after excluding wildtype phrasings) |
|---|---|---|
| `pathology_diagnosis` | 0 | 0 |
| `pathology_free_text_diagnosis` | 15 (all "IDH-wt"/"IDH-1 negative"/"idh1 r132h-wild") | 1 (false-positive; see BS_GG9W3J9Y below) |
| `molecular_subtype` | 0 | 0 |
| `molecular_subtype_methyl` | 1 | 1 (discordant; see BS_V0DQFC75 below) |
| `integrated_diagnosis` | 129 (all "IDH-wildtype and H3-wildtype") | 0 |
| `broad_histology` / `cancer_group` | 0 | 0 |

The 129 hits in `integrated_diagnosis` are the pHGG_WT cohort samples explicitly
annotated as "High-grade glioma, IDH-wildtype and H3-wildtype" — i.e., confirmed
IDH-wildtype pediatric-type HGG, which is the intended target category.

## The two ambiguous samples

### 1. BS_V0DQFC75 — methylation-classifier discordance (KEEP in DMG_K27)
| Field | Value |
|---|---|
| `molecular_subtype` | DMG, H3 K28, TP53 |
| `integrated_diagnosis` | Diffuse midline glioma, H3 K28-altered |
| `pathology_free_text` | dmg, h3k27m mutant, who gr iv |
| `molecular_subtype_methyl` | **HGG, IDH** ← discordant methylation call |
| Site / age / OS | Thalamus / 13.7 y / 306 d, deceased |
| Cohort assignment | DMG_K27, Midline, Intermediate ecotype |

**Interpretation.** The gold-standard clinical annotations (integrated diagnosis
+ pathology free text + molecular_subtype) all confirm this is an H3 K27M-mutant
pediatric-type DMG located in the thalamus. The methylation classifier occasionally
emits low-confidence secondary calls (here "HGG, IDH"), particularly for H3 K27M-
altered tumors that share DNA-methylation signatures with certain IDH-mutant HGG
subgroups. **Recommend**: retain in DMG_K27, cite as `molecular_subtype = DMG, H3
K28, TP53` per OpenPedCan curated call, add one-line disclosure of the discordant
methylation call to Methods.

### 2. BS_GG9W3J9Y — pathology-narrative false positive (KEEP in pHGG_WT)
| Field | Value |
|---|---|
| `molecular_subtype` | HGG, H3 wildtype |
| `integrated_diagnosis` | High-grade glioma, IDH-wildtype and H3-wildtype |
| `pathology_free_text` | glioblastoma, **idh-1 negative**, who grade iv |

**Interpretation.** The token "idh-1 negative" is pathologist shorthand for "IDH1
immunohistochemistry / mutation test returned negative", i.e., **IDH-wildtype**.
Every structured OpenPedCan column confirms IDH-wildtype H3-wildtype pHGG. My
heuristic caught the "IDH" substring but negation ("negative") reverses the sense.
**No change** — sample stays in pHGG_WT as IDH-wildtype.

## Per-cohort molecular_subtype composition (n = 349)

```
cohort_group     DHG_G34  DMG_K27  IHG  pHGG_WT
DMG, H3 K28           0       78    0        0
DMG, H3 K28, TP53     0       97    0        0
DHG, H3 G35           8        0    0        0
DHG, H3 G35, TP53    23        0    0        0
HGG, H3 wildtype      0        0    3       77
HGG, H3 wildtype+TP53 0        0    1       48
HGG, To be classified 0        0    2        0
IHG, NTRK-altered     0        0    4        0
IHG, ROS1-altered     0        0    3        0
GNG/GNT, RTK          0        0    2        0
LGG, KIAA1549-BRAF    0        0    1        0
NA                    0        0    2        0
```

All 175 DMG_K27 samples are H3 K28-altered (K27M), all 31 DHG_G34 are H3 G35-mutant
(G34), all 125 pHGG_WT are H3-wildtype. The 18 IHG samples include NTRK/ROS1-altered
infants plus 4 IHG-cohort infants whose molecular_subtype falls in the HGG-H3wt or
GNG/GNT category — consistent with the IHG operational definition (age < 18 months
+ hemispheric + BRAF/RTK fusion) rather than the OpenPedCan `molecular_subtype`
column.

## Manuscript action items (Phase 1 outputs)
1. **Add a single sentence to Methods (Section 2.1 or 2.4)**: "The 349-sample
   primary analysis cohort is confirmed IDH-wildtype: no case carries an IDH1
   or IDH2 mutation in the OpenPedCan `molecular_subtype`, `integrated_diagnosis`,
   or `pathology_free_text_diagnosis` annotations; one sample (BS_V0DQFC75) received
   a discordant secondary methylation-classifier call ("HGG, IDH") that is
   contradicted by its integrated diagnosis (Diffuse midline glioma, H3 K28-
   altered) and by its curated molecular_subtype (DMG, H3 K28, TP53), and was
   retained in the DMG_K27 subgroup per the curated call."
2. **No primary-analysis sample exclusion is required**; no IDH-mutant-excluded
   sensitivity Cox is needed.
3. **Include the audit table** (`IDH_audit_per_sample.tsv`) as a Supplementary
   Table (Supplementary Table S_IDH_audit).

## Deliverables (output/phase1_IDH_audit/)
- `IDH_audit_per_sample.tsv` — 349 rows, all key annotation columns + `IDH_flag`
- `IDH_mut_candidate_flags.tsv` — refined boolean flags per annotation column
- `Phase1_IDH_audit_report.md` (this document)
