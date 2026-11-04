# Figure 5A target correction

The curated signature has AUROC **0.925** for the annotated CD4/CD8 T-cell compartment versus myeloid cells; the four-specimen range is **0.927–0.976**. NK cells are not separately resolved in this reference.

The original notebook excludes unclear cells and retains only Myeloid, CD4 and CD8 annotations. The retained H5AD confirms 13,920 Myeloid, 2,385 CD4 and 2,314 CD8 cells (18,619 total); 93 unclear cells are excluded. Consequently the revised binary labels are exactly the complement of the original labels. For unchanged scores, AUROC of the complementary target is exactly one minus the original AUROC, including ties. The stored source table records three decimals, so no extra precision or unavailable per-cell scores are inferred. The endpoints of the per-specimen range reverse under this transformation.

Only the curated signature target changes. The other 15 directional AUROCs and all genes/scores are unchanged. Current targets give 16/16 lineage-concordant signatures and a median AUROC of 0.879. The original author/year attribution is not reinstated. Figure 5A has no numbers inside the bars.
