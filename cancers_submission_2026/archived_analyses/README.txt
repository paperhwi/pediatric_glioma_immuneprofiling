DEPOSITED ANALYSES — not submitted with the Cancers manuscript
https://github.com/paperhwi/pediatric_glioma_immuneprofiling

These analyses were part of the earlier revision package for the previous journal. They
are retained here, with the code that produces them, so that the record is complete and
so that a reader who wants them can find them. They are not cited in the manuscript.

figures/
  FigureS1, FigureS2, FigureS7A-C     same-resource pooled n = 702 cross-method
                                      concordance, the cross-method consensus heatmap,
                                      and the pooled re-clustering sensitivity analysis.
                                      Internal sensitivity analyses on the same resource,
                                      not independent external validation.
  FigureS14A-C                        same-resource BRAF/RTK-altered cohort-family
                                      projection and pathway-score sensitivity.
  FigureS15                           immune-checkpoint-related transcript expression by
                                      ecotype. Descriptive; does not establish
                                      checkpoint-inhibitor sensitivity.
  FigureS16                           CAR-T target-related transcript expression by
                                      ecotype (B4GALNT1, ST8SIA1, IL13RA2, ERBB2, CD276).
                                      Transcript abundance does not measure cell-surface
                                      antigen density or treatment response.
  FigureS17                           audit copy of the molecular-group-adjusted
                                      differential-expression analysis that became main
                                      Figure 6A-C. Retained for traceability only.
  FigureS18_small_strata_KM           Kaplan-Meier curves for the DHG, H3 G34-mutant
                                      (n = 15) and infant-type hemispheric (n = 14)
                                      strata. Descriptive only; the nominal DHG log-rank
                                      P = 0.0035 is not evidence of subtype
                                      generalizability.
  FigureS10A                          subtype-stratified Kaplan-Meier curves including
                                      the very small strata; superseded by the submitted
                                      Supplementary Figure S11.
  FigureS19_HSPC_programs             hematopoietic stem and progenitor programs across
                                      the ecotypes (Azimuth 2023 bone-marrow reference,
                                      gseapy ssGSEA).

tables/
  T4_HSPC_ssGSEA_scores.tsv           per-sample scores for the progenitor modules
  T4_KW_by_ecotype.tsv                Kruskal-Wallis tests of those modules by ecotype
  T4_Dunn_posthoc.tsv                 Dunn post-hoc comparisons
  T4_residualised_specificity.tsv     effect sizes before and after residualizing on
                                      mature CD14 monocyte and macrophage content

The principal finding of the progenitor analysis is a qualified negative: committed
granulocyte-monocyte progenitor programs vary across ecotypes, but the signal falls from
epsilon-squared 0.205 to 0.007 after residualizing on mature myeloid content, and
uncommitted hematopoietic stem cell programs do not distinguish Lymphocyte-inflamed from
Myeloid-dominant (Dunn z = 0.29, P = 0.775). A distinct myeloid-progenitor or stem-like
subset can be neither supported nor excluded from bulk data.

The code and clean, output-stripped notebooks for these analyses are deposited in the adjacent
analysis_code/revision_scripts and analysis_code/notebooks directories.
