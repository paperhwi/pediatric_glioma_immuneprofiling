Figure legends

Figure 1 Study design and cohort

A, analysis workflow from clinical annotation, expression and immune features to consensus clustering. B, derivation of the final cohort: 343 pediatric-type diffuse high-grade gliomas, one RNA specimen per participant, after exclusion of six non-eligible specimens. The three modal immune states contain 122 Immune-high, 144 Immune-intermediate and 77 Immune-low tumours.

Figure 2 LM22 quality and cross-method comparison

A, distributions of LM22 deconvolution quality measures; 88 of 343 tumours (25.7%) have permutation P < 0.05. B, cross-method Spearman concordance for matched immune-cell axes; four of 12 comparisons pass Benjamini–Hochberg q < 0.05. C, method-specific cell estimates across the three immune states. The methods estimate related but non-identical quantities.

Figure 3 Immune-state cluster selection

A–E, PAC, silhouette, consensus CDF, gap statistic and minimum cluster size across candidate k values on the 46-feature matrix. PAC and silhouette favor k = 2; the one-standard-error gap rule selects k = 3. F, projection of samples; G, discriminating immune features. Final three-state labels are the modal labels across seeds 42–66; median pairwise adjusted Rand index is 0.94.

Figure 4 Immune states by location and molecular group

A, immune-state fractions by anatomical class (n = 326; Monte Carlo chi-square P = 0.0063). B, immune states by integrated molecular group (n = 343; P = 0.0001). C, midline versus hemispheric location within H3/IDH-wildtype pHGG (n = 86; P = 0.5061). D, sequential PERMANOVA on the immune-feature matrix (n = 253); partial R-squared values depend on term order.

Figure 5 Single-cell and spatial support

A, directional immune-signature discrimination in annotated single cells. B, immune-cell expression assignment for the top 300 Immune-high versus Immune-low genes (294/300, 98.0%; matched-random mean 69.7%; permutation P = 0.005). C–E, Visium spots from three diffuse midline glioma sections (577, 856 and 217 spots) mapped to bulk immune-state centroids and evaluated for spatial and cell-mixture patterns. Spot assignments are mixed-cell estimates.

Figure 6 Adjusted differential expression

A–C, molecular-group-adjusted gene-wise comparisons of log2(TPM + 1) across immune states (18,564 genes). Colored genes satisfy BH q < 0.05 and adjusted absolute log2 fold change > 1. D, overlap of genes up-regulated versus Immune-low: 404 shared genes, comprising 93.7% of the Immune-intermediate set. E–F, effect-size relationship across contrasts (orthogonal slope 1.90; Pearson r = 0.78; median shared-gene ratio 1.51).

Figure 7 Overall survival

A, Kaplan–Meier overall survival by immune state for 245 tumours (206 deaths), with numbers at risk; log-rank P = 0.058. B, Cox model adjusted for molecular group, age at diagnosis and sex. Immune-high versus Immune-low HR 0.62 (95% CI 0.42–0.91; P = 0.016); Immune-intermediate versus Immune-low HR 0.84 (0.57–1.22; P = 0.35). Sensitivity analyses are in the supplementary tables and figure.

Figure S1 survival analyses

Sensitivity analyses of the prespecified survival model, within-subtype associations and incremental value. The primary set includes 245 tumours and 206 deaths. Excluding autopsy specimens gives n = 189, 150 deaths and Immune-high versus Immune-low HR 0.71 (95% CI 0.44–1.14; P = 0.16). Interpret these as observational associations.

Figure S2 immune-state structure

Robustness analyses of the immune-state assignment, including within-group reclustering, k = 2 versus k = 3, RNA-library effects and comparisons with published tumour-microenvironment frameworks. The final labels are modal assignments over seeds 42–66 in the 343-tumour cohort.

Figure S3 immune features

Immune-feature and programme associations with survival and immune state. Analyses include adjustment for molecular group, age and sex, with additional immune-axis models where indicated. Reported q values use Benjamini–Hochberg correction; tumour-cell programme findings are exploratory.

Figure S4 ScPCA axis linkage

Single-cell donor-level comparison with the n = 343 bulk immune axis. The harmonised classifier has cross-validated accuracy 299/343 (87.2%); projected donor states are four Immune-low, seven Immune-intermediate and six Immune-high. The bulk axis correlates with single-cell immune-cell fraction (Spearman rho = 0.78; permutation P = 0.0003). Cell-specific secondary tests remain exploratory. Ordered donor-state tests use Immune-low < Immune-intermediate < Immune-high.

Figure S5 Within cohort standardization and assignment consistency

Transition counts between primary retained-input labels (rows) and labels after feature re-standardization within the 343 eligible participants (columns). Both partitions use k = 3, 25 seeds 42–66, 500 subsamples per seed, 80% item resampling, KMeans n_init = 10, average linkage of 1 minus consensus, T-cell-cytotoxicity-based names, and per-participant modal aggregation. In total, 334/343 labels were unchanged (97.4%; ARI = 0.915). Primary counts are 122/144/77 and re-standardized counts are 116/149/78 for Immune-high/Immune-intermediate/Immune-low. Counts and ARI are descriptive cohort measures; full seed and assignment records are in Supplementary Table S6.

Figure S6 Continuous immune axis survival association with RNA library control

Unpenalized Cox HRs per standard-deviation increase in the fixed n = 343 immune axis. All models adjust for molecular group, age, and sex. Base, additional categorical RNA-library adjustment (four levels, poly-A stranded reference), and RNA-library-stratified models are shown for all survival-evaluable patients (n = 245, 206 deaths) and after autopsy exclusion (n = 189, 150 deaths). Points show HRs and whiskers are Wald 95% confidence intervals; P values are two-sided. The vertical line denotes HR = 1. Additional library adjustment attenuated the full-cohort axis association to HR 0.861 (95% CI, 0.726–1.022; P = 0.088), and combined autopsy exclusion and library adjustment gave HR 0.942 (95% CI, 0.766–1.160; P = 0.575). Exact estimates and diagnostics are in Supplementary Table S6.
