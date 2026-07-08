# Re-annotate ecotypes based on observed data pattern
# Key finding: lymphocyte-inflamed and myeloid signals co-occur in pediatric glioma
# Re-label cluster 2 (intermediate) and cluster 3 (pan-immune inflamed)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
})

ROOT <- "C:/Users/scott/Downloads/Open PBTA"
OUT  <- file.path(ROOT, "output")
FIG  <- file.path(OUT, "figs_clustering")

# -------- Load k=3 cluster assignments --------
ccp3 <- readRDS(file.path(OUT, "ccp_main_k3.rds"))
clus <- ccp3$class
cm   <- ccp3$consensus

clust_means <- fread(file.path(OUT, "ecotype_k3_feature_means.tsv"))

cat("Cluster sizes:\n")
print(table(clus))

cat("\nKey feature means per cluster:\n")
key_feats <- c("SS_Microglia_Klemm2020", "SS_MDM_Klemm2020",
               "SS_T_Cell_Cytotoxicity", "SS_NK_Cell_Activity",
               "SS_M1_Macrophage", "SS_M2_Macrophage",
               "SS_MHC_Class_II", "SS_IFN_Gamma_Response",
               "SS_TGFb_Immunosuppression", "SS_MAPK_Activity")
print(clust_means[, c("cluster", intersect(key_feats, colnames(clust_means))), with = FALSE])

# -------- Manual annotation based on biology --------
# Cluster 1 (n=79):  uniformly LOW everything           -> Immune-desert
# Cluster 2 (n=157): intermediate, slightly desert-like -> Intermediate
# Cluster 3 (n=113): uniformly HIGH T-cell + macrophage + microglia + MHC + IFN -> Inflamed (mixed lymph+myeloid)
ecotype_map <- c("1" = "Immune-desert",
                 "2" = "Intermediate",
                 "3" = "Inflamed")

clus_df <- data.frame(
  Kids_First_Biospecimen_ID = names(clus),
  cluster = clus,
  ecotype = ecotype_map[as.character(clus)]
)
cat("\nFinal ecotype distribution:\n")
print(table(clus_df$ecotype))

# -------- Reload metadata --------
ann <- fread(file.path(OUT, "cohort_for_deconvolution.tsv"))
meta <- ann[match(clus_df$Kids_First_Biospecimen_ID, ann$Kids_First_Biospecimen_ID), ]
clus_df$cohort_group   <- meta$cohort_group
clus_df$location_class <- meta$location_class
clus_df$age_dev_group  <- meta$age_dev_group
clus_df$age_years      <- as.numeric(meta$age_years)
clus_df$reported_gender<- meta$reported_gender
clus_df$OS_days        <- as.numeric(meta$OS_days)
clus_df$OS_status      <- meta$OS_status

fwrite(clus_df, file.path(OUT, "ecotype_assignment_k3_annotated.tsv"))

# -------- Re-run association tests --------
cat("\n=== Ecotype x cohort_group ===\n")
tab1 <- table(clus_df$ecotype, clus_df$cohort_group)
print(tab1)
fisher_p1 <- tryCatch(
  fisher.test(tab1, simulate.p.value = TRUE, B = 10000)$p.value,
  error = function(e) NA)
cat(sprintf("Fisher p = %.2e\n", fisher_p1))

cat("\n=== Ecotype x location_class (excl. Ambiguous/Other) ===\n")
loc_ok <- clus_df[clus_df$location_class %in% c("Midline", "Hemispheric", "Posterior_fossa"), ]
tab2 <- table(loc_ok$ecotype, loc_ok$location_class)
print(tab2)
fisher_p2 <- tryCatch(
  fisher.test(tab2, simulate.p.value = TRUE, B = 10000)$p.value,
  error = function(e) NA)
cat(sprintf("Fisher p = %.2e\n", fisher_p2))

cat("\n=== Ecotype x age_dev_group ===\n")
tab3 <- table(clus_df$ecotype, clus_df$age_dev_group)
print(tab3)
fisher_p3 <- tryCatch(
  fisher.test(tab3, simulate.p.value = TRUE, B = 10000)$p.value,
  error = function(e) NA)
cat(sprintf("Fisher p = %.2e\n", fisher_p3))

# -------- Microglia/MDM ratio analysis per ecotype --------
cat("\n=== Microglia vs MDM ratio per ecotype ===\n")
ss <- fread(file.path(OUT, "ssGSEA_brain_immune_scores.tsv"))
mg_row  <- as.numeric(ss[ss$Signature == "Microglia_Klemm2020", -1, with = FALSE])
mdm_row <- as.numeric(ss[ss$Signature == "MDM_Klemm2020", -1, with = FALSE])
samples <- colnames(ss)[-1]
ratio_df <- data.frame(sample = samples,
                       microglia = mg_row, mdm = mdm_row,
                       mg_mdm_ratio = mg_row - mdm_row)
ratio_df$ecotype <- clus_df$ecotype[match(ratio_df$sample, clus_df$Kids_First_Biospecimen_ID)]
ratio_df$cohort_group <- meta$cohort_group[match(ratio_df$sample, clus_df$Kids_First_Biospecimen_ID)]
ratio_df <- ratio_df[!is.na(ratio_df$ecotype), ]

cat("\nMean microglia-MDM z-score difference (positive = microglia-tilted):\n")
print(aggregate(mg_mdm_ratio ~ ecotype, data = ratio_df, FUN = function(x) round(mean(x), 3)))

cat("\nMean by cohort_group:\n")
print(aggregate(mg_mdm_ratio ~ cohort_group + ecotype, data = ratio_df,
                FUN = function(x) round(mean(x), 3)))

fwrite(ratio_df, file.path(OUT, "microglia_mdm_ratio_by_ecotype.tsv"))

# -------- Stats summary --------
stats_df <- data.frame(
  comparison = c("ecotype_x_cohort_group", "ecotype_x_location_class", "ecotype_x_age_dev_group"),
  fisher_p_simulated = c(fisher_p1, fisher_p2, fisher_p3),
  N = c(nrow(clus_df), nrow(loc_ok), nrow(clus_df))
)
fwrite(stats_df, file.path(OUT, "ecotype_association_pvalues.tsv"))

# -------- Re-generate visualizations with correct annotation --------
cat("\n=== Re-generating figures with correct ecotype labels ===\n")

ec_colors <- c("Inflamed" = "#3B82F6",
               "Intermediate" = "#F59E0B",
               "Immune-desert" = "#9CA3AF")
cg_colors <- c("DMG_K27" = "#1F77B4", "DHG_G34" = "#FF7F0E",
               "pHGG_WT" = "#2CA02C", "IHG" = "#9467BD")
loc_colors <- c("Midline" = "#D62728", "Hemispheric" = "#17BECF",
                "Posterior_fossa" = "#8C564B", "Ambiguous" = "#BCBD22",
                "Other/NA" = "#7F7F7F")

# Reload feature matrix
feat <- fread(file.path(OUT, "clustering_feature_matrix_z.tsv"))
feat_z <- as.matrix(feat[, -1, with = FALSE])
rownames(feat_z) <- feat$feature
mat_main <- feat_z[, clus_df$Kids_First_Biospecimen_ID]

# Heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
col_anno <- HeatmapAnnotation(
  Ecotype  = clus_df$ecotype,
  Cohort   = clus_df$cohort_group,
  Location = clus_df$location_class,
  col = list(Ecotype = ec_colors, Cohort = cg_colors, Location = loc_colors),
  show_legend = TRUE, annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9)
)

png(file.path(FIG, "ecotype_heatmap_main.png"), width = 14, height = 10, units = "in", res = 200)
ht <- Heatmap(
  mat_main, name = "z-score", col = col_fun,
  top_annotation = col_anno,
  show_column_names = FALSE,
  column_split = factor(clus_df$ecotype, levels = c("Inflamed","Intermediate","Immune-desert")),
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  row_names_gp = gpar(fontsize = 8),
  row_split = ifelse(grepl("^QTS", rownames(mat_main)), "quanTIseq", "ssGSEA"),
  column_title_gp = gpar(fontsize = 11, fontface = "bold")
)
draw(ht, merge_legend = TRUE)
dev.off()
cat("  Wrote updated ecotype_heatmap_main.png\n")

# PCA
pca_df <- fread(file.path(OUT, "ecotype_pca_main.tsv"))
pca_df$ecotype <- clus_df$ecotype[match(pca_df$Kids_First_Biospecimen_ID, clus_df$Kids_First_Biospecimen_ID)]
fwrite(pca_df, file.path(OUT, "ecotype_pca_main.tsv"))
png(file.path(FIG, "ecotype_pca_main.png"), width = 7, height = 6, units = "in", res = 200)
plot(pca_df$PC1, pca_df$PC2, col = ec_colors[pca_df$ecotype], pch = 19, cex = 0.9,
     xlab = "PC1", ylab = "PC2",
     main = "PCA — Main cohort by ecotype (n=349)")
legend("topright", names(ec_colors), col = ec_colors, pch = 19, bty = "n", cex = 0.9)
dev.off()
cat("  Wrote updated ecotype_pca_main.png\n")

# UMAP
umap_df <- fread(file.path(OUT, "ecotype_umap_main.tsv"))
umap_df$ecotype <- clus_df$ecotype[match(umap_df$Kids_First_Biospecimen_ID, clus_df$Kids_First_Biospecimen_ID)]
fwrite(umap_df, file.path(OUT, "ecotype_umap_main.tsv"))
png(file.path(FIG, "ecotype_umap_main.png"), width = 7, height = 6, units = "in", res = 200)
plot(umap_df$UMAP1, umap_df$UMAP2, col = ec_colors[umap_df$ecotype], pch = 19, cex = 0.9,
     xlab = "UMAP1", ylab = "UMAP2",
     main = "UMAP — Main cohort by ecotype (n=349)")
legend("topright", names(ec_colors), col = ec_colors, pch = 19, bty = "n", cex = 0.9)
dev.off()
cat("  Wrote updated ecotype_umap_main.png\n")

# Stacked barplot — ecotype by cohort_group
png(file.path(FIG, "ecotype_distribution_by_cohort.png"), width = 8, height = 5, units = "in", res = 200)
prop_tab <- prop.table(table(clus_df$cohort_group, clus_df$ecotype), margin = 1)
ec_order <- c("Inflamed", "Intermediate", "Immune-desert")
prop_tab <- prop_tab[, ec_order]
par(mar = c(5, 4, 4, 8), xpd = TRUE)
barplot(t(prop_tab), col = ec_colors[colnames(prop_tab)],
        las = 1, ylab = "Proportion",
        main = "Ecotype distribution across pediatric glioma subgroups")
legend("topright", inset = c(-0.25, 0), legend = colnames(prop_tab),
       fill = ec_colors[colnames(prop_tab)], bty = "n", cex = 0.85)
dev.off()
cat("  Wrote updated ecotype_distribution_by_cohort.png\n")

# Ecotype by location
png(file.path(FIG, "ecotype_distribution_by_location.png"), width = 8, height = 5, units = "in", res = 200)
prop_tab2 <- prop.table(table(loc_ok$location_class, loc_ok$ecotype), margin = 1)
prop_tab2 <- prop_tab2[, intersect(ec_order, colnames(prop_tab2))]
par(mar = c(5, 4, 4, 8), xpd = TRUE)
barplot(t(prop_tab2), col = ec_colors[colnames(prop_tab2)],
        las = 1, ylab = "Proportion",
        main = "Ecotype distribution by anatomical location")
legend("topright", inset = c(-0.25, 0), legend = colnames(prop_tab2),
       fill = ec_colors[colnames(prop_tab2)], bty = "n", cex = 0.85)
dev.off()
cat("  Wrote updated ecotype_distribution_by_location.png\n")

# Microglia-MDM ratio boxplot by cohort within ecotype
png(file.path(FIG, "microglia_mdm_ratio_by_ecotype.png"), width = 9, height = 5, units = "in", res = 200)
boxplot(mg_mdm_ratio ~ cohort_group + ecotype, data = ratio_df,
        las = 2, col = "lightblue", outline = TRUE,
        main = "Microglia − MDM ssGSEA score difference",
        ylab = "Microglia (Klemm) − MDM (Klemm) z-score",
        xlab = "")
abline(h = 0, col = "red", lty = 2)
dev.off()
cat("  Wrote microglia_mdm_ratio_by_ecotype.png\n")

cat("\nDone.\n")
