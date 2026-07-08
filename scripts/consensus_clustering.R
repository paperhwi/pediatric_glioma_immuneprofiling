# Consensus clustering on quanTIseq + brain-tuned ssGSEA features
# Main analysis: independent-primary-plus cohort (n=349) — Main 4 groups only
# Secondary: include BRAF_ALT (702 total) as supplementary projection

suppressPackageStartupMessages({
  library(ConsensusClusterPlus)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(uwot)
  library(cluster)
  library(ComplexHeatmap)
  library(circlize)
})

ROOT <- "C:/Users/scott/Downloads/Open PBTA"
OUT  <- file.path(ROOT, "output")
FIG  <- file.path(OUT, "figs_clustering")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# -------- 1. Load cohort + features --------
cat("=== [1] Load cohort and feature matrices ===\n")
ann <- fread(file.path(OUT, "cohort_for_deconvolution.tsv"))
main_ids <- fread(file.path(OUT, "cohort_main_final.tsv"))$Kids_First_Biospecimen_ID
ann$is_main <- ann$Kids_First_Biospecimen_ID %in% main_ids
cat(sprintf("  Main cohort (n=%d), All (n=%d)\n", sum(ann$is_main), nrow(ann)))

qts <- fread(file.path(OUT, "immunedeconv_quantiseq.tsv"))
ss  <- fread(file.path(OUT, "ssGSEA_brain_immune_scores.tsv"))

# -------- 2. Feature selection --------
cat("\n=== [2] Feature selection ===\n")
# quanTIseq: drop "uncharacterized cell" (residual fraction)
qts_feat <- qts[qts$cell_type != "uncharacterized cell", ]
qts_mat  <- as.matrix(qts_feat[, -1, with = FALSE])
rownames(qts_mat) <- paste0("QTS_", gsub("[^A-Za-z0-9]+", "_", qts_feat$cell_type))
cat(sprintf("  quanTIseq features: %d (immune cells, excl. uncharacterized)\n", nrow(qts_mat)))

# ssGSEA: keep immune-relevant signatures (drop Stemness, Cell_Cycle for clustering, keep separately)
ss_keep <- c(
  "Microglia_Klemm2020", "MDM_Klemm2020",
  "MgTAM_Antunes2021", "MoTAM_Antunes2021",
  "T_Cell_Cytotoxicity", "T_Cell_Exhaustion",
  "NK_Cell_Activity", "Tregs_Friebel2020",
  "Dendritic_Cell_Activation", "Neutrophil_Activation",
  "MHC_Class_I", "MHC_Class_II",
  "IFN_Gamma_Response", "IFN_Alpha_Response",
  "Chemokine_T_Cell_Recruitment",
  "M1_Macrophage", "M2_Macrophage",
  "TGFb_Immunosuppression",
  "Glioma_Inflammatory_Wang2017"
)
ss_feat <- ss[ss$Signature %in% ss_keep, ]
ss_mat  <- as.matrix(ss_feat[, -1, with = FALSE])
rownames(ss_mat) <- paste0("SS_", ss_feat$Signature)
cat(sprintf("  ssGSEA features: %d (immune-relevant)\n", nrow(ss_mat)))

# Align column order
common_samples <- intersect(colnames(qts_mat), colnames(ss_mat))
qts_mat <- qts_mat[, common_samples]
ss_mat  <- ss_mat[, common_samples]
cat(sprintf("  Shared samples: %d\n", length(common_samples)))

# Concatenate and z-score per feature (row)
feat_raw <- rbind(qts_mat, ss_mat)
feat_z   <- t(scale(t(feat_raw)))
storage.mode(feat_z) <- "double"
cat(sprintf("  Final feature matrix: %d features x %d samples (z-scored)\n",
            nrow(feat_z), ncol(feat_z)))

# Save feature matrix
fwrite(as.data.frame(feat_z) |> tibble::rownames_to_column("feature"),
       file.path(OUT, "clustering_feature_matrix_z.tsv"))

# -------- 3. Run on MAIN cohort only --------
cat("\n=== [3] Consensus clustering — MAIN cohort (n=349) ===\n")
main_ids_avail <- intersect(main_ids, colnames(feat_z))
cat(sprintf("  Available main cohort samples: %d / %d\n", length(main_ids_avail), length(main_ids)))

mat_main <- feat_z[, main_ids_avail]

# ConsensusClusterPlus expects features x samples, will use 1-pearson distance + km
ccp_main <- ConsensusClusterPlus::ConsensusClusterPlus(
  d = mat_main,
  maxK = 6,
  reps = 500,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "km",
  distance = "euclidean",
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2",
  seed = 42,
  plot = "png",
  title = file.path(FIG, "ConsensusCluster_main")
)

# Save consensus matrices for k=2..6
for (k in 2:6) {
  cm <- ccp_main[[k]]$consensusMatrix
  cls <- ccp_main[[k]]$consensusClass
  saveRDS(list(consensus = cm, class = cls),
          file.path(OUT, sprintf("ccp_main_k%d.rds", k)))
}

# PAC score for each k
cat("\n  PAC (Proportion of Ambiguous Clustering) by k:\n")
pac <- sapply(2:6, function(k) {
  cm <- ccp_main[[k]]$consensusMatrix
  x  <- cm[upper.tri(cm)]
  mean(x > 0.1 & x < 0.9)
})
pac_df <- data.frame(k = 2:6, PAC = pac)
print(pac_df)
fwrite(pac_df, file.path(OUT, "consensus_clustering_PAC_main.tsv"))

# Silhouette for each k
cat("\n  Silhouette width by k:\n")
sil <- sapply(2:6, function(k) {
  cls <- ccp_main[[k]]$consensusClass
  d <- as.dist(1 - ccp_main[[k]]$consensusMatrix)
  s <- cluster::silhouette(cls, d)
  mean(s[, "sil_width"])
})
sil_df <- data.frame(k = 2:6, silhouette = sil)
print(sil_df)
fwrite(sil_df, file.path(OUT, "consensus_clustering_silhouette_main.tsv"))

# Choose optimal k (data-driven): lowest PAC ∩ high silhouette
optimal_k <- 3  # default hypothesis; will report selection rationale
best_pac_k <- pac_df$k[which.min(pac_df$PAC)]
best_sil_k <- sil_df$k[which.max(sil_df$silhouette)]
cat(sprintf("\n  Best PAC at k = %d; Best silhouette at k = %d\n", best_pac_k, best_sil_k))
cat(sprintf("  Pre-specified hypothesis: k = 3\n"))

# Use k = 3 as primary as per pre-specified hypothesis; also save k=2 and k=4 as sensitivity
for (k in c(2, 3, 4)) {
  ec <- data.frame(
    Kids_First_Biospecimen_ID = names(ccp_main[[k]]$consensusClass),
    cluster_k = paste0("E", ccp_main[[k]]$consensusClass)
  )
  fwrite(ec, file.path(OUT, sprintf("ecotype_assignment_k%d_main.tsv", k)))
}

# -------- 4. Annotate ecotypes (k=3) by dominant features --------
cat("\n=== [4] Annotate ecotypes (k=3) ===\n")
clus <- ccp_main[[3]]$consensusClass
clus_df <- data.frame(
  Kids_First_Biospecimen_ID = names(clus),
  cluster = clus
)

# Compute per-cluster mean of each feature
feat_long <- as.data.frame(t(mat_main))
feat_long$cluster <- clus[rownames(feat_long)]
clust_means <- aggregate(. ~ cluster, data = feat_long, FUN = mean)
fwrite(clust_means, file.path(OUT, "ecotype_k3_feature_means.tsv"))

cat("\n  Per-cluster mean (z-score) of key features:\n")
print(round(clust_means, 2))

# Auto-annotation logic:
# Lymphocyte-inflamed: high T-cell cytotoxicity, IFN, MHC, CD8
# Myeloid-dominant:    high Microglia/MDM/M2 + neutrophil
# Immune-desert:       low across the board
get_top_cluster <- function(feature_name) {
  if (!feature_name %in% colnames(clust_means)) return(NA)
  v <- clust_means[[feature_name]]
  clust_means$cluster[which.max(v)]
}
get_bottom_cluster <- function(feature_name) {
  if (!feature_name %in% colnames(clust_means)) return(NA)
  v <- clust_means[[feature_name]]
  clust_means$cluster[which.min(v)]
}

lym_cluster <- get_top_cluster("SS_T_Cell_Cytotoxicity")
mye_cluster <- get_top_cluster("SS_MoTAM_Antunes2021")
# Desert = the remaining cluster
desert_cluster <- setdiff(1:3, c(lym_cluster, mye_cluster))[1]

ecotype_map <- c()
ecotype_map[as.character(lym_cluster)]    <- "Lymphocyte-inflamed"
ecotype_map[as.character(mye_cluster)]    <- "Myeloid-dominant"
ecotype_map[as.character(desert_cluster)] <- "Immune-desert"

cat(sprintf("\n  Cluster %d -> Lymphocyte-inflamed\n", lym_cluster))
cat(sprintf("  Cluster %d -> Myeloid-dominant\n", mye_cluster))
cat(sprintf("  Cluster %d -> Immune-desert\n", desert_cluster))

clus_df$ecotype <- ecotype_map[as.character(clus_df$cluster)]
fwrite(clus_df, file.path(OUT, "ecotype_assignment_k3_annotated.tsv"))

# -------- 5. PCA + UMAP visualization --------
cat("\n=== [5] PCA + UMAP ===\n")
# PCA
pc <- prcomp(t(mat_main), center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  Kids_First_Biospecimen_ID = rownames(pc$x),
  PC1 = pc$x[, 1], PC2 = pc$x[, 2],
  ecotype = clus_df$ecotype[match(rownames(pc$x), clus_df$Kids_First_Biospecimen_ID)]
)
fwrite(pca_df, file.path(OUT, "ecotype_pca_main.tsv"))

# UMAP
um <- uwot::umap(t(mat_main), n_neighbors = 15, min_dist = 0.3, n_components = 2, seed = 42)
umap_df <- data.frame(
  Kids_First_Biospecimen_ID = colnames(mat_main),
  UMAP1 = um[, 1], UMAP2 = um[, 2],
  ecotype = clus_df$ecotype[match(colnames(mat_main), clus_df$Kids_First_Biospecimen_ID)]
)
fwrite(umap_df, file.path(OUT, "ecotype_umap_main.tsv"))

# -------- 6. Ecotype x cohort_group / location / age contingency --------
cat("\n=== [6] Association tests ===\n")
meta <- ann[match(clus_df$Kids_First_Biospecimen_ID, ann$Kids_First_Biospecimen_ID), ]
clus_df$cohort_group   <- meta$cohort_group
clus_df$location_class <- meta$location_class
clus_df$age_dev_group  <- meta$age_dev_group
clus_df$age_years      <- as.numeric(meta$age_years)
clus_df$reported_gender<- meta$reported_gender
clus_df$OS_days        <- as.numeric(meta$OS_days)
clus_df$OS_status      <- meta$OS_status

# ecotype x cohort_group
cat("\n  Ecotype x cohort_group:\n")
tab1 <- table(clus_df$ecotype, clus_df$cohort_group)
print(tab1)
fisher_p1 <- tryCatch(
  fisher.test(tab1, simulate.p.value = TRUE, B = 10000)$p.value,
  error = function(e) NA
)
cat(sprintf("  Fisher's exact (simulated) p = %.2e\n", fisher_p1))

# ecotype x location_class (exclude Ambiguous / Other from association as planned)
loc_ok <- clus_df[clus_df$location_class %in% c("Midline", "Hemispheric", "Posterior_fossa"), ]
cat("\n  Ecotype x location_class (excluding Ambiguous/Other):\n")
tab2 <- table(loc_ok$ecotype, loc_ok$location_class)
print(tab2)
fisher_p2 <- tryCatch(
  fisher.test(tab2, simulate.p.value = TRUE, B = 10000)$p.value,
  error = function(e) NA
)
cat(sprintf("  Fisher's exact (simulated) p = %.2e\n", fisher_p2))

# ecotype x age_dev_group
cat("\n  Ecotype x age_dev_group:\n")
tab3 <- table(clus_df$ecotype, clus_df$age_dev_group)
print(tab3)
fisher_p3 <- tryCatch(
  fisher.test(tab3, simulate.p.value = TRUE, B = 10000)$p.value,
  error = function(e) NA
)
cat(sprintf("  Fisher's exact (simulated) p = %.2e\n", fisher_p3))

# Save tables
write_xtab <- function(t, fname) {
  df <- as.data.frame.matrix(t) |> tibble::rownames_to_column("ecotype")
  fwrite(df, file.path(OUT, fname))
}
write_xtab(tab1, "ecotype_x_cohort_group.tsv")
write_xtab(tab2, "ecotype_x_location.tsv")
write_xtab(tab3, "ecotype_x_age_dev.tsv")

stats_df <- data.frame(
  comparison = c("ecotype_x_cohort_group", "ecotype_x_location_class", "ecotype_x_age_dev_group"),
  fisher_p_simulated = c(fisher_p1, fisher_p2, fisher_p3),
  N = c(nrow(clus_df), nrow(loc_ok), nrow(clus_df))
)
fwrite(stats_df, file.path(OUT, "ecotype_association_pvalues.tsv"))

# -------- 7. Heatmap (ComplexHeatmap) --------
cat("\n=== [7] Heatmap ===\n")
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
ec_colors <- c("Lymphocyte-inflamed" = "#3B82F6",
               "Myeloid-dominant"    = "#EF4444",
               "Immune-desert"       = "#9CA3AF")
cg_colors <- c("DMG_K27" = "#1F77B4", "DHG_G34" = "#FF7F0E",
               "pHGG_WT" = "#2CA02C", "IHG" = "#9467BD")
loc_colors <- c("Midline" = "#D62728", "Hemispheric" = "#17BECF",
                "Posterior_fossa" = "#8C564B", "Ambiguous" = "#BCBD22",
                "Other/NA" = "#7F7F7F")

# Column annotation
col_anno <- HeatmapAnnotation(
  Ecotype  = clus_df$ecotype[match(colnames(mat_main), clus_df$Kids_First_Biospecimen_ID)],
  Cohort   = clus_df$cohort_group[match(colnames(mat_main), clus_df$Kids_First_Biospecimen_ID)],
  Location = clus_df$location_class[match(colnames(mat_main), clus_df$Kids_First_Biospecimen_ID)],
  col = list(Ecotype = ec_colors, Cohort = cg_colors, Location = loc_colors),
  show_legend = TRUE, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 9)
)

png(file.path(FIG, "ecotype_heatmap_main.png"), width = 14, height = 10, units = "in", res = 200)
ht <- Heatmap(
  mat_main,
  name = "z-score",
  col = col_fun,
  top_annotation = col_anno,
  show_column_names = FALSE,
  cluster_columns = TRUE,
  column_split = factor(clus_df$ecotype[match(colnames(mat_main), clus_df$Kids_First_Biospecimen_ID)],
                        levels = c("Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert")),
  cluster_rows = TRUE,
  row_names_gp = gpar(fontsize = 8),
  row_split = ifelse(grepl("^QTS", rownames(mat_main)), "quanTIseq", "ssGSEA"),
  column_title_gp = gpar(fontsize = 11, fontface = "bold")
)
draw(ht, merge_legend = TRUE)
dev.off()
cat("  Wrote ecotype_heatmap_main.png\n")

# PCA plot
png(file.path(FIG, "ecotype_pca_main.png"), width = 7, height = 6, units = "in", res = 200)
plot(pca_df$PC1, pca_df$PC2, col = ec_colors[pca_df$ecotype], pch = 19, cex = 0.9,
     xlab = sprintf("PC1 (%.1f%%)", 100 * pc$sdev[1]^2 / sum(pc$sdev^2)),
     ylab = sprintf("PC2 (%.1f%%)", 100 * pc$sdev[2]^2 / sum(pc$sdev^2)),
     main = "PCA — Main cohort by ecotype (n=349)")
legend("topright", names(ec_colors), col = ec_colors, pch = 19, bty = "n", cex = 0.9)
dev.off()
cat("  Wrote ecotype_pca_main.png\n")

# UMAP plot
png(file.path(FIG, "ecotype_umap_main.png"), width = 7, height = 6, units = "in", res = 200)
plot(umap_df$UMAP1, umap_df$UMAP2, col = ec_colors[umap_df$ecotype], pch = 19, cex = 0.9,
     xlab = "UMAP1", ylab = "UMAP2",
     main = "UMAP — Main cohort by ecotype (n=349)")
legend("topright", names(ec_colors), col = ec_colors, pch = 19, bty = "n", cex = 0.9)
dev.off()
cat("  Wrote ecotype_umap_main.png\n")

# Stacked barplot: ecotype proportion by cohort_group
png(file.path(FIG, "ecotype_distribution_by_cohort.png"), width = 8, height = 5, units = "in", res = 200)
prop_tab <- prop.table(table(clus_df$cohort_group, clus_df$ecotype), margin = 1)
par(mar = c(5, 4, 4, 8), xpd = TRUE)
barplot(t(prop_tab), col = ec_colors[colnames(prop_tab)],
        las = 1, ylab = "Proportion",
        main = "Ecotype distribution across pediatric glioma subgroups")
legend("topright", inset = c(-0.25, 0), legend = rownames(t(prop_tab)),
       fill = ec_colors[colnames(prop_tab)], bty = "n", cex = 0.85)
dev.off()
cat("  Wrote ecotype_distribution_by_cohort.png\n")

# Ecotype distribution by location_class (excluding Ambiguous)
png(file.path(FIG, "ecotype_distribution_by_location.png"), width = 8, height = 5, units = "in", res = 200)
prop_tab2 <- prop.table(table(loc_ok$location_class, loc_ok$ecotype), margin = 1)
par(mar = c(5, 4, 4, 8), xpd = TRUE)
barplot(t(prop_tab2), col = ec_colors[colnames(prop_tab2)],
        las = 1, ylab = "Proportion",
        main = "Ecotype distribution by anatomical location")
legend("topright", inset = c(-0.25, 0), legend = rownames(t(prop_tab2)),
       fill = ec_colors[colnames(prop_tab2)], bty = "n", cex = 0.85)
dev.off()
cat("  Wrote ecotype_distribution_by_location.png\n")

cat("\nAll outputs written to:\n", OUT, "\n", sep = "")
cat("Figures to:\n", FIG, "\n", sep = "")
cat("\nDone.\n")
