# ssGSEA pathway comparison across ecotypes
# Kruskal-Wallis + Dunn post-hoc + boxplots

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(FSA)
  library(rstatix)
})

ROOT <- "C:/Users/scott/Downloads/Open PBTA"
OUT  <- file.path(ROOT, "output")
FIG  <- file.path(OUT, "figs_pathway")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# -------- Load --------
clus <- fread(file.path(OUT, "ecotype_assignment_k3_annotated.tsv"))
ss   <- fread(file.path(OUT, "ssGSEA_brain_immune_scores.tsv"))

# Wide → long
ss_long <- pivot_longer(ss, -Signature,
                        names_to = "Kids_First_Biospecimen_ID",
                        values_to = "score")

merged <- inner_join(ss_long, clus[, c("Kids_First_Biospecimen_ID","ecotype","cohort_group")],
                     by = "Kids_First_Biospecimen_ID")
merged$ecotype <- factor(merged$ecotype,
                         levels = c("Immune-desert","Intermediate","Inflamed"))
cat(sprintf("Merged: %d (signatures x samples)\n", nrow(merged)))

ec_colors <- c("Inflamed" = "#3B82F6",
               "Intermediate" = "#F59E0B",
               "Immune-desert" = "#9CA3AF")

# -------- 1. Kruskal-Wallis per signature --------
cat("\n=== [1] Kruskal-Wallis test per signature ===\n")
sigs <- unique(merged$Signature)
kw_res <- data.frame()
for (s in sigs) {
  d <- merged[merged$Signature == s, ]
  kw <- kruskal.test(score ~ ecotype, data = d)
  kw_res <- rbind(kw_res,
                  data.frame(Signature = s,
                             KW_chi2 = kw$statistic,
                             KW_df = kw$parameter,
                             KW_p = kw$p.value))
}
kw_res$KW_p_BH <- p.adjust(kw_res$KW_p, method = "BH")
kw_res <- kw_res[order(kw_res$KW_p), ]
print(kw_res)
fwrite(kw_res, file.path(OUT, "ssGSEA_KruskalWallis_by_ecotype.tsv"))

# -------- 2. Dunn post-hoc per signature --------
cat("\n=== [2] Dunn post-hoc per signature ===\n")
dunn_all <- data.frame()
for (s in sigs) {
  d <- merged[merged$Signature == s, ]
  dt <- tryCatch(
    dunnTest(score ~ ecotype, data = d, method = "bh")$res,
    error = function(e) NULL
  )
  if (!is.null(dt)) {
    dt$Signature <- s
    dunn_all <- rbind(dunn_all, dt)
  }
}
dunn_all <- dunn_all[, c("Signature", "Comparison", "Z", "P.unadj", "P.adj")]
print(head(dunn_all, 30))
fwrite(dunn_all, file.path(OUT, "ssGSEA_Dunn_posthoc_by_ecotype.tsv"))

# -------- 3. Per-signature boxplots --------
cat("\n=== [3] Boxplots per signature ===\n")
plot_box <- function(sig, file) {
  d <- merged[merged$Signature == sig, ]
  kw_p <- kw_res$KW_p[kw_res$Signature == sig]
  kw_p_bh <- kw_res$KW_p_BH[kw_res$Signature == sig]
  p <- ggplot(d, aes(x = ecotype, y = score, fill = ecotype)) +
    geom_boxplot(outlier.size = 0.6, width = 0.7) +
    scale_fill_manual(values = ec_colors) +
    labs(title = sig,
         subtitle = sprintf("Kruskal-Wallis p = %.2e (BH-adj. p = %.2e)", kw_p, kw_p_bh),
         x = "Ecotype", y = "ssGSEA score (z)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
  ggsave(file, p, width = 5, height = 4, dpi = 200)
}

# Save key brain/immune signatures
key_sigs <- c(
  "Microglia_Klemm2020", "MDM_Klemm2020",
  "MgTAM_Antunes2021", "MoTAM_Antunes2021",
  "T_Cell_Cytotoxicity", "NK_Cell_Activity",
  "MHC_Class_I", "MHC_Class_II",
  "IFN_Gamma_Response", "IFN_Alpha_Response",
  "MAPK_Activity", "TGFb_Immunosuppression",
  "M1_Macrophage", "M2_Macrophage",
  "Chemokine_T_Cell_Recruitment", "T_Cell_Exhaustion",
  "Dendritic_Cell_Activation", "Neutrophil_Activation",
  "Tregs_Friebel2020", "Glioma_Inflammatory_Wang2017",
  "DAM_KerenShaul2017", "Microglia_Core_Homeostatic",
  "Cell_Cycle_Proliferation", "Stemness_Brain_Tumor"
)
for (s in key_sigs) {
  plot_box(s, file.path(FIG, sprintf("box_%s.png", s)))
}
cat(sprintf("  Wrote %d signature boxplots\n", length(key_sigs)))

# -------- 4. Multi-panel composite figure --------
cat("\n=== [4] Composite figure (key 12 signatures) ===\n")
hilite_sigs <- c(
  "Microglia_Klemm2020", "MDM_Klemm2020",
  "MgTAM_Antunes2021", "MoTAM_Antunes2021",
  "T_Cell_Cytotoxicity", "NK_Cell_Activity",
  "MHC_Class_II", "IFN_Gamma_Response",
  "MAPK_Activity", "TGFb_Immunosuppression",
  "Chemokine_T_Cell_Recruitment", "DAM_KerenShaul2017"
)
sub_long <- merged[merged$Signature %in% hilite_sigs, ]
sub_long$Signature <- factor(sub_long$Signature, levels = hilite_sigs)
sub_long$kw_p <- kw_res$KW_p_BH[match(sub_long$Signature, kw_res$Signature)]
sub_long$facet_label <- sprintf("%s\n(KW p_BH = %.1e)",
                                 sub_long$Signature, sub_long$kw_p)
sub_long$facet_label <- factor(sub_long$facet_label, levels = unique(sub_long$facet_label[order(sub_long$Signature)]))

p_comp <- ggplot(sub_long, aes(x = ecotype, y = score, fill = ecotype)) +
  geom_boxplot(outlier.size = 0.4, width = 0.7) +
  facet_wrap(~ facet_label, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = ec_colors) +
  labs(x = NULL, y = "ssGSEA score (z)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(size = 8, face = "bold"))
ggsave(file.path(FIG, "Pathway_composite_12sig.png"), p_comp,
       width = 13, height = 9, dpi = 200)
cat("  Wrote Pathway_composite_12sig.png\n")

# -------- 5. Per-cohort_group stratified boxplots for key signatures --------
cat("\n=== [5] Cohort-stratified boxplots ===\n")
cohort_sigs <- c("Microglia_Klemm2020", "MDM_Klemm2020",
                 "T_Cell_Cytotoxicity", "MAPK_Activity",
                 "MHC_Class_II", "IFN_Gamma_Response")
for (s in cohort_sigs) {
  d <- merged[merged$Signature == s, ]
  d$cohort_group <- factor(d$cohort_group,
                            levels = c("DMG_K27","DHG_G34","pHGG_WT","IHG","BRAF_ALT"))
  d <- d[!is.na(d$cohort_group), ]
  p <- ggplot(d, aes(x = cohort_group, y = score, fill = ecotype)) +
    geom_boxplot(outlier.size = 0.4, position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = ec_colors) +
    labs(title = s, x = NULL, y = "ssGSEA score (z)",
         fill = "Ecotype") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG, sprintf("cohort_strat_%s.png", s)), p,
         width = 8, height = 4.5, dpi = 200)
}
cat(sprintf("  Wrote %d cohort-stratified plots\n", length(cohort_sigs)))

# -------- 6. Per-ecotype effect size summary (rank-biserial) --------
cat("\n=== [6] Effect size table ===\n")
es_all <- data.frame()
for (s in sigs) {
  d <- merged[merged$Signature == s, c("score","ecotype")]
  comps <- list(c("Inflamed","Immune-desert"),
                c("Inflamed","Intermediate"),
                c("Intermediate","Immune-desert"))
  for (cp in comps) {
    sub <- d[d$ecotype %in% cp, ]
    if (nrow(sub) < 4) next
    sub$ecotype <- droplevels(sub$ecotype)
    w <- tryCatch(wilcox.test(score ~ ecotype, data = sub, conf.int = TRUE),
                  error = function(e) NULL)
    if (is.null(w)) next
    # rank-biserial r effect size
    n1 <- sum(sub$ecotype == cp[1]); n2 <- sum(sub$ecotype == cp[2])
    U  <- w$statistic
    r  <- (2 * U) / (n1 * n2) - 1
    es_all <- rbind(es_all, data.frame(
      Signature = s,
      Comparison = sprintf("%s vs %s", cp[1], cp[2]),
      W = unname(U), wilcox_p = w$p.value,
      rank_biserial_r = round(r, 3)
    ))
  }
}
es_all$wilcox_p_BH <- p.adjust(es_all$wilcox_p, method = "BH")
fwrite(es_all, file.path(OUT, "ssGSEA_pairwise_effect_size.tsv"))
print(head(es_all[order(abs(es_all$rank_biserial_r), decreasing = TRUE), ], 15))

cat("\nDone.\n")
