#!/usr/bin/env Rscript
# ==============================================================================
# DESeq2 Differential Expression: Inflamed vs Immune-desert
# Pediatric Glioma Immune Ecotype Project
# ==============================================================================
#
# USAGE:
#   1. Set working directory to your project folder
#   2. Ensure gene-counts-rsem-expected_count-collapsed.rds is in data/
#   3. Run: Rscript DESeq2_Inflamed_vs_Desert.R
#   OR source in RStudio
#
# REQUIREMENTS: DESeq2, ggplot2, ggrepel, EnhancedVolcano (optional)
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})

# ---- 0. Paths (edit if needed) -----------------------------------------------
BASE_DIR   <- getwd()  # should be "Open PBTA" folder
setwd("C:/Users/scott/Downloads/Open PBTA")
getwd()
COUNT_RDS  <- file.path(BASE_DIR, "data",
                        "gene-counts-rsem-expected_count-collapsed.rds")
ECOTYPE_F  <- file.path(BASE_DIR, "output",
                        "ecotype_assignment_k3_annotated.tsv")
OUT_DIR    <- file.path(BASE_DIR, "output")

cat("=== DESeq2 Inflamed vs Immune-desert ===\n")
cat("Count RDS:", COUNT_RDS, "\n")
cat("Ecotype file:", ECOTYPE_F, "\n\n")

# ---- 1. Load count matrix ----------------------------------------------------
cat("[1/6] Loading count matrix...\n")
counts_full <- readRDS(COUNT_RDS)
cat("  Full matrix:", nrow(counts_full), "genes x", ncol(counts_full), "samples\n")

# ---- 2. Load ecotype assignments ---------------------------------------------
cat("[2/6] Loading ecotype assignments...\n")
eco <- read.csv(ECOTYPE_F, stringsAsFactors = FALSE)
cat("  Ecotype distribution:\n")
print(table(eco$ecotype))

# Subset to Inflamed vs Immune-desert only
eco_sub <- eco[eco$ecotype %in% c("Inflamed", "Immune-desert"), ]
cat("\n  Inflamed vs Desert:", nrow(eco_sub), "samples\n")

# ---- 3. Subset & prepare count matrix ----------------------------------------
cat("[3/6] Preparing count matrix...\n")
shared <- intersect(eco_sub$Kids_First_Biospecimen_ID, colnames(counts_full))
cat("  Matched samples:", length(shared), "\n")

counts_sub <- counts_full[, shared]

# RSEM expected_count can be fractional — round to integers for DESeq2
counts_sub <- round(as.matrix(counts_sub))

# Build colData
coldata <- eco_sub[match(shared, eco_sub$Kids_First_Biospecimen_ID), ]
rownames(coldata) <- coldata$Kids_First_Biospecimen_ID

# Set reference level: Immune-desert as baseline → positive LFC = up in Inflamed
coldata$ecotype <- factor(coldata$ecotype, levels = c("Immune-desert", "Inflamed"))
coldata$cohort_group <- factor(coldata$cohort_group)

cat("  Final design:\n")
print(table(coldata$ecotype, coldata$cohort_group))

# ---- 4. DESeq2 ---------------------------------------------------------------
cat("\n[4/6] Running DESeq2 (design: ~ cohort_group + ecotype)...\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData   = coldata,
  design    = ~ cohort_group + ecotype
)

# Pre-filter: keep genes with ≥10 total counts across samples
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("  Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)

res <- results(dds, contrast = c("ecotype", "Inflamed", "Immune-desert"),
               alpha = 0.05)
cat("  DESeq2 complete.\n")
cat("  Summary:\n")
summary(res)

# ---- 5. Results table ---------------------------------------------------------
cat("[5/6] Extracting results...\n")

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

# DEG thresholds: |log2FC| >= 1 and padj < 0.05
res_df$significant <- ifelse(
  !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1,
  ifelse(res_df$log2FoldChange > 0, "Up_Inflamed", "Up_Desert"),
  "NS"
)

cat("\n  DEG counts (|log2FC| >= 1, padj < 0.05):\n")
print(table(res_df$significant))

# Save full results
out_csv <- file.path(OUT_DIR, "DESeq2_Inflamed_vs_Desert_results.csv")
write.csv(res_df, out_csv, row.names = FALSE)
cat("  Full results saved:", out_csv, "\n")

# Save significant DEGs only
deg_sig <- res_df[res_df$significant != "NS", ]
deg_sig <- deg_sig[order(deg_sig$log2FoldChange, decreasing = TRUE), ]
out_deg <- file.path(OUT_DIR, "DESeq2_significant_DEGs.csv")
write.csv(deg_sig, out_deg, row.names = FALSE)
cat("  Significant DEGs saved:", out_deg, "\n")
cat("  Up in Inflamed:", sum(deg_sig$significant == "Up_Inflamed"), "\n")
cat("  Up in Desert:", sum(deg_sig$significant == "Up_Desert"), "\n")

# Top 20 by fold-change
cat("\n  Top 20 DEGs (by |log2FC|):\n")
top20 <- head(deg_sig[order(-abs(deg_sig$log2FoldChange)), ], 20)
print(top20[, c("gene", "log2FoldChange", "padj", "significant")])

# ---- 6. Volcano plot ---------------------------------------------------------
cat("\n[6/6] Generating volcano plot...\n")

# Color palette matching manuscript ecotype colors
res_df$color <- ifelse(res_df$significant == "Up_Inflamed", "#E64B35",
                ifelse(res_df$significant == "Up_Desert", "#91D1C2", "grey80"))

# Label top genes
top_label <- rbind(
  head(res_df[res_df$significant == "Up_Inflamed", ], 15),
  head(res_df[res_df$significant == "Up_Desert", ], 15)
)

# Immune-relevant genes to always label if significant
immune_genes <- c("HAVCR2", "LAG3", "PDCD1", "CD274", "CTLA4", "TIGIT",
                  "IDO1", "SIGLEC15", "CD276", "B4GALNT1", "IL13RA2",
                  "ERBB2", "CXCL9", "CXCL10", "CCL5", "IFNG", "GZMB",
                  "PRF1", "CD8A", "CD4", "FOXP3", "CD163", "CSF1R",
                  "TGFB1", "IL10", "VEGFA", "HLA-A", "HLA-DRA")
immune_label <- res_df[res_df$gene %in% immune_genes & res_df$significant != "NS", ]
top_label <- rbind(top_label, immune_label)
top_label <- top_label[!duplicated(top_label$gene), ]

p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), size = 0.8, alpha = 0.6) +
  scale_color_manual(
    values = c("Up_Inflamed" = "#E64B35", "Up_Desert" = "#91D1C2", "NS" = "grey80"),
    labels = c("Up_Inflamed" = "Up in Inflamed",
               "Up_Desert"   = "Up in Immune-desert",
               "NS"          = "Not significant"),
    name = ""
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_text_repel(
    data = top_label,
    aes(label = gene),
    size = 2.5, max.overlaps = 30, segment.size = 0.2,
    fontface = "italic", color = "black"
  ) +
  labs(
    title = "Differential Expression: Inflamed vs Immune-desert",
    subtitle = paste0("DESeq2 | design: ~ cohort_group + ecotype | ",
                      sum(res_df$significant == "Up_Inflamed"), " up / ",
                      sum(res_df$significant == "Up_Desert"), " down"),
    x = expression(log[2]~"Fold Change (Inflamed / Immune-desert)"),
    y = expression(-log[10]~"(p-value)")
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

# Save
out_png <- file.path(OUT_DIR, "figs_300dpi",
                     "DESeq2_volcano_Inflamed_vs_Desert.png")
out_pdf <- file.path(OUT_DIR, "figs_300dpi",
                     "DESeq2_volcano_Inflamed_vs_Desert.pdf")

ggsave(out_png, p, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 8, height = 6, bg = "white")
cat("  Volcano plot saved:\n    ", out_png, "\n    ", out_pdf, "\n")

# ---- Summary -----------------------------------------------------------------
cat("\n=== DONE ===\n")
cat("Outputs in:", OUT_DIR, "\n")
cat("  1. DESeq2_Inflamed_vs_Desert_results.csv  (all genes)\n")
cat("  2. DESeq2_significant_DEGs.csv             (|log2FC|≥1, padj<0.05)\n")
cat("  3. figs_300dpi/DESeq2_volcano_*.png/pdf    (volcano plot, 300dpi)\n")

sessionInfo()
