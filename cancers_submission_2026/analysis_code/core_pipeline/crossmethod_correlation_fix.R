# Post-hoc cross-method myeloid correlation + summary stats
# Reads saved immunedeconv + ssGSEA outputs and computes Spearman correlation

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

OUT <- "C:/Users/scott/Downloads/Open PBTA/output"

# -------- Load all results --------
mcp  <- fread(file.path(OUT, "immunedeconv_mcp_counter.tsv"))
qts  <- fread(file.path(OUT, "immunedeconv_quantiseq.tsv"))
epic <- fread(file.path(OUT, "immunedeconv_epic.tsv"))
xcl  <- fread(file.path(OUT, "immunedeconv_xcell.tsv"))
ss   <- fread(file.path(OUT, "ssGSEA_brain_immune_scores.tsv"))

cat("Cell types per method:\n")
cat("  MCP-counter:\n"); print(mcp$cell_type)
cat("  quanTIseq:\n");   print(qts$cell_type)
cat("  EPIC:\n");        print(epic$cell_type)
cat("  xCell (first 30):\n"); print(head(xcl$cell_type, 30))
cat("  ssGSEA signatures:\n"); print(ss$Signature)

# -------- Helper to grab per-sample score from a method --------
get_score <- function(df, pattern, label) {
  match_rows <- grep(pattern, df$cell_type, value = TRUE, ignore.case = TRUE)
  if (length(match_rows) == 0) {
    cat(sprintf("    [%s] no match for pattern '%s'\n", label, pattern))
    return(NULL)
  }
  cat(sprintf("    [%s] matched: %s\n", label, paste(match_rows, collapse = ", ")))
  sub <- df[df$cell_type %in% match_rows, , drop = FALSE]
  vals <- colSums(as.matrix(sub[, -1, with = FALSE]))
  data.frame(sample = names(vals), score = unname(vals))
}

get_ssgsea <- function(df, sig, label) {
  if (!sig %in% df$Signature) {
    cat(sprintf("    [%s] signature '%s' not found\n", label, sig))
    return(NULL)
  }
  vals <- as.numeric(df[df$Signature == sig, -1, with = FALSE])
  data.frame(sample = colnames(df)[-1], score = vals)
}

# -------- Myeloid / Macrophage cross-method --------
cat("\n=== Myeloid signals across methods ===\n")
mcp_mac  <- get_score(mcp,  "Macrophage/Monocyte|Macrophage|Monocyte", "MCP")
qts_mac  <- get_score(qts,  "Macrophage", "QTS")
epic_mac <- get_score(epic, "Macrophage", "EPIC")
xcl_mac  <- get_score(xcl,  "^Macrophage|^Monocyte", "xCell")
mg_ss    <- get_ssgsea(ss, "Microglia_Klemm2020", "Microglia(Klemm)")
mgcore   <- get_ssgsea(ss, "Microglia_Core_Homeostatic", "MicrogliaCore(Bohlen)")
mdm_ss   <- get_ssgsea(ss, "MDM_Klemm2020", "MDM(Klemm)")
mgtam    <- get_ssgsea(ss, "MgTAM_Antunes2021", "MgTAM(Antunes)")
motam    <- get_ssgsea(ss, "MoTAM_Antunes2021", "MoTAM(Antunes)")

merge_named <- function(lst) {
  lst <- lst[!sapply(lst, is.null)]
  Reduce(function(a, b) merge(a, b, by = "sample"),
         lapply(seq_along(lst), function(i) {
           x <- lst[[i]]; colnames(x) <- c("sample", names(lst)[i]); x
         }))
}

mye <- merge_named(list(
  MCP_macro_mono       = mcp_mac,
  QTS_macro            = qts_mac,
  EPIC_macro           = epic_mac,
  xCell_macro_mono     = xcl_mac,
  Microglia_Klemm      = mg_ss,
  MicrogliaCore_Bohlen = mgcore,
  MDM_Klemm            = mdm_ss,
  MgTAM_Antunes        = mgtam,
  MoTAM_Antunes        = motam
))

cor_mye <- cor(mye[, -1], method = "spearman", use = "pairwise.complete.obs")
cat("\n  Spearman correlation matrix (myeloid):\n")
print(round(cor_mye, 3))

# -------- T-cell / CD8 cross-method --------
cat("\n=== T-cell / CD8 signals across methods ===\n")
mcp_cd8  <- get_score(mcp,  "T cell CD8|CD8", "MCP")
qts_cd8  <- get_score(qts,  "T cell CD8", "QTS")
epic_cd8 <- get_score(epic, "CD8", "EPIC")
xcl_cd8  <- get_score(xcl,  "^CD8", "xCell")
tcyt_ss  <- get_ssgsea(ss, "T_Cell_Cytotoxicity", "Cytolytic")

tcd8 <- merge_named(list(
  MCP_CD8     = mcp_cd8,
  QTS_CD8     = qts_cd8,
  EPIC_CD8    = epic_cd8,
  xCell_CD8   = xcl_cd8,
  Cytolytic_ssGSEA = tcyt_ss
))

cor_cd8 <- cor(tcd8[, -1], method = "spearman", use = "pairwise.complete.obs")
cat("\n  Spearman correlation matrix (CD8/cytolytic):\n")
print(round(cor_cd8, 3))

# -------- NK cross-method --------
cat("\n=== NK cell signals ===\n")
mcp_nk  <- get_score(mcp,  "NK", "MCP")
qts_nk  <- get_score(qts,  "NK cell", "QTS")
epic_nk <- get_score(epic, "NK cell", "EPIC")
xcl_nk  <- get_score(xcl,  "^NK", "xCell")
nk_ss   <- get_ssgsea(ss, "NK_Cell_Activity", "NK_activity")

nk_all <- merge_named(list(
  MCP_NK      = mcp_nk,
  QTS_NK      = qts_nk,
  EPIC_NK     = epic_nk,
  xCell_NK    = xcl_nk,
  NK_ssGSEA   = nk_ss
))
cor_nk <- cor(nk_all[, -1], method = "spearman", use = "pairwise.complete.obs")
cat("\n  Spearman correlation matrix (NK):\n")
print(round(cor_nk, 3))

# -------- Save outputs --------
fwrite(mye,  file.path(OUT, "myeloid_signals_merged.tsv"))
fwrite(tcd8, file.path(OUT, "cd8_signals_merged.tsv"))
fwrite(nk_all, file.path(OUT, "nk_signals_merged.tsv"))

fwrite(as.data.frame(cor_mye) |> tibble::rownames_to_column("method"),
       file.path(OUT, "spearman_myeloid_crossmethod.tsv"))
fwrite(as.data.frame(cor_cd8) |> tibble::rownames_to_column("method"),
       file.path(OUT, "spearman_cd8_crossmethod.tsv"))
fwrite(as.data.frame(cor_nk) |> tibble::rownames_to_column("method"),
       file.path(OUT, "spearman_nk_crossmethod.tsv"))

# -------- Long-format ssGSEA --------
ss_long <- pivot_longer(ss, -Signature, names_to = "Kids_First_Biospecimen_ID", values_to = "score")
fwrite(ss_long, file.path(OUT, "ssGSEA_brain_immune_scores_long.tsv"))

cat("\nDone.\n")
cat("Files written:\n")
cat("  - myeloid_signals_merged.tsv  (per-sample Macrophage/Microglia/MDM signals from 4+ssGSEA)\n")
cat("  - cd8_signals_merged.tsv\n")
cat("  - nk_signals_merged.tsv\n")
cat("  - spearman_*_crossmethod.tsv  (cross-method correlation matrices)\n")
cat("  - ssGSEA_brain_immune_scores_long.tsv  (long format for plotting)\n")
