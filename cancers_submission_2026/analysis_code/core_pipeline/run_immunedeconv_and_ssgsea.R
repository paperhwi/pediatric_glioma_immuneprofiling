# Run immunedeconv (xCell, quanTIseq, MCP-counter, EPIC) + brain-immune ssGSEA
# Cross-validation deconvolution + brain-tuned signature scoring for pediatric glioma TME

suppressPackageStartupMessages({
  library(immunedeconv)
  library(GSVA)
  library(GSEABase)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(readr)
})

ROOT <- "C:/Users/scott/Downloads/Open PBTA"
OUT  <- file.path(ROOT, "output")
dir.create(OUT, showWarnings = FALSE)

cat("\n=== [1] Load TPM matrix ===\n")
tpm_path <- file.path(OUT, "tpm_for_cibersortx.tsv")
tpm <- fread(tpm_path, sep = "\t", data.table = FALSE)
rownames(tpm) <- tpm$GeneSymbol
tpm$GeneSymbol <- NULL
tpm_mat <- as.matrix(tpm)
cat(sprintf("  TPM matrix: %d genes x %d samples\n", nrow(tpm_mat), ncol(tpm_mat)))

# Sanity: tibble/data.frame to matrix
storage.mode(tpm_mat) <- "double"

# -------- 2. immunedeconv cross-validation --------
cat("\n=== [2] immunedeconv methods ===\n")
deconv_results <- list()

cat("  [2a] MCP-counter (fastest, ~1-2 min) ...\n")
t0 <- Sys.time()
deconv_results$mcp_counter <- deconvolute(tpm_mat, "mcp_counter")
cat(sprintf("       done in %.1f sec\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cat("  [2b] quanTIseq (~10-20 min) ...\n")
t0 <- Sys.time()
deconv_results$quantiseq <- deconvolute(tpm_mat, "quantiseq", tumor = TRUE)
cat(sprintf("       done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cat("  [2c] EPIC (~5-10 min) ...\n")
t0 <- Sys.time()
deconv_results$epic <- deconvolute(tpm_mat, "epic", tumor = TRUE)
cat(sprintf("       done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cat("  [2d] xCell (~20-40 min, slowest) ...\n")
t0 <- Sys.time()
deconv_results$xcell <- deconvolute(tpm_mat, "xcell")
cat(sprintf("       done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# -------- 3. Save each immunedeconv result --------
cat("\n=== [3] Save immunedeconv outputs ===\n")
for (nm in names(deconv_results)) {
  out_file <- file.path(OUT, sprintf("immunedeconv_%s.tsv", nm))
  fwrite(deconv_results[[nm]], out_file, sep = "\t")
  cat(sprintf("  Wrote %s  [%d cell types x %d samples]\n",
              basename(out_file),
              nrow(deconv_results[[nm]]),
              ncol(deconv_results[[nm]]) - 1))
}

# -------- 4. Brain-immune ssGSEA via GSVA --------
cat("\n=== [4] Brain-immune ssGSEA scoring ===\n")
gmt_path <- file.path(ROOT, "data", "brain_immune_signatures.gmt")
sig_lines <- readLines(gmt_path)
sig_list <- lapply(sig_lines, function(l) {
  parts <- strsplit(l, "\t", fixed = TRUE)[[1]]
  list(name = parts[1], desc = parts[2], genes = parts[-(1:2)])
})
gene_sets <- setNames(lapply(sig_list, function(s) s$genes),
                       sapply(sig_list, function(s) s$name))
cat(sprintf("  Loaded %d gene signatures from %s\n", length(gene_sets), basename(gmt_path)))
cat("  Signatures:\n")
for (n in names(gene_sets)) {
  in_mat <- sum(gene_sets[[n]] %in% rownames(tpm_mat))
  cat(sprintf("    %-30s  %d/%d genes present in TPM matrix\n", n, in_mat, length(gene_sets[[n]])))
}

cat("\n  Running ssGSEA (GSVA::gsva, method='ssgsea') ...\n")
t0 <- Sys.time()
# Use new ssgseaParam API if available (GSVA >= 1.50), else legacy
ssgsea_score <- tryCatch({
  param <- GSVA::ssgseaParam(exprData = tpm_mat,
                              geneSets = gene_sets,
                              normalize = TRUE)
  GSVA::gsva(param)
}, error = function(e) {
  cat("  Fallback to legacy GSVA API\n")
  GSVA::gsva(tpm_mat, gene_sets, method = "ssgsea", ssgsea.norm = TRUE)
})
cat(sprintf("  done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(sprintf("  ssGSEA score matrix: %d signatures x %d samples\n",
            nrow(ssgsea_score), ncol(ssgsea_score)))

# Save: signatures x samples (rows = signatures)
ssgsea_df <- as.data.frame(ssgsea_score) |>
  tibble::rownames_to_column("Signature")
out_ssgsea <- file.path(OUT, "ssGSEA_brain_immune_scores.tsv")
fwrite(ssgsea_df, out_ssgsea, sep = "\t")
cat(sprintf("  Wrote %s\n", basename(out_ssgsea)))

# Also save long-format for plotting
ssgsea_long <- as.data.frame(ssgsea_score) |>
  tibble::rownames_to_column("Signature") |>
  tidyr::pivot_longer(-Signature, names_to = "Kids_First_Biospecimen_ID", values_to = "score")
fwrite(ssgsea_long, file.path(OUT, "ssGSEA_brain_immune_scores_long.tsv"), sep = "\t")

# -------- 5. Cross-method correlation (myeloid focus) --------
cat("\n=== [5] Cross-method correlation: macrophage/myeloid signal ===\n")
# Pull macrophage/myeloid columns from each method
get_score <- function(df, pattern) {
  # df has 'cell_type' first col and one col per sample
  match_rows <- grep(pattern, df$cell_type, value = TRUE, ignore.case = TRUE)
  if (length(match_rows) == 0) return(NULL)
  cat(sprintf("    matching rows: %s\n", paste(match_rows, collapse=", ")))
  sub <- df[df$cell_type %in% match_rows, , drop = FALSE]
  # Sum across matching rows
  vals <- colSums(as.matrix(sub[, -1]))
  data.frame(sample = names(vals), score = unname(vals))
}

cat("  MCP-counter macrophages/monocytes:\n")
mcp_mac <- get_score(deconv_results$mcp_counter, "Macrophage|Monocytic")
cat("  quanTIseq macrophages:\n")
qts_mac <- get_score(deconv_results$quantiseq, "Macrophage")
cat("  xCell macrophages:\n")
xc_mac  <- get_score(deconv_results$xcell, "^Macrophages|^Monocytes")
cat("  Klemm 2020 microglia ssGSEA:\n")
mg_ss   <- data.frame(sample = colnames(ssgsea_score),
                       score = ssgsea_score["Microglia_Klemm2020", ])
cat("  MDM (Klemm) ssGSEA:\n")
mdm_ss  <- data.frame(sample = colnames(ssgsea_score),
                       score = ssgsea_score["MDM_Klemm2020", ])

# Merge and compute pairwise Spearman
merged <- Reduce(function(a, b) merge(a, b, by = "sample"),
                 list(rename_with(mcp_mac, ~ c("sample", "MCP_macro_mono")),
                      rename_with(qts_mac, ~ c("sample", "QTS_macro")),
                      rename_with(xc_mac,  ~ c("sample", "xCell_macro_mono")),
                      rename_with(mg_ss,   ~ c("sample", "Microglia_Klemm_ssGSEA")),
                      rename_with(mdm_ss,  ~ c("sample", "MDM_Klemm_ssGSEA"))))

cor_mat <- cor(merged[, -1], method = "spearman")
cat("\n  Spearman correlation matrix (myeloid signals):\n")
print(round(cor_mat, 3))

fwrite(as.data.frame(cor_mat) |> tibble::rownames_to_column("method"),
       file.path(OUT, "myeloid_crossmethod_spearman.tsv"), sep = "\t")
fwrite(merged, file.path(OUT, "myeloid_signals_merged.tsv"), sep = "\t")

# -------- 6. Quick summary report --------
cat("\n=== [6] Summary ===\n")
cat(sprintf("  immunedeconv methods run: %s\n",
            paste(names(deconv_results), collapse = ", ")))
cat(sprintf("  ssGSEA signatures scored: %d (output: ssGSEA_brain_immune_scores.tsv)\n",
            nrow(ssgsea_score)))
cat("\nAll outputs written to:\n", OUT, "\n", sep = "")

cat("\nDone.\n")
