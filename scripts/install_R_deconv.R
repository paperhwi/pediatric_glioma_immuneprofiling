# Install immunedeconv + GSVA stack
# Run once, takes ~20-40 minutes for first-time installation

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = parallel::detectCores() - 1)

cat("=== Installing core dependencies ===\n")
core_cran <- c("remotes", "BiocManager", "data.table", "tidyverse",
               "matrixStats", "Matrix", "Rcpp", "RcppArmadillo")
to_install <- setdiff(core_cran, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

cat("\n=== Installing Bioconductor packages (GSVA, GSEABase, SummarizedExperiment) ===\n")
bioc_pkgs <- c("GSVA", "GSEABase", "SummarizedExperiment", "preprocessCore", "limma")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}

cat("\n=== Installing immunedeconv (includes xCell, quanTIseq, MCP-counter) ===\n")
if (!requireNamespace("immunedeconv", quietly = TRUE)) {
  remotes::install_github("omnideconv/immunedeconv",
                          dependencies = TRUE, upgrade = "never", quiet = FALSE)
}

cat("\n=== Final verification ===\n")
for (p in c("GSVA", "GSEABase", "immunedeconv", "xCell", "MCPcounter",
            "quantiseqr", "EPIC", "preprocessCore", "limma")) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("%-20s %s\n", p, ifelse(ok, "OK", "MISSING")))
}
