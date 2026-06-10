# ============================================================================
# install.R — reproduce the package environment.
# Pinned to R 4.3.x + Bioconductor 3.18 (the release that matches R 4.3.1).
# Tested on R 4.3.1 / Bioconductor 3.18 / Windows 11 WITHOUT Rtools or JAGS:
#   the CORE stack installs from CRAN/Bioconductor binaries (no compiler needed)
#   and the whole pipeline runs via the built-in fallbacks.
#
# The three "primary" tools (inferCNV, BayesPrism, spacexr) are OPTIONAL: they
# need extra system tooling (JAGS for inferCNV; a C/C++ toolchain such as Rtools
# for the GitHub packages). They are attempted only if you opt in:
#     Sys.setenv(INSTALL_OPTIONAL = "1")    # before sourcing this file
# If they are absent the pipeline still runs and clearly logs "[FALLBACK]".
# ============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# 1. Bioconductor manager, pinned release ------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (as.character(BiocManager::version()) != "3.18")
  message("NOTE: expected Bioconductor 3.18 (R 4.3.x). Found ",
          BiocManager::version(), " — versions in renv.lock may differ.")

# 2. CORE packages (binaries on Windows; no compiler required) ----------------
cran_pkgs <- c("yaml", "jsonlite", "data.table", "Matrix", "ggplot2",
               "pheatmap", "gridExtra", "simsurv", "survival", "nnls", "renv")
bioc_pkgs <- c("splatter", "SingleCellExperiment", "SummarizedExperiment",
               "scater", "scran", "bluster", "igraph", "matrixStats")

need <- function(p) !requireNamespace(p, quietly = TRUE)
for (p in cran_pkgs) if (need(p)) install.packages(p, type = "binary")
for (p in bioc_pkgs) if (need(p)) BiocManager::install(p, ask = FALSE, update = FALSE)

# 3. OPTIONAL primary tools (opt in via INSTALL_OPTIONAL=1) -------------------
# This is the TESTED recipe (R 4.3.1 + Rtools43 + JAGS 4.3.2 on Windows).
if (Sys.getenv("INSTALL_OPTIONAL") == "1") {
  message(">> Installing optional primary tools (need Rtools + JAGS)...")
  # Seurat & inferCNV require Matrix >= 1.6.4 (newer than the one shipped with R 4.3.1)
  if (utils::packageVersion("Matrix") < "1.6.4") install.packages("Matrix", type = "binary")
  if (need("Seurat")) install.packages("Seurat", type = "binary")   # inferCNV import
  # inferCNV (+ rjags, which links to JAGS at load time)
  BiocManager::install(c("rjags", "infercnv"), ask = FALSE, update = FALSE)
  # BayesPrism's fiddliest deps first, then the two GitHub packages
  for (p in c("snowfall", "gplots")) if (need(p)) install.packages(p, type = "binary")
  for (p in c("NMF", "BiocParallel")) if (need(p)) BiocManager::install(p, ask = FALSE, update = FALSE)
  if (need("remotes")) install.packages("remotes")
  try(remotes::install_github("Danko-Lab/BayesPrism/BayesPrism", upgrade = "never"))
  try(remotes::install_github("dmcable/spacexr", upgrade = "never"))
  message("NOTE: inferCNV needs JAGS at RUNTIME. If JAGS is a user-local install ",
          "(AppData) with no registry entry, set JAGS_HOME to the JAGS version dir, ",
          "e.g. .../JAGS-4.3.2 (NOT the x64 subdir). The pipeline auto-detects this ",
          "via ensure_jags() in R/cnv_fallback.R, and .Renviron sets it for `cd`-based runs.")
} else {
  message(">> Skipping optional tools (set INSTALL_OPTIONAL=1 to attempt them). ",
          "Fallbacks will be used.")
}

# 4. report -------------------------------------------------------------------
core_ok <- all(vapply(c(cran_pkgs, bioc_pkgs), requireNamespace,
                      logical(1), quietly = TRUE))
cat("\nCORE packages installed:", core_ok, "\n")
for (p in c("infercnv", "BayesPrism", "spacexr"))
  cat(sprintf("  optional %-11s: %s\n", p,
      if (requireNamespace(p, quietly = TRUE)) "present (PRIMARY path)"
      else "absent -> [FALLBACK]"))
if (!core_ok) stop("Some core packages failed to install — see messages above.")
cat("\nEnvironment ready. Run:  Rscript run_all.R\n")
