#!/usr/bin/env Rscript
# =============================================================================
# 05_signature_vs_celltype.R
# Score each TCGA bulk sample by the cisplatin-response gene signatures
# (ssGSEA + GSVA), then correlate every signature score against every estimated
# cell-type proportion from Stage 3. Output = signature x cell-type heatmap.
#
# Run on rbgo-server2:
#   source /userhome/conda/etc/profile.d/conda.sh && conda activate hgsoc_pipeline
#   pwd                          # confirm you are in the repo root
#   Rscript 05_signature_vs_celltype.R
#
# The machinery here is the same one validated on simulated data (planted
# score<->proportion links recovered by both methods; flat negative control
# stayed null). This script swaps the sandbox reimplementation for the real
# Bioconductor GSVA/ssGSEA and points it at real inputs.
# =============================================================================

suppressPackageStartupMessages({
  library(GSVA)
  library(jsonlite)
})
source("R/adapters.R")                           # harmonize_bulk_to_symbols()
set.seed(1000)                                   # pipeline seed convention

# ------------------------------------------------------------------- CONFIG
# Edit these paths to your real objects. Defaults follow the repo layout.
CFG <- list(
  # bulk expression matrix: genes (rows) x samples (cols). RDS of a matrix,
  # or a SummarizedExperiment (an assay is pulled), or a TSV/TXT table whose
  # first column holds the gene IDs. Gene IDs must be SYMBOLS to match the
  # signatures; versioned Ensembl IDs from a TSV are harmonized on load.
  bulk_path       = "../tcga_bulk_data/gene_count_matrix.tsv",
  bulk_assay      = "logcpm",                    # assay name if a SummarizedExperiment; else ignored
  bulk_is_counts  = TRUE,                         # TRUE -> raw integer counts (kcdf="Poisson")

  # Stage 3 deconvolution output: cell-type proportions.
  # CSV or RDS; samples x celltypes OR celltypes x samples (auto-oriented).
  prop_path       = "results/03_estimated_fractions.csv",
  prop_engine     = "NNLS",                       # label only, for the figure caption (e.g. "BayesPrism")

  # gene sets produced earlier
  genesets_path   = "cisplatin_response_genes_HGSOC.json",
  genesets_use    = NULL,                          # NULL = all sets in the JSON; or c("resistance_up_core","matrix_CAF",...)

  # scoring
  min_set_size    = 2,
  ssgsea_alpha    = 0.25,
  gsva_max_diff   = TRUE,

  # stats
  cor_method      = "spearman",
  clr_pseudocount = 1e-6,

  # output
  out_dir         = "results",
  fig_dir         = "figs",
  tag             = "tcga_ov"
)

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$fig_dir, showWarnings = FALSE, recursive = TRUE)
logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

# ------------------------------------------------------------------- load bulk
# TSV/TXT branch: the real TCGA-OV matrix ships as gene_count_matrix.tsv with
# versioned Ensembl rownames (ENSG...\.N), so it is routed through the SAME
# harmonize_bulk_to_symbols() that 03_deconvolution.R uses -- version strip,
# org.Hs.eg.db ENSEMBL->SYMBOL, duplicate symbols collapsed BY SUM. Matching
# Stage 3 exactly matters: these scores are correlated against Stage 3's
# proportions, so the two must sit in the same gene space.
load_bulk <- function(path, assay, is_counts) {
  stopifnot(file.exists(path))
  if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      df <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE)
    } else {
      df <- read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    }
    gene_ids <- as.character(df[[1]])
    mat <- as.matrix(df[, -1, drop = FALSE])
    rownames(mat) <- gene_ids
    storage.mode(mat) <- "double"
    mat <- harmonize_bulk_to_symbols(mat)
    storage.mode(mat) <- "double"
    return(mat)
  }
  obj <- readRDS(path)
  if (inherits(obj, "SummarizedExperiment")) {
    stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))
    a <- if (assay %in% SummarizedExperiment::assayNames(obj)) assay else 1L
    mat <- as.matrix(SummarizedExperiment::assay(obj, a))
  } else {
    mat <- as.matrix(obj)
  }
  storage.mode(mat) <- "double"
  mat
}

bulk <- load_bulk(CFG$bulk_path, CFG$bulk_assay, CFG$bulk_is_counts)
logmsg(sprintf("bulk: %d genes x %d samples", nrow(bulk), ncol(bulk)))

# ------------------------------------------------------------------- load proportions
load_props <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    p <- as.matrix(readRDS(path))
  } else {
    df <- read.csv(path, row.names = 1, check.names = FALSE)
    p <- as.matrix(df)
  }
  storage.mode(p) <- "double"
  p
}
props <- load_props(CFG$prop_path)

# orient so rows = samples, cols = celltypes, using overlap with bulk sample IDs
bulk_samples <- colnames(bulk)
if (length(intersect(rownames(props), bulk_samples)) <
    length(intersect(colnames(props), bulk_samples))) {
  props <- t(props)
}
logmsg(sprintf("proportions: %d samples x %d cell types (%s)",
               nrow(props), ncol(props), paste(colnames(props), collapse = ", ")))

# ------------------------------------------------------------------- align samples
common <- intersect(bulk_samples, rownames(props))
if (length(common) < 3) stop("Fewer than 3 samples shared between bulk and proportions. Check sample IDs.")
if (length(common) < 0.5 * min(ncol(bulk), nrow(props)))
  warning(sprintf("Only %d shared samples; check ID formatting (e.g. TCGA barcodes vs aliquot IDs).", length(common)))
bulk  <- bulk[, common, drop = FALSE]
props <- props[common, , drop = FALSE]
logmsg(sprintf("aligned on %d shared samples", length(common)))

# ------------------------------------------------------------------- gene sets
gj <- fromJSON(CFG$genesets_path, simplifyVector = TRUE)
genesets <- gj$genesets
if (!is.null(CFG$genesets_use)) genesets <- genesets[CFG$genesets_use]
genesets <- lapply(genesets, function(x) unique(as.character(x)))

# gene-ID sanity: signatures are HGNC symbols
present <- rownames(bulk)
overlap_frac <- mean(unlist(genesets) %in% present)
if (overlap_frac < 0.3)
  warning(sprintf("Only %.0f%% of signature genes found in bulk rownames. Are bulk IDs Ensembl rather than symbols? Map first.",
                  100 * overlap_frac))

# restrict to present genes and drop sets below min size
genesets <- lapply(genesets, function(g) intersect(g, present))
kept <- names(genesets)[vapply(genesets, length, 1L) >= CFG$min_set_size]
dropped <- setdiff(names(genesets), kept)
if (length(dropped))
  logmsg(sprintf("dropped %d set(s) below min size %d: %s",
                 length(dropped), CFG$min_set_size, paste(dropped, collapse = ", ")))
genesets <- genesets[kept]
logmsg(sprintf("scoring %d gene sets", length(genesets)))

# ------------------------------------------------------------------- score (GSVA 1.50 param API, legacy fallback)
kcdf <- if (CFG$bulk_is_counts) "Poisson" else "Gaussian"

run_scores <- function(mat, gs, kind) {
  # kind %in% c("gsva","ssgsea"); returns sets x samples matrix
  out <- tryCatch({
    if (kind == "gsva") {
      gp <- GSVA::gsvaParam(exprData = mat, geneSets = gs, kcdf = kcdf,
                            maxDiff = CFG$gsva_max_diff, minSize = CFG$min_set_size)
    } else {
      gp <- GSVA::ssgseaParam(exprData = mat, geneSets = gs,
                              alpha = CFG$ssgsea_alpha, normalize = TRUE,
                              minSize = CFG$min_set_size)
    }
    GSVA::gsva(gp, verbose = FALSE)
  }, error = function(e) {
    logmsg(sprintf("param API failed for %s (%s); using legacy gsva() signature", kind, conditionMessage(e)))
    if (kind == "gsva") {
      GSVA::gsva(mat, gs, method = "gsva", kcdf = kcdf,
                 mx.diff = CFG$gsva_max_diff, min.sz = CFG$min_set_size, verbose = FALSE)
    } else {
      GSVA::gsva(mat, gs, method = "ssgsea", ssgsea.norm = TRUE,
                 min.sz = CFG$min_set_size, verbose = FALSE)
    }
  })
  as.matrix(out)
}

logmsg("running GSVA ...");   sc_gsva   <- run_scores(bulk, genesets, "gsva")
logmsg("running ssGSEA ..."); sc_ssgsea <- run_scores(bulk, genesets, "ssgsea")

saveRDS(list(gsva = sc_gsva, ssgsea = sc_ssgsea),
        file.path(CFG$out_dir, sprintf("%s_signature_scores.rds", CFG$tag)))
write.csv(t(sc_gsva),   file.path(CFG$out_dir, sprintf("%s_scores_GSVA.csv",   CFG$tag)))
write.csv(t(sc_ssgsea), file.path(CFG$out_dir, sprintf("%s_scores_ssGSEA.csv", CFG$tag)))

# ------------------------------------------------------------------- CLR proportions (compositional-safe)
clr <- function(P, eps) {
  L <- log(P + eps)
  L - rowMeans(L)
}
props_clr <- clr(props, CFG$clr_pseudocount)
# drop zero-variance cell types (e.g. an undetectable adipocyte blind spot)
nzv <- apply(props_clr, 2, function(x) sd(x) > 1e-8)
if (any(!nzv))
  logmsg(sprintf("dropping zero-variance cell type(s) from correlation: %s",
                 paste(colnames(props_clr)[!nzv], collapse = ", ")))
props_clr <- props_clr[, nzv, drop = FALSE]

# ------------------------------------------------------------------- correlation + BH
corr_block <- function(scores) {           # scores: sets x samples
  S <- t(scores)                            # samples x sets
  sets <- colnames(S); cts <- colnames(props_clr)
  rho <- matrix(NA_real_, length(sets), length(cts), dimnames = list(sets, cts))
  pv  <- rho
  for (s in sets) for (c in cts) {
    ct <- suppressWarnings(cor.test(S[, s], props_clr[, c], method = CFG$cor_method))
    rho[s, c] <- unname(ct$estimate); pv[s, c] <- ct$p.value
  }
  padj <- matrix(p.adjust(pv, "BH"), nrow(pv), ncol(pv), dimnames = dimnames(pv))
  list(rho = rho, pv = pv, padj = padj)
}
res_gsva   <- corr_block(sc_gsva)
res_ssgsea <- corr_block(sc_ssgsea)

# tidy long output
tidy_long <- function(res, method) {
  df <- expand.grid(signature = rownames(res$rho), celltype = colnames(res$rho),
                    stringsAsFactors = FALSE)
  df$spearman_rho <- mapply(function(s, c) res$rho[s, c], df$signature, df$celltype)
  df$p_value      <- mapply(function(s, c) res$pv[s, c],  df$signature, df$celltype)
  df$padj_BH      <- mapply(function(s, c) res$padj[s, c], df$signature, df$celltype)
  df$method <- method
  df
}
long <- rbind(tidy_long(res_gsva, "GSVA"), tidy_long(res_ssgsea, "ssGSEA"))
write.csv(long, file.path(CFG$out_dir, sprintf("%s_signature_celltype_correlations.csv", CFG$tag)), row.names = FALSE)
write.csv(res_gsva$rho,   file.path(CFG$out_dir, sprintf("%s_corr_rho_GSVA.csv",   CFG$tag)))
write.csv(res_ssgsea$rho, file.path(CFG$out_dir, sprintf("%s_corr_rho_ssGSEA.csv", CFG$tag)))

# ------------------------------------------------------------------- heatmap
stars <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", "")))

draw_heatmap <- function(res, method, file) {
  disp <- matrix(sprintf("%.2f%s", res$rho, stars(res$padj)),
                 nrow(res$rho), dimnames = dimnames(res$rho))
  cap <- sprintf("%s: signature score vs cell-type proportion (%s)  |  Spearman on CLR-proportions, BH-adjusted  |  n=%d TCGA-OV samples",
                 method, CFG$prop_engine, length(common))
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pal <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
    pheatmap::pheatmap(res$rho, color = pal, breaks = seq(-1, 1, length.out = 102),
                       display_numbers = disp, number_color = "black",
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       main = cap, filename = file, width = 8, height = 7)
  } else {                                   # base-R fallback
    png(file, width = 1100, height = 900, res = 130)
    par(mar = c(9, 12, 4, 2))
    image(seq_len(ncol(res$rho)), seq_len(nrow(res$rho)), t(res$rho[nrow(res$rho):1, , drop = FALSE]),
          axes = FALSE, xlab = "", ylab = "", zlim = c(-1, 1),
          col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), main = cap)
    axis(1, seq_len(ncol(res$rho)), colnames(res$rho), las = 2)
    axis(2, seq_len(nrow(res$rho)), rev(rownames(res$rho)), las = 2)
    dev.off()
  }
}
draw_heatmap(res_gsva,   "GSVA",   file.path(CFG$fig_dir, sprintf("%s_heatmap_GSVA.png",   CFG$tag)))
draw_heatmap(res_ssgsea, "ssGSEA", file.path(CFG$fig_dir, sprintf("%s_heatmap_ssGSEA.png", CFG$tag)))

logmsg("done.")
logmsg("outputs:")
logmsg(sprintf("  %s/%s_signature_celltype_correlations.csv (tidy long)", CFG$out_dir, CFG$tag))
logmsg(sprintf("  %s/%s_corr_rho_{GSVA,ssGSEA}.csv", CFG$out_dir, CFG$tag))
logmsg(sprintf("  %s/%s_scores_{GSVA,ssGSEA}.csv", CFG$out_dir, CFG$tag))
logmsg(sprintf("  %s/%s_heatmap_{GSVA,ssGSEA}.png", CFG$fig_dir, CFG$tag))

# sessionInfo for reproducibility
writeLines(capture.output(sessionInfo()),
           file.path(CFG$out_dir, sprintf("%s_sessionInfo_stage05.txt", CFG$tag)))
