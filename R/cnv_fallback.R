# ============================================================================
# R/cnv_fallback.R — separate malignant (aneuploid) from stroma (diploid).
# PRIMARY: inferCNV (config tools$use_infercnv = true; needs JAGS + infercnv).
# FALLBACK: the SAME core idea inferCNV uses, minus the HMM:
#   1. log-normalise; 2. subtract the mean of a known DIPLOID reference
#      (stroma/immune) gene-by-gene; 3. smooth along genome order within each
#      chromosome (moving average = "are whole regions shifted up/down?");
#   4. per-cell aneuploidy score = mean |smoothed deviation|; 5. call malignant
#      if the score exceeds the reference cells' null distribution.
# REFERENCE-BASELINE RISK (noted in the dissertation): if the reference is
# mislabelled or stressed, every call is suspect. We check that reference cells
# themselves score low and warn loudly otherwise.
# ============================================================================

# edge-safe running mean (averages whatever part of the window is available)
.runmean <- function(x, w) {
  n <- length(x)
  if (w < 2 || n < 2) return(x)
  half <- (w - 1) %/% 2
  cs <- c(0, cumsum(x))
  out <- numeric(n)
  for (i in seq_len(n)) {
    lo <- max(1, i - half); hi <- min(n, i + half)
    out[i] <- (cs[hi + 1] - cs[lo]) / (hi - lo + 1)
  }
  out
}

# numeric chromosome index for ordering ("chr1".."chr22","chrX" -> 1..23)
.chr_index <- function(chr) {
  v <- sub("^chr", "", chr)
  v[v == "X"] <- "23"; v[v == "Y"] <- "24"
  as.integer(v)
}

# build the per-cell CNV signal matrix (genes x cells, smoothed, ref-centred)
cnv_signal <- function(sce, ref_cells, window) {
  logc <- get_logcounts(sce)
  go <- gene_order_df(sce)
  ord <- order(.chr_index(go$chr), go$start)        # genome order
  logc <- logc[ord, , drop = FALSE]; go <- go[ord, ]
  ref_mean <- rowMeans(logc[, ref_cells, drop = FALSE])
  centred <- logc - ref_mean                        # diploid-relative
  smoothed <- centred
  for (ch in unique(go$chr)) {                       # smooth within chromosome
    rows <- which(go$chr == ch)
    if (length(rows) >= 2)
      smoothed[rows, ] <- apply(centred[rows, , drop = FALSE], 2,
                                .runmean, w = min(window, length(rows)))
  }
  list(smoothed = smoothed, gene_order = go)
}

# FALLBACK caller -------------------------------------------------------------
call_malignant_fallback <- function(sce, ref_types, cfg) {
  labels <- as.character(SummarizedExperiment::colData(sce)$celltype)
  ref_cells <- which(labels %in% ref_types)
  G <- nrow(sce)
  window <- max(5, round(G / 40)); if (window %% 2 == 0) window <- window + 1

  sig <- cnv_signal(sce, ref_cells, window)
  score <- colMeans(abs(sig$smoothed))               # per-cell aneuploidy

  # null from the diploid reference; threshold = robust upper tail
  ref_score <- score[ref_cells]
  thr <- max(stats::quantile(ref_score, 0.99, names = FALSE),
             mean(ref_score) + 3 * stats::sd(ref_score))
  call <- ifelse(score > thr, "malignant", "diploid")

  # reference-baseline sanity: what fraction of the reference itself trips?
  ref_fpr <- mean(score[ref_cells] > thr)

  # arm-level mean signal in called-malignant cells (to compare vs planted map)
  mal_cells <- which(call == "malignant")
  arm_profile <- tapply(seq_len(nrow(sig$smoothed)), sig$gene_order$arm,
                        function(rows) mean(sig$smoothed[rows, mal_cells]))
  list(call = call, score = score, threshold = thr, ref_fpr = ref_fpr,
       arm_profile = arm_profile, method = "windowed-expression [FALLBACK]")
}

# inferCNV imports rjags, which needs JAGS. JAGS may be a user-local install with
# no registry entry, and rjags wants JAGS_HOME = the version dir (it appends
# /x64/bin itself). Auto-detect and set it if not already valid.
ensure_jags <- function() {
  jh <- Sys.getenv("JAGS_HOME")
  if (nzchar(jh) && dir.exists(file.path(jh, "x64", "bin"))) return(invisible(jh))
  roots <- c(file.path(Sys.getenv("LOCALAPPDATA"), "Programs", "JAGS"),
             "C:/Program Files/JAGS", "C:/Program Files (x86)/JAGS")
  for (r in roots) if (dir.exists(r)) {
    vers <- list.dirs(r, recursive = FALSE)
    vers <- vers[dir.exists(file.path(vers, "x64", "bin"))]
    if (length(vers)) { jh <- sort(vers, decreasing = TRUE)[1]
                        Sys.setenv(JAGS_HOME = jh); return(invisible(jh)) }
  }
  invisible("")
}

# PRIMARY inferCNV wrapper (guarded; reuses the scoring on inferCNV's denoised
# matrix). Returns same structure, or NULL if infercnv is unavailable/fails.
call_malignant_infercnv <- function(sce, ref_types, cfg, out_dir = "results/infercnv") {
  ensure_jags()
  if (!requireNamespace("infercnv", quietly = TRUE)) return(NULL)
  files <- write_infercnv_files(sce, out_dir, annot_col = "celltype")
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = files$counts,
    annotations_file  = files$annotations,
    gene_order_file   = files$gene_order,
    ref_group_names   = ref_types)
  obj <- infercnv::run(obj, cutoff = 0.1, out_dir = out_dir,
                       cluster_by_groups = TRUE, denoise = TRUE,
                       HMM = FALSE, num_threads = cfg$threads,
                       no_plot = TRUE)               # skip plots: faster, fewer failure points
  # Score inferCNV's denoised matrix exactly like the fallback. IMPORTANT:
  # inferCNV may reorder/drop cells, so index reference cells BY NAME and realign
  # the result to the original SCE cell order (cells inferCNV dropped -> diploid).
  expr <- obj@expr.data
  labels <- as.character(SummarizedExperiment::colData(sce)$celltype)
  ref_names <- colnames(sce)[labels %in% ref_types]
  ref_in_expr <- which(colnames(expr) %in% ref_names)
  s <- colMeans(abs(expr - rowMeans(expr[, ref_in_expr, drop = FALSE])))
  thr <- mean(s[ref_in_expr]) + 3 * stats::sd(s[ref_in_expr])
  score <- setNames(rep(min(s), ncol(sce)), colnames(sce))   # dropped cells -> low
  score[names(s)] <- s
  call <- ifelse(score > thr, "malignant", "diploid")
  ref_cells_all <- which(labels %in% ref_types)
  list(call = call, score = score, threshold = thr,
       ref_fpr = mean(score[ref_cells_all] > thr),
       arm_profile = NULL, method = "inferCNV")
}

# dispatcher
run_cnv <- function(sce, ref_types, cfg) {
  if (isTRUE(cfg$tools$use_infercnv)) ensure_jags()   # set JAGS_HOME before loading
  if (isTRUE(cfg$tools$use_infercnv) &&
      requireNamespace("infercnv", quietly = TRUE)) {
    log_msg("CNV calling: inferCNV (PRIMARY)")
    res <- try(call_malignant_infercnv(sce, ref_types, cfg), silent = TRUE)
    if (!inherits(res, "try-error") && !is.null(res)) return(res)
    log_msg("inferCNV unavailable/failed -> using fallback")
  }
  log_msg("[FALLBACK] CNV calling: windowed expression vs diploid reference")
  call_malignant_fallback(sce, ref_types, cfg)
}
