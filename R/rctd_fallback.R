# ============================================================================
# R/rctd_fallback.R — deconvolve spatial spots into cell-type fractions.
# PRIMARY: RCTD from spacexr (config tools$use_rctd = true; needs the GitHub pkg).
# FALLBACK: run the same NNLS deconvolution (R/deconvolve.R) independently on
#           each spot's summed counts. A spot is just a tiny bulk sample, so the
#           reference-based mixing logic is identical to Aim 3.
# ============================================================================

# PRIMARY RCTD wrapper (guarded). spot_counts: genes x spots; coords: x,y df.
# Permissive UMI/counts thresholds so spots are not silently dropped, and we
# realign RCTD's output to the ORIGINAL spot order (RCTD can reorder/filter).
spatial_deconv_rctd <- function(signatures, spot_counts, sce_ref, coords, cfg) {
  if (!requireNamespace("spacexr", quietly = TRUE)) return(NULL)
  sp <- as_spatialrna(coords, spot_counts, sce_ref)
  ref <- spacexr::Reference(sp$ref_counts, sp$ref_types)
  puck <- spacexr::SpatialRNA(coords = sp$coords, counts = sp$spot_counts,
                              use_fake_coords = FALSE)
  rctd <- spacexr::create.RCTD(puck, ref, max_cores = cfg$threads,
                               CELL_MIN_INSTANCE = 5, UMI_min = 1,
                               counts_MIN = 1, UMI_min_sigma = 1)
  rctd <- spacexr::run.RCTD(rctd, doublet_mode = "full")
  w <- as.matrix(rctd@results$weights)
  w <- sweep(w, 1, pmax(rowSums(w), 1e-8), "/")     # normalise to fractions
  out <- matrix(NA_real_, ncol(spot_counts), ncol(w),
                dimnames = list(colnames(spot_counts), colnames(w)))
  common <- intersect(rownames(out), rownames(w))
  out[common, ] <- w[common, ]                      # align to input spot order
  out
}

# dispatcher -> list(props = spots x types, method)
run_spatial_deconv <- function(signatures, spot_counts, sce_ref, coords, cfg,
                               markers = NULL) {
  if (isTRUE(cfg$tools$use_rctd) &&
      requireNamespace("spacexr", quietly = TRUE)) {
    log_msg("spatial deconvolution: RCTD (PRIMARY)")
    res <- try(spatial_deconv_rctd(signatures, spot_counts, sce_ref, coords, cfg),
               silent = TRUE)
    if (!inherits(res, "try-error") && !is.null(res)) {
      common <- intersect(colnames(signatures), colnames(res))
      return(list(props = res[, common, drop = FALSE], method = "RCTD"))
    }
    log_msg("RCTD unavailable/failed -> using fallback")
  }
  log_msg("[FALLBACK] spatial deconvolution: per-spot NNLS")
  list(props = deconvolve_nnls(signatures, spot_counts, markers = markers),
       method = "NNLS [FALLBACK]")
}
