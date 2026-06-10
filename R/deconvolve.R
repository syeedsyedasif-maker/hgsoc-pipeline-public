# ============================================================================
# R/deconvolve.R — bulk deconvolution.
# PRIMARY: BayesPrism (config tools$use_bayesprism = true; needs the GitHub pkg).
# FALLBACK: non-negative least squares (NNLS) on per-cell-type signatures, then
#           scale to sum to 1. This is the same idea every reference-based
#           deconvolution method uses (express bulk as a non-negative mix of
#           cell-type profiles); it just skips BayesPrism's Bayesian machinery.
# Every fallback log line is tagged [FALLBACK] so results are never silently
# attributed to the named tool.
# ============================================================================

# --- pure-R NNLS via projected gradient descent (no compiled dependency) -----
# Solve  min ||S w - b||^2  s.t. w >= 0 , then normalise w to sum 1 = fractions.
.nnls_pg <- function(AtA, Atb, L, n, iters = 800) {
  w <- rep(1 / n, n)
  for (i in seq_len(iters)) {
    grad <- AtA %*% w - Atb
    w <- w - grad / L
    w[w < 0] <- 0
  }
  s <- sum(w)
  if (s > 0) w <- w / s
  as.numeric(w)
}

# Pick the most cell-type-specific genes from the reference signatures. For each
# type we keep the genes where it is highest above the next-best type (good
# markers). This conditions the NNLS so similar stroma/immune types don't get
# merged, exactly as marker-based methods (CIBERSORT etc.) do.
select_markers <- function(signatures, n_per_type = 50) {
  L <- log1p(as.matrix(signatures))
  markers <- character(0)
  for (t in colnames(L)) {
    others <- L[, setdiff(colnames(L), t), drop = FALSE]
    score <- L[, t] - apply(others, 1, max)       # enrichment over next-best type
    top <- names(sort(score, decreasing = TRUE))[seq_len(min(n_per_type, nrow(L)))]
    markers <- union(markers, top)
  }
  markers
}

# signatures: genes x types (mean per-cell profile). bulk: genes x samples.
# We FIRST normalise each signature and each bulk sample to sum 1 over ALL shared
# genes -- this turns them into expression-FRACTION profiles, so the model
# becomes a simplex mixture  b = M %*% fractions  and NNLS recovers fractions.
# Only THEN do we subset to marker genes (without renormalising, which would
# break the identity): markers just down-weight uninformative, confound-noisy
# genes. Recovered weights are renormalised to sum to 1. This is the same logic
# reference-based methods (CIBERSORT, BayesPrism) use.
deconvolve_nnls <- function(signatures, bulk, markers = NULL) {
  genes <- intersect(rownames(signatures), rownames(bulk))
  S <- as.matrix(signatures[genes, , drop = FALSE])
  B <- as.matrix(bulk[genes, , drop = FALSE])
  S <- sweep(S, 2, pmax(colSums(S), 1e-8), "/")    # profiles over ALL genes
  B <- sweep(B, 2, pmax(colSums(B), 1e-8), "/")
  if (!is.null(markers)) {                         # subset AFTER normalising
    mk <- intersect(markers, genes)
    S <- S[mk, , drop = FALSE]; B <- B[mk, , drop = FALSE]
  }
  use_pkg <- requireNamespace("nnls", quietly = TRUE)   # proper active-set NNLS
  AtA <- crossprod(S)                                   # (only needed by fallback)
  L <- max(eigen(AtA, only.values = TRUE)$values)
  if (!is.finite(L) || L <= 0) L <- 1
  out <- matrix(0, ncol(B), ncol(S),
                dimnames = list(colnames(B), colnames(S)))
  for (j in seq_len(ncol(B))) {
    w <- if (use_pkg) nnls::nnls(S, B[, j])$x
         else .nnls_pg(AtA, crossprod(S, B[, j]), L, ncol(S))
    s <- sum(w); if (s > 0) w <- w / s
    out[j, ] <- w
  }
  out
}

# --- BayesPrism primary path (guarded; runs only if package present) ---------
deconvolve_bayesprism <- function(sce_ref, bulk_mat, type_col = "celltype",
                                  n_cores = 1) {
  if (!requireNamespace("BayesPrism", quietly = TRUE)) return(NULL)
  bp <- as_bayesprism(sce_ref, bulk_mat, type_col = type_col)
  prism <- BayesPrism::new.prism(
    reference        = round(bp$reference),       # BayesPrism expects counts
    mixture          = round(bp$bulk),
    input.type       = "count.matrix",
    cell.type.labels = bp$cell.type.labels,
    cell.state.labels = bp$cell.type.labels,
    key              = NULL)
  res <- BayesPrism::run.prism(prism = prism, n.cores = n_cores)
  theta <- BayesPrism::get.fraction(bp = res, which.theta = "final",
                                    state.or.type = "type")
  theta  # samples x types
}

# --- InstaPrism: fast BayesPrism-equivalent (used as the iteration engine) ----
# bulk: genes x samples. Returns samples x types, or NULL if unavailable.
deconvolve_instaprism <- function(sce_ref, bulk, type_col = "celltype", n_cores = 1) {
  if (!requireNamespace("InstaPrism", quietly = TRUE)) return(NULL)
  # InstaPrism_legacy() calls magrittr's %>% and pbapply's pboptions() without
  # importing them, so attach those packages to make the symbols resolvable.
  for (dep in c("magrittr", "pbapply"))
    if (requireNamespace(dep, quietly = TRUE))
      suppressWarnings(suppressMessages(require(dep, character.only = TRUE, quietly = TRUE)))
  ref <- round(get_counts(sce_ref))                       # genes x cells
  labels <- as.character(SummarizedExperiment::colData(sce_ref)[[type_col]])
  common <- intersect(rownames(ref), rownames(bulk))
  # current InstaPrism() wants a prebuilt reference; InstaPrism_legacy() takes
  # raw single-cell data directly, which is what we have.
  res <- InstaPrism::InstaPrism_legacy(input_type = "raw",
           sc_Expr   = ref[common, , drop = FALSE],
           bulk_Expr = round(as.matrix(bulk[common, , drop = FALSE])),
           cell.type.labels = labels, cell.state.labels = labels,
           n.core = n_cores)
  t(res@Post.ini.ct@theta)                                # samples x types
}

# --- MuSiC: naive cross-subject comparator. bulk: genes x samples ------------
# Needs subject IDs in the reference to model cross-subject variance; we assign
# pseudo-subjects. Returns samples x types, or NULL if unavailable.
deconvolve_music <- function(sce_ref, bulk, type_col = "celltype", n_subjects = 5) {
  if (!requireNamespace("MuSiC", quietly = TRUE)) return(NULL)
  common <- intersect(rownames(sce_ref), rownames(bulk))
  ref <- sce_ref[common, ]
  set.seed(1)
  SummarizedExperiment::colData(ref)$sampleID <-
    paste0("subj", sample(seq_len(n_subjects), ncol(ref), replace = TRUE))
  types <- sort(unique(as.character(SummarizedExperiment::colData(ref)[[type_col]])))
  est <- MuSiC::music_prop(
    bulk.mtx = round(as.matrix(bulk[common, , drop = FALSE])),
    sc.sce = ref, clusters = type_col, samples = "sampleID",
    select.ct = types, verbose = FALSE)
  as.matrix(est$Est.prop.weighted)                        # samples x types
}

# --- dispatcher: pick primary if requested+available, else fallback ----------
# Returns list(props = samples x types, method = character).
run_deconvolution <- function(signatures, bulk, sce_ref, bulk_mat, cfg,
                              markers = NULL) {
  if (isTRUE(cfg$tools$use_bayesprism) &&
      requireNamespace("BayesPrism", quietly = TRUE)) {
    log_msg("deconvolution: BayesPrism (PRIMARY)")
    res <- try(deconvolve_bayesprism(sce_ref, bulk_mat,
                                     n_cores = cfg$threads), silent = TRUE)
    if (!inherits(res, "try-error") && !is.null(res)) {
      common <- intersect(colnames(signatures), colnames(res))
      return(list(props = res[, common, drop = FALSE], method = "BayesPrism"))
    }
    log_msg("BayesPrism unavailable/failed -> using fallback")
  }
  log_msg("[FALLBACK] deconvolution: NNLS on cell-type marker signatures")
  list(props = deconvolve_nnls(signatures, bulk, markers = markers),
       method = "NNLS [FALLBACK]")
}
