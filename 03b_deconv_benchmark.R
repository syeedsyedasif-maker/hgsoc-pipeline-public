# ============================================================================
# 03b_deconv_benchmark.R : COMPARATIVE deconvolution benchmark (Aim 3 extra)
# ----------------------------------------------------------------------------
# WHAT:  Deconvolve TWO kinds of bulk with several methods and compare accuracy:
#          (1) "pseudobulk"      : bulk built from the SAME population as the
#                                  single-cell reference (the easy case)
#          (2) "independent draw": bulk built from a SEPARATE splatter batch
#                                  (its own per-gene technical profile), with a
#                                  DIFFERENT composition and extra cell types
#                                  dropped from the reference (the hard case)
#        Methods: NNLS (naive), MuSiC (naive), InstaPrism (fast BayesPrism), and
#        optionally true BayesPrism once for confirmation.
# WHY:   Tests the Hippen-2023 claim that naive methods look great on pseudobulk but
#        DEGRADE under reference<->bulk mismatch, while BayesPrism stays robust.
#        This is a HYPOTHESIS TEST, not a PASS/FAIL gate: a naive method scoring
#        worse on the independent draw is the EXPECTED, desired result.
# INPUT: config.yml (benchmark: block). Generates its own data, does NOT touch
#        00-04 or data/sce.rds.
# OUTPUT: results/03b_benchmark.csv (+ per-type detail), figs/03b_benchmark.png,
#         truth/benchmark_notes.md
# HOW TO READ: the table shows each method's mean Pearson on pseudobulk vs
#        independent draw. The "drop" column (pseudobulk - independent) is the
#        degradation. Expect naive methods to drop MORE than InstaPrism/BayesPrism.
# NOTE:  We realise the "independent draw" with splatter's BATCH model rather than
#        a fresh different-seed simulation ON PURPOSE: a fresh sim re-randomises
#        which gene means which cell type, so the reference would be meaningless
#        and EVERY method would score ~0. Batches keep gene identity fixed while
#        giving the target its own technical profile (= different patient/platform).
# ============================================================================

suppressPackageStartupMessages({
  library(splatter); library(SingleCellExperiment)
  library(SummarizedExperiment); library(ggplot2)
})
source("R/utils.R"); source("R/adapters.R"); source("R/metrics.R"); source("R/deconvolve.R")
cfg <- load_config()
bm  <- cfg$benchmark
init_script("03b_deconv_benchmark.R : comparative deconvolution", cfg)

# ---- 1. simulate ONE population with TWO batches ----------------------------
gnames <- group_levels(cfg)
gp <- unlist(cfg$group_prob[gnames]); gp <- gp / sum(gp)
de_prob   <- setNames(rep(cfg$splatter$de_prob,   length(gnames)), gnames)
de_facLoc <- setNames(rep(cfg$splatter$de_facLoc, length(gnames)), gnames)
de_prob["malignant_omental"]   <- cfg$splatter$omental_de_prob
de_facLoc["malignant_omental"] <- cfg$splatter$omental_de_facLoc
N <- cfg$sizes$n_cells; n1 <- round(N / 2)
log_msg("simulating", cfg$sizes$n_genes, "genes x", N,
        "cells in 2 batches (batch1=reference source, batch2=independent target)")
params <- setParams(newSplatParams(),
  nGenes = cfg$sizes$n_genes, batchCells = c(n1, N - n1),
  batch.facLoc = bm$independent_batch_facLoc,
  batch.facScale = bm$independent_batch_facScale,
  group.prob = as.numeric(gp), de.prob = as.numeric(de_prob),
  de.facLoc = as.numeric(de_facLoc), de.facScale = cfg$splatter$de_facScale,
  seed = cfg$seed)
sim <- splatSimulate(params, method = "groups", verbose = FALSE)
grp_idx <- as.integer(sub("Group", "", as.character(sim$Group)))
sim$celltype <- celltype_of_group(gnames[grp_idx])
counts <- get_counts(sim)
ct <- as.character(sim$celltype); types_all <- sort(unique(ct))
b1 <- sim$Batch == "Batch1"; b2 <- !b1
pool_easy <- split(which(b1), ct[b1])     # reference + pseudobulk source
pool_hard <- split(which(b2), ct[b2])     # independent-draw source

# ---- 2. reference (drop types, downsample) + signatures ---------------------
ref_types <- setdiff(types_all, bm$drop_from_reference)
ncpt <- bm$ref_cells_per_type
ref_cells <- unlist(lapply(ref_types, function(t) {
  idx <- pool_easy[[t]]; if (length(idx) > ncpt) sample(idx, ncpt) else idx }))
sce_ref <- sim[, ref_cells]
signatures <- build_reference_signatures(sce_ref, "celltype")
markers <- select_markers(signatures, bm$markers_per_type)
log_msg("reference types:", paste(ref_types, collapse = ", "),
        "| dropped from ref (present in bulk):", paste(bm$drop_from_reference, collapse = ", "))

# ---- 3. build the two bulks -------------------------------------------------
collapse_ct <- function(gprob) tapply(unlist(gprob), celltype_of_group(names(gprob)), sum)
easy_props <- collapse_ct(cfg$group_prob)[types_all]
hard_props <- collapse_ct(bm$target_group_prob)[types_all]
draw_dirichlet <- function(a) { g <- rgamma(length(a), shape = a, rate = 1); g / sum(g) }

build_bulk <- function(pool, props, n_samples, n_cells) {
  alpha <- cfg$survival$dirichlet_alpha * as.numeric(props)
  true <- matrix(0, n_samples, length(types_all), dimnames = list(paste0("s", 1:n_samples), types_all))
  bulk <- matrix(0, nrow(counts), n_samples,
                 dimnames = list(rownames(counts), paste0("s", 1:n_samples)))
  for (i in seq_len(n_samples)) {
    n_k <- as.integer(rmultinom(1, n_cells, draw_dirichlet(alpha)))
    b <- numeric(nrow(counts)); contrib <- setNames(numeric(length(types_all)), types_all)
    for (j in seq_along(types_all)) {
      idx <- pool[[types_all[j]]]
      if (n_k[j] == 0 || is.null(idx)) next
      cl <- sample(idx, n_k[j], replace = TRUE)
      s <- rowSums(counts[, cl, drop = FALSE]); b <- b + s; contrib[j] <- sum(s)
    }
    bulk[, i] <- b; true[i, ] <- contrib / sum(contrib)
  }
  list(bulk = bulk, true = true)
}
apply_confounds <- function(bulk) {
  if (isTRUE(cfg$confounds$batch_shift$enabled))
    bulk <- bulk * exp(rnorm(nrow(bulk), 0, cfg$confounds$batch_shift$sd))
  if (isTRUE(cfg$confounds$mrna_enrichment$enabled)) {
    k <- round(cfg$confounds$mrna_enrichment$subset_frac * nrow(bulk))
    cap <- sample(nrow(bulk), k); bulk[cap, ] <- bulk[cap, ] * cfg$confounds$mrna_enrichment$factor
  }
  bulk
}
easy <- build_bulk(pool_easy, easy_props, bm$n_bulk, bm$cells_per_bulk)
hard <- build_bulk(pool_hard, hard_props, bm$n_bulk, bm$cells_per_bulk)
easy$bulk <- apply_confounds(easy$bulk); hard$bulk <- apply_confounds(hard$bulk)
log_msg(sprintf("built %d pseudobulk + %d independent-draw samples", bm$n_bulk, bm$n_bulk))

# ---- 4. run each method on BOTH bulks, with a time budget -------------------
budget <- bm$time_budget_sec
run_method <- function(fn) {
  t0 <- Sys.time()
  res <- tryCatch({ setTimeLimit(elapsed = budget, transient = TRUE); v <- fn(); setTimeLimit(); v },
                  error = function(e) { setTimeLimit(); structure(list(msg = conditionMessage(e)), class = "bmerr") })
  list(res = res, sec = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
score <- function(props, true) {
  if (is.null(props) || inherits(props, "bmerr"))
    return(list(mp = NA_real_, rmse = NA_real_, acc = NULL))
  props <- props[rownames(true), , drop = FALSE]
  acc <- per_type_accuracy(props, true)
  list(mp = mean(ifelse(is.na(acc$pearson), 0, acc$pearson)),
       rmse = mean(acc$rmse, na.rm = TRUE), acc = acc)
}

methods <- list()
if (isTRUE(bm$methods$nnls))
  methods[["NNLS"]] <- function(bulk) deconvolve_nnls(signatures, bulk, markers)
if (isTRUE(bm$methods$music))
  methods[["MuSiC"]] <- function(bulk) deconvolve_music(sce_ref, bulk)
if (isTRUE(bm$methods$instaprism))
  methods[["InstaPrism"]] <- function(bulk) deconvolve_instaprism(sce_ref, bulk, n_cores = cfg$threads)
if (isTRUE(bm$methods$bayesprism_confirm))
  methods[["BayesPrism"]] <- function(bulk) deconvolve_bayesprism(sce_ref, t(bulk), n_cores = cfg$threads)

conditions <- list(pseudobulk = easy, `independent-draw` = hard)
rows <- list(); detail <- list()
for (mname in names(methods)) {
  for (cname in names(conditions)) {
    cond <- conditions[[cname]]
    r <- run_method(function() methods[[mname]](cond$bulk))
    if (inherits(r$res, "bmerr") || is.null(r$res)) {
      status <- if (is.null(r$res)) "unavailable (not installed)"
                else if (grepl("elapsed time", r$res$msg)) sprintf("did not complete (> %ds)", budget)
                else paste0("error: ", substr(r$res$msg, 1, 60))
      sc <- list(mp = NA, rmse = NA, acc = NULL)
    } else { status <- "ok"; sc <- score(r$res, cond$true) }
    rows[[length(rows) + 1]] <- data.frame(method = mname, bulk = cname,
      mean_pearson = round(sc$mp, 3), mean_rmse = round(sc$rmse, 3),
      runtime_sec = round(r$sec, 1), status = status)
    if (!is.null(sc$acc)) { sc$acc$method <- mname; sc$acc$bulk <- cname; detail[[length(detail) + 1]] <- sc$acc }
    log_msg(sprintf("%-11s | %-16s | meanR=%s rmse=%s (%.0fs) [%s]",
            mname, cname, ifelse(is.na(sc$mp), "NA", sprintf("%.3f", sc$mp)),
            ifelse(is.na(sc$rmse), "NA", sprintf("%.3f", sc$rmse)), r$sec, status))
  }
}
tab <- do.call(rbind, rows)

# ---- 5. comparison table (wide) + degradation ------------------------------
wide <- reshape(tab[, c("method", "bulk", "mean_pearson")], idvar = "method",
                timevar = "bulk", direction = "wide")
names(wide) <- sub("mean_pearson\\.", "", names(wide))
wide$degradation <- round(wide$pseudobulk - wide$`independent-draw`, 3)
cat("\n================ COMPARISON: mean Pearson (per recoverable cell type) ============\n")
print(wide, row.names = FALSE)
cat("\n(per condition runtimes / status)\n"); print(tab, row.names = FALSE)
write.csv(tab, "results/03b_benchmark.csv", row.names = FALSE)
if (length(detail)) write.csv(do.call(rbind, detail), "results/03b_benchmark_per_type.csv", row.names = FALSE)

# ---- 6. figure --------------------------------------------------------------
p <- ggplot(tab[!is.na(tab$mean_pearson), ],
            aes(method, mean_pearson, fill = bulk)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", mean_pearson)),
            position = position_dodge(0.8), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(pseudobulk = "#4C9F70", `independent-draw` = "#C0504D")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Comparative deconvolution: pseudobulk (easy) vs independent draw (hard)",
       subtitle = "Hippen-2023 hypothesis: naive methods degrade more than BayesPrism-class methods",
       x = NULL, y = "mean Pearson vs true fractions", fill = "bulk type") +
  theme_bw()
ggsave("figs/03b_benchmark.png", p, width = 8, height = 5, dpi = 120)
log_msg("figure -> figs/03b_benchmark.png")

# ---- 7. verdict + notes -----------------------------------------------------
naive <- intersect(c("NNLS", "MuSiC"), wide$method)
robust <- intersect(c("InstaPrism", "BayesPrism"), wide$method)
naive_drop  <- mean(wide$degradation[wide$method %in% naive], na.rm = TRUE)
robust_drop <- mean(wide$degradation[wide$method %in% robust], na.rm = TRUE)
verdict <- if (is.finite(naive_drop) && is.finite(robust_drop)) {
  if (naive_drop > robust_drop + 0.02)
    sprintf("EXPECTED pattern seen: naive methods degraded more (mean drop %.3f) than BayesPrism-class (%.3f).", naive_drop, robust_drop)
  else sprintf("Pattern NOT clearly seen at this scale: naive drop %.3f vs BayesPrism-class %.3f (report as-is).", naive_drop, robust_drop)
} else "Not enough methods completed to judge the pattern."
cat("\nVERDICT:", verdict, "\n")

bp_rows <- tab[tab$method == "BayesPrism", ]
bp_line <- if (nrow(bp_rows) == 0)
  "- True BayesPrism: not run this time. Set bayesprism_confirm to attempt it."
else if (all(grepl("did not complete", bp_rows$status)))
  paste0("- True BayesPrism: attempted once for confirmation but did not complete within the ",
         budget, "s per-condition budget at this scale. Every Bayesian number above is therefore ",
         "InstaPrism, not BayesPrism.")
else "- True BayesPrism: completed. See the table above."

notes <- c(
  "# Comparative deconvolution benchmark (03b): notes",
  "",
  "## Hypothesis (Hippen 2023)",
  "Naive reference-based methods (MuSiC, with NNLS as a simple baseline) score well on pseudobulk,",
  "where the bulk is built from the same cells as the reference. They are expected to degrade when the",
  "bulk is an independent draw: different patients or platform, a different composition, and cell types",
  "present in the bulk but missing from the reference. BayesPrism-family methods should hold up better.",
  "",
  "## What was observed",
  paste0("Profile and scale: ", cfg$.profile, " (", cfg$sizes$n_genes, " genes, ", N,
         " cells in 2 batches, ", bm$n_bulk, " bulk samples per condition)."),
  "Mean Pearson per recoverable cell type, pseudobulk then independent-draw, with the change:",
  paste(capture.output(print(wide, row.names = FALSE)), collapse = "\n"),
  "",
  paste("VERDICT:", verdict),
  "",
  "## Platform-mismatch setting (stated for transparency)",
  paste0("The independent-draw bulk is splatter Batch 2 with a per-gene batch effect set to ",
         "independent_batch_facLoc = ", bm$independent_batch_facLoc,
         " and independent_batch_facScale = ", bm$independent_batch_facScale,
         " (the benchmark block in config.yml)."),
  "This is a substantial shift, chosen to represent a realistic scRNA-versus-bulk platform gap. The",
  "separation between methods depends on this choice: with a weak mismatch the naive methods do not",
  "degrade and all methods look similar. The result is therefore conditional on a substantial",
  "mismatch, recorded here rather than hidden.",
  "",
  "## Engines and scale",
  paste0("- Reference downsampled to ", ncpt, " cells per type; ", length(markers), " marker genes."),
  "- InstaPrism is the engine behind every Bayesian number here. It is a documented fast approximation",
  "  of BayesPrism.",
  bp_line,
  paste0("- Per-method time budget: ", budget, "s. Anything over is logged as 'did not complete'."),
  "Per-condition status and runtime:",
  paste(capture.output(print(tab, row.names = FALSE)), collapse = "\n"),
  "",
  "## Design choice: batches, not a fresh different-seed simulation",
  "The independent target is splatter Batch 2. It is a separate sample with its own per-gene technical",
  "profile (the batch effect, standing in for a different patient or platform), plus a different",
  "composition and extra reference dropouts. A fresh simulation with a different seed would",
  "re-randomise which gene means which cell type. The reference would then be meaningless and every",
  "method would score near zero, which is not a useful comparison. Batches keep gene identity fixed.",
  "",
  "## Caveat",
  "This validates machinery on simulated data, not biology. It reproduces the shape of the Hippen 2023",
  "robustness contrast under controlled conditions. It is not evidence about real tissue."
)
writeLines(notes, "truth/benchmark_notes.md")
log_msg("notes -> truth/benchmark_notes.md")
log_msg("DONE. This is a hypothesis test, not a build gate (no PASS/FAIL emitted).")
