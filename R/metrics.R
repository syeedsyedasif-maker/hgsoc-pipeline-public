# ============================================================================
# R/metrics.R — scoring functions used to compare recovered signal vs truth.
# All base R (no extra packages) so validation never fails to load.
# Each function is small and self-explanatory so you can audit the maths.
# ============================================================================

# Adjusted Rand Index: agreement between two clusterings (1 = identical,
# ~0 = random). Used to score recovered cell groups vs true cell types.
adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  comb2 <- function(x) sum(x * (x - 1) / 2)
  sum_c <- comb2(rowSums(tab))
  sum_k <- comb2(colSums(tab))
  sum_t <- comb2(as.vector(tab))
  expected <- sum_c * sum_k / (n * (n - 1) / 2)
  max_idx  <- (sum_c + sum_k) / 2
  if (max_idx == expected) return(1)            # degenerate (all identical)
  (sum_t - expected) / (max_idx - expected)
}

# Confusion matrix + accuracy for a binary/multiclass call vs truth.
confusion <- function(truth, predicted) {
  cm <- table(truth = truth, predicted = predicted)
  acc <- sum(diag(cm)) / sum(cm)
  list(matrix = cm, accuracy = acc)
}

# Precision / recall / F1 for SET recovery (e.g. recovered DE genes vs planted).
precision_recall_f1 <- function(predicted, truth) {
  predicted <- unique(predicted); truth <- unique(truth)
  tp <- length(intersect(predicted, truth))
  fp <- length(setdiff(predicted, truth))
  fn <- length(setdiff(truth, predicted))
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  recall    <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  list(tp = tp, fp = fp, fn = fn,
       precision = precision, recall = recall, f1 = f1)
}

rmse    <- function(est, truth) sqrt(mean((est - truth)^2))
pearson <- function(est, truth) suppressWarnings(stats::cor(est, truth))

# Per-cell-type deconvolution accuracy: Pearson + RMSE for each type across
# patients/spots. `est` and `truth` are matrices [samples x celltypes].
per_type_accuracy <- function(est, truth) {
  types <- intersect(colnames(est), colnames(truth))
  out <- data.frame(celltype = types,
                    pearson = NA_real_, rmse = NA_real_,
                    mean_true = NA_real_, mean_est = NA_real_)
  for (i in seq_along(types)) {
    t <- types[i]
    out$pearson[i]   <- pearson(est[, t], truth[, t])
    out$rmse[i]      <- rmse(est[, t], truth[, t])
    out$mean_true[i] <- mean(truth[, t])
    out$mean_est[i]  <- mean(est[, t])
  }
  out
}

# Global Pearson colocalisation (context only): note this is the WRONG statistic
# when colocalisation is regional (e.g. CAF only in the core while malignant
# spans core + margin) -- it can go negative. Reported, not gated.
spatial_coloc_test <- function(x, y, n_perm = 999, seed = 1) {
  set.seed(seed)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  obs <- suppressWarnings(stats::cor(x, y))
  if (is.na(obs)) return(list(r = NA, p = NA, n_perm = n_perm))
  null <- replicate(n_perm, suppressWarnings(stats::cor(x, sample(y))))
  p <- (1 + sum(null >= obs)) / (1 + n_perm)
  list(r = obs, p = p, n_perm = n_perm)
}

# Permutation test for a difference in means between two spatial groups.
# Used for the desmoplastic signal: among malignant-rich (tumour-bed) spots, is
# inferred CAF higher in the core than the margin? `group` is TRUE for the core.
# Null shuffles the core/margin labels, so it tests spatial concentration without
# assuming any distribution.
permutation_diff_test <- function(x, group, n_perm = 999, seed = 1) {
  set.seed(seed)
  g <- as.logical(group)
  obs <- mean(x[g]) - mean(x[!g])
  null <- replicate(n_perm, { gp <- sample(g); mean(x[gp]) - mean(x[!gp]) })
  p <- (1 + sum(null >= obs)) / (1 + n_perm)
  list(diff = obs, p = p, n_perm = n_perm)
}
