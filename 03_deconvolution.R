# ============================================================================
# 03_deconvolution.R  —  AIM 3: deconvolve bulk + link to survival
# ----------------------------------------------------------------------------
# WHAT:  Build ~120 patient pseudobulk samples by summing single cells in KNOWN,
#        varying proportions (with the realistic confounds applied), give each
#        patient a survival outcome whose hazard depends on malignant + CAF
#        fraction, then deconvolve the bulk back into cell-type fractions
#        (BayesPrism, or the NNLS fallback) and test the fraction<->survival link.
# WHY:   Aim 3 — can we recover composition from bulk under reference/bulk
#        mismatch, and does composition carry prognostic signal?
# INPUT: data/sce.rds (+ truth/cell_metadata.csv)
# OUTPUT: results/03_*.csv, figs/03_*.png
# HOW TO READ THE OUTPUT:
#   * per-type scatter: estimated vs true fraction should hug the diagonal;
#     we report Pearson + RMSE PER cell type (not just overall).
#   * adipocyte BLIND SPOT: adipocyte is in the bulk but NOT in the reference,
#     so it cannot be detected — we print the missed fraction. This is a
#     RESULT (it mirrors dissociation loss), not a bug.
#   * survival: Cox + KM by tertile; higher malignant/CAF fraction => worse
#     survival, matching the planted hazard. 120 patients is a DEMO, not powered.
# CONFOUNDS applied to the bulk only (see config): reference<->bulk batch shift
#   and mRNA-enrichment, so the reference and bulk are NOT biologically identical.
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(survival)
})
source("R/utils.R"); source("R/adapters.R"); source("R/metrics.R")
source("R/deconvolve.R")
cfg <- load_config()
init_script("03_deconvolution.R — Aim 3 (deconvolution + survival)", cfg)

sce <- readRDS("data/sce.rds")
counts <- get_counts(sce)
celltype <- as.character(sce$celltype)
types_all <- sort(unique(celltype))
P  <- cfg$sizes$n_patients
nc <- cfg$sizes$cells_per_patient

# pools of cell indices per type (to sample from)
pool <- split(seq_len(ncol(sce)), celltype)

# ---- 1. per-patient KNOWN proportions (Dirichlet around real composition) ---
base_comp <- table(celltype)[types_all] / length(celltype)
alpha <- cfg$survival$dirichlet_alpha * as.numeric(base_comp)
draw_dirichlet <- function(a) { g <- rgamma(length(a), shape = a, rate = 1); g / sum(g) }

# We track TWO known truths per patient:
#   true_cellfrac  = fraction of CELLS of each type (the biology; drives survival)
#   true_countfrac = fraction of EXPRESSION (counts) each type contributes to the
#                    bulk -- this is what reference-based deconvolution actually
#                    estimates, so it is the fair target for accuracy scoring.
# They differ because cell types carry different amounts of mRNA (e.g. malignant
# cells with CNV gains express more). We document this distinction in the README.
true_cellfrac  <- matrix(0, P, length(types_all), dimnames = list(NULL, types_all))
true_countfrac <- matrix(0, P, length(types_all), dimnames = list(NULL, types_all))
bulk <- matrix(0, nrow(sce), P, dimnames = list(rownames(sce), paste0("patient", 1:P)))
for (i in seq_len(P)) {
  p_i <- draw_dirichlet(alpha)
  n_k <- as.integer(rmultinom(1, nc, p_i))          # integer cells per type
  true_cellfrac[i, ] <- n_k / sum(n_k)
  b <- numeric(nrow(sce)); contrib <- setNames(numeric(length(types_all)), types_all)
  for (j in seq_along(types_all)) {
    if (n_k[j] == 0) next
    cells_k <- sample(pool[[types_all[j]]], n_k[j], replace = TRUE)
    sub <- rowSums(counts[, cells_k, drop = FALSE])
    b <- b + sub; contrib[j] <- sum(sub)
  }
  bulk[, i] <- b
  true_countfrac[i, ] <- contrib / sum(contrib)
}
log_msg(sprintf("built %d patient pseudobulks (%d cells each)", P, nc))

# ---- 2. CONFOUNDS on the bulk only (reference stays clean) ------------------
if (isTRUE(cfg$confounds$batch_shift$enabled)) {
  bfac <- exp(rnorm(nrow(bulk), 0, cfg$confounds$batch_shift$sd))   # per-gene
  bulk <- bulk * bfac
  log_msg("confound: reference<->bulk batch shift applied (sd =",
          cfg$confounds$batch_shift$sd, ")")
}
if (isTRUE(cfg$confounds$mrna_enrichment$enabled)) {
  k <- round(cfg$confounds$mrna_enrichment$subset_frac * nrow(bulk))
  capt <- sample(nrow(bulk), k)
  bulk[capt, ] <- bulk[capt, ] * cfg$confounds$mrna_enrichment$factor
  log_msg(sprintf("confound: mRNA-enrichment boosted %d genes x%.1f",
                  k, cfg$confounds$mrna_enrichment$factor))
}

# ---- 3. reference signatures (adipocyte OMITTED = blind spot) ---------------
blind <- isTRUE(cfg$confounds$adipocyte_blindspot$enabled)
ref_types <- if (blind) setdiff(types_all, "adipocyte") else types_all
sce_ref <- sce[, celltype %in% ref_types]
signatures <- build_reference_signatures(sce_ref, type_col = "celltype")
markers <- select_markers(signatures, cfg$deconvolution$markers_per_type)
log_msg("reference cell types:", paste(colnames(signatures), collapse = ", "),
        if (blind) "(adipocyte deliberately excluded)" else "")
log_msg("selected", length(markers), "marker genes for deconvolution")

# ---- 4. deconvolve ----------------------------------------------------------
dec <- run_deconvolution(signatures, bulk, sce_ref, t(bulk), cfg, markers = markers)
est <- dec$props
log_msg("deconvolution method:", dec$method)
write.csv(data.frame(patient = rownames(est), est),
          "results/03_estimated_fractions.csv", row.names = FALSE)
write.csv(data.frame(patient = paste0("patient", 1:P), true_countfrac),
          "results/03_true_fractions.csv", row.names = FALSE)

# ---- 5. survival outcome: hazard ~ malignant + CAF fraction -----------------
# Driven by the true CELL fractions (biology); we then test whether the
# DECONVOLVED fractions recover that prognostic signal.
mal_true <- true_cellfrac[, "malignant"]; caf_true <- true_cellfrac[, "fibroblast"]
sv <- cfg$survival
if (isTRUE(sv$use_simsurv) && requireNamespace("simsurv", quietly = TRUE)) {
  log_msg("survival: simsurv (PRIMARY)")
  xdat <- data.frame(id = 1:P, mal = mal_true, caf = caf_true)
  st <- simsurv::simsurv(dist = "exponential", lambdas = sv$baseline_hazard,
          betas = c(mal = sv$beta_malignant, caf = sv$beta_caf),
          x = xdat, maxt = sv$max_followup)
  time <- st$eventtime; status <- st$status
} else {
  log_msg("[FALLBACK] survival: exponential inverse-CDF hazard model")
  rate <- sv$baseline_hazard * exp(sv$beta_malignant * mal_true + sv$beta_caf * caf_true)
  t_evt <- rexp(P, rate)
  time <- pmin(t_evt, sv$max_followup); status <- as.integer(t_evt <= sv$max_followup)
}
log_msg(sprintf("events: %d / %d (%.0f%% censored)", sum(status), P,
                100 * mean(status == 0)))

# ============================  VALIDATION (Aim 3)  ===========================
# (a) per-cell-type accuracy vs the EXPRESSION-fraction truth (incl adipocyte)
acc <- per_type_accuracy(est, true_countfrac)
cat("\n  Per-cell-type deconvolution accuracy (vs expression fractions):\n")
print(acc, row.names = FALSE)
write.csv(acc, "results/03_per_type_accuracy.csv", row.names = FALSE)
# honest mean: a type that is NEVER detected (constant 0 -> NA correlation)
# counts as 0, not silently dropped.
mean_pearson <- mean(ifelse(is.na(acc$pearson), 0, acc$pearson))

# (b) adipocyte blind spot
adipo_true <- mean(true_countfrac[, "adipocyte"])
adipo_detected <- "adipocyte" %in% colnames(est)
log_msg(sprintf("ADIPOCYTE BLIND SPOT: true mean fraction = %.1f%%, detected by reference = %s",
                100 * adipo_true, if (adipo_detected) "YES" else "NO (missed, as expected)"))

# (c) fraction <-> survival: Cox on ESTIMATED fractions + KM by tertile
est_mal <- est[, "malignant"]; est_caf <- est[, "fibroblast"]
cox <- coxph(Surv(time, status) ~ est_mal + est_caf)
cox_s <- summary(cox)
hr_mal <- cox_s$coefficients["est_mal", "exp(coef)"]
p_mal  <- cox_s$coefficients["est_mal", "Pr(>|z|)"]
log_msg(sprintf("Cox (estimated fractions): malignant HR=%.2f p=%.3g | CAF HR=%.2f",
                hr_mal, p_mal, cox_s$coefficients["est_caf", "exp(coef)"]))

tert <- cut(est_mal, stats::quantile(est_mal, c(0, 1/3, 2/3, 1)),
            labels = c("low", "mid", "high"), include.lowest = TRUE)
sd_km <- survdiff(Surv(time, status) ~ tert)
p_lr <- 1 - pchisq(sd_km$chisq, length(sd_km$n) - 1)
log_msg(sprintf("KM by estimated-malignant tertile: log-rank p = %.3g", p_lr))

# figures
png("figs/03_fraction_scatter.png", width = 900, height = 600, res = 110)
op <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
for (t in colnames(est)) {
  plot(true_countfrac[, t], est[, t], pch = 19, cex = 0.5, col = "steelblue",
       xlab = "true", ylab = "estimated",
       main = sprintf("%s (r=%.2f)", t, acc$pearson[acc$celltype == t]))
  abline(0, 1, col = "grey40", lty = 2)
}
par(op); dev.off()
log_msg("figure -> figs/03_fraction_scatter.png")

png("figs/03_km_by_malignant_tertile.png", width = 700, height = 500, res = 110)
plot(survfit(Surv(time, status) ~ tert), col = c("forestgreen", "orange", "firebrick"),
     lwd = 2, xlab = "time (months)", ylab = "survival",
     main = sprintf("Aim 3: KM by estimated malignant tertile (log-rank p=%.3g)", p_lr))
legend("bottomleft", legend = levels(tert),
       col = c("forestgreen", "orange", "firebrick"), lwd = 2, bty = "n")
dev.off()
log_msg("figure -> figs/03_km_by_malignant_tertile.png")

write.csv(data.frame(
  metric = c("mean_pearson", "adipocyte_true_frac", "adipocyte_detected",
             "cox_malignant_HR", "cox_malignant_p", "km_logrank_p"),
  value  = c(mean_pearson, adipo_true, as.numeric(adipo_detected),
             hr_mal, p_mal, p_lr)),
  "results/03_metrics.csv", row.names = FALSE)

log_msg("NOTE: 120 patients is a DEMO of the machinery, not a powered survival study.")

# sub-results (each logged), then the stage headline
pass_fail("03_deconvolution", "adipocyte blind spot demonstrated",
          (adipo_true > 0) && !adipo_detected,
          sprintf("true=%.1f%% detected=%s", 100 * adipo_true, adipo_detected))
pass_fail("03_deconvolution", "fraction-survival association recovered (direction)",
          hr_mal > 1,
          sprintf("malignant HR=%.2f (planted: higher->worse) KM p=%.3g", hr_mal, p_lr))
pass_fail("03_deconvolution", "per-cell-type fractions recovered under confounds",
          mean_pearson >= cfg$thresholds$decon_pearson_min,
          sprintf("mean Pearson=%.3f (>=%.2f)", mean_pearson,
                  cfg$thresholds$decon_pearson_min))
