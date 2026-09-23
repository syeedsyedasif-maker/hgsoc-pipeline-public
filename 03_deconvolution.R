# ============================================================================
# 03_deconvolution.R  —  AIM 3: deconvolve real TCGA-OV bulk + link to survival
# ----------------------------------------------------------------------------
# WHAT:  Load real TCGA-OV bulk RNA-seq (tcga_bulk_data/gene_count_matrix.tsv),
#        harmonize its gene IDs to the GSE165897 reference's symbol space,
#        deconvolve into cell-type fractions (BayesPrism, or the NNLS
#        fallback) against the real GSE165897-derived reference, and test the
#        deconvolved fraction<->survival link using REAL patient survival
#        (vital_status/days_to_death/days_to_last_follow_up, fetched by
#        tcga_bulk_data/fetch_tcga.R).
# WHY:   Aim 3 — does malignant/CAF composition estimated from real bulk carry
#        real prognostic signal?
# INPUT: data/sce.rds, tcga_bulk_data/gene_count_matrix.tsv,
#        tcga_bulk_data/sample_metadata.tsv (must include vital_status/
#        os_time_days/os_status -- see fetch_tcga.R)
# OUTPUT: results/03_estimated_fractions.csv, results/03_metrics.csv,
#        figs/03_km_by_malignant_tertile.png
# HOW TO READ THE OUTPUT: there is no ground-truth composition for real bulk
#        (unlike the synthetic-data version of this script), so there is no
#        fraction-accuracy PASS/FAIL and no truth-vs-estimate scatter figure
#        here. The adipocyte blind spot is reported as information, not
#        graded (GSE165897 has zero adipocytes, so it structurally can't be
#        detected -- see 00_load_real_data.R). The Cox/KM survival test on
#        REAL patient outcomes is this stage's headline real result.
# CONFOUNDS: batch_shift / mrna_enrichment stay available (still gated by
#        config.yml, still applied to `bulk` if enabled) but are OFF in the
#        `real` config profile -- GSE165897 (scRNA) vs TCGA-OV (bulk) is
#        already a genuine, uncontrolled reference<->bulk platform gap;
#        stacking an uncalibrated synthetic multiplier on top of a real one
#        adds noise, not realism.
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(survival)
  library(data.table)
})
source("R/utils.R"); source("R/adapters.R"); source("R/metrics.R")
source("R/deconvolve.R")
cfg <- load_config()
init_script("03_deconvolution.R — Aim 3 (deconvolution + survival)", cfg)

sce <- readRDS("data/sce.rds")
celltype <- as.character(sce$celltype)
types_all <- sort(unique(celltype))

# ---- 1. load real TCGA-OV bulk + harmonize gene IDs to the sc reference's ---
# symbol space (TCGA is versioned Ensembl; GSE165897-derived sce is HGNC
# symbols -- see R/adapters.R::harmonize_bulk_to_symbols()).
TCGA_DIR <- Sys.getenv("TCGA_DIR", "../tcga_bulk_data")
bulk_raw <- data.table::fread(file.path(TCGA_DIR, "gene_count_matrix.tsv"))
meta_tcga <- read.delim(file.path(TCGA_DIR, "sample_metadata.tsv"), stringsAsFactors = FALSE)

bulk_mat <- as.matrix(bulk_raw[, -1, with = FALSE])
rownames(bulk_mat) <- bulk_raw[[1]]
rm(bulk_raw); gc()
bulk_mat <- harmonize_bulk_to_symbols(bulk_mat)

# De-duplicate patients: some TCGA-OV patients have both a -01 (primary) and
# -02 (recurrent) sample. Keep -01 only -- survival time/status is attached
# at the CASE level, so keeping both would silently duplicate one observation
# in the Cox model (pseudo-replication).
sample_type <- sub("^TCGA-[^-]+-[^-]+-([0-9]+).*$", "\\1", colnames(bulk_mat))
n_before_dedup <- ncol(bulk_mat)
bulk_mat <- bulk_mat[, sample_type == "01", drop = FALSE]
log_msg(sprintf("kept %d of %d TCGA-OV samples (sample type -01/primary tumor only)",
                ncol(bulk_mat), n_before_dedup))

meta_tcga <- meta_tcga[match(colnames(bulk_mat), meta_tcga$sample_id), ]
# fetch_tcga.R already writes only survival-complete rows, but re-check here
# defensively in case sample_metadata.tsv was regenerated some other way.
has_surv <- !is.na(meta_tcga$sample_id) & !is.na(meta_tcga$os_time_days) &
            meta_tcga$os_time_days > 0
bulk_mat  <- bulk_mat[, has_surv, drop = FALSE]
meta_tcga <- meta_tcga[has_surv, ]

P <- ncol(bulk_mat)
bulk <- bulk_mat
log_msg(sprintf("loaded %d TCGA-OV bulk samples with usable survival data", P))

# ---- 2. CONFOUNDS on the bulk only (reference stays clean) ------------------
# Off by default in the `real` profile (config.yml) -- see header note.
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

# ---- 3. reference signatures (adipocyte structurally absent, not toggled) ---
# GSE165897 has zero adipocytes at all (see 00_load_real_data.R), so
# adipocyte_blindspot's toggle is vacuous here -- ref_types == types_all
# regardless, and "adipocyte" simply never appears in types_all.
blind <- isTRUE(cfg$confounds$adipocyte_blindspot$enabled)
ref_types <- if (blind) setdiff(types_all, "adipocyte") else types_all
sce_ref <- sce[, celltype %in% ref_types]
signatures <- build_reference_signatures(sce_ref, type_col = "celltype")
markers <- select_markers(signatures, cfg$deconvolution$markers_per_type)
log_msg("reference cell types:", paste(colnames(signatures), collapse = ", "))
log_msg("selected", length(markers), "marker genes for deconvolution")

# ---- 4. deconvolve ----------------------------------------------------------
dec <- run_deconvolution(signatures, bulk, sce_ref, t(bulk), cfg, markers = markers)
est <- dec$props
log_msg("deconvolution method:", dec$method)
write.csv(data.frame(sample_id = rownames(est), est),
          "results/03_estimated_fractions.csv", row.names = FALSE)

# ---- 5. real survival (from tcga_bulk_data/sample_metadata.tsv) ------------
time   <- meta_tcga$os_time_days
status <- meta_tcga$os_status
log_msg(sprintf("events: %d / %d (%.0f%% censored)", sum(status), P,
                100 * mean(status == 0)))

# ============================  Aim 3 output  =================================
# (a) No ground-truth composition exists for real TCGA bulk, so there is no
# per-cell-type accuracy check here (unlike the synthetic-data version).

# (b) adipocyte blind spot -- informational only, nothing to grade against
# since GSE165897 has no adipocyte reference cells at all.
adipo_detected <- "adipocyte" %in% colnames(est)
log_msg(sprintf("ADIPOCYTE BLIND SPOT: reference has zero adipocyte cells (GSE165897); detected = %s",
                if (adipo_detected) "YES (unexpected)" else "NO (expected -- no reference to detect from)"))

# (c) fraction <-> survival: Cox on ESTIMATED fractions + KM by tertile
est_mal <- est[, "malignant"]; est_caf <- est[, "fibroblast"]
cox <- coxph(Surv(time, status) ~ est_mal + est_caf)
cox_s <- summary(cox)
hr_mal <- cox_s$coefficients["est_mal", "exp(coef)"]
p_mal  <- cox_s$coefficients["est_mal", "Pr(>|z|)"]
log_msg(sprintf("Cox (estimated fractions, real survival): malignant HR=%.2f p=%.3g | CAF HR=%.2f",
                hr_mal, p_mal, cox_s$coefficients["est_caf", "exp(coef)"]))

tert <- cut(est_mal, stats::quantile(est_mal, c(0, 1/3, 2/3, 1)),
            labels = c("low", "mid", "high"), include.lowest = TRUE)
sd_km <- survdiff(Surv(time, status) ~ tert)
p_lr <- 1 - pchisq(sd_km$chisq, length(sd_km$n) - 1)
log_msg(sprintf("KM by estimated-malignant tertile: log-rank p = %.3g", p_lr))

png("figs/03_km_by_malignant_tertile.png", width = 700, height = 500, res = 110)
plot(survfit(Surv(time, status) ~ tert), col = c("forestgreen", "orange", "firebrick"),
     lwd = 2, xlab = "time (days)", ylab = "survival",
     main = sprintf("Aim 3: KM by estimated malignant tertile (log-rank p=%.3g)", p_lr))
legend("bottomleft", legend = levels(tert),
       col = c("forestgreen", "orange", "firebrick"), lwd = 2, bty = "n")
dev.off()
log_msg("figure -> figs/03_km_by_malignant_tertile.png")

write.csv(data.frame(
  metric = c("n_patients", "adipocyte_detected", "cox_malignant_HR",
             "cox_malignant_p", "km_logrank_p"),
  value  = c(P, as.numeric(adipo_detected), hr_mal, p_mal, p_lr)),
  "results/03_metrics.csv", row.names = FALSE)

# Pipeline-integrity check, not a biological-direction gate: the earlier
# hr_mal > 1 gate was calibrated to the deleted PLANTED synthetic hazard and
# has no basis on real data. Confirm the Cox model converged to finite
# estimates; the actual HR/CI/p and log-rank p are logged above in full for
# write-up regardless of which direction they point.
pass_fail("03_deconvolution", "Cox model converged on real TCGA-OV survival data",
          is.finite(hr_mal) && is.finite(p_mal),
          sprintf("malignant HR=%.2f p=%.3g | KM log-rank p=%.3g", hr_mal, p_mal, p_lr))
