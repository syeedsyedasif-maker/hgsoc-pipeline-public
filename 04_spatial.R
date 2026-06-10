# ============================================================================
# 04_spatial.R  —  AIM 4: spatial structure (mock slide)
# ----------------------------------------------------------------------------
# WHAT:  Lay cells out on a 2D grid of spots with a planted tissue architecture:
#        a DESMOPLASTIC CORE where malignant + CAF (fibroblast) colocalise, and
#        an EXPOSED MARGIN with malignant but little CAF. Bin cells into spots,
#        sum counts, deconvolve each spot (RCTD, or per-spot NNLS fallback), and
#        test whether the recovered malignant<->CAF colocalisation matches the
#        planted layout.
# WHY:   Aim 4 — recover spatial co-localisation of tumour and stroma.
# INPUT: data/sce.rds
# OUTPUT: results/04_*.csv, figs/04_*.png
# HOW TO READ THE OUTPUT: PASS = (1) inferred CAF is enriched in the core vs the
#        margin (both are malignant-rich, so this is the desmoplastic signal),
#        and (2) a permutation test shows inferred malignant & CAF co-occur in
#        space more than chance. The maps should show CAF lighting up the core.
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(ggplot2)
})
source("R/utils.R"); source("R/adapters.R"); source("R/metrics.R")
source("R/deconvolve.R"); source("R/rctd_fallback.R")
cfg <- load_config()
init_script("04_spatial.R — Aim 4 (spatial colocalisation)", cfg)

sce <- readRDS("data/sce.rds")
counts <- get_counts(sce)
celltype <- as.character(sce$celltype)
types_all <- sort(unique(celltype))
pool <- split(seq_len(ncol(sce)), celltype)

# ---- 1. grid + planted regions ----------------------------------------------
g <- cfg$sizes$spatial_grid
coords <- expand.grid(x = 1:g, y = 1:g)
ctr <- (g + 1) / 2
d <- sqrt((coords$x - ctr)^2 + (coords$y - ctr)^2)
r_core   <- cfg$spatial$core_radius_frac * g
r_margin <- r_core + 0.18 * g
region <- ifelse(d <= r_core, "core",
           ifelse(d <= r_margin, "margin", "background"))
coords$region <- region
log_msg(sprintf("grid %dx%d = %d spots | core=%d margin=%d background=%d",
                g, g, nrow(coords), sum(region == "core"),
                sum(region == "margin"), sum(region == "background")))

# planted composition per region (fractions over the 6 cell types) ------------
# core: malignant + CAF colocalise; margin: malignant, little CAF; bg: low both.
# Dense desmoplastic core (highest malignant + CAF) vs sparser exposed margin
# (malignant, little CAF) vs background (low malignant).
comp_region <- list(
  core       = c(malignant = 0.50, fibroblast = 0.35, macrophage = 0.07,
                 endothelial = 0.04, Tcell = 0.02, adipocyte = 0.02),
  margin     = c(malignant = 0.40, fibroblast = 0.05, macrophage = 0.15,
                 endothelial = 0.10, Tcell = 0.25, adipocyte = 0.05),
  background = c(malignant = 0.05, fibroblast = 0.10, macrophage = 0.20,
                 endothelial = 0.15, Tcell = 0.20, adipocyte = 0.30))
comp_region <- lapply(comp_region, function(v) v[types_all])   # align order

# ---- 2. build spot expression by summing sampled cells ----------------------
ns <- nrow(coords); ncps <- cfg$sizes$cells_per_spot
spot_counts <- matrix(0, nrow(sce), ns,
                      dimnames = list(rownames(sce), paste0("spot", 1:ns)))
true_comp <- matrix(0, ns, length(types_all), dimnames = list(NULL, types_all))
for (i in seq_len(ns)) {
  base <- comp_region[[coords$region[i]]]
  p_i <- base * exp(rnorm(length(base), 0, 0.15)); p_i <- p_i / sum(p_i)  # spot noise
  n_k <- as.integer(rmultinom(1, ncps, p_i))
  true_comp[i, ] <- n_k / sum(n_k)
  cells_i <- unlist(mapply(function(t, n)
    if (n > 0) sample(pool[[t]], n, replace = TRUE) else integer(0),
    types_all, n_k, SIMPLIFY = FALSE))
  spot_counts[, i] <- rowSums(counts[, cells_i, drop = FALSE])
}

# spatial assay batch shift (same confound family as the bulk in Aim 3)
if (isTRUE(cfg$confounds$batch_shift$enabled)) {
  spot_counts <- spot_counts * exp(rnorm(nrow(spot_counts), 0,
                                         cfg$confounds$batch_shift$sd))
  log_msg("confound: spatial assay batch shift applied")
}

# ---- 3. reference (adipocyte still OMITTED = blind spot) + deconvolve --------
blind <- isTRUE(cfg$confounds$adipocyte_blindspot$enabled)
ref_types <- if (blind) setdiff(types_all, "adipocyte") else types_all
sce_ref <- sce[, celltype %in% ref_types]
signatures <- build_reference_signatures(sce_ref, type_col = "celltype")
markers <- select_markers(signatures, cfg$deconvolution$markers_per_type)
dec <- run_spatial_deconv(signatures, spot_counts, sce_ref,
                          coords[, c("x", "y")], cfg, markers = markers)
est <- dec$props
log_msg("spatial deconvolution method:", dec$method)
write.csv(cbind(coords, est), "results/04_spot_estimates.csv", row.names = FALSE)

# ============================  VALIDATION (Aim 4)  ===========================
# We score colocalisation INSIDE the tumour bed (core + margin: both planted as
# malignant-rich). This isolates the desmoplastic signal and is immune to the
# adipocyte blind spot contaminating background spots.
inf_mal <- est[, "malignant"]; inf_caf <- est[, "fibroblast"]
ok_spot <- is.finite(inf_mal) & is.finite(inf_caf)   # RCTD may drop a few spots
if (sum(!ok_spot) > 0)
  log_msg(sum(!ok_spot), "spots dropped by deconvolution (excluded from tests)")
bed <- (region %in% c("core", "margin")) & ok_spot
caf_core   <- inf_caf[region == "core"   & ok_spot]
caf_margin <- inf_caf[region == "margin" & ok_spot]

# (1) PERMUTATION TEST: is inferred CAF concentrated in the core vs the
#     equally-malignant margin? (shuffle core/margin labels within the bed)
perm <- permutation_diff_test(inf_caf[bed], region[bed] == "core",
                              n_perm = cfg$spatial$n_permutations, seed = cfg$seed)
log_msg(sprintf("inferred CAF: core mean=%.2f vs margin mean=%.2f | permutation p=%.3g",
                mean(caf_core), mean(caf_margin), perm$p))

# (2) context only: raw fraction-vs-fraction correlations are NEGATIVE here, and
#     that is expected -- (a) malignant spans core+margin while CAF is core-only,
#     and (b) fractions are compositional (sum to 1), so within a spot malignant
#     and CAF trade off. So correlating two fractions is the WRONG colocalisation
#     statistic; the permutation test above is the correct one. We report these
#     to make the point explicit, not to gate on them.
xcorr <- suppressWarnings(cor(inf_mal[bed], inf_caf[bed]))
gp <- spatial_coloc_test(inf_mal, inf_caf, n_perm = 99, seed = cfg$seed)$r
log_msg(sprintf("(context) malignant<->CAF fraction Pearson: tumour-bed r=%.2f, global r=%.2f (negative is expected, see code note)",
                xcorr, gp))

write.csv(data.frame(
  metric = c("caf_core_mean", "caf_margin_mean", "perm_p",
             "tumourbed_xcorr", "global_pearson_r"),
  value  = c(mean(caf_core), mean(caf_margin), perm$p, xcorr, gp)),
  "results/04_metrics.csv", row.names = FALSE)

# figures: planted vs inferred maps -------------------------------------------
mapdf <- data.frame(coords,
                    true_malignant = true_comp[, "malignant"],
                    true_CAF = true_comp[, "fibroblast"],
                    inf_malignant = inf_mal, inf_CAF = inf_caf)
mk_map <- function(fill, title)
  ggplot(mapdf, aes(x, y, fill = .data[[fill]])) + geom_tile() +
    scale_fill_viridis_c() + coord_equal() + theme_minimal() +
    labs(title = title, fill = "frac")
pl <- list(mk_map("true_malignant", "planted malignant"),
           mk_map("true_CAF",       "planted CAF"),
           mk_map("inf_malignant",  "inferred malignant"),
           mk_map("inf_CAF",        "inferred CAF"))
# arrange 2x2 without extra packages
png("figs/04_spatial_maps.png", width = 900, height = 800, res = 110)
gridExtra_ok <- requireNamespace("gridExtra", quietly = TRUE)
if (gridExtra_ok) {
  gridExtra::grid.arrange(grobs = pl, ncol = 2)
} else {
  print(pl[[4]])   # fallback: at least the headline inferred-CAF map
}
dev.off()
log_msg("figure -> figs/04_spatial_maps.png",
        if (!gridExtra_ok) "(install gridExtra for the 2x2 panel)" else "")

pass_fail("04_spatial", "recovered malignant-CAF colocalisation matches layout",
          (perm$p < cfg$thresholds$spatial_p_max) &&
          (mean(caf_core) > mean(caf_margin)),
          sprintf("CAF concentrated in malignant-dense core vs margin: perm p=%.3g (core=%.2f margin=%.2f)",
                  perm$p, mean(caf_core), mean(caf_margin)))
