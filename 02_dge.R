# ============================================================================
# 02_dge.R  —  AIM 2: chemoresistance / chemo-stress DE programme
# ----------------------------------------------------------------------------
# WHAT:  Paired pseudobulk differential expression between post-NACT and
#        treatment-naive malignant (EOC) cells, one profile per patient x
#        phase (NOT per-cell -- see note below), CPM + log2, paired t-test
#        across patients + BH FDR.
# WHY:   Aim 2 — recover a treatment-associated transcriptional programme in
#        malignant cells. On real data (GSE165897) the axis with a direct,
#        literal chemoresistance framing is treatment_phase (treatment-naive
#        vs post-neoadjuvant-chemotherapy), not anatomic site.
# STATS NOTE: testing per-cell (as an earlier synthetic-data version of this
#        script did) would be pseudoreplication -- cells from the same patient
#        are correlated, so treating ~8,800 malignant cells as independent
#        replicates massively inflates false discoveries (Squair et al. 2021,
#        "Confronting false discoveries in single-cell differential
#        expression"). This version aggregates to one raw-count profile per
#        patient x treatment_phase first, then tests across patients (n=11,
#        or fewer if any group is dropped for too few cells) -- underpowered
#        by design, which is an honest property of n=11 real patients, not a
#        bug to compensate for by reverting to per-cell testing.
# BIOLOGICAL CAVEAT: Zhang et al. 2022 (the source study) frame the post-NACT
#        vs treatment-naive contrast as capturing BOTH subclonal selection
#        under chemotherapy pressure AND an acute, microenvironment-induced
#        transcriptional stress response -- NOT acquired resistance alone.
#        See results/02_README.md.
# INPUT: data/sce.rds
# OUTPUT: results/02_de_results.csv, results/02_plausibility_check.csv,
#        results/02_README.md, figs/02_volcano.png
# HOW TO READ THE OUTPUT: there is no planted ground truth on real data, so
#        there is no PASS/FAIL gate here (unlike the synthetic-data version).
#        results/02_plausibility_check.csv reports (does not grade) overlap
#        with a small set of established HGSOC platinum-resistance genes.
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment)
  library(matrixStats); library(ggplot2)
})
source("R/utils.R"); source("R/adapters.R"); source("R/metrics.R")
cfg <- load_config()
init_script("02_dge.R — Aim 2 (chemoresistance/stress-response DE)", cfg)

sce <- readRDS("data/sce.rds")
sce_mal <- sce[, sce$celltype == "malignant"]

# ---- 1. pseudobulk aggregation: one profile per patient x treatment_phase ---
MIN_CELLS <- 20L
grp_key <- paste(sce_mal$patient_id, sce_mal$treatment_phase, sep = "||")
counts  <- get_counts(sce_mal)

n_per_grp <- table(grp_key)
keep_grp  <- names(n_per_grp)[n_per_grp >= MIN_CELLS]
dropped   <- setdiff(names(n_per_grp), keep_grp)
if (length(dropped) > 0)
  log_msg("dropping", length(dropped), "patient x phase group(s) with <", MIN_CELLS,
          "malignant cells:", paste(dropped, collapse = "; "))

pb <- sapply(keep_grp, function(g) rowSums(counts[, grp_key == g, drop = FALSE]))
rownames(pb) <- rownames(sce_mal)
pb_meta <- as.data.frame(do.call(rbind, strsplit(colnames(pb), "\\|\\|")),
                         stringsAsFactors = FALSE)
colnames(pb_meta) <- c("patient_id", "treatment_phase")

# Keep only patients with BOTH phases present after the min-cell filter --
# required for the paired design.
paired_patients <- intersect(pb_meta$patient_id[pb_meta$treatment_phase == "post-NACT"],
                              pb_meta$patient_id[pb_meta$treatment_phase == "treatment-naive"])
log_msg(length(paired_patients), "of", length(unique(sce_mal$patient_id)),
        "patients have both phases with >=", MIN_CELLS, "malignant cells (paired n)")
stopifnot("Fewer than 3 paired patients -- not enough for a paired t-test" =
            length(paired_patients) >= 3)

post_idx  <- which(pb_meta$treatment_phase == "post-NACT"       & pb_meta$patient_id %in% paired_patients)
naive_idx <- which(pb_meta$treatment_phase == "treatment-naive" & pb_meta$patient_id %in% paired_patients)
post  <- pb[, post_idx[order(pb_meta$patient_id[post_idx])],   drop = FALSE]
naive <- pb[, naive_idx[order(pb_meta$patient_id[naive_idx])], drop = FALSE]
stopifnot(identical(sort(pb_meta$patient_id[post_idx]), sort(pb_meta$patient_id[naive_idx])))
n <- ncol(post)

# ---- 2. CPM normalise + log2, then paired t-test (vectorised) ---------------
cpm <- function(m) sweep(m, 2, pmax(colSums(m), 1), "/") * 1e6
log2cpm_post  <- log2(cpm(post)  + 1)
log2cpm_naive <- log2(cpm(naive) + 1)

d    <- log2cpm_post - log2cpm_naive           # per-patient (post - naive), paired
m_d  <- rowMeans(d)
se   <- rowSds(d) / sqrt(n)
tstat <- m_d / se
df <- n - 1
pval <- 2 * pt(-abs(tstat), df)
valid <- is.finite(pval) & se > 0
pval[!valid] <- 1
fdr <- p.adjust(pval, "BH")

de <- data.frame(gene = rownames(sce_mal), log2fc = m_d,
                 pval = pval, fdr = fdr,
                 direction = ifelse(m_d > 0, "up_in_postNACT", "down_in_postNACT"))
de$significant <- de$fdr < 0.05
write.csv(de[order(de$fdr), ], "results/02_de_results.csv", row.names = FALSE)
save_peek(de[order(de$fdr), ], "results/peek_02_top_de.csv")
log_msg(sprintf("paired pseudobulk DE: n=%d patient pairs, %d genes significant at FDR<0.05",
                n, sum(de$significant)))

# ============================  Aim 2 output (no planted truth)  =============
# thr is kept only as a visual reference line on the volcano plot below (the
# same effect-size convention used in the synthetic pipeline), NOT as a
# recovery/PASS-FAIL gate -- there is no real ground truth to gate against.
thr <- cfg$truth$resistance_log2fc_threshold

lit_genes <- c("ABCB1", "ABCC2", "ERCC1", "BRCA1", "BRCA2", "TP53",
               "NFE2L2", "GSTP1", "AKR1C1", "AKR1C2")
de_ranked <- de[order(de$fdr), ]
lit_hit <- de_ranked[de_ranked$gene %in% lit_genes, c("gene", "log2fc", "pval", "fdr", "significant")]
lit_hit$rank_overall <- match(lit_hit$gene, de_ranked$gene)
write.csv(lit_hit, "results/02_plausibility_check.csv", row.names = FALSE)
log_msg(sprintf("literature-plausibility check: %d/%d genes found in reference; %d significant at FDR<0.05",
                nrow(lit_hit), length(lit_genes), sum(lit_hit$significant, na.rm = TRUE)))

writeLines(c(
  "Aim 2 plausibility check -- read before citing this gene set",
  "==============================================================",
  "",
  sprintf("results/02_de_results.csv reports paired pseudobulk DE between post-NACT and"),
  sprintf("treatment-naive malignant (EOC) cells (n=%d patient pairs%s). There is no", n,
          if (length(dropped) > 0) ", after dropping thin patient x phase groups -- see the run log" else ""),
  "planted ground truth for real data, so results/02_plausibility_check.csv reports",
  "(does not grade/gate) overlap with a small set of established HGSOC platinum-",
  "resistance genes.",
  "",
  "CAVEAT: Zhang et al. 2022 (the source study, Science Advances) frame the post-NACT",
  "vs treatment-naive contrast as capturing BOTH subclonal selection under",
  "chemotherapy pressure AND an acute, microenvironment-induced transcriptional",
  "stress response -- these two processes are NOT separated by this DE axis. Do not",
  "present the recovered gene set as a clean, isolated acquired-resistance",
  "programme without this caveat."
), "results/02_README.md")

# volcano: literature resistance genes highlighted (reported, not gated)
de$is_lit <- de$gene %in% lit_genes
p <- ggplot(de, aes(log2fc, -log10(pmax(pval, 1e-300)))) +
  geom_point(aes(colour = is_lit, size = is_lit), alpha = 0.6) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                      labels = c("other gene", "literature resistance gene")) +
  scale_size_manual(values = c(`FALSE` = 0.6, `TRUE` = 1.6), guide = "none") +
  geom_vline(xintercept = c(-thr, thr), linetype = "dashed") +
  labs(title = "Aim 2: post-NACT vs treatment-naive malignant DE (paired pseudobulk)",
       subtitle = sprintf("n=%d patient pairs; red = literature chemoresistance genes (reported, not gated)", n),
       x = "log2 fold change (post-NACT / treatment-naive)", y = "-log10 p", colour = NULL) +
  theme_bw()
ggsave("figs/02_volcano.png", p, width = 7, height = 5, dpi = 120)
log_msg("figure -> figs/02_volcano.png")

log_msg("Aim 2 complete -- no PASS/FAIL gate on real data (see results/02_README.md)")
