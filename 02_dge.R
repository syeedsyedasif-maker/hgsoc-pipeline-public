# ============================================================================
# 02_dge.R  —  AIM 2: the resistance DE programme
# ----------------------------------------------------------------------------
# WHAT:  Differential expression between malignant_omental and malignant_ascites
#        cells (a per-gene Welch t-test on log-normalised counts + BH FDR).
# WHY:   Aim 2 — recover the cisplatin-resistance signal we planted as DE
#        up in the omental tumour cells.
# INPUT: data/sce.rds (+ truth/resistance_genes.csv for scoring)
# OUTPUT: results/02_de_results.csv, figs/02_volcano.png
# HOW TO READ THE OUTPUT: PASS = the genes we call DE match the planted
#        resistance set well (F1 >= threshold). The volcano highlights the
#        planted genes; they should sit in the significant tails. Both malignant
#        groups carry the SAME CNV, so copy-number cancels out of this contrast.
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment)
  library(matrixStats); library(ggplot2)
})
source("R/utils.R"); source("R/adapters.R"); source("R/metrics.R")
cfg <- load_config()
init_script("02_dge.R — Aim 2 (resistance DE)", cfg)

sce <- readRDS("data/sce.rds")
sce_mal <- sce[, sce$celltype == "malignant"]
grp <- as.character(sce_mal$group)
om  <- which(grp == "malignant_omental")
asc <- which(grp == "malignant_ascites")
log_msg(sprintf("malignant cells: %d omental vs %d ascites", length(om), length(asc)))

# ---- per-gene Welch t-test (vectorised) -------------------------------------
# logNormCounts is base-2, so (mean_omental - mean_ascites) ~ log2 fold change.
X <- get_logcounts(sce_mal)
m_om  <- rowMeans(X[, om]);  m_asc <- rowMeans(X[, asc])
v_om  <- rowVars(X[, om]);   v_asc <- rowVars(X[, asc])
n_om  <- length(om);         n_asc <- length(asc)
se <- sqrt(v_om / n_om + v_asc / n_asc)
log2fc <- m_om - m_asc
tstat <- log2fc / se
df <- (v_om / n_om + v_asc / n_asc)^2 /
      ((v_om / n_om)^2 / (n_om - 1) + (v_asc / n_asc)^2 / (n_asc - 1))
pval <- 2 * pt(-abs(tstat), df)
valid <- is.finite(pval) & se > 0
pval[!valid] <- 1
fdr <- p.adjust(pval, "BH")

de <- data.frame(gene = rownames(sce_mal), log2fc = log2fc,
                 pval = pval, fdr = fdr,
                 direction = ifelse(log2fc > 0, "up_in_omental", "down_in_omental"))
thr <- cfg$truth$resistance_log2fc_threshold
de$recovered <- de$fdr < 0.05 & abs(de$log2fc) >= thr
write.csv(de[order(de$fdr), ], "results/02_de_results.csv", row.names = FALSE)
save_peek(de[order(de$fdr), ], "results/peek_02_top_de.csv")

# ============================  VALIDATION (Aim 2)  ===========================
truth <- read.csv("truth/resistance_genes.csv", stringsAsFactors = FALSE)
truth_genes <- truth$gene
recovered_genes <- de$gene[de$recovered]
prf <- precision_recall_f1(recovered_genes, truth_genes)
de$is_truth <- de$gene %in% truth_genes

log_msg(sprintf("planted=%d  recovered=%d  TP=%d FP=%d FN=%d",
                length(truth_genes), length(recovered_genes), prf$tp, prf$fp, prf$fn))
log_msg(sprintf("precision=%.3f  recall=%.3f  F1=%.3f",
                prf$precision, prf$recall, prf$f1))
write.csv(data.frame(metric = c("precision", "recall", "f1", "tp", "fp", "fn"),
                     value = c(prf$precision, prf$recall, prf$f1,
                               prf$tp, prf$fp, prf$fn)),
          "results/02_metrics.csv", row.names = FALSE)

# volcano: planted resistance genes highlighted
p <- ggplot(de, aes(log2fc, -log10(pmax(pval, 1e-300)))) +
  geom_point(aes(colour = is_truth, size = is_truth), alpha = 0.6) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                      labels = c("other gene", "planted resistance gene")) +
  scale_size_manual(values = c(`FALSE` = 0.6, `TRUE` = 1.6), guide = "none") +
  geom_vline(xintercept = c(-thr, thr), linetype = "dashed") +
  labs(title = "Aim 2: omental vs ascites malignant DE (volcano)",
       subtitle = "red = genes we planted as resistance programme",
       x = "log2 fold change (omental / ascites)", y = "-log10 p", colour = NULL) +
  theme_bw()
ggsave("figs/02_volcano.png", p, width = 7, height = 5, dpi = 120)
log_msg("figure -> figs/02_volcano.png")

pass_fail("02_dge", "recovered resistance genes match planted set",
          prf$f1 >= cfg$thresholds$f1_min,
          sprintf("F1=%.3f (>=%.2f) P=%.3f R=%.3f",
                  prf$f1, cfg$thresholds$f1_min, prf$precision, prf$recall))
