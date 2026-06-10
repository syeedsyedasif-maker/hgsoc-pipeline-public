# ============================================================================
# 01_qc_cluster.R  —  AIM 1: find the malignant cells
# ----------------------------------------------------------------------------
# WHAT:  Quality-control + normalise + cluster the cells, then use a CNV signal
#        (inferCNV, or the windowed-expression fallback) with the stroma as a
#        DIPLOID reference to separate malignant (aneuploid) from stroma.
# WHY:   Aim 1 — identify the tumour compartment before studying it.
# INPUT: data/sce.rds  (+ truth/cell_metadata.csv for scoring)
# OUTPUT: results/01_*.csv, figs/01_*.png, data/sce_clustered.rds
# HOW TO READ THE OUTPUT: PASS = the malignant/diploid calls match the planted
#        labels well (Adjusted Rand Index >= threshold). We also print how often
#        the reference cells themselves get mis-called (ref_fpr) — if that is
#        high the reference is unreliable and ALL calls are suspect (the
#        reference-baseline risk).
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment)
  library(scater); library(ggplot2)
})
source("R/utils.R"); source("R/adapters.R"); source("R/metrics.R")
source("R/cnv_fallback.R")
cfg <- load_config()
init_script("01_qc_cluster.R — Aim 1 (malignant vs stroma)", cfg)

sce <- readRDS("data/sce.rds")
log_msg("loaded", ncol(sce), "cells x", nrow(sce), "genes")

# ---- 1. QC: drop obvious low-quality cells (gentle on simulated data) -------
qc <- perCellQCMetrics(sce)
discard <- isOutlier(qc$sum, log = TRUE, type = "lower") |
           isOutlier(qc$detected, log = TRUE, type = "lower")
log_msg("QC discards", sum(discard), "of", ncol(sce), "cells")
sce <- sce[, !discard]

# ---- 2. normalise + PCA -----------------------------------------------------
sce <- logNormCounts(sce)
npc <- min(20L, ncol(sce) - 1L, nrow(sce) - 1L)
sce <- runPCA(sce, ncomponents = npc)

# ---- 3. cluster (graph-based via scran/bluster; kmeans fallback) ------------
clusters <- NULL
if (requireNamespace("scran", quietly = TRUE) &&
    requireNamespace("bluster", quietly = TRUE)) {
  log_msg("clustering: SNN graph + Louvain (scran/bluster)")
  clusters <- scran::clusterCells(sce, use.dimred = "PCA",
                BLUSPARAM = bluster::NNGraphParam(cluster.fun = "louvain"))
} else {
  log_msg("[FALLBACK] clustering: kmeans on PCs")
  k <- length(unique(sce$celltype))
  clusters <- factor(kmeans(reducedDim(sce, "PCA"), centers = k,
                            nstart = 10)$cluster)
}
sce$cluster <- clusters
log_msg("found", nlevels(clusters), "clusters (true cell types:",
        length(unique(sce$celltype)), ")")

# ---- 4. CNV calling: malignant vs diploid -----------------------------------
ref_types <- setdiff(unique(sce$celltype), "malignant")   # known diploid stroma
cnv <- run_cnv(sce, ref_types, cfg)
sce$cnv_score <- cnv$score
sce$cnv_call  <- cnv$call
log_msg(sprintf("CNV method: %s | reference false-call rate = %.1f%%",
                cnv$method, 100 * cnv$ref_fpr))
if (cnv$ref_fpr > 0.10)
  log_msg("WARNING: reference cells frequently mis-called -> reference baseline",
          "may be unreliable (see header note).")

saveRDS(sce, "data/sce_clustered.rds")

# ============================  VALIDATION (Aim 1)  ===========================
# (delete this block when you swap in real data)
truth_bin <- ifelse(sce$celltype == "malignant", "malignant", "diploid")
cm  <- confusion(truth_bin, sce$cnv_call)
ari_call    <- adjusted_rand_index(truth_bin, sce$cnv_call)   # malignant vs not
ari_cluster <- adjusted_rand_index(sce$cluster, sce$celltype) # clustering quality

cat("\n  Confusion matrix (rows=truth, cols=call):\n")
print(cm$matrix)
log_msg(sprintf("call accuracy = %.3f | ARI(call vs truth) = %.3f | ARI(cluster vs celltype) = %.3f",
                cm$accuracy, ari_call, ari_cluster))

write.csv(as.data.frame(cm$matrix), "results/01_confusion_malignant.csv",
          row.names = FALSE)
write.csv(data.frame(metric = c("call_accuracy", "ari_call_vs_truth",
                                "ari_cluster_vs_celltype", "ref_false_call_rate"),
                     value = c(cm$accuracy, ari_call, ari_cluster, cnv$ref_fpr)),
          "results/01_metrics.csv", row.names = FALSE)

# figure: CNV score by true cell type (malignant should stand out)
p <- ggplot(data.frame(celltype = sce$celltype, score = sce$cnv_score),
            aes(celltype, score, fill = celltype)) +
  geom_boxplot(outlier.size = 0.4) +
  geom_hline(yintercept = cnv$threshold, linetype = "dashed") +
  labs(title = "Aim 1: CNV aneuploidy score by cell type",
       subtitle = "dashed line = malignant-call threshold",
       y = "mean |smoothed CNV signal|") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                     legend.position = "none")
ggsave("figs/01_cnv_score_by_celltype.png", p, width = 7, height = 4, dpi = 120)
log_msg("figure -> figs/01_cnv_score_by_celltype.png")

pass_fail("01_qc_cluster", "malignant calls recover true labels",
          ari_call >= cfg$thresholds$ari_min,
          sprintf("ARI=%.3f (>=%.2f) acc=%.3f ref_fpr=%.2f",
                  ari_call, cfg$thresholds$ari_min, cm$accuracy, cnv$ref_fpr))
