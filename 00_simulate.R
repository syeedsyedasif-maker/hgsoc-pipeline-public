# ============================================================================
# 00_simulate.R  —  GENERATE THE DATA *AND* THE KNOWN TRUTH
# ----------------------------------------------------------------------------
# WHAT:  Build one simulated HGSOC single-cell dataset with splatter, give the
#        genes pseudo-genomic coordinates, overlay clonal CNVs on the malignant
#        cells, and WRITE DOWN the ground truth (resistance genes, CNV map,
#        per-cell labels) so later scripts can be scored against it.
# WHY:   This is the foundation for all four aims. Because we plant the answer
#        here, every downstream "discovery" can be checked for correctness.
# INPUT: config.yml
# OUTPUT: data/sce.rds  (the central SingleCellExperiment) and truth/*.csv
# HOW TO READ THE OUTPUT: the PASS line means the malignant cells really do show
#        higher expression on "gain" arms and lower on "loss" arms than stroma
#        (i.e. the CNV signal is present and strong), and a non-empty resistance
#        gene set was planted. A FAIL means the planted signal didn't take.
# NOTE:  This is the ONLY script you delete when you swap in real data.
# ============================================================================

suppressPackageStartupMessages({
  library(splatter); library(scater)
  library(SingleCellExperiment); library(SummarizedExperiment)
})
source("R/utils.R"); source("R/adapters.R")
cfg <- load_config()
init_script("00_simulate.R — generate truth", cfg)

G <- cfg$sizes$n_genes
N <- cfg$sizes$n_cells
gnames <- group_levels(cfg)                       # fixed group order
gp <- unlist(cfg$group_prob[gnames]); gp <- gp / sum(gp)

# ---- 1. per-group DE settings -----------------------------------------------
# Baseline DE makes every cell type distinct (identity). malignant_omental gets
# stronger/extra DE => that becomes the recoverable RESISTANCE programme (Aim 2).
de_prob   <- setNames(rep(cfg$splatter$de_prob,   length(gnames)), gnames)
de_facLoc <- setNames(rep(cfg$splatter$de_facLoc, length(gnames)), gnames)
de_prob["malignant_omental"]   <- cfg$splatter$omental_de_prob
de_facLoc["malignant_omental"] <- cfg$splatter$omental_de_facLoc

# ---- 2. simulate counts with splatter ---------------------------------------
log_msg("simulating", G, "genes x", N, "cells across", length(gnames), "groups")
params <- newSplatParams()
params <- setParams(params,
  nGenes = G, batchCells = N, group.prob = as.numeric(gp),
  de.prob = as.numeric(de_prob), de.facLoc = as.numeric(de_facLoc),
  de.facScale = cfg$splatter$de_facScale, seed = cfg$seed)
sim <- splatSimulate(params, method = "groups", verbose = FALSE)

# relabel Group1..K -> our names; add coarse celltype (collapses the malignants)
grp_idx <- as.integer(sub("Group", "", as.character(sim$Group)))
sim$group    <- factor(gnames[grp_idx], levels = gnames)
sim$celltype <- celltype_of_group(as.character(sim$group))

# ---- 3. site as a property of EVERY cell ------------------------------------
# malignant groups are site-bound; shared types drawn to match per-site abundance.
site <- character(N)
site[sim$group == "malignant_ascites"] <- "ascites"
site[sim$group == "malignant_omental"] <- "omental"
for (ct in names(cfg$site_omental_prob)) {
  idx <- which(sim$group == ct)
  p_om <- cfg$site_omental_prob[[ct]]
  site[idx] <- sample(c("omental", "ascites"), length(idx), replace = TRUE,
                      prob = c(p_om, 1 - p_om))
}
sim$site <- site

# ---- 4. map genes -> chromosome arms (pseudo-genome, in gene order) ---------
chrs <- c(as.character(1:22), "X")
arms <- unlist(lapply(chrs, function(c) c(paste0(c, "p"), paste0(c, "q"))))
arm_of_gene <- arms[as.integer(cut(seq_len(G), breaks = length(arms),
                                   labels = FALSE))]
chr_of_gene <- paste0("chr", sub("[pq]$", "", arm_of_gene))
start <- integer(G); stop <- integer(G)
for (ch in unique(chr_of_gene)) {                 # increasing position within chr
  idx <- which(chr_of_gene == ch)
  start[idx] <- seq_along(idx) * 1000L; stop[idx] <- start[idx] + 800L
}
rowData(sim)$chr <- chr_of_gene; rowData(sim)$arm <- arm_of_gene
rowData(sim)$start <- start;     rowData(sim)$stop <- stop

# ---- 5. overlay CNV on malignant cells (clonal: both malignant groups) ------
arm_factor <- setNames(rep(1, length(arms)), arms)
arm_factor[cfg$cnv$gains]  <- cfg$cnv$gain_factor
arm_factor[cfg$cnv$losses] <- cfg$cnv$loss_factor
gene_factor <- arm_factor[arm_of_gene]            # length G
mal_cells <- which(sim$celltype == "malignant")
cts <- counts(sim)
cts[, mal_cells] <- round(sweep(cts[, mal_cells, drop = FALSE], 1,
                                gene_factor, "*"))
counts(sim) <- cts
rowData(sim)$cnv_factor <- gene_factor
log_msg("CNV overlaid on", length(mal_cells), "malignant cells;",
        sum(arm_factor > 1), "gain arms,", sum(arm_factor < 1), "loss arms")

# ---- 6. normalise (so the central object always carries logcounts) ----------
sim <- logNormCounts(sim)

# ---- 7. extract & write the planted truth -----------------------------------
gi_om  <- which(gnames == "malignant_omental")
gi_asc <- which(gnames == "malignant_ascites")
defac_om  <- rowData(sim)[[paste0("DEFacGroup", gi_om)]]
defac_asc <- rowData(sim)[[paste0("DEFacGroup", gi_asc)]]
log2fc <- log2(defac_om / defac_asc)
is_resist <- abs(log2fc) >= cfg$truth$resistance_log2fc_threshold
resist_tbl <- data.frame(gene = rownames(sim)[is_resist],
                         log2fc = log2fc[is_resist],
                         direction = ifelse(log2fc[is_resist] > 0,
                                            "up_in_omental", "down_in_omental"))
write.csv(resist_tbl, "truth/resistance_genes.csv", row.names = FALSE)

cnv_map <- data.frame(arm = names(arm_factor), factor = as.numeric(arm_factor),
                      state = ifelse(arm_factor > 1, "gain",
                                ifelse(arm_factor < 1, "loss", "neutral")))
write.csv(cnv_map, "truth/cnv_map.csv", row.names = FALSE)

cell_meta <- data.frame(cell = colnames(sim),
                        celltype = sim$celltype, group = as.character(sim$group),
                        site = sim$site,
                        aneuploid = sim$celltype == "malignant")
write.csv(cell_meta, "truth/cell_metadata.csv", row.names = FALSE)

go <- gene_order_df(sim)
write.table(go[, c("gene", "chr", "start", "stop")], "truth/gene_order.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

prop_tbl <- as.data.frame(table(site = sim$site, celltype = sim$celltype))
write.csv(prop_tbl, "truth/composition_by_site.csv", row.names = FALSE)

saveRDS(sim, "data/sce.rds")
save_peek(resist_tbl, "results/peek_resistance_truth.csv")
save_peek(cnv_map[cnv_map$state != "neutral", ], "results/peek_cnv_truth.csv")

# ---- 8. validate the planted signal really took -----------------------------
logc <- get_logcounts(sim)
stroma <- which(sim$celltype != "malignant")
gain_genes <- which(gene_factor > 1); loss_genes <- which(gene_factor < 1)
gain_effect <- mean(logc[gain_genes, mal_cells]) - mean(logc[gain_genes, stroma])
loss_effect <- mean(logc[loss_genes, mal_cells]) - mean(logc[loss_genes, stroma])
log_msg(sprintf("gain-arm effect (mal-stroma logcounts) = %+.2f (want > 0)", gain_effect))
log_msg(sprintf("loss-arm effect (mal-stroma logcounts) = %+.2f (want < 0)", loss_effect))
log_msg("resistance genes planted:", nrow(resist_tbl))

ok <- file.exists("data/sce.rds") &&
      all(file.exists(c("truth/resistance_genes.csv", "truth/cnv_map.csv",
                        "truth/cell_metadata.csv", "truth/gene_order.tsv"))) &&
      nrow(resist_tbl) > 0 && gain_effect > 0 && loss_effect < 0
pass_fail("00_simulate", "truth planted (CNV signal + resistance set)", ok,
          sprintf("n_resist=%d gain=%+.2f loss=%+.2f", nrow(resist_tbl),
                  gain_effect, loss_effect))
