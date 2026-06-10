# ============================================================================
# R/adapters.R — the ONLY place the central SingleCellExperiment is reshaped
# into the formats each downstream tool wants. Keeps conversions in one spot
# so they can't drift out of sync between scripts.
#
#   get_counts / get_logcounts ........ pull matrices out of the SCE
#   gene_order_df ..................... gene -> chr/start/stop (+ arm) table
#   write_infercnv_files .............. inferCNV counts + annotation + gene order
#   build_reference_signatures ........ genes x celltypes mean profile (for NNLS)
#   as_bayesprism ..................... matrices in BayesPrism's expected shape
#   as_spatialrna ..................... coords + counts in spacexr's expected shape
# Tool OBJECTS (CreateInfercnvObject, SpatialRNA, ...) are built inside the tool
# wrappers, guarded by requireNamespace, so this file needs no heavy packages.
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

get_counts    <- function(sce) as.matrix(SummarizedExperiment::assay(sce, "counts"))
get_logcounts <- function(sce) as.matrix(SummarizedExperiment::assay(sce, "logcounts"))

# rowData(sce) carries chr/start/stop/arm planted in 00_simulate.R.
gene_order_df <- function(sce) {
  rd <- as.data.frame(SummarizedExperiment::rowData(sce))
  data.frame(gene = rownames(sce),
             chr = rd$chr, start = rd$start, stop = rd$stop, arm = rd$arm,
             stringsAsFactors = FALSE)
}

# Write the three files inferCNV expects.
#   - counts matrix (genes x cells), tab-separated
#   - annotation file: cell <tab> group   (group used to pick the diploid ref)
#   - gene ordering : gene <tab> chr <tab> start <tab> stop   (no header)
# `annot_col` is the colData column holding the group label (e.g. celltype).
write_infercnv_files <- function(sce, out_dir, annot_col = "celltype") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  counts <- get_counts(sce)
  f_counts <- file.path(out_dir, "infercnv_counts.tsv")
  f_annot  <- file.path(out_dir, "infercnv_annotations.tsv")
  f_order  <- file.path(out_dir, "infercnv_gene_order.tsv")

  utils::write.table(counts, f_counts, sep = "\t", quote = FALSE,
                     col.names = NA)                       # leading tab for rownames
  annot <- data.frame(cell = colnames(sce),
                      group = SummarizedExperiment::colData(sce)[[annot_col]])
  utils::write.table(annot, f_annot, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  go <- gene_order_df(sce)[, c("gene", "chr", "start", "stop")]
  utils::write.table(go, f_order, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  list(counts = f_counts, annotations = f_annot, gene_order = f_order)
}

# genes x celltypes matrix of MEAN counts per cell — the reference "signature"
# the NNLS / BayesPrism deconvolution leans on. `type_col` collapses groups
# (e.g. both malignant_* -> "malignant") when set to a colData column.
build_reference_signatures <- function(sce, type_col = "celltype") {
  counts <- get_counts(sce)
  labels <- as.character(SummarizedExperiment::colData(sce)[[type_col]])
  types  <- sort(unique(labels))
  sig <- sapply(types, function(t) rowMeans(counts[, labels == t, drop = FALSE]))
  colnames(sig) <- types
  sig  # genes x types
}

# Shape inputs the way BayesPrism's run.prism() wants them.
#   reference: cells x genes ; bulk: samples x genes ; cell.type.labels: per cell
as_bayesprism <- function(sce_ref, bulk_mat, type_col = "celltype") {
  ref <- t(get_counts(sce_ref))                       # cells x genes
  labels <- as.character(SummarizedExperiment::colData(sce_ref)[[type_col]])
  common <- intersect(colnames(ref), colnames(bulk_mat))
  list(reference = ref[, common, drop = FALSE],
       bulk = bulk_mat[, common, drop = FALSE],
       cell.type.labels = labels)
}

# Shape inputs for spacexr: coords (data.frame x,y with spot rownames), spot
# counts (genes x spots, integer), and a single-cell reference (counts genes x
# cells + a NAMED cell-type factor whose names match the reference cell columns).
as_spatialrna <- function(coords, spot_counts, sce_ref, type_col = "celltype") {
  ref_counts <- round(get_counts(sce_ref))
  ref_types <- factor(as.character(SummarizedExperiment::colData(sce_ref)[[type_col]]))
  names(ref_types) <- colnames(sce_ref)                 # spacexr needs names
  cd <- as.data.frame(coords); rownames(cd) <- colnames(spot_counts)
  list(coords = cd,
       spot_counts = round(spot_counts),                # spacexr expects counts
       ref_counts = ref_counts,
       ref_types = ref_types)
}
