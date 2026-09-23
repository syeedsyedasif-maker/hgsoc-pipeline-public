# ============================================================================
# 00_load_real_data_site.R  —  GSE180661 (MSK SPECTRUM) site-comparison loader
# ----------------------------------------------------------------------------
# WHAT:  Build data/sce_site_comparison.rds from GSE180661 (Vazquez-Garcia et
#        al. 2022, Nature, "Ovarian cancer mutational processes drive site-
#        specific immune evasion"): 41 treatment-naive HGSOC patients, 276
#        samples, profiled across Adnexa/Omentum/Ascites/Bowel/Peritoneum/UQ.
#        Produces the same colData/rowData contract as 00_load_real_data.R
#        (celltype, site, cell_subtype, rowData chr/start/stop/arm, counts +
#        logcounts) so the same downstream idioms apply.
# WHY:   A SECOND, SEPARATE real-data input for an ascites-vs-omentum cell-type
#        COMPOSITION question. This is NOT the resistance pipeline: it does not
#        read or write data/sce.rds, and 01-04 do not consume this object.
#        GSE165897 (the resistance dataset) has NO true ascites samples -- its
#        loader had to map Peritoneum -> "ascites" as a clinical stand-in, with
#        an explicit CAVEAT. GSE180661 has 108,287 genuine Ascites cells
#        collected alongside 208,894 Omentum cells, so no stand-in is needed
#        and none is invented here.
# INPUT: raw_data/GSE180661/GSE180661_matrix.h5        (32,555,276,423 bytes)
#        raw_data/GSE180661/GSE180661_GEO_cells.tsv.gz (42,728,783 bytes)
# OUTPUT: data/sce_site_comparison.rds
# HOW TO READ: table(sce$site) gives the ascites/omental/other split;
#        table(sce$patient_id, sce$site) shows the paired structure. ALWAYS
#        read composition alongside sce$sort_parameters -- see CAVEAT 1.
#
# CONFIRMED (checked directly against both files, this session -- an earlier
# draft of a GSE180661 loader in this project was written blind against assumed
# column names and was wrong, so nothing below is assumed):
#   cells.tsv.gz: 929,690 cells x 18 columns -- cell_id, sample, cell_type,
#     percent.mt, nCount_RNA, nFeature_RNA, umap50_1, umap50_2, cluster_label,
#     cluster_label_sub, cell_type_super, patient_id, tumor_subsite, tumor_site,
#     tumor_supersite, sort_parameters, therapy, surgery.
#     tumor_supersite in {Adnexa 322651, Omentum 208894, Ascites 108287,
#     Bowel 105500, Peritoneum 90852, UQ 52094, Other 41412}.
#     therapy is "pre-Rx" and surgery is "S1" for ALL cells -- there is no
#     treatment axis in this cohort, which is exactly why it is orthogonal to
#     the GSE165897 resistance question.
#     41 patients; 25 have ascites, 35 have omentum, 21 have BOTH (paired).
#   matrix.h5: NOT a 10x CellRanger HDF5 -- see the read_csr_rows() note below.
#     AnnData/h5sparse layout: /X/{data,indices,indptr} with attributes
#     h5sparse_format="csr", h5sparse_shape=[1376121, 32223]; /obs is a COMPOUND
#     holding ONLY a barcode `index` column (no metadata -- all of it comes from
#     the TSV); /var is a COMPOUND with index (HGNC symbol), gene_ids (Ensembl),
#     feature_types. All 32,223 gene symbols are unique (no dedup needed).
#     /X/data values are integer-valued raw UMI counts. /X/indices are 0-based.
#     Barcode join is exact: 929,690/929,690 (100.00%) of TSV cell_ids are found
#     in /obs/index. The h5 carries 446,431 extra cells the published annotation
#     dropped at QC; those are not recoverable as labelled cells and are ignored.
#
# CAVEAT 1 -- SORTING CONFOUND, the single most important thing about this
#   dataset for a COMPOSITION question. Samples were FACS-sorted into CD45+ and
#   CD45- fractions, and the fraction mix is NOT balanced across the two sites:
#       CD45-   Ascites 29.7%   Omentum 47.0%
#       CD45+   Ascites 70.3%   Omentum 50.9%
#       unsorted Ascites 0.0%   Omentum  2.2%
#   Raw per-site cell-type proportions therefore measure the sorting protocol as
#   much as the biology (e.g. T cells are 46% of ascites vs 28% of omentum, in
#   the direction the CD45+ over-representation alone would predict). Any
#   composition comparison MUST condition on sort_parameters (compare within
#   CD45- and within CD45+ separately) or reweight. sort_parameters is retained
#   in colData precisely so this cannot be forgotten downstream. Do not compare
#   raw table(celltype, site) proportions and call the difference biological.
# CAVEAT 2 -- the design is paired but not balanced: 21/41 patients contributed
#   both ascites and omentum, 4 are ascites-only and 14 omentum-only. A paired
#   within-patient test on the 21 is the well-controlled analysis; the unpaired
#   remainder is at best a sensitivity check.
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(S4Vectors)
  library(data.table); library(Matrix); library(scater); library(biomaRt)
  library(rhdf5)
})
source("R/utils.R")
cfg <- load_config()
init_script("00_load_real_data_site.R — GSE180661 (Vazquez-Garcia et al. 2022)", cfg)

# NOTE: no pass_fail() anywhere in this script on purpose. This is a separate
# side input, not a run_all.R stage; writing to results/validation_log.tsv would
# corrupt run_all.R's "n/n checks PASSED" summary for the main pipeline.

RAW_DIR   <- Sys.getenv("GSE180661_DIR", "raw_data/GSE180661")
h5_path   <- file.path(RAW_DIR, "GSE180661_matrix.h5")
cells_path <- file.path(RAW_DIR, "GSE180661_GEO_cells.tsv.gz")
stopifnot(
  "Set GSE180661_DIR or place GSE180661_matrix.h5 + GSE180661_GEO_cells.tsv.gz in raw_data/GSE180661 (see header)" =
    file.exists(h5_path) && file.exists(cells_path)
)

# ---- 1. cell metadata ---------------------------------------------------------
log_msg("reading cell metadata from", cells_path)
meta <- data.table::fread(cells_path, sep = "\t", header = TRUE)
log_msg("cellInfo:", nrow(meta), "annotated cells x", ncol(meta), "columns")

# ---- 2. map GEO labels -> pipeline vocabulary --------------------------------
# Mapping is done on cluster_label (63 fine-grained values), NOT on the coarse
# cell_type column, for the same reason 00_load_real_data.R maps on cell_subtype
# rather than cell_type: the coarse column silently mixes categories that must
# not share a reference signature. Confirmed by cross-tabulation:
#   cell_type "Myeloid.cell"        contains cDC1/cDC2/mDC/pDC (~17.9k dendritic cells)
#   cell_type "T.cell"              contains NK.CD56/NK.cytotoxic/Cycling.T.NK (~59.8k NK)
#   cell_type "Fibroblast"          contains Pericyte (10.2k mural cells)
#   cell_type "Ovarian.cancer.cell" contains Ciliated.cell.1/2 (4.4k normal epithelium)
# Collapsing those into macrophage/Tcell/fibroblast/malignant would bias every
# one of those signatures, so they are DROPPED rather than mis-bucketed --
# the same call the GSE165897 loader made for NK/B/Mast/DC/Plasma/ILC.
cluster_map <- c(
  setNames(rep("malignant", 8),
           c(paste0("Cancer.cell.", 1:6), paste0("Cycling.cancer.cell.", 1:2))),
  setNames(rep("fibroblast", 8),
           c("Activated.CAF.TGFb", "Mesothelial.CAF.IL1", "Activated.CAF.IGF1",
             "Angiogenic.CAF", "Activated.CAF.ISG", "Cycling.CAF",
             "Early.CAF.1", "Early.CAF.2")),
  setNames(rep("macrophage", 8),
           c("M2.SELENOP", "M2.MARCO", "M1.S100A8", "M2.CXCL10", "Cycling.M",
             "Clearing.M", "M2.COL1A1", "M2.MMP9")),
  setNames(rep("endothelial", 10),
           c("Venous.EC", "Arterial.EC", "Quiescent.EC", "Vascular.EC",
             "Sprouting.EC", "Microvascular.EC", "Cycling.EC", "CD36.EC",
             "Lymphatic.EC", "Migrating.EC")),
  setNames(rep("Tcell", 7),
           c("CD8.T.cytotoxic", "CD4.T.naive", "CD4.T.reg", "CD8.T.ISG",
             "CD8.T.CXCL13", "CD4.T.CXCL13", "CD4.Th17"))
  # DROPPED (no home in the pipeline's 5-type vocabulary, and mis-bucketing any
  # of them would bias a reference signature):
  #   NK.CD56, NK.cytotoxic, Cycling.T.NK   -- NK / ambiguous T-NK
  #   Pericyte                              -- mural, not a CAF
  #   Ciliated.cell.1, Ciliated.cell.2      -- normal tubal epithelium, not malignant
  #   B.*, Plasma.*, cDC1, cDC2, mDC, pDC, Mast.cell
)
meta[, celltype := unname(cluster_map[cluster_label])]

n_drop <- sum(is.na(meta$celltype))
log_msg("cell-type mapping: kept", sum(!is.na(meta$celltype)), "cells, DROPPED",
        n_drop, sprintf("(%.1f%%)", 100 * n_drop / nrow(meta)),
        "with no mapped cell type")
log_msg("dropped cluster_label breakdown:\n",
        paste(capture.output(sort(table(meta$cluster_label[is.na(meta$celltype)]),
                                  decreasing = TRUE)), collapse = "\n "))

# site: this dataset has REAL ascites, so it is used directly -- no clinical
# stand-in mapping of the kind GSE165897 needed (Peritoneum -> "ascites").
# Every other supersite (Adnexa/Bowel/Peritoneum/UQ/Other) becomes "other"
# rather than being force-mapped into one of the two compared compartments or
# silently dropped: keeping them makes the cohort auditable and lets a later
# analysis subset explicitly, and nothing here needs them excluded.
site_map <- c("Ascites" = "ascites", "Omentum" = "omental")
meta[, site := unname(site_map[tumor_supersite])]
meta[is.na(site), site := "other"]

meta <- meta[!is.na(celltype)]
setnames(meta, "cluster_label", "cell_subtype")          # 63 fine-grained labels
setnames(meta, "cluster_label_sub", "cell_subtype_fine") # finer still, ~76% NA
setnames(meta, "cell_type", "cell_type_coarse")          # the paper's coarse call
# cell_subtype is kept even though celltype derives from it -- same reasoning as
# the GSE165897 loader: finer CAF / T-cell / macrophage resolution may be wanted
# later, and collapsing it away here would make that unrecoverable.

# ---- 3. counts: read only the retained cells out of the CSR matrix -----------
# DropletUtils::read10xCounts() is NOT usable here and is deliberately not used.
# It expects the 10x CellRanger layout (matrix/{data,indices,indptr,shape},
# matrix/features/{id,name}, matrix/barcodes). h5ls() on this file shows an
# AnnData/h5sparse layout instead: /X/{data,indices,indptr} + /obs + /var, with
# the dimensions carried in an h5sparse_shape ATTRIBUTE rather than a shape
# dataset. Calling read10xCounts() on it fails; this reader replaces it.
#
# The second reason not to slurp the whole matrix: total nnz is 2,700,985,203,
# which exceeds .Machine$integer.max (2,147,483,647). A dgCMatrix stores its row
# indices in an integer slot, so the FULL matrix cannot be represented in R at
# all. Restricting to the cells with a mapped cell type brings nnz to
# ~1.76e9, which does fit -- with ~18% headroom. The guard below is not
# hypothetical; widening the cell-type map could push it over.
h5_shape <- h5readAttributes(h5_path, "X")$h5sparse_shape
h5_fmt   <- h5readAttributes(h5_path, "X")$h5sparse_format
stopifnot("GSE180661 matrix.h5 is not CSR -- the row-slice reader below assumes it" =
            identical(as.character(h5_fmt), "csr"))
n_cells_h5 <- as.integer(h5_shape[1]); n_genes <- as.integer(h5_shape[2])
log_msg("h5 matrix:", n_cells_h5, "cells x", n_genes, "genes (CSR, cells as rows)")

var_tbl <- as.data.frame(h5read(h5_path, "var"))
obs_ids <- as.data.frame(h5read(h5_path, "obs"))$index
stopifnot("h5 /var row count disagrees with h5sparse_shape[2]" = nrow(var_tbl) == n_genes,
          "h5 /obs row count disagrees with h5sparse_shape[1]" = length(obs_ids) == n_cells_h5)

# indptr is int64 and its last entry overflows int32 -- read as double or the
# tail silently becomes NA.
indptr <- as.numeric(h5read(h5_path, "X/indptr", bit64conversion = "double"))

row_of <- match(meta$cell_id, obs_ids)
match_rate <- mean(!is.na(row_of))
log_msg(sprintf("barcode join: %d/%d retained cells matched an h5 barcode (%.2f%%)",
                sum(!is.na(row_of)), length(row_of), 100 * match_rate))
if (match_rate < 1)
  stop("Barcode join is not exact (", sum(is.na(row_of)), " unmatched). Direct ",
       "inspection this session showed 929,690/929,690 TSV cell_ids present in ",
       "/obs/index -- investigate before trusting this object.", call. = FALSE)

# Keep h5 row order: CSR slicing reads contiguous spans, so sorted rows make the
# reads sequential rather than random-seeking across a 30GB file.
ord    <- order(row_of)
meta   <- meta[ord]
row_of <- row_of[ord]

lens <- indptr[row_of + 1] - indptr[row_of]
nnz  <- sum(lens)
log_msg("retained nnz:", format(nnz, big.mark = ","),
        sprintf("(%.1f%% of int32 limit)", 100 * nnz / .Machine$integer.max))
if (nnz >= .Machine$integer.max)
  stop("Retained non-zeros (", format(nnz, big.mark = ","), ") exceed R's int32 ",
       "limit, so a dgCMatrix cannot hold them. Narrow the cell set (e.g. restrict ",
       "to site %in% c('ascites','omental')) or switch to an on-disk ",
       "HDF5Array/DelayedArray backend.", call. = FALSE)

# Read the CSR rows in blocks. Per block we read ONE contiguous span covering
# the block's first..last retained row (which also pulls in the unretained rows
# interleaved between them) and then subset to the wanted positions -- far fewer,
# far larger reads than one h5read() per cell.
read_csr_rows <- function(block_size = 20000L) {
  x_all <- numeric(nnz); i_all <- integer(nnz)
  out_pos <- 0
  blocks <- split(seq_along(row_of), ceiling(seq_along(row_of) / block_size))
  for (bi in seq_along(blocks)) {
    b     <- blocks[[bi]]
    rows  <- row_of[b]
    lo    <- indptr[min(rows)]           # 0-based start of the span
    hi    <- indptr[max(rows) + 1]       # 0-based end (exclusive)
    span  <- h5read(h5_path, "X/data",    index = list((lo + 1):hi))
    spani <- h5read(h5_path, "X/indices", index = list((lo + 1):hi))
    # positions of the wanted rows within the span, 1-based
    sel <- sequence(lens[b], from = indptr[rows] - lo + 1)
    n   <- length(sel)
    x_all[(out_pos + 1):(out_pos + n)] <- as.numeric(span[sel])
    i_all[(out_pos + 1):(out_pos + n)] <- as.integer(spani[sel])
    out_pos <- out_pos + n
    rm(span, spani, sel); gc(verbose = FALSE)
    if (bi %% 5 == 0 || bi == length(blocks))
      log_msg(sprintf("  ... block %d/%d (%.1f%% of nnz read)",
                      bi, length(blocks), 100 * out_pos / nnz))
  }
  stopifnot("CSR read did not fill the expected number of non-zeros" = out_pos == nnz)
  list(x = x_all, i = i_all)
}
log_msg("reading counts for", length(row_of), "cells from the CSR matrix ...")
csr <- read_csr_rows()

# A CSR (cells x genes) triple IS a CSC (genes x cells) triple -- the transpose
# is free. So this builds the genes-x-cells orientation the SCE wants directly,
# with no explicit t() on a 1.76e9-element matrix.
p_vec <- c(0, cumsum(lens))
counts <- Matrix::sparseMatrix(
  i = csr$i, p = as.integer(p_vec), x = csr$x,
  dims = c(n_genes, length(row_of)), index1 = FALSE, repr = "C")
rownames(counts) <- var_tbl$index
colnames(counts) <- meta$cell_id
rm(csr, indptr); gc(verbose = FALSE)
log_msg("counts matrix:", nrow(counts), "genes x", ncol(counts), "cells")
stopifnot("counts are not integer-valued -- expected raw UMI counts" =
            all(counts@x == round(counts@x)))

meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta_df$cell_id
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = DataFrame(meta_df))
rowData(sce)$gene_id <- var_tbl$gene_ids   # Ensembl IDs, kept as a passthrough
rm(counts); gc(verbose = FALSE)

log_msg("composition by celltype:\n",
        paste(capture.output(table(sce$celltype)), collapse = "\n "))
log_msg("composition by site:\n",
        paste(capture.output(table(sce$site)), collapse = "\n "))
log_msg("sort_parameters x site (SEE CAVEAT 1 -- this is a confound):\n",
        paste(capture.output(table(sce$sort_parameters, sce$site)), collapse = "\n "))

# ---- 4. gene coordinates (rowData chr/start/stop/arm) ------------------------
# Same approach as 00_load_real_data.R: query by HGNC symbol, fail fast on a
# connectivity probe, and spot-check known genes before trusting the dedup.
log_msg("querying Ensembl (biomaRt) for gene coordinates ...")
mart <- tryCatch(
  biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl"),
  error = function(e) {
    stop("Could not reach Ensembl via biomaRt: ", conditionMessage(e),
         "\nCheck the server's network/proxy/firewall settings, or try ",
         "biomaRt::useEnsembl(..., mirror = \"useast\"/\"asia\") as a fallback. ",
         "This must be resolved before the gene-coordinate lookup can proceed.",
         call. = FALSE)
  })
probe <- biomaRt::getBM(
  attributes = c("external_gene_name", "chromosome_name"),
  filters = "external_gene_name", values = c("TP53", "BRCA1", "MYC"), mart = mart)
if (nrow(probe) < 3)
  stop("biomaRt connectivity probe failed (expected TP53/BRCA1/MYC, got ",
       nrow(probe), " rows) -- Ensembl is unreachable from this environment.",
       call. = FALSE)
log_msg("biomaRt connectivity OK")

coords <- biomaRt::getBM(
  attributes = c("external_gene_name", "chromosome_name", "start_position",
                 "end_position", "band"),
  filters = "external_gene_name", values = rownames(sce), mart = mart)

std_chr <- c(as.character(1:22), "X")
coords  <- coords[coords$chromosome_name %in% std_chr & nzchar(coords$band), ]
# getBM() can return several rows per symbol (a scaffold/alt-contig hit next to
# the primary chromosome); !duplicated() keeps whichever came first, which is
# arbitrary -- hence the spot-check below.
coords  <- coords[!duplicated(coords$external_gene_name), ]
rownames(coords) <- coords$external_gene_name

known_genes <- c(TP53 = "17", BRCA1 = "17", MYC = "8", EPCAM = "2", PTPRC = "1")
present <- intersect(names(known_genes), rownames(coords))
mismatched <- present[coords[present, "chromosome_name"] != known_genes[present]]
if (length(mismatched) > 0)
  stop("Gene-coordinate spot-check failed for: ", paste(mismatched, collapse = ", "),
       " -- dedup picked a wrong chromosome (likely a scaffold/alt-contig row). ",
       "Fix the coordinate dedup logic before trusting these coordinates.", call. = FALSE)
log_msg("gene-coordinate spot-check PASSED (", length(present), "/", length(known_genes),
        "genes checked:", paste(present, collapse = ", "), ")")

n_genes_before <- nrow(sce)
mapped <- intersect(rownames(sce), rownames(coords))
log_msg(length(mapped), "of", n_genes_before, "genes mapped to a chromosome + band",
        sprintf("(%.1f%% dropped;", 100 * (1 - length(mapped) / n_genes_before)),
        "expect a real drop -- clone-based/novel/pseudogene symbols won't resolve)")
sce <- sce[mapped, ]
coords <- coords[mapped, ]

rowData(sce)$chr   <- paste0("chr", coords$chromosome_name)
rowData(sce)$start <- coords$start_position
rowData(sce)$stop  <- coords$end_position
rowData(sce)$arm   <- paste0(coords$chromosome_name, substr(coords$band, 1, 1))

# ---- 5. normalise ------------------------------------------------------------
log_msg("normalising (scater::logNormCounts) ...")
sce <- scater::logNormCounts(sce)

# ---- 6. write ----------------------------------------------------------------
if (!dir.exists("data")) dir.create("data")
# Compressed on purpose. An earlier version of this script used compress = FALSE
# on the assumption that gzip would cost more wall-clock than the disk it saved;
# that was wrong by a wide margin. The uncompressed write produced a 147GB file
# from a 54GB in-memory object, and a 600MB slice sampled from the middle of it
# compresses 146:1 -- i.e. the payload is almost entirely redundant bytes, so
# gzip is both fast here and worth roughly an order of magnitude of disk.
saveRDS(sce, "data/sce_site_comparison.rds")
log_msg("wrote data/sce_site_comparison.rds:", nrow(sce), "genes x", ncol(sce), "cells,",
        length(unique(sce$patient_id)), "patients")

log_msg("patient x site:\n",
        paste(capture.output(table(sce$patient_id, sce$site)), collapse = "\n "))
tab <- table(sce$patient_id, sce$site)
both <- rownames(tab)[tab[, "ascites"] > 0 & tab[, "omental"] > 0]
log_msg("patients with BOTH ascites and omentum (paired design):", length(both))

# Deliberately no truth/ files and no pass_fail(): real data has no planted
# ground truth, and this object is not a run_all.R stage.
