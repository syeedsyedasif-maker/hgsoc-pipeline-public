# ============================================================================
# 00_load_real_data.R  —  REPLACES 00_simulate.R with real HGSOC scRNA-seq
# ----------------------------------------------------------------------------
# WHAT:  Build data/sce.rds from GSE165897 (Zhang et al. 2022, Science Advances,
#        "Longitudinal single-cell RNA-seq analysis reveals stress-promoted
#        chemoresistance in metastatic ovarian cancer"): 11 HGSOC patients
#        profiled treatment-naive and post-neoadjuvant-chemotherapy (post-NACT),
#        across omentum/peritoneum/mesentery/tumor sites. Produces the same
#        contract 01-04 expect: colData$celltype, colData$site, rowData
#        chr/start/stop/arm, assays "counts" + "logcounts" -- plus two NEW
#        colData columns that earlier real-data drafts of this pipeline didn't
#        have: treatment_phase (drives Aim 2's DE contrast) and cell_subtype
#        (kept as a passthrough for finer-grained stroma analysis later, e.g.
#        CAF-1/2/3 resolution for Aim 4).
# WHY:   Swaps the synthetic splatter data for a real labelled reference. This
#        REPLACES an earlier draft of this script that targeted GSE180661 (MSK
#        SPECTRUM) -- that dataset's 30GB matrix.h5 was never actually
#        downloaded/verified; GSE165897 is ~105MB, already on disk, and its
#        schema has been fully confirmed by direct inspection.
# INPUT: GSE165897_cellInfo_HGSOC.tsv.gz, GSE165897_UMIcounts_HGSOC.tsv.gz
#        (already present one directory above the repo root; see RAW_DIR below)
# OUTPUT: data/sce.rds
#
# CONFIRMED (checked directly against both files, this session):
#   cellInfo: 51,786 cells x columns cell, sample, patient_id, treatment_phase,
#     anatomical_location, cell_type, cell_subtype, nCount_RNA, nFeature_RNA,
#     percent.mt. cell_type in {Immune, EOC, Stromal}; cell_subtype has 27
#     values (CAF-1/2/3, Macrophages, T-cells, Endothelial, ..., EOC_C1..C12);
#     NO adipocyte category exists in this cohort at all. anatomical_location
#     in {Omentum, Peritoneum, Mesentery, Tumor} -- no literal "ascites" label.
#     treatment_phase in {treatment-naive, post-NACT}. 11 unique patient_id,
#     and malignant (EOC) cells are present in BOTH phases for all 11 patients
#     (relevant for 02_dge.R's paired pseudobulk design).
#   UMIcounts: genes x cells, 32,847 genes x 51,786 cells, HGNC gene SYMBOLS
#     (not Ensembl IDs, incl. legacy/clone-based names like RP11-34P13.7).
#     Column names match cellInfo$cell 1:1 exactly (verified: sorted diff is
#     empty) -- no barcode cleanup needed.
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(S4Vectors)
  library(data.table); library(Matrix); library(scater); library(biomaRt)
})
source("R/utils.R")
cfg <- load_config()
init_script("00_load_real_data.R — GSE165897 (Zhang et al. 2022)", cfg)

RAW_DIR       <- Sys.getenv("GSE165897_DIR", "..")   # files sit one dir above the repo root
cellinfo_path <- file.path(RAW_DIR, "GSE165897_cellInfo_HGSOC.tsv.gz")
counts_path   <- file.path(RAW_DIR, "GSE165897_UMIcounts_HGSOC.tsv.gz")
stopifnot(
  "Set GSE165897_DIR or place the two GSE165897_*.tsv.gz files one directory above the repo root (see header)" =
    file.exists(cellinfo_path) && file.exists(counts_path)
)

# ---- 1. counts ---------------------------------------------------------------
# fread, not read.delim: 32,847 x 51,786 is too big for base read.delim to be
# practical. UMI counts are >90% zero, so coerce to sparse immediately -- a
# dense matrix this size would waste tens of GB of RAM.
log_msg("reading counts from", counts_path, "(large file; this can take a while)")
dt <- data.table::fread(counts_path, sep = "\t", header = TRUE)
gene_ids <- dt[[1]]
counts <- as.matrix(dt[, -1, with = FALSE])
rownames(counts) <- gene_ids
counts <- Matrix::Matrix(counts, sparse = TRUE)
rm(dt); gc()
log_msg("matrix:", nrow(counts), "genes x", ncol(counts), "cells")

# ---- 2. cell metadata: join on cell / column name -----------------------------
meta <- read.delim(cellinfo_path, stringsAsFactors = FALSE)
rownames(meta) <- meta$cell
# Barcode match was verified 1:1 (sorted diff empty) by direct inspection this
# session -- no cleanup/prefix-stripping needed, just a cheap regression guard.
stopifnot("GSE165897 barcode mismatch: matrix columns and cellInfo$cell no longer line up" =
            setequal(colnames(counts), rownames(meta)))
meta <- meta[colnames(counts), ]

sce <- SingleCellExperiment(assays = list(counts = counts), colData = DataFrame(meta))
rm(counts); gc()

# ---- 3. map GEO labels -> pipeline vocabulary (celltype, site, treatment_phase)
# cell_type alone is too coarse (only Immune/EOC/Stromal); use cell_subtype for
# everything except malignant.
subtype_map <- c(
  "CAF-1"       = "fibroblast",
  "CAF-2"       = "fibroblast",
  "CAF-3"       = "fibroblast",
  "Macrophages" = "macrophage",
  "T-cells"     = "Tcell",
  "Endothelial" = "endothelial"
  # DROPPED (no home in the pipeline's fixed 4-stroma-type vocabulary):
  # Mesothelial, B-cells, NK, Mast-cells, DC-1, DC-2, Plasma-cells, pDC, ILC --
  # real immune/stromal subtypes, but lumping e.g. NK into Tcell would bias
  # that reference signature, so they're dropped rather than mis-bucketed.
  # No adipocyte category exists in this cohort at all -- expected (matches
  # the structural adipocyte blind spot already modelled in config.yml), not
  # a gap to fill.
)
celltype <- ifelse(sce$cell_type == "EOC", "malignant",
                    unname(subtype_map[sce$cell_subtype]))

# site: Omentum -> "omental", Peritoneum -> "ascites" (closest clinical analog:
# peritoneal carcinomatosis and malignant ascites are the same intraperitoneal
# spread compartment in HGSOC). CAVEAT: this merges solid peritoneal tumour
# deposits with free-floating malignant ascites, which the literature usually
# treats as distinct compartments -- acceptable ONLY because no downstream
# script (03_deconvolution.R, 04_spatial.R) currently branches on sce$site.
# Revisit this mapping before adding any site-stratified analysis.
# Mesentery (2/11 patients, treatment-naive only) and Tumor (1/11 patients,
# post-NACT only) get "other" rather than being force-mapped or dropped --
# neither is well-powered enough to serve as a real site anyway, and nothing
# downstream needs them excluded.
site_map <- c("Omentum" = "omental", "Peritoneum" = "ascites")
site <- unname(site_map[sce$anatomical_location])
site[is.na(site)] <- "other"

keep <- !is.na(celltype)
log_msg("keeping", sum(keep), "of", length(keep),
        "cells with a mapped cell type (dropped: unmapped immune/stromal subtypes)")
sce <- sce[, keep]
sce$celltype   <- celltype[keep]
sce$site       <- site[keep]
# treatment_phase, patient_id, cell_subtype pass through unchanged (cell_subtype
# is kept even though celltype already derives from it -- may be needed for
# finer CAF-1/2/3 resolution later, e.g. Aim 4).
log_msg("composition by celltype:\n",
        paste(capture.output(table(sce$celltype)), collapse = "\n "))
log_msg("composition by cell_subtype:\n",
        paste(capture.output(table(sce$cell_subtype)), collapse = "\n "))
log_msg("composition by site:\n", paste(capture.output(table(sce$site)), collapse = "\n "))
log_msg("composition by treatment_phase:\n",
        paste(capture.output(table(sce$treatment_phase)), collapse = "\n "))

# ---- 4. gene coordinates for Aim 1 CNV windowing (rowData chr/start/stop/arm)
# GSE165897 has zero Ensembl IDs (all HGNC symbols), so query by symbol only --
# no auto-detect branch needed here (unlike the earlier GSE180661 draft).
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
# Connectivity smoke-test on a handful of well-known genes before the full
# ~30k-gene query, so a network problem fails fast and clearly instead of
# hanging or silently returning an empty/partial result for the full list.
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
# getBM() can return multiple rows per symbol (e.g. a scaffold/alt-contig hit
# alongside the primary chromosome); !duplicated() keeps whichever came back
# first, which is arbitrary. Spot-check known genes below to catch a bad pick.
coords  <- coords[!duplicated(coords$external_gene_name), ]
rownames(coords) <- coords$external_gene_name

known_genes <- c(TP53 = "17", BRCA1 = "17", MYC = "8", EPCAM = "2", PTPRC = "1")
present <- intersect(names(known_genes), rownames(coords))
mismatched <- present[coords[present, "chromosome_name"] != known_genes[present]]
if (length(mismatched) > 0)
  stop("Gene-coordinate spot-check failed for: ", paste(mismatched, collapse = ", "),
       " -- dedup picked a wrong chromosome (likely a scaffold/alt-contig row). ",
       "Fix the coordinate dedup logic before trusting CNV windowing.", call. = FALSE)
log_msg("gene-coordinate spot-check passed (", length(present), "/", length(known_genes),
        "genes checked:", paste(present, collapse = ", "), ")")

mapped <- intersect(rownames(sce), rownames(coords))
log_msg(length(mapped), "of", nrow(sce), "genes mapped to a chromosome + band",
        "(unmapped genes dropped -- CNV windowing needs coordinates; expect a real",
        "drop here since legacy/clone-based/pseudogene symbols won't resolve)")
sce <- sce[mapped, ]
coords <- coords[mapped, ]

rowData(sce)$chr   <- paste0("chr", coords$chromosome_name)
rowData(sce)$start <- coords$start_position
rowData(sce)$stop  <- coords$end_position
rowData(sce)$arm   <- paste0(coords$chromosome_name, substr(coords$band, 1, 1))  # e.g. "8q"

# ---- 5. normalise (central object always carries logcounts, per adapters.R) --
sce <- scater::logNormCounts(sce)

# ---- 6. write ----------------------------------------------------------------
if (!dir.exists("data")) dir.create("data")
saveRDS(sce, "data/sce.rds")
log_msg("wrote data/sce.rds:", nrow(sce), "genes x", ncol(sce), "cells,",
        length(unique(sce$patient_id)), "patients")

# No truth/ files: real data has no planted ground truth. 01_qc_cluster.R's
# validation stays meaningful as-is (sce$celltype is a real published
# annotation, not fabricated) -- see README "Swapping in real data" for what
# changed in 02_dge.R/03_deconvolution.R/04_spatial.R instead of blanket
# deleting their VALIDATION blocks.
