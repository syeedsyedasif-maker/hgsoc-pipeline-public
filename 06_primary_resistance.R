# ============================================================================
# 06_primary_resistance.R  —  primary resistant-vs-sensitive analysis (GSE165897)
# ----------------------------------------------------------------------------
# WHAT:  Implements docs/primary_spec.md: pre-specified pseudobulk gene-set tests
#        P1 (CAF, cheon2014_collagen) and P2 (malignant, resistance_up_intrinsic)
#        in treatment-naive samples, PFI resistant vs sensitive, limma::fry,
#        Holm across the two. Plus the stress positive control, secondary
#        analyses S1-S5 (plan s6.1-6.5) and the sensitivity analyses added in
#        DEVIATIONS.md D1, D4, D5.
# WHY:   docs/ANALYSIS_PLAN.md v1.1 s3-s6 and DEVIATIONS.md D1-D5. Every
#        parameter is in config.yml -> real: primary:.
# INPUT: data/sce_clustered.rds (author celltype + Stage 01 cnv_call),
#        metadata/gse165897_labels.csv, metadata/gse165897_sample_flags.csv,
#        cisplatin_response_genes_HGSOC.json
# OUTPUT (results/06_primary/, or $HGSOC_OUT_DIR):
#        06_pseudobulk_inventory.csv  every pseudobulk, included or why not
#        06_geneset_tests.csv         one row per gene-set test (+ adjusted p)
#        06_genelevel_<celltype>.csv  S5 gene-level tables (exploratory)
#        06_sessionInfo.txt
# HOW TO READ THE OUTPUT: the primary result is the `positive` column of the
#        P1/P2 rows of 06_geneset_tests.csv (Holm p < alpha AND Direction Up).
#        The console prints counts only, never test statistics. Rows with
#        estimable = FALSE carry the reason (plan s3 minimums).
# TEST MODE: HGSOC_SHUFFLE_LABELS=<seed> permutes pfi_group and crs_group across
#        patients before any analysis (to check the code without seeing the
#        true comparison). It refuses to write into config's out_dir, and every
#        output row is stamped with run_mode.
# CAVEAT: n = 11 patients (7 resistant / 4 sensitive by PFI). Results are
#        hypothesis-generating (plan s8). Open ambiguities, each with the
#        default implemented here, are listed in docs/primary_spec.md s4.
# Run:   SIM_PROFILE=real Rscript 06_primary_resistance.R
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment)
  library(edgeR); library(limma); library(data.table); library(jsonlite); library(Matrix)
})
source("R/utils.R"); source("R/adapters.R")
options(warn = 1)                               # print warnings as they happen
cfg <- load_config()
if (cfg$.profile != "real" || is.null(cfg$primary))
  stop("06_primary_resistance.R must run with SIM_PROFILE=real (needs config real: primary:)")
init_script("06_primary_resistance.R — primary resistant-vs-sensitive (plan v1.1 + D1-D5)", cfg)
pc <- cfg$primary

# ---- 0. output location + test mode ------------------------------------------
abs_path <- function(p) { p <- path.expand(p); if (!startsWith(p, "/")) p <- file.path(getwd(), p)
                          gsub("/+", "/", sub("/+$", "", p)) }
out_dir <- Sys.getenv("HGSOC_OUT_DIR", pc$out_dir)
shuffle <- Sys.getenv("HGSOC_SHUFFLE_LABELS", "")
if (nzchar(shuffle) && startsWith(paste0(abs_path(out_dir), "/"), paste0(abs_path(pc$out_dir), "/")))
  stop("Shuffled-label runs must set HGSOC_OUT_DIR outside ", pc$out_dir)
run_mode <- if (nzchar(shuffle)) paste0("SHUFFLED_LABELS_seed", shuffle) else "TRUE_LABELS"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_msg("run mode:", run_mode, "| outputs ->", out_dir)

# ---- 1. inputs ----------------------------------------------------------------
sce <- readRDS(pc$input_sce)
cd  <- as.data.table(as.data.frame(colData(sce)))
labels <- fread(pc$labels, colClasses = "character")
stopifnot("labels must cover exactly the SCE patients" = setequal(labels$patient_id, unique(cd$patient_id)),
          all(labels$pfi_group %in% c("resistant", "sensitive")),
          all(labels$crs_group %in% c("poor_response", "response")))
if (nzchar(shuffle)) {
  set.seed(as.integer(shuffle))
  labels[, pfi_group := sample(pfi_group)]
  labels[, crs_group := sample(crs_group)]
  set.seed(cfg$seed)
  log_msg("LABEL-SHUFFLE TEST MODE: pfi_group and crs_group permuted -- NOT the true comparison")
}
flags <- fread(pc$sample_flags)
stopifnot("sample flags must match SCE samples" = setequal(flags$sample, unique(cd$sample)))
fw <- dcast(flags, patient_id ~ timepoint, value.var = c("site", "frozen"))
stopifnot("D4 patient list disagrees with sample flags" =
            setequal(fw[frozen_TN != frozen_PN, patient_id], pc$freezing_discordant),
          "D5(a) patient list disagrees with sample flags" =
            setequal(fw[site_TN == site_PN, patient_id], pc$site_matched))
baseline_site <- setNames(ifelse(fw$site_TN == "Peritoneum", "peritoneum", "other"), fw$patient_id)

gs <- jsonlite::fromJSON(pc$genesets)$genesets
needed <- c(pc$p1_set, pc$p2_set, pc$pc_set, pc$s3_sets, pc$genelevel_only_sets)
stopifnot("gene set missing from JSON" = all(needed %in% names(gs)))

PHASES <- c("treatment-naive", "post-NACT")
celltypes <- sort(unique(cd$celltype))
pops <- c(setNames(lapply(celltypes, function(ct) which(cd$celltype == ct)), celltypes),
          list(malignant_both = which(cd$celltype == "malignant" & cd$cnv_call == "malignant")))  # D1
PB <- lapply(pops, function(ix) pseudobulk_counts(sce, ix))
rm(sce); invisible(gc())

# ---- 2. one analysis = select pseudobulks -> check minimums -> voom -> fry ----
inventory <- list(); tests <- list()
LEVELS <- list(pfi_group = c(ref = "sensitive", test = "resistant"),
               crs_group = c(ref = "response",  test = "poor_response"))

run_analysis <- function(id, family, population, design, label = "pfi_group",
                         phases = "treatment-naive", sets = character(), expected = NA_character_,
                         restrict = NULL, exclude = character(), site_covariate = FALSE,
                         exempt = FALSE, gene_level = FALSE, note = "") {
  pb <- PB[[population]]
  m <- as.data.table(pb$meta)[, col := rownames(pb$meta)][treatment_phase %in% phases]
  m[, group := labels[[label]][match(patient_id, labels$patient_id)]]
  m[, reason := NA_character_]
  m[patient_id %in% exclude, reason := "excluded by analysis definition"]
  if (!is.null(restrict)) m[is.na(reason) & !patient_id %in% restrict, reason := "outside restricted patient set"]
  m[is.na(reason) & n_cells < pc$min_cells, reason := sprintf("fewer than %d cells", pc$min_cells)]
  if (design %in% c("paired", "interaction")) {
    paired_ok <- m[is.na(reason), .N, by = patient_id][N == length(phases), patient_id]
    m[is.na(reason) & !patient_id %in% paired_ok, reason := "other timepoint excluded (unpaired)"]
  }
  m[, included := is.na(reason)]
  inventory[[id]] <<- m[, .(run_mode, analysis_id = id, population, patient_id, treatment_phase,
                            group = if (design == "paired") NA_character_ else group,
                            n_cells, included, reason)]
  k <- m[included == TRUE]
  lv <- LEVELS[[label]]
  if (design == "paired") {
    n_test <- uniqueN(k$patient_id); n_ref <- NA_integer_
    enough <- n_test >= pc$min_patients_per_group
    grp_txt <- sprintf("pairs=%d", n_test)
  } else {
    pts <- unique(k[, .(patient_id, group)])
    n_test <- sum(pts$group == lv[["test"]]); n_ref <- sum(pts$group == lv[["ref"]])
    enough <- min(n_test, n_ref) >= pc$min_patients_per_group
    grp_txt <- sprintf("%s=%d %s=%d", lv[["test"]], n_test, lv[["ref"]], n_ref)
  }
  log_msg(sprintf("%-14s %-15s %-11s pseudobulks=%2d | %s | cells/pseudobulk %s",
                  id, population, design, nrow(k), grp_txt,
                  if (nrow(k)) sprintf("%d-%d", min(k$n_cells), max(k$n_cells)) else "-"))

  base <- data.table(run_mode, analysis_id = id, family, population, design,
                     grouping = if (design == "paired") "treatment_phase" else label,
                     test_level = if (design == "paired") "post-NACT" else lv[["test"]],
                     ref_level  = if (design == "paired") "treatment-naive" else lv[["ref"]],
                     n_patients_test = n_test, n_patients_ref = n_ref, n_pseudobulk = nrow(k),
                     site_covariate, expected = expected, note = note)
  set_names <- if (length(sets)) sets else NA_character_
  skip <- function(why) {
    log_msg(sprintf("%-14s NOT ESTIMABLE: %s", id, why))
    tests[[id]] <<- cbind(base, gene_set = set_names, estimable = FALSE, tested = FALSE, reason = why)
    invisible(NULL)
  }
  if (!enough) return(skip(sprintf("fewer than %d patients per group", pc$min_patients_per_group)))

  # design matrices (spec s2; A3 fixed patient effect; A11 site covariate)
  if (design == "group") {
    g <- factor(k$group, levels = unname(lv))
    X <- if (site_covariate) {
      site <- factor(baseline_site[k$patient_id], levels = c("other", "peritoneum"))
      model.matrix(~ site + g)
    } else model.matrix(~ g)
    contr <- as.numeric(colnames(X) == paste0("g", lv[["test"]]))
    filt_group <- g
  } else if (design == "paired") {
    tp <- factor(k$treatment_phase, levels = PHASES)
    X <- model.matrix(~ factor(patient_id) + tp, data = k)
    contr <- as.numeric(colnames(X) == "tppost-NACT")
    filt_group <- tp
  } else {                                           # interaction: ~ patient + group:timepoint
    post <- k$treatment_phase == "post-NACT"
    X <- cbind(model.matrix(~ 0 + factor(patient_id), data = k),
               test_PN = as.numeric(post & k$group == lv[["test"]]),
               ref_PN  = as.numeric(post & k$group == lv[["ref"]]))
    contr <- c(rep(0, ncol(X) - 2), 1, -1)
    filt_group <- interaction(k$group, k$treatment_phase, drop = TRUE)
  }
  if (qr(X)$rank < ncol(X)) return(skip("design matrix rank-deficient"))
  if (nrow(X) - ncol(X) < 1) return(skip("no residual degrees of freedom"))

  y <- DGEList(pb$counts[, k$col, drop = FALSE])
  keep <- filterByExpr(y, group = filt_group)                      # spec A2
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(y, X)
  cf <- contrasts.fit(lmFit(v, X), contr)

  rows <- lapply(sets, function(s) {
    surv <- intersect(gs[[s]], rownames(v))
    r <- cbind(base, gene_set = s, set_size = length(gs[[s]]), genes_passing_filter = length(surv),
               estimable = TRUE)
    if (length(surv) == 0 || (length(surv) < pc$min_set_genes && !exempt))
      return(cbind(r, tested = FALSE,
                   reason = sprintf("fewer than %d set genes pass filterByExpr", pc$min_set_genes)))
    fr <- fry(v, index = list(surv), design = X, contrast = contr)
    cbind(r, tested = TRUE, reason = NA_character_, NGenes = fr$NGenes, Direction = as.character(fr$Direction),
          PValue = fr$PValue, mean_log2FC = mean(cf$coefficients[surv, 1]))  # spec A8
  })
  if (length(rows)) tests[[id]] <<- rbindlist(rows, fill = TRUE)

  if (gene_level) {                                                # S5, exploratory
    tt <- as.data.table(topTable(eBayes(cf), number = Inf, sort.by = "P"), keep.rownames = "gene")
    for (s in names(gs)) tt[, (paste0("in_", s)) := gene %in% gs[[s]]]
    tt[, `:=`(run_mode = run_mode, analysis_id = id, n_patients_test = n_test, n_patients_ref = n_ref)]
    fwrite(tt, file.path(out_dir, sprintf("06_genelevel_%s.csv", population)))
  }
  invisible(NULL)
}

# ---- 3. the analyses (docs/primary_spec.md s2) ----------------------------------
log_msg("--- primary (plan s4) ---")
run_analysis("P1", "primary", "fibroblast", "group", sets = pc$p1_set, expected = "Up")
run_analysis("P2", "primary", "malignant",  "group", sets = pc$p2_set, expected = "Up")
log_msg("--- positive control (plan s5) ---")
run_analysis("PC", "positive_control", "malignant", "paired", phases = PHASES,
             sets = pc$pc_set, expected = "Up", exempt = TRUE)
log_msg("--- secondary (plan s6.1-6.5) ---")
run_analysis("S1-P1", "S1", "fibroblast", "group", label = "crs_group", sets = pc$p1_set, expected = "Up",
             note = "D2: expected not estimable")
run_analysis("S1-P2", "S1", "malignant",  "group", label = "crs_group", sets = pc$p2_set, expected = "Up")
run_analysis("S2-P2", "S2", "malignant_both", "group", sets = pc$p2_set, expected = "Up",
             note = "D1: author EOC AND inferCNV-called malignant")
for (p in pc$s3_populations) run_analysis(paste0("S3-", p), "S3", p, "group", sets = pc$s3_sets)
run_analysis("S4-CAF", "S4", "fibroblast", "interaction", phases = PHASES, sets = pc$p1_set)
run_analysis("S4-mal", "S4", "malignant",  "interaction", phases = PHASES, sets = pc$p2_set)
for (ct in celltypes) run_analysis(paste0("S5-", ct), "S5", ct, "group", gene_level = TRUE)
log_msg("--- sensitivity (D4, D5) ---")
run_analysis("D4-PC", "sensitivity", "malignant", "paired", phases = PHASES, sets = pc$pc_set,
             expected = "Up", exempt = TRUE, exclude = pc$freezing_discordant, note = "D4")
run_analysis("D5a-PC", "sensitivity_direction_only", "malignant", "paired", phases = PHASES,
             sets = pc$pc_set, expected = "Up", exempt = TRUE, restrict = pc$site_matched,
             note = "D5(a): direction check only")
run_analysis("D5a-S4-CAF", "sensitivity_direction_only", "fibroblast", "interaction", phases = PHASES,
             sets = pc$p1_set, restrict = pc$site_matched, note = "D5(a): direction check only")
run_analysis("D5a-S4-mal", "sensitivity_direction_only", "malignant", "interaction", phases = PHASES,
             sets = pc$p2_set, restrict = pc$site_matched, note = "D5(a): direction check only")
run_analysis("D5b-P1", "sensitivity", "fibroblast", "group", sets = pc$p1_set, expected = "Up",
             site_covariate = TRUE, note = "D5(b): baseline site covariate")
run_analysis("D5b-P2", "sensitivity", "malignant", "group", sets = pc$p2_set, expected = "Up",
             site_covariate = TRUE, note = "D5(b): baseline site covariate")

# ---- 4. multiplicity (spec s2, A5, A15) + decision columns -----------------------
res <- rbindlist(tests, fill = TRUE)
for (col in c("PValue", "mean_log2FC")) if (!col %in% names(res)) res[, (col) := NA_real_]
if (!"Direction" %in% names(res)) res[, Direction := NA_character_]
res[, `:=`(p_adj = NA_real_, adj_method = "none")]
res[family == "primary" & tested == TRUE,
    `:=`(p_adj = p.adjust(PValue, "holm", n = 2L), adj_method = "holm (n=2)")]
for (f in c("S1", "S2", "S3", "S4"))
  res[family == f & tested == TRUE, `:=`(p_adj = p.adjust(PValue, "BH"), adj_method = "BH")]
res[, positive := NA]
res[family == "primary", positive := tested == TRUE & p_adj < pc$alpha & Direction == expected]
res[, check_passed := NA]
res[analysis_id == "PC", check_passed := tested == TRUE & PValue < pc$alpha & Direction == expected]  # A4
fwrite(res, file.path(out_dir, "06_geneset_tests.csv"))   # S5 rows appear only when not estimable
fwrite(rbindlist(inventory), file.path(out_dir, "06_pseudobulk_inventory.csv"))
save_session_info(file.path(out_dir, "06_sessionInfo.txt"))

ne <- res[estimable == FALSE, unique(paste0(analysis_id, " (", reason, ")"))]
log_msg(sprintf("done: %d set-test rows written; %d analyses not estimable%s", nrow(res),
                length(ne), if (length(ne)) paste0(": ", paste(ne, collapse = "; ")) else ""))
log_msg("test statistics are in", file.path(out_dir, "06_geneset_tests.csv"), "(not printed)")
