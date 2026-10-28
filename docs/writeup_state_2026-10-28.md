# WRITEUP_STATE — evidence pack

Generated 2026-10-28 15:13 GMT on rbgo-server2 by a read-only collector. Nothing was edited, rerun, committed or deleted.
Every number below is copied from a file on disk; the source file is named next to it. Where a planned output does not exist, it says **NOT RUN**.

## 1. Repo state

- Repo: `hgsoc-pipeline-public/`
- Branch: `real-data-migration`  (HEAD `70a0ecf`, 2026-09-23 13:29:27 +0100)

`git log --oneline --since=2026-09-01`
```
70a0ecf Loader: store logcounts as sparse dgCMatrix
f71ac90 Export hgsoc_ccc conda env (liana-py 1.10.0, Python 3.13.2)
5569990 Pin liana-py environment for cell-cell communication (60 packages)
43db3c7 Analysis plan v1.1: gene-set amendments before results
6a55e73 Pre-specified analysis plan v1.0
ea9a9a9 Migrate pipeline to real data (GSE165897 + TCGA-OV) and add site/signature stages
```

`git status --short` (captured before this pack was written)
```
 M site_load_log.txt
?? WRITEUP_STATE.md
?? writeup_bundle/
```

The one modified file, `site_load_log.txt`, is a re-run of `00_load_real_data_site.R` (GSE180661) that **failed** at the biomaRt step ("Could not reach Ensembl via biomaRt", exit status 1). The committed version of the log records the earlier successful run.

## 2. ANALYSIS_PLAN.md (verbatim)

- Path: `docs/ANALYSIS_PLAN.md`
- Version (from file header): **1.1**
- Last commit: `43db3c7 2026-09-23 13:08:17 +0100 Analysis plan v1.1: gene-set amendments before results`
- All commits touching it:
```
43db3c7 2026-09-23 13:08:17 +0100 Analysis plan v1.1: gene-set amendments before results
6a55e73 2026-09-23 12:38:54 +0100 Pre-specified analysis plan v1.0
```

````markdown
# Pre-specified analysis plan

**Project:** Patient-level single-cell analysis of platinum resistance in HGSOC (GSE165897), with validation in TCGA-OV
**Author:** Syeed Syed Asif
**Version:** 1.1, 23 September 2026. Changes from 1.0 are listed in DEVIATIONS.md; all were made before any resistant-versus-sensitive result was seen.
**Status:** Committed before any resistant-versus-sensitive comparison was run. Earlier Stage 2 runs compared post-NACT with treatment-naive samples only.

## 1. Question

Do malignant cells and cancer-associated fibroblasts (CAFs) differ between platinum-resistant and platinum-sensitive patients before treatment, and does a signature derived from that difference associate with platinum response in TCGA-OV?

## 2. Data and labels

- **Single-cell data:** GSE165897 count data (GEO), 11 patients with paired treatment-naive (TN) and post-NACT (PN) samples.
- **Primary label (PFI):** resistant if platinum-free interval < 183 days; sensitive if >= 183 days. Patients without a recorded PFI are excluded from primary analyses. Values are taken from Zhang et al. (2022), Table 1, and stored in `metadata/gse165897_labels.csv` with the source noted per row.
- **Secondary label (CRS):** chemotherapy response score (Böhm et al., 2015). CRS 3 = responder; CRS 1 or 2 = non-responder.
- **Cell types:** author annotations. CAFs = the author fibroblast annotation. Malignant cells = author epithelial ovarian cancer (EOC) annotation for primary analyses; inferCNV-confirmed malignant cells as a sensitivity analysis.

## 3. Pseudobulk construction

- Raw counts summed per patient x timepoint x cell type.
- A pseudobulk sample is included only if it contains >= 20 cells.
- A contrast is run only if >= 3 patients per group remain. Otherwise it is reported as not estimable.
- Genes filtered with `edgeR::filterByExpr`; TMM normalisation; `limma::voom`.
- A gene set is tested only if >= 5 of its genes survive filtering (the positive control is exempt).

## 4. Primary analyses

Two tests, both in TN samples, PFI labels, using `limma::fry` (directional), Holm-adjusted across the two:

| ID | Cell type | Gene set | Expected direction in resistant |
|---|---|---|---|
| P1 | CAFs | Cheon et al. (2014) collagen-remodelling signature (`cheon2014_collagen`, 10 genes) | Up |
| P2 | Malignant | `resistance_up_intrinsic` (14 genes; `resistance_up_core` without the stress and extrinsic genes) | Up |

A primary result is positive only if Holm-adjusted p < 0.05 in the expected direction. Effect sizes are reported as the mean log2 fold change across set genes, with the number of patients per group.

## 5. Positive control

`stress_IEG_state_GSE165897` (FOS, JUN, NR4A1) should be higher in PN than TN malignant cells (paired, patient as blocking factor). If this is not reproduced, the pipeline is checked before any other result is interpreted.

## 6. Secondary analyses (exploratory; Benjamini-Hochberg within each family)

1. P1 and P2 repeated with CRS labels.
2. P2 repeated with inferCNV-confirmed malignant cells.
3. `drug_efflux` and `extrinsic_TME` in malignant cells, CAFs and macrophages, where >= 5 genes survive filtering. Sets below the 5-gene minimum (`matrix_CAF`, `sensitivity_up`, `partial_EMT_Xu2022`) are reported only within the gene-level results.
4. Resistance x timepoint interaction: paired design, `~ patient + group:timepoint`.
5. Gene-level differential expression per cell type. Reported as exploratory only; no gene-level result is claimed as a finding.
6. Cell-type composition: `speckle::propeller` on per-sample proportions (logit transform).
7. Cell-cell communication in TN samples, inferred per patient. Pre-specified axes: malignant to CAF, CAF to malignant, macrophage to malignant. Patient-level interaction scores compared between groups (Wilcoxon). Tool: liana-py (rank_aggregate by sample).

## 7. TCGA-OV validation (run once)

- **Signature rule:** genes from the primary gene set with log2 fold change in the expected direction and nominal p < 0.05 in that primary contrast. If the primary test is not significant, or fewer than 5 genes qualify, the full pre-specified gene set is used and reported as a literature-set test, not a single-cell-derived validation.
- **Scoring:** GSVA on TCGA-OV primary-tumour counts (Poisson kernel), gene IDs harmonised with `harmonize_bulk_to_symbols()`.
- **Outcome 1:** platinum status (resistant versus sensitive; "too early" and unknown excluded), logistic regression.
- **Outcome 2:** progression-free survival, Cox proportional hazards.
- **Covariates (both):** age, stage, residual disease.
- **Reporting:** odds ratio or hazard ratio per 1 SD of score, with 95% CI; Holm-adjusted across the two outcomes.
- **Sensitivity:** both models repeated with the ESTIMATE stromal score as an added covariate.
- **Label source:** a curated reannotation of TCGA-OV platinum status if obtainable; otherwise TCGA clinical platinum status. The choice is recorded in `DEVIATIONS.md` before the model is run.
- **Not permitted:** re-selecting genes, re-weighting or tuning on TCGA; machine-learning models.

## 8. Claims

Results are hypothesis-generating given n <= 11. No claim of clinical utility is made.

## 9. Deviations

Any change after this commit is logged in `DEVIATIONS.md` with the date, the reason, and whether it was made before or after seeing the affected result.
````

`DEVIATIONS.md` (referenced by the plan), verbatim:

````markdown
# Deviations from the pre-specified analysis plan

1. **23 September 2026, before results.** P1 gene set added: `cheon2014_collagen` (10 genes; Cheon et al. 2014, Clin Cancer Res 20(3):711-723, doi 10.1158/1078-0432.CCR-13-1256). The Cheon signature named in plan v1.0 was not previously in the project gene-set JSON.
2. **23 September 2026, before results.** P2 gene set changed from `resistance_up_core` to `resistance_up_intrinsic` (14 genes: `resistance_up_core` without FOS, JUN, NR4A1, COL11A1, LOX, FABP4, IL6, STAT3, STAB1), to avoid overlap with the stress positive control (`stress_IEG_state_GSE165897`) and with stromal genes.
3. **23 September 2026, before results.** Gene sets with fewer than 5 genes (`matrix_CAF`, `sensitivity_up`, `partial_EMT_Xu2022`) moved from set-level testing to gene-level reporting. The 5-gene minimum is unchanged.
4. **23 September 2026, before results.** liana-py (rank_aggregate by sample) chosen as the cell-cell communication tool for secondary analysis 7.
````

## 3. gse165897_labels.csv (verbatim)

- Path: `metadata/gse165897_labels.csv` — 1 line(s), 51 bytes
- Last commit: `6a55e73 2026-09-23 12:38:54 +0100 Pre-specified analysis plan v1.0`

```
patient_id,pfi_days,crs,pfi_group,crs_group,source
```

**The file contains only the header row. No patient has a PFI or CRS label.**

Verification against Zhang et al. 2022 Table 1: **none found.**
- `git log --all -i --grep=verif --grep=Zhang --grep='Table 1'`: 0 matching commits.
- The only text mentioning Table 1 is the plan itself ("Values are taken from Zhang et al. (2022), Table 1"), which describes intent, not a completed check. DEVIATIONS.md does not mention the labels.

## 4. Stage 2 code check

**No script in the repo runs a resistant-vs-sensitive DE.** A search of all `.R`/`.py`/`.ipynb` files for `PFI`, `crs`, `resistant`, `fry`, `voom`, `edgeR`, `limma`, `propeller`, `liana` returns no analysis code (only `coxph` in `03_deconvolution.R`).

The only Stage 2 DE script is `02_dge.R` (last modified 2026-07-15 16:12:03), which tests **post-NACT vs treatment-naive** malignant cells. Its aggregation and test lines:

```r
 47  # ---- 1. pseudobulk aggregation: one profile per patient x treatment_phase ---
 48  MIN_CELLS <- 20L
 49  grp_key <- paste(sce_mal$patient_id, sce_mal$treatment_phase, sep = "||")
 50  counts  <- get_counts(sce_mal)
 51  
 52  n_per_grp <- table(grp_key)
 53  keep_grp  <- names(n_per_grp)[n_per_grp >= MIN_CELLS]
 54  dropped   <- setdiff(names(n_per_grp), keep_grp)
 55  if (length(dropped) > 0)
 56    log_msg("dropping", length(dropped), "patient x phase group(s) with <", MIN_CELLS,
 57            "malignant cells:", paste(dropped, collapse = "; "))
 58  
 59  pb <- sapply(keep_grp, function(g) rowSums(counts[, grp_key == g, drop = FALSE]))
 60  rownames(pb) <- rownames(sce_mal)
   ...
 81  # ---- 2. CPM normalise + log2, then paired t-test (vectorised) ---------------
 82  cpm <- function(m) sweep(m, 2, pmax(colSums(m), 1), "/") * 1e6
 83  log2cpm_post  <- log2(cpm(post)  + 1)
 84  log2cpm_naive <- log2(cpm(naive) + 1)
 85  
 86  d    <- log2cpm_post - log2cpm_naive           # per-patient (post - naive), paired
 87  m_d  <- rowMeans(d)
 88  se   <- rowSds(d) / sqrt(n)
 89  tstat <- m_d / se
 90  df <- n - 1
 91  pval <- 2 * pt(-abs(tstat), df)
 92  valid <- is.finite(pval) & se > 0
 93  pval[!valid] <- 1
 94  fdr <- p.adjust(pval, "BH")
```

- **Pseudobulk by patient × treatment_phase** (raw counts summed with `rowSums`, groups < 20 cells dropped), **not** a per-cell test.
- Test: CPM + log2, then a **paired t-test written in base R** (`pt()`, `matrixStats::rowSds`), BH FDR via `p.adjust`. **No DE package is used** (no edgeR, limma/voom or limma::fry, which the plan specifies).
- Cell type: author EOC annotation (`celltype == "malignant"`); CAFs are not tested.

## 5. Inventory: results/, figs/, data/

Labels: **REAL** = produced from GSE165897, GSE180661 or TCGA-OV. **SIM** = synthetic/simulation output. **REAL cells, SIM layout** = real GSE165897 cells placed on a simulated slide with a planted architecture.

| File | Size | Modified | Origin | Produced by |
|---|---:|---|---|---|
| `data/sce_clustered.rds` | 355MB | 2026-08-05 16:19 | REAL | 01_qc_cluster.R |
| `data/sce.rds` | 348MB | 2026-07-21 15:06 | REAL | 00_load_real_data.R (GSE165897) |
| `data/sce_site_comparison.rds` | 137GB | 2026-08-07 14:37 | REAL | 00_load_real_data_site.R (GSE180661); predates commit 70a0ecf |
| `figs/01_cnv_score_by_celltype.png` | 36KB | 2026-08-05 16:19 | REAL | 01_qc_cluster.R (log says profile: default) |
| `figs/02_volcano.png` | 97KB | 2026-08-06 14:54 | REAL | 02_dge.R (log says profile: default) |
| `figs/03b_benchmark.png` | 34KB | 2026-07-15 13:11 | SIM | initial commit 9fb71b7 (synthetic pipeline) |
| `figs/03_fraction_scatter.png` | 15KB | 2026-07-15 13:11 | SIM | initial commit 9fb71b7 (synthetic pipeline) |
| `figs/03_km_by_malignant_tertile.png` | 35KB | 2026-08-06 16:45 | REAL | 03_deconvolution.R (profile: real) |
| `figs/04_spatial_maps.png` | 41KB | 2026-08-06 16:49 | REAL cells, SIM layout | 04_spatial.R (mock slide) |
| `figs/tcga_celltype_kde.png` | 147KB | 2026-08-06 17:49 | REAL | **no script in repo** (proportion_survival_plot.R is 0 bytes) |
| `figs/tcga_celltype_pie.png` | 140KB | 2026-08-06 17:33 | REAL | pie_chart_plot.R |
| `figs/tcga_celltype_score_distributions.png` | 413KB | 2026-08-06 17:40 | REAL | distribution_plot.R |
| `figs/tcga_celltype_vs_survival_scatter.png` | 530KB | 2026-08-06 17:45 | REAL | **no script in repo** (proportion_survival_plot.R is 0 bytes) |
| `figs/tcga_ov_heatmap_GSVA.png` | 232KB | 2026-08-19 15:23 | REAL | 05_signature_vs_celltype.R |
| `figs/tcga_ov_heatmap_ssGSEA.png` | 233KB | 2026-08-19 15:23 | REAL | 05_signature_vs_celltype.R |
| `results/01_confusion_malignant.csv` | 135B | 2026-08-05 16:19 | REAL | 01_qc_cluster.R (log says profile: default) |
| `results/01_metrics.csv` | 174B | 2026-08-05 16:19 | REAL | 01_qc_cluster.R (log says profile: default) |
| `results/01_run_log.txt` | 35KB | 2026-08-05 16:19 | REAL | 01_qc_cluster.R (log says profile: default) |
| `results/02_de_results.csv` | 1.8MB | 2026-08-06 14:54 | REAL | 02_dge.R (log says profile: default) |
| `results/02_plausibility_check.csv` | 803B | 2026-08-06 14:54 | REAL | 02_dge.R (log says profile: default) |
| `results/02_README.md` | 931B | 2026-08-06 14:54 | REAL | 02_dge.R (log says profile: default) |
| `results/02_run_log.txt` | 968B | 2026-08-06 14:54 | REAL | 02_dge.R (log says profile: default) |
| `results/03_estimated_fractions.csv` | 30KB | 2026-08-06 16:45 | REAL | 03_deconvolution.R (profile: real) |
| `results/03_metrics.csv` | 163B | 2026-08-06 16:45 | REAL | 03_deconvolution.R (profile: real) |
| `results/03_run_log.txt` | 1.7KB | 2026-08-06 16:45 | REAL | 03_deconvolution.R (profile: real) |
| `results/04_metrics.csv` | 178B | 2026-08-06 16:49 | REAL cells, SIM layout | 04_spatial.R (mock slide) |
| `results/04_run_log.txt` | 1.3KB | 2026-08-06 16:49 | REAL cells, SIM layout | 04_spatial.R (mock slide) |
| `results/04_spot_estimates.csv` | 71KB | 2026-08-06 16:49 | REAL cells, SIM layout | 04_spatial.R (mock slide) |
| `results/05_run_log.txt` | 2.0KB | 2026-08-19 15:23 | REAL | 05_signature_vs_celltype.R |
| `results/infercnv/01_incoming_data.infercnv_obj` | 259MB | 2026-08-05 15:22 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/02_reduced_by_cutoff.infercnv_obj` | 166MB | 2026-08-05 15:23 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/03_normalized_by_depth.infercnv_obj` | 195MB | 2026-08-05 15:25 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/04_logtransformed.infercnv_obj` | 199MB | 2026-08-05 15:26 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/08_remove_ref_avg_from_obs_logFC.infercnv_obj` | 2.0GB | 2026-08-05 15:29 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/09_apply_max_centered_expr_threshold.infercnv_obj` | 2.0GB | 2026-08-05 15:31 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/10_smoothed_by_chr.infercnv_obj` | 2.1GB | 2026-08-05 15:46 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/11_recentered_cells_by_chr.infercnv_obj` | 2.1GB | 2026-08-05 15:49 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/12_remove_ref_avg_from_obs_adjust.infercnv_obj` | 1.7GB | 2026-08-05 15:52 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/14_invert_log_transform.infercnv_obj` | 1.6GB | 2026-08-05 15:55 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/15_tumor_subclusters.leiden.infercnv_obj` | 1.6GB | 2026-08-05 16:03 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/22_denoise.leiden.NF_NA.SD_1.5.NL_FALSE.infercnv_obj` | 422MB | 2026-08-05 16:16 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/infercnv_annotations.tsv` | 1.6MB | 2026-08-05 15:15 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/infercnv_counts.tsv` | 1.7GB | 2026-08-05 15:15 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/infercnv_gene_order.tsv` | 649KB | 2026-08-05 15:15 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/infercnv_subclusters.heatmap_thresholds.txt` | 280B | 2026-08-05 16:08 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/infercnv_subclusters.observation_groupings.txt` | 514KB | 2026-08-05 16:08 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/infercnv_subclusters.png` | 1.3MB | 2026-08-05 16:12 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/preliminary.infercnv_obj` | 1.6GB | 2026-08-05 16:15 | REAL | inferCNV via 01_qc_cluster.R |
| `results/infercnv/run.final.infercnv_obj` | 422MB | 2026-08-05 16:17 | REAL | inferCNV via 01_qc_cluster.R |
| `results/peek_02_top_de.csv` | 832B | 2026-08-06 14:54 | REAL | 02_dge.R (log says profile: default) |
| `results/tcga_celltype_scores_long.csv` | 47KB | 2026-08-06 17:40 | REAL | distribution_plot.R |
| `results/tcga_celltype_survival_scatter_data.csv` | 71KB | 2026-08-06 17:45 | REAL | **no script in repo** (proportion_survival_plot.R is 0 bytes) |
| `results/tcga_mean_celltype_proportions.csv` | 169B | 2026-08-06 17:33 | REAL | pie_chart_plot.R |
| `results/tcga_ov_corr_rho_GSVA.csv` | 1.3KB | 2026-08-19 15:23 | REAL | 05_signature_vs_celltype.R |
| `results/tcga_ov_corr_rho_ssGSEA.csv` | 1.3KB | 2026-08-19 15:23 | REAL | 05_signature_vs_celltype.R |
| `results/tcga_ov_scores_GSVA.csv` | 94KB | 2026-08-19 15:23 | REAL | 05_signature_vs_celltype.R |
| `results/tcga_ov_scores_ssGSEA.csv` | 90KB | 2026-08-19 15:23 | REAL | 05_signature_vs_celltype.R |
| `results/tcga_ov_sessionInfo_stage05.txt` | 2.9KB | 2026-08-19 15:23 | REAL | 05_signature_vs_celltype.R |
| `results/tcga_ov_signature_celltype_correlations.csv` | 11KB | 2026-08-19 15:23 | REAL | 05_signature_vs_celltype.R |
| `results/tcga_ov_signature_scores.rds` | 74KB | 2026-08-19 15:23 | REAL | 05_signature_vs_celltype.R |
| `results/validation_log.tsv` | 1013B | 2026-08-06 16:49 | MIXED | appended by every stage; includes 2026-07-21 rows of unclear profile |

`truth/` (not requested, listed for completeness) holds the synthetic planted truth: **SIM**.

## 6. Real-data results

Analyses in the plan that have **no output on disk**: every PFI/CRS contrast, P1, P2, the positive control, limma::fry gene-set tests, speckle::propeller, liana-py, and the plan's TCGA-OV validation. All results on disk date from 2026-08-05 to 2026-08-19, **before** the plan was committed (2026-09-23).

### Aim 1 — inferCNV on GSE165897

Run: 2026-08-05 (`results/01_run_log.txt`). CNV method: **inferCNV (PRIMARY)**, not the fallback. The log header says `profile: default`, though the input is the real `data/sce.rds` (42,126 cells × 21,544 genes).

**Cells per patient before and after QC.** "Raw" counts come from `../GSE165897_cellInfo_HGSOC.tsv.gz`. "In SCE" counts come from `colData(data/sce_clustered.rds)`, after `00_load_real_data.R` dropped cell types outside the 6-type vocabulary. QC removed 0 cells (log: "QC discards 0 of 42126 cells"; inferCNV: "removing 0 % of cells").

| patient_id | raw cells (GSE165897) | cells in SCE (6 types) | cells after QC |
|---|---:|---:|---:|
| EOC1005 | 4359 | 3824 | 3824 |
| EOC136 | 6103 | 5742 | 5742 |
| EOC153 | 3781 | 2812 | 2812 |
| EOC227 | 4193 | 2802 | 2802 |
| EOC3 | 5645 | 4180 | 4180 |
| EOC349 | 2954 | 2336 | 2336 |
| EOC372 | 5382 | 4898 | 4898 |
| EOC443 | 6585 | 5353 | 5353 |
| EOC540 | 3964 | 3357 | 3357 |
| EOC733 | 5819 | 4518 | 4518 |
| EOC87 | 3001 | 2304 | 2304 |
| **TOTAL** | **51786** | **42126** | **42126** |

**inferCNV-called malignant cells per patient and timepoint** (`cnv_call == "malignant"` in `data/sce_clustered.rds`; total 5,807):

| patient_id | post-NACT | treatment-naive |
|---|---:|---:|
| EOC1005 | 330 | 432 |
| EOC136 | 2064 | 175 |
| EOC153 | 57 | 119 |
| EOC227 | 180 | 64 |
| EOC3 | 61 | 19 |
| EOC349 | 56 | 20 |
| EOC372 | 100 | 162 |
| EOC443 | 78 | 771 |
| EOC540 | 8 | 169 |
| EOC733 | 133 | 629 |
| EOC87 | 93 | 87 |

**Agreement with the author EOC annotation** (`results/01_confusion_malignant.csv`, `results/01_metrics.csv`):

| | called diploid | called malignant |
|---|---:|---:|
| author non-EOC (reference) | 33142 | 178 |
| author EOC | 3177 | 5629 |

- **Reference false-positive rate: 0.00534** (178 / 33,320 reference cells; `ref_false_call_rate` = 0.0053421368547419). `validation_log.tsv` prints the same value rounded to `ref_fpr=0.01`.
- Call accuracy 0.920; ARI(call vs author label) 0.643; ARI(cluster vs cell type) 0.285.
- Gate: **[FAIL]**, ARI 0.643 < `ari_min` 0.70. `config.yml` says the `real` profile deliberately does not override `ari_min` until a justified decision is made.
- Reference cells called malignant, by type (sce_clustered.rds): fibroblast 11, macrophage 17, Tcell 150, endothelial 0.

### Aim 2 — pseudobulk DE

**Resistant vs sensitive (PFI): NOT RUN.** `metadata/gse165897_labels.csv` has no rows.
**Resistant vs sensitive (CRS): NOT RUN.** Same reason.

| Item | PFI | CRS |
|---|---|---|
| Patients per group | NOT RUN (no labels) | NOT RUN (no labels) |
| Cells per pseudobulk | NOT RUN | NOT RUN |
| Genes at FDR < 0.05 | NOT RUN | NOT RUN |
| Top 20 genes | NOT RUN | NOT RUN |
| P1 (CAF, `cheon2014_collagen`, limma::fry): raw p / Holm p | NOT RUN | NOT RUN |
| P2 (malignant, `resistance_up_intrinsic`, limma::fry): raw p / Holm p | NOT RUN | NOT RUN |
| Positive control (`stress_IEG_state_GSE165897`, PN vs TN) | NOT RUN | — |

**What exists instead (not a plan analysis):** `02_dge.R`, post-NACT vs treatment-naive in author-annotated malignant cells. Run 2026-08-06 14:54. Log header says `profile: default`.

- Paired n = **10 patients** (`02_run_log.txt`): EOC540 post-NACT was dropped with < 20 malignant cells.
- Cells per pseudobulk (author EOC cells per patient × phase from `colData(data/sce_clustered.rds)`; 02 uses the same `celltype` from `data/sce.rds`):

| patient_id | post-NACT | treatment-naive |
|---|---:|---:|
| EOC1005 | 418 | 487 |
| EOC136 | 2427 | 205 |
| EOC153 | 116 | 330 |
| EOC227 | 174 | 57 |
| EOC3 | 223 | 203 |
| EOC349 | 66 | 30 |
| EOC372 | 131 | 281 |
| EOC443 | 103 | 1086 |
| EOC540 | 18 (dropped, < 20) | 202 |
| EOC733 | 372 | 1550 |
| EOC87 | 191 | 136 |

- Genes tested: 21,544. **Genes at FDR < 0.05: 0** (log line, and 0 rows with `fdr < 0.05` in `results/02_de_results.csv`). Smallest FDR = 0.121.
- Top 20 genes by FDR (`results/02_de_results.csv`, file order; log2FC = post-NACT − treatment-naive, values rounded):

| # | gene | log2FC | p | FDR |
|---:|---|---:|---:|---:|
| 1 | S100A1 | 1.269 | 1.31e-05 | 0.121 |
| 2 | CD14 | 1.369 | 2.43e-05 | 0.121 |
| 3 | R3HCC1 | -0.451 | 2.42e-05 | 0.121 |
| 4 | KLF6 | 1.043 | 1.78e-05 | 0.121 |
| 5 | AXIN1 | -0.649 | 2.81e-05 | 0.121 |
| 6 | EIF4A2 | 0.529 | 4.39e-05 | 0.158 |
| 7 | AUNIP | -0.780 | 1.68e-04 | 0.321 |
| 8 | KLHL23 | -0.767 | 3.65e-04 | 0.321 |
| 9 | EXOG | -0.705 | 1.99e-04 | 0.321 |
| 10 | IGFBP7 | 1.748 | 3.95e-04 | 0.321 |
| 11 | ANKRD37 | 0.890 | 3.62e-04 | 0.321 |
| 12 | APBB3 | 0.819 | 1.83e-04 | 0.321 |
| 13 | HLA-DPB1 | 1.661 | 2.26e-04 | 0.321 |
| 14 | PNRC1 | 1.167 | 1.96e-04 | 0.321 |
| 15 | WASF1 | -0.749 | 3.74e-04 | 0.321 |
| 16 | RP9 | -0.685 | 4.17e-04 | 0.321 |
| 17 | PDK4 | 1.509 | 3.08e-04 | 0.321 |
| 18 | RAB11FIP1 | 0.620 | 3.11e-04 | 0.321 |
| 19 | DUSP5 | 1.417 | 4.13e-04 | 0.321 |
| 20 | EHF | 1.047 | 2.77e-04 | 0.321 |

- Literature plausibility check (`results/02_plausibility_check.csv`): 10/10 genes present, 0 significant (lowest FDR: BRCA2, 0.582).
- No gene-set test was run in this script.

### Aim 3 — module scoring and liana-py (GSE165897)

| Item | Status |
|---|---|
| Module/gene-set score summaries by resistant/sensitive group (single cell) | **NOT RUN** (no labels, no scoring script for GSE165897) |
| liana-py rank_aggregate by sample, TN samples | **NOT RUN** (env installed, no script and no output on disk) |
| Top 20 malignant/fibroblast interactions | **NOT RUN** |
| Patient-level interaction Wilcoxon tests | **NOT RUN** |

The repo's own "Aim 3" is a different analysis: TCGA-OV deconvolution plus survival, listed under Aim 4 below as related output.

### Aim 4 — TCGA-OV validation

**The plan's validation (Section 7: GSVA score of the single-cell-derived signature, logistic model of platinum status plus Cox model of PFS, adjusted for age, stage and residual disease) is NOT RUN.**

| Item | Planned validation |
|---|---|
| Engine | NOT RUN |
| n | NOT RUN |
| Covariates | NOT RUN |
| OR / HR with 95% CI, p | NOT RUN |
| Run date | NOT RUN |
| Label-source choice recorded in DEVIATIONS.md | Not recorded |

**Related TCGA-OV outputs that do exist (not the plan analysis):**

*Stage 03, deconvolution and overall survival* (`results/03_metrics.csv`, `results/03_run_log.txt`, validation log timestamp 2026-08-06 16:45:01, profile `real`):
- Deconvolution engine: **NNLS [FALLBACK]** on 250 marker genes. BayesPrism was not used.
- Survival engine: `survival::coxph`, outcome **overall survival** (`os_time_days` / `os_status`), not PFS or platinum status.
- n = **425** samples, 265 events (38% censored). These are 425 samples from **424 unique patients**: TCGA-23-1023 contributes both `-01A` and `-01R` (`results/03_estimated_fractions.csv`).
- Covariates: estimated malignant fraction + estimated fibroblast fraction only. **No clinical covariates.**
- Malignant fraction **HR = 2.16, p = 0.0728** (`cox_malignant_HR` 2.16089254713295, `cox_malignant_p` 0.0728459435051877). CAF HR = 1.88 (log only). **95% CI: not written to any file.**
- KM by malignant tertile: log-rank p = 0.502.
- Gate: [PASS] "Cox model converged". Two earlier runs that day (15:32 and 16:00) logged [FAIL] with HR = NA.

*Stage 05, GSVA/ssGSEA signature scores vs deconvolved fractions* (`results/05_run_log.txt`, 2026-08-19): GSVA 1.50.0 (Poisson kernel on counts) and ssGSEA. 11 gene sets were scored on 425 aligned samples. Outputs are Spearman correlations (`results/tcga_ov_corr_rho_{GSVA,ssGSEA}.csv`); no outcome model was fitted. The gene sets are the JSON literature sets, not a single-cell-derived signature, and they predate the v1.1 sets `cheon2014_collagen` and `resistance_up_intrinsic`.

*Stage 04, spatial colocalisation (repo "Aim 4")*: a mock slide of real GSE165897 cells on a 30×30 grid with a **planted** core/margin layout. Engine: per-spot NNLS [FALLBACK]. CAF core mean 0.347 vs margin 0.082, permutation p = 0.001; [PASS]. This shows the planted layout was recovered. It is not an observation about real tissue.

## 7. PASS / FAIL / FALLBACK lines and run-log tails

### 7a. results/validation_log.tsv (verbatim)
```
time	stage	label	status	detail
2026-07-21 15:40:11.596008	01_qc_cluster	malignant calls recover true labels	FAIL	ARI=0.153 (>=0.70) acc=0.815 ref_fpr=0.01
2026-07-21 15:54:29.196845	01_qc_cluster	malignant calls recover true labels	FAIL	ARI=0.153 (>=0.70) acc=0.815 ref_fpr=0.01
2026-08-05 16:19:18.502926	01_qc_cluster	malignant calls recover true labels	FAIL	ARI=0.643 (>=0.70) acc=0.920 ref_fpr=0.01
2026-08-06 15:32:40.532807	03_deconvolution	Cox model converged on real TCGA-OV survival data	FAIL	malignant HR=NA p=NA | KM log-rank p=1
2026-08-06 16:00:21.335939	03_deconvolution	Cox model converged on real TCGA-OV survival data	FAIL	malignant HR=NA p=NA | KM log-rank p=1
2026-08-06 16:45:01.016419	03_deconvolution	Cox model converged on real TCGA-OV survival data	PASS	malignant HR=2.16 p=0.0728 | KM log-rank p=0.502
2026-08-06 16:49:56.325072	04_spatial	recovered malignant-CAF colocalisation matches layout	PASS	CAF concentrated in malignant-dense core vs margin: perm p=0.001 (core=0.35 margin=0.08)
```

### 7b. Every PASS/FAIL/FALLBACK line in any log
```
results/01_run_log.txt:497:[FAIL] 01_qc_cluster :: malignant calls recover true labels (ARI=0.643 (>=0.70) acc=0.920 ref_fpr=0.01)
results/03_run_log.txt:20:  - [FALLBACK] deconvolution: NNLS on cell-type marker signatures 
results/03_run_log.txt:21:  - deconvolution method: NNLS [FALLBACK] 
results/03_run_log.txt:31:[PASS] 03_deconvolution :: Cox model converged on real TCGA-OV survival data (malignant HR=2.16 p=0.0728 | KM log-rank p=0.502)
results/04_run_log.txt:15:  - [FALLBACK] spatial deconvolution: per-spot NNLS 
results/04_run_log.txt:16:  - spatial deconvolution method: NNLS [FALLBACK] 
results/04_run_log.txt:24:[PASS] 04_spatial :: recovered malignant-CAF colocalisation matches layout (CAF concentrated in malignant-dense core vs margin: perm p=0.001 (core=0.35 margin=0.08))
verify_site_log.txt:146:spot-check on saved object: PASSED 
install_log.txt:654:  optional infercnv   : absent -> [FALLBACK]
install_log.txt:655:  optional BayesPrism : absent -> [FALLBACK]
install_log.txt:656:  optional spacexr    : absent -> [FALLBACK]
```

### 7c. Last 30 lines: `results/01_run_log.txt` (2026-08-05 16:19:18)
```
INFO [2026-08-05 16:08:57] Colors for breaks:  #00008B,#24249B,#4848AB,#6D6DBC,#9191CC,#B6B6DD,#DADAEE,#FFFFFF,#EEDADA,#DDB6B6,#CC9191,#BC6D6D,#AB4848,#9B2424,#8B0000
INFO [2026-08-05 16:08:57] Quantiles of plotted data range: 0.898859657273169,0.979334336148186,1,1.02292202234853,1.10465459673775
INFO [2026-08-05 16:09:06] plot_cnv_references:Start
INFO [2026-08-05 16:09:06] Reference data size: Cells= 33320 Genes= 6543
INFO [2026-08-05 16:11:19] plot_cnv_references:Number reference groups= 288
INFO [2026-08-05 16:11:22] plot_cnv_references:Plotting heatmap.
INFO [2026-08-05 16:12:19] Colors for breaks:  #00008B,#24249B,#4848AB,#6D6DBC,#9191CC,#B6B6DD,#DADAEE,#FFFFFF,#EEDADA,#DDB6B6,#CC9191,#BC6D6D,#AB4848,#9B2424,#8B0000
INFO [2026-08-05 16:12:19] Quantiles of plotted data range: 0.898859657273169,0.987094088175052,1,1.01056385378286,1.10465459673775
INFO [2026-08-05 16:15:10] 

	STEP 22: Denoising

INFO [2026-08-05 16:15:10] ::process_data:Remove noise, noise threshold defined via ref mean sd_amplifier:  1.5
INFO [2026-08-05 16:15:10] denoising using mean(normal) +- sd_amplifier * sd(normal) per gene per cell across all data
INFO [2026-08-05 16:15:19] :: **** clear_noise_via_ref_quantiles **** : removing noise between bounds:  0.962784069622794 - 1.03840830047824
Warning messages:
1: In asMethod(object) :
  sparse->dense coercion: allocating vector of size 6.8 GiB
2: package ‘future’ was built under R version 4.3.3 
  - CNV method: inferCNV | reference false-call rate = 0.5% 

  Confusion matrix (rows=truth, cols=call):
           predicted
truth       diploid malignant
  diploid     33142       178
  malignant    3177      5629
  - call accuracy = 0.920 | ARI(call vs truth) = 0.643 | ARI(cluster vs celltype) = 0.285 
  - figure -> figs/01_cnv_score_by_celltype.png 

[FAIL] 01_qc_cluster :: malignant calls recover true labels (ARI=0.643 (>=0.70) acc=0.920 ref_fpr=0.01)
```

### 7c. Last 30 lines: `results/02_run_log.txt` (2026-08-06 14:54:46)
```
There were 11 warnings (use warnings() to see them)
Warning message:
package ‘yaml’ was built under R version 4.3.3 
==========================================================================
  02_dge.R — Aim 2 (chemoresistance/stress-response DE)   [profile: default | seed: 1000 | threads: 1]
==========================================================================
Warning message:
In asMethod(object) :
  sparse->dense coercion: allocating vector of size 1.4 GiB
  - dropping 1 patient x phase group(s) with < 20 malignant cells: EOC540||post-NACT 
  - 10 of 11 patients have both phases with >= 20 malignant cells (paired n) 
  - peek -> results/peek_02_top_de.csv 
  - paired pseudobulk DE: n=10 patient pairs, 0 genes significant at FDR<0.05 
  - literature-plausibility check: 10/10 genes found in reference; 0 significant at FDR<0.05 
  - figure -> figs/02_volcano.png 
  - Aim 2 complete -- no PASS/FAIL gate on real data (see results/02_README.md) 
```

### 7c. Last 30 lines: `results/03_run_log.txt` (2026-08-06 16:45:01)
```
There were 12 warnings (use warnings() to see them)
Warning message:
package ‘yaml’ was built under R version 4.3.3 
==========================================================================
  03_deconvolution.R — Aim 3 (deconvolution + survival)   [profile: real | seed: 1000 | threads: 1]
==========================================================================
            used   (Mb) gc trigger   (Mb)  max used   (Mb)
Ncells   6312044  337.1   10241934  547.0   8145122  435.0
Vcells 324755339 2477.7  450596366 3437.8 356611808 2720.8

  harmonize_bulk_to_symbols: 24432 of 60660 TCGA genes dropped (no ENSEMBL->SYMBOL mapping)
  - kept 427 of 434 TCGA-OV samples (sample type -01/primary tumor only) 
  - loaded 425 TCGA-OV bulk samples with usable survival data 
Warning message:
In asMethod(object) :
  sparse->dense coercion: allocating vector of size 6.8 GiB
  - reference cell types: endothelial, fibroblast, macrophage, malignant, Tcell 
  - selected 250 marker genes for deconvolution 
  - [FALLBACK] deconvolution: NNLS on cell-type marker signatures 
  - deconvolution method: NNLS [FALLBACK] 
  - events: 265 / 425 (38% censored) 
  - ADIPOCYTE BLIND SPOT: reference has zero adipocyte cells (GSE165897); detected = NO (expected -- no reference to detect from) 
  - Cox (estimated fractions, real survival): malignant HR=2.16 p=0.0728 | CAF HR=1.88 
  - KM by estimated-malignant tertile: log-rank p = 0.502 
Fontconfig warning: using without calling FcInit()
null device 
          1 
  - figure -> figs/03_km_by_malignant_tertile.png 

[PASS] 03_deconvolution :: Cox model converged on real TCGA-OV survival data (malignant HR=2.16 p=0.0728 | KM log-rank p=0.502)
```

### 7c. Last 30 lines: `results/04_run_log.txt` (2026-08-06 16:49:56)
```
nohup: ignoring input
There were 11 warnings (use warnings() to see them)
Warning message:
package ‘yaml’ was built under R version 4.3.3 
==========================================================================
  04_spatial.R — Aim 4 (spatial colocalisation)   [profile: real | seed: 1000 | threads: 1]
==========================================================================
Warning message:
In asMethod(object) :
  sparse->dense coercion: allocating vector of size 6.8 GiB
  - grid 30x30 = 900 spots | core=256 margin=392 background=252 
Warning message:
In asMethod(object) :
  sparse->dense coercion: allocating vector of size 6.8 GiB
  - [FALLBACK] spatial deconvolution: per-spot NNLS 
  - spatial deconvolution method: NNLS [FALLBACK] 
  - inferred CAF: core mean=0.35 vs margin mean=0.08 | permutation p=0.001 
  - (context) malignant<->CAF fraction Pearson: tumour-bed r=-0.29, global r=-0.32 (negative is expected, see code note) 
Fontconfig warning: using without calling FcInit()
null device 
          1 
  - figure -> figs/04_spatial_maps.png  

[PASS] 04_spatial :: recovered malignant-CAF colocalisation matches layout (CAF concentrated in malignant-dense core vs margin: perm p=0.001 (core=0.35 margin=0.08))
```

### 7c. Last 30 lines: `results/05_run_log.txt` (2026-08-19 15:23:59)
```
8: package ‘IRanges’ was built under R version 4.3.3 
9: package ‘GenomeInfoDb’ was built under R version 4.3.2 
10: package ‘Biobase’ was built under R version 4.3.3 

  harmonize_bulk_to_symbols: 24432 of 60660 TCGA genes dropped (no ENSEMBL->SYMBOL mapping)
[15:00:53]  bulk: 36152 genes x 434 samples 
[15:00:53]  proportions: 425 samples x 5 cell types (endothelial, fibroblast, macrophage, malignant, Tcell) 
[15:00:53]  aligned on 425 shared samples 
[15:00:53]  dropped 1 set(s) below min size 2: prognostic_only 
[15:00:53]  scoring 11 gene sets 
[15:00:53]  running GSVA ... 
Warning messages:
1: In .filterFeatures_newAPI(dataMatrix) :
  1272 rows with constant values throughout the samples.
2: In .filterFeatures_newAPI(dataMatrix) :
  Rows with constant values are discarded.
[15:23:41]  running ssGSEA ... 
[1] "Calculating ranks..."
[1] "Calculating absolute values from ranks..."
[1] "Normalizing..."
Warning message:
In .filterFeatures_newAPI(dataMatrix, dropConstantRows = FALSE) :
  1272 rows with constant values throughout the samples.
Fontconfig warning: using without calling FcInit()
[15:23:59]  done. 
[15:23:59]  outputs: 
[15:23:59]    results/tcga_ov_signature_celltype_correlations.csv (tidy long) 
[15:23:59]    results/tcga_ov_corr_rho_{GSVA,ssGSEA}.csv 
[15:23:59]    results/tcga_ov_scores_{GSVA,ssGSEA}.csv 
[15:23:59]    figs/tcga_ov_heatmap_{GSVA,ssGSEA}.png 
```

### 7c. Last 30 lines: `site_load_log.txt` (2026-09-23 13:34:44)
```
   singlet, live, CD45+   48059   68121 202401
   singlet, live, U           0    3462  16512 
  - querying Ensembl (biomaRt) for gene coordinates ... 
Error: Could not reach Ensembl via biomaRt: invalid 'file' argument
Check the server's network/proxy/firewall settings, or try biomaRt::useEnsembl(..., mirror = "useast"/"asia") as a fallback. This must be resolved before the gene-coordinate lookup can proceed.
Execution halted
Command exited with non-zero status 1
	Command being timed: "Rscript 00_load_real_data_site.R"
	User time (seconds): 3629.17
	System time (seconds): 3466.25
	Percent of CPU this job got: 100%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:58:04
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 70091652
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 106
	Minor (reclaiming a frame) page faults: 1725448282
	Voluntary context switches: 9725
	Involuntary context switches: 39188
	Swaps: 0
	File system inputs: 63714984
	File system outputs: 519320
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 1
```

### 7c. Last 30 lines: `verify_site_log.txt` (2026-08-07 15:04:13)
```
  SPECTRUM-OV-112    1948    3480  9058
  SPECTRUM-OV-115       0    6317 16921
  SPECTRUM-OV-116    4989    6098 10382
  SPECTRUM-OV-118       0       0 22624

patients total: 41 | with ascites: 25 | with omentum: 35 | PAIRED (both): 21 
paired IDs: SPECTRUM-OV-007, SPECTRUM-OV-009, SPECTRUM-OV-024, SPECTRUM-OV-037, SPECTRUM-OV-041, SPECTRUM-OV-042, SPECTRUM-OV-050, SPECTRUM-OV-051, SPECTRUM-OV-054, SPECTRUM-OV-065, SPECTRUM-OV-068, SPECTRUM-OV-070, SPECTRUM-OV-075, SPECTRUM-OV-077, SPECTRUM-OV-080, SPECTRUM-OV-081, SPECTRUM-OV-083, SPECTRUM-OV-105, SPECTRUM-OV-107, SPECTRUM-OV-112, SPECTRUM-OV-116 

-- gene/cell counts --
genes after coordinate mapping: 21650 
cells: 700227 
rowData cols: gene_id, chr, start, stop, arm 
any missing chr/start/stop/arm? FALSE 

-- coordinate spot-check (re-verified on the SAVED object) --
   gene expected   got   ok
1  TP53    chr17 chr17 TRUE
2 BRCA1    chr17 chr17 TRUE
3   MYC     chr8  chr8 TRUE
4 EPCAM     chr2  chr2 TRUE
5 PTPRC     chr1  chr1 TRUE
spot-check on saved object: PASSED 

-- assay sanity --
counts integer-valued: TRUE 
cells with zero total counts: 0 
logcounts range: 0.3128746 10.51706 
any non-finite in logcounts (sampled): FALSE 

done in 11.8 min
```

### 7c. Last 30 lines: `raw_data/GSE180661/download_log.txt` (2026-08-07 12:24:36)
```
+ curl -sSL --retry 5 --retry-delay 10 -C - -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180661/suppl/GSE180661_GEO_cells.tsv.gz
+ echo 'CELLS_DONE rc=0'
CELLS_DONE rc=0
+ curl -sSL --retry 5 --retry-delay 10 -C - -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180661/suppl/GSE180661_matrix.h5
+ echo 'MATRIX_DONE rc=0'
MATRIX_DONE rc=0
+ ls -l
total 31834000
-rw-rw-r-- 1 syeed syeed         365 Aug  7 12:24 download_log.txt
-rw-rw-r-- 1 syeed syeed    42728783 Aug  7 11:46 GSE180661_GEO_cells.tsv.gz
-rw-rw-r-- 1 syeed syeed 32555276423 Aug  7 12:24 GSE180661_matrix.h5
```

### 7c. Last 30 lines: `install_log.txt` (2026-07-21 14:35:43)
```
Updating HTML index of packages in '.Library'
Making 'packages.html' ... done
'getOption("repos")' replaces Bioconductor standard repositories, see
'help("repositories", package = "BiocManager")' for details.
Replacement repositories:
    CRAN: https://cloud.r-project.org
Bioconductor version 3.18 (BiocManager 1.30.27), R 4.3.1 (2023-06-16)
Installing package(s) 'igraph'
Warning: dependency ‘Matrix’ is not available
trying URL 'https://cloud.r-project.org/src/contrib/igraph_2.3.3.tar.gz'
Content type 'application/x-gzip' length 5215250 bytes (5.0 MB)
==================================================
downloaded 5.0 MB

ERROR: dependency ‘Matrix’ is not available for package ‘igraph’
* removing ‘/userhome/conda/envs/hgsoc_pipeline/lib/R/library/igraph’

The downloaded source packages are in
	‘/tmp/RtmpVRcRDp/downloaded_packages’
Updating HTML index of packages in '.Library'
Making 'packages.html' ... done
There were 50 or more warnings (use warnings() to see the first 50)
>> Skipping optional tools (set INSTALL_OPTIONAL=1 to attempt them). Fallbacks will be used.

CORE packages installed: FALSE 
  optional infercnv   : absent -> [FALLBACK]
  optional BayesPrism : absent -> [FALLBACK]
  optional spacexr    : absent -> [FALLBACK]
Error: Some core packages failed to install — see messages above.
Execution halted
```

## 8. Environment

### sessionInfo files on disk

| File | Next to | R | Platform | Key packages |
|---|---|---|---|---|
| `sessionInfo.txt` | . | R version 4.3.1 (2023-06-16 ucrt) | x86_64-w64-mingw32/x64 (64-bit) | data.table_1.17.0 edgeR_4.0.16 GenomeInfoDb_1.38.8 ggplot2_4.0.3 infercnv_1.18.1 jsonlite_2.0.0 limma_3.58.1 Matrix_1.6-5 matrixStats_1.5.0 scater_1.30.1 scuttle_1.12.0 Seurat_5.2.1 SingleCellExperiment_1.24.0 splatter_1.26.0 SummarizedExperiment_1.32.0 survival_3.5-5 yaml_2.3.10  |
| `results/tcga_ov_sessionInfo_stage05.txt` | results | R version 4.3.1 (2023-06-16) | x86_64-conda-linux-gnu (64-bit) | data.table_1.17.8 GenomeInfoDb_1.38.1 GSVA_1.50.0 jsonlite_2.0.0 Matrix_1.6-5 matrixStats_1.5.0 SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0  |

- Root `sessionInfo.txt` (modified 2026-07-15 13:11:27) is from **Windows 11, R 4.3.1 ucrt, time zone Asia/Dubai**, with `splatter` attached: it belongs to the synthetic pipeline, not to any real-data run on this server.
- `results/tcga_ov_sessionInfo_stage05.txt` is the only sessionInfo for a real-data run (Stage 05, Linux conda, R 4.3.1).
- **No sessionInfo exists for real-data stages 00, 01, 02, 03, 04.**
- `renv.lock`: R 4.3.1, Bioconductor 3.18.

### Current hgsoc_pipeline env
```
R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
SingleCellExperiment   1.24.0
infercnv               1.20.0
GSVA                   1.50.0
survival               3.8.3
edgeR                  4.0.16
limma                  3.58.1
scater                 1.30.1
scran                  1.30.0
matrixStats            1.5.0
BayesPrism             NOT INSTALLED
spacexr                NOT INSTALLED
speckle                NOT INSTALLED
```

### hgsoc_ccc env (`/userhome/syeed/envs/hgsoc_ccc`)
```
Python 3.13.2
liana      1.10.0
scanpy     1.12.4
anndata    0.13.4
mudata     0.4.1
pandas     2.3.3
numpy      2.5.3
```
Matches `env_hgsoc_ccc.yml` / `requirements_ccc.txt` (python 3.13.2, liana 1.10.0). The env exists but **no liana output exists anywhere on disk**.

## 9. Inconsistencies noticed

**Plan vs what was run**
1. **No result on disk was produced under the analysis plan.** The plan (v1.0/v1.1) was committed 2026-09-23; every results file is dated 2026-08-05 to 2026-08-19. There are no commits since 2026-09-23 13:29, so no analysis has followed the plan yet.
2. **`metadata/gse165897_labels.csv` is header-only.** No PFI or CRS values have been entered, and no commit or note records a check against Zhang et al. 2022 Table 1. Every primary and secondary contrast depends on this file.
3. **Stage 2 does not match the plan.** `02_dge.R` compares post-NACT with treatment-naive samples, not resistant with sensitive. It uses CPM plus a hand-written paired t-test, not `filterByExpr`/TMM/`voom`/`fry`. It tests malignant cells only, not CAFs, and runs no gene-set tests. Its pseudobulk is by patient × phase, not patient × timepoint × cell type.
4. **The aim numbering differs between your brief and the repo.** In the repo, Aim 3 = TCGA deconvolution plus OS survival and Aim 4 = mock-slide spatial. Your brief has Aim 3 = module scoring/liana and Aim 4 = TCGA validation. The plan itself numbers analyses by section (4 primary, 6 secondary, 7 TCGA), not by aim.
5. **TCGA outcome and model differ from the plan.** Stage 03 models overall survival on deconvolved fractions with no clinical covariates and records no 95% CI. The plan specifies PFS plus platinum status, a GSVA signature score, and adjustment for age, stage and residual disease.
6. **Two of the three named tools fell back.** BayesPrism (Stage 03) and RCTD (Stage 04) both logged `[FALLBACK]` to NNLS. Only inferCNV ran as the primary tool.

**Run provenance**
7. **The profile labels conflict.** The `01_run_log.txt` and `02_run_log.txt` headers say `profile: default`, but both ran on the real `data/sce.rds`. `03` and `04` say `profile: real`. The `default` and `real` profiles differ in `config.yml` (confounds, thresholds), so the exact config used for 01/02 cannot be read from the logs.
8. **`validation_log.tsv` mixes runs.** It includes two 2026-07-21 Stage 01 FAIL rows (ARI 0.153) and two failed 2026-08-06 Stage 03 runs (HR = NA). The files do not say which data or profile the 07-21 rows used.
9. **Stage 01 fails its own gate on real data** (ARI 0.643 < 0.70) and is still the input to 02–04. `config.yml` says the `real` threshold awaits a justified decision, and none is recorded.
10. **The rounded ref_fpr differs between files.** `validation_log.tsv` shows `ref_fpr=0.01` (rounded with `%.2f`); `01_metrics.csv` has 0.00534.
11. **One TCGA patient is duplicated.** TCGA-23-1023 enters the Stage 03 Cox model twice (`-01A`, `-01R`), so n = 425 samples covers 424 patients. Stage 03 says 427 primary samples kept, then 425 with survival; Stage 05 reads 434 bulk samples and aligns 425.
12. **The site-comparison object is out of date.** `data/sce_site_comparison.rds` (147 GB, 2026-08-07) predates the "store logcounts as sparse" loader change (commit 70a0ecf). The re-run after that commit failed at biomaRt, and that is the uncommitted `site_load_log.txt`. GSE180661 does not appear in the analysis plan.

**Figures and results**
13. **Three outputs have no generating script.** `figs/tcga_celltype_vs_survival_scatter.png`, `figs/tcga_celltype_kde.png` and `results/tcga_celltype_survival_scatter_data.csv` cannot be traced to code: `proportion_survival_plot.R` is 0 bytes, and no script mentions `kde`.
14. **Two figures are synthetic leftovers.** `figs/03b_benchmark.png` and `figs/03_fraction_scatter.png` come from the initial synthetic commit (9fb71b7) and have no matching file in `results/`. Neither should be used as a real-data figure.
15. **The Stage 04 PASS is not a biological finding.** It shows recovery of a planted layout on a mock slide.
16. **The Stage 05 gene sets are older than the plan.** It scores `resistance_up_core`, not `resistance_up_intrinsic`, and has no `cheon2014_collagen`.

**Environment**
17. **No sessionInfo exists for real-data stages 00–04.** The root `sessionInfo.txt` comes from a Windows machine (R 4.3.1 ucrt, Asia/Dubai, `splatter` attached), so it describes the synthetic pipeline, not this server.
18. **The R build versions are mixed.** Real runs use R 4.3.1 (conda, matches `renv.lock`), but every run log warns that packages were built under R 4.3.2/4.3.3 (yaml, GSVA, jsonlite, future, SingleCellExperiment, matrixStats, …).
19. **GenomeInfoDb versions differ.** It is 1.38.8 in the Windows sessionInfo and 1.38.1 in the server Stage 05 sessionInfo.
20. **The `hgsoc_ccc` env is in an unexpected place.** It lives at `/userhome/syeed/envs/hgsoc_ccc`, not under `/userhome/conda/envs/` with `hgsoc_pipeline`. Scripts that assume the conda root path will not find it.
21. **Package versions differ between the Windows sessionInfo and the current server env:** infercnv 1.18.1 vs **1.20.0**, survival 3.5-5 vs **3.8.3**, data.table 1.17.0 vs 1.17.8. The real Stage 01 inferCNV run has no sessionInfo, so its exact inferCNV version is not recorded. The current env has 1.20.0.
22. **Tools are missing from the current env.** BayesPrism and spacexr are **not installed**, which explains the Stage 03/04 `[FALLBACK]`s. `speckle` (needed for plan secondary analysis 6, propeller) is **not installed** either.
