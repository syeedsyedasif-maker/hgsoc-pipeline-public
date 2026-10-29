# Primary analysis specification: stage 06 (GSE165897)

**Status:** draft, written 29 October 2026, before any resistant-versus-sensitive analysis was run. It describes what `06_primary_resistance.R` implements. The sources are `docs/ANALYSIS_PLAN.md` v1.1 ("Plan") and `DEVIATIONS.md` items 1–4 and D1–D5 ("Dev").
**Not yet decided:** the items under "Ambiguities" below. Each has a stated default that the script implements. Change the default here and in the script before the true-label run if you disagree.

## 0. Inputs

| Input | File | Basis |
|---|---|---|
| Counts and cell labels | `data/sce_clustered.rds`: raw `counts`, author labels in `celltype` (from `cell_subtype`), our inferCNV call in `cnv_call` | Plan §2: "Cell types: author annotations"; Dev D1 needs `cnv_call`. See A12. |
| Patient labels | `metadata/gse165897_labels.csv` (`pfi_group`, `crs_group`) | Plan §2; Dev D2 |
| Sample flags | `metadata/gse165897_sample_flags.csv` (`site`, `frozen`) | Dev D4, D5 |
| Gene sets | `cisplatin_response_genes_HGSOC.json` → `genesets` | Plan §4, §6 |
| Parameters | `config.yml` → `real: primary:` | CLAUDE.md: every parameter in config |

Cell populations (Plan §2: "CAFs = the author fibroblast annotation. Malignant cells = author epithelial ovarian cancer (EOC) annotation for primary analyses"):
- **CAF** = `celltype == "fibroblast"` (author CAF-1, CAF-2, CAF-3)
- **Malignant** = `celltype == "malignant"` (author EOC)
- **Malignant, both** = author EOC **and** `cnv_call == "malignant"`. Dev D1: "'inferCNV-confirmed malignant cells' means cells called malignant by both our inferCNV and the authors."
- **Macrophage** = `celltype == "macrophage"`. Plan §6.3 names macrophages.

Groups:
- **PFI** = `pfi_group`: resistant (test) vs sensitive (reference). Plan §2: "resistant if platinum-free interval < 183 days".
- **CRS** = `crs_group`: poor_response (test) vs response (reference). Dev D2: "CRS 1 (poor response) vs CRS 2-3."

## 1. Machinery shared by every analysis

| Step | Implementation | Quote |
|---|---|---|
| Aggregation | Sum raw counts per patient × timepoint × population. Each patient × timepoint is one GEO sample. | Plan §3: "Raw counts summed per patient x timepoint x cell type." |
| Minimum cells | Drop a pseudobulk with < 20 cells (`primary.min_cells`) | Plan §3: "included only if it contains >= 20 cells." |
| Minimum patients | If any group has < 3 patients (`primary.min_patients_per_group`), write the analysis as `estimable = FALSE` with the reason. No test is run. | Plan §3: "A contrast is run only if >= 3 patients per group remain. Otherwise it is reported as not estimable." |
| Gene filter | `edgeR::filterByExpr(y, group = <see A2>)`, default thresholds | Plan §3: "Genes filtered with `edgeR::filterByExpr`" |
| Normalisation | `edgeR::calcNormFactors(method = "TMM")`, then `limma::voom(y, design)` | Plan §3: "TMM normalisation; `limma::voom`." |
| Set test | `limma::fry(v, index, design, contrast)`. Report `NGenes`, `Direction`, `PValue`. | Plan §4: "using `limma::fry` (directional)" |
| Set size rule | Test a set only if ≥ 5 of its genes pass the filter (`primary.min_set_genes`). The positive control is exempt. | Plan §3: "A gene set is tested only if >= 5 of its genes survive filtering (the positive control is exempt)." |
| Effect size | Mean of the voom/`lmFit` contrast log2FC over the set genes that pass the filter, plus n patients per group | Plan §4: "mean log2 fold change across set genes, with the number of patients per group." |

## 2. Analyses

| ID | Population | Timepoint | Grouping / design | Gene set (JSON key) | Expected | Multiplicity | Quote |
|---|---|---|---|---|---|---|---|
| **P1** | CAF | TN | PFI; `~ group` | `cheon2014_collagen` (10) | Up in resistant | Holm across P1, P2 | Plan §4 table; Dev item 1 |
| **P2** | Malignant | TN | PFI; `~ group` | `resistance_up_intrinsic` (14) | Up in resistant | Holm across P1, P2 | Plan §4 table; Dev item 2 |
| **PC** | Malignant | TN + PN, complete pairs | `~ patient + timepoint` (PN vs TN) | `stress_IEG_state_GSE165897` (3), exempt from size rule | Up in PN | none (pipeline check) | Plan §5 |
| **S1-P1** | CAF | TN | CRS; `~ group` | `cheon2014_collagen` | Up in poor_response | BH, family S1 | Plan §6.1. Dev D2: "the CAF contrast is reported as not estimable" |
| **S1-P2** | Malignant | TN | CRS; `~ group` | `resistance_up_intrinsic` | Up in poor_response | BH, family S1 | Plan §6.1; Dev D2 |
| **S2-P2** | Malignant, both | TN | PFI; `~ group` | `resistance_up_intrinsic` | Up in resistant | BH, family S2 (1 test) | Plan §6.2; Dev D1 |
| **S3** ×6 | Malignant, CAF, macrophage | TN | PFI; `~ group` | `drug_efflux` (5), `extrinsic_TME` (7) | not specified (A6) | BH, family S3 | Plan §6.3: "where >= 5 genes survive filtering" |
| **S4-CAF**, **S4-mal** | CAF / malignant | TN + PN, complete pairs | PFI; `~ patient + group:timepoint`; contrast (resistant PN−TN) − (sensitive PN−TN) | P1 set / P2 set (A6) | not specified | BH, family S4 | Plan §6.4: "paired design, `~ patient + group:timepoint`" |
| **S5** | each author cell type where estimable (A7) | TN | PFI; `~ group`; `voom` + `eBayes`, all genes | gene level; columns flag membership of `matrix_CAF`, `sensitivity_up`, `partial_EMT_Xu2022` and all other sets | — | BH across genes, within each cell type | Plan §6.5: "Reported as exploratory only"; Plan §6.3 / Dev item 3: below-minimum sets "reported only within the gene-level results" |
| **D1** | = S2-P2 | | | | | | Dev D1 adds no new analysis; it defines S2's population |
| **D4-PC** | Malignant | TN + PN, pairs, excluding EOC1005 and EOC3 | as PC | as PC | Up in PN | none (sensitivity) | Dev D4 |
| **D5a-PC**, **D5a-S4-CAF**, **D5a-S4-mal** | as PC / S4 | as PC / S4, restricted to EOC153, EOC227, EOC372, EOC443, EOC540 | as PC / S4 | as PC / S4 | direction only | none; labelled "direction check only" | Dev D5(a) |
| **D5b-P1**, **D5b-P2** | as P1 / P2 | TN | PFI; `~ baseline_site + group`, site = Peritoneum vs other | as P1 / P2 | Up in resistant | none (sensitivity) | Dev D5(b) |

The patient lists for D4 and D5a are stored in `config.yml`. The script also derives them from `gse165897_sample_flags.csv` (frozen status differs within a pair; TN site equals PN site) and stops if the two disagree.

**Primary decision rule** (Plan §4): "A primary result is positive only if Holm-adjusted p < 0.05 in the expected direction." Implemented as `p_holm < 0.05` **and** `Direction == "Up"`. See A1.

**Out of scope for stage 06**, with its own reason:
- Plan §6.6 composition (`speckle::propeller`): not a P1/P2 analysis. It needs per-sample proportions, not pseudobulk.
- Plan §6.7 / Dev D3 liana-py: "run only if the primary tests and the TCGA validation are complete by 24 Nov 2026".
- Plan §7 TCGA-OV validation: "run once". It needs the P1/P2 result first.
- GSE300897: Dev D3, "not analysed".

## 3. Outputs

All outputs go to `results/06_primary/` (`primary.out_dir`, overridable with `HGSOC_OUT_DIR`):

| File | Contents |
|---|---|
| `06_pseudobulk_inventory.csv` | Every pseudobulk considered: analysis ID, patient, timepoint, group, n_cells, included or reason excluded |
| `06_geneset_tests.csv` | One row per set test: ID, family, population, design, n patients per group, set size, genes passing the filter, `estimable`, reason, `Direction`, `PValue`, adjustment method, adjusted p, mean log2FC, expected direction, `positive` (P1/P2 only), note |
| `06_genelevel_<celltype>.csv` | S5 `topTable` for each estimable cell type, plus set-membership columns |
| `06_sessionInfo.txt` | `sessionInfo()` at the end of the run |

The console log reports only counts, inclusions and warnings. Test statistics go only to the CSVs.

**Label-shuffle test mode:** with `HGSOC_SHUFFLE_LABELS=<seed>`, the script permutes `pfi_group` and `crs_group` independently across patients, using that seed, before any analysis. It refuses to write into `results/06_primary/` in this mode, so `HGSOC_OUT_DIR` must point elsewhere.

## 4. Ambiguities (default implemented → for you to confirm)

- **A1. Directional p for fry.** The plan says "directional" and "in the expected direction" but does not say one-sided or two-sided. **Default:** fry's two-sided directional `PValue`, with `Direction` required to match. The p is not halved.
- **A2. Grouping passed to `filterByExpr`.** Not specified. **Default:** the contrast's own grouping. That is `pfi_group` or `crs_group` for group designs, timepoint for PC, and group × timepoint for S4. Filtering is done separately in each analysis's own samples.
- **A3. "Patient as blocking factor" (PC).** It could mean a fixed patient effect or `duplicateCorrelation(block = patient)`. **Default:** fixed effect, `~ patient + timepoint`. This matches the S4 formula in Plan §6.4.
- **A4. When is the positive control "reproduced"?** Plan §5 gives no threshold. **Default:** `PValue < 0.05` and `Direction == "Up"` (PN > TN), with no adjustment.
- **A5. Definition of the BH "families"** (Plan §6: "Benjamini-Hochberg within each family"). **Default:** each numbered item of Plan §6 is one family. S1 = {S1-P1, S1-P2}, S2 = {S2-P2}, S3 = the 6 tests, S4 = the 2 tests, S5 = genes within each cell type. Tests that are not estimable or have too few genes are left out of their family's adjustment.
- **A6. S3 and S4 details.** Plan §6.3 gives no expected direction. Plan §6.4 names no population, gene set or label. **Default:** no expected direction for either; Direction and p are reported. S4 uses PFI, the P1 set in CAFs and the P2 set in malignant cells.
- **A7. Which cell types S5 covers** ("Gene-level differential expression per cell type"). **Default:** every author cell type in the SCE (malignant, fibroblast, macrophage, Tcell, endothelial) that meets the 20-cell and 3-patient rules. That is decided from cell counts alone.
- **A8. Effect-size genes.** **Default:** only set genes that pass `filterByExpr`, unweighted mean of the contrast coefficient.
- **A9. Which analyses D5(a) covers** ("before-after analyses"). **Default:** PC and S4. From cell counts alone:
  - D5a-PC has 4 complete pairs (EOC540 post-NACT has 18 malignant cells).
  - D5a-S4-mal has 1 resistant patient with complete pairs, so it is not estimable.
  - D5a-S4-CAF has 2 sensitive patients with complete pairs (EOC227 has under 20 CAFs), so it is not estimable.
- **A10. Minimum for paired analyses.** The "≥ 3 patients per group" rule is written for group contrasts. **Default:** paired analyses (PC, D4, D5a-PC) need ≥ 3 complete pairs. A patient enters a paired analysis only if both pseudobulks have ≥ 20 cells.
- **A11. D5(b) site coding and multiplicity.** **Default:** baseline site = TN `site` from `gse165897_sample_flags.csv`, coded Peritoneum vs other (Omentum, Mesentery). p is reported unadjusted and kept out of the Holm pair. If site and group make the design rank-deficient, the analysis is marked not estimable.
- **A12. Input object.** The plan does not name one. **Default:** `data/sce_clustered.rds`. It holds the same 42,126 cells and counts as `data/sce.rds`, plus `cnv_call`, which S2 needs.
- **A13. S2 cell minimum.** **Default:** the 20-cell rule is applied after intersecting author EOC with `cnv_call`. Some treatment-naive samples may drop (for example, EOC3 has 19 inferCNV-called malignant cells in TN).
- **A14. Plan §7 signature rule.** It needs "nominal p < 0.05" genes from the primary contrast. Stage 06 writes per-gene contrast statistics for the P1/P2 populations in the S5 gene-level tables. No gene selection is done here.
- **A15. Holm when one primary is not estimable.** The plan says "Holm-adjusted across the two". **Default:** `p.adjust(..., "holm", n = 2)` is always used. If only one of P1/P2 can be tested, its p is still multiplied by 2, so the result is never less conservative than planned.
- **A16. Gene-set test rows for S5.** S5 is gene-level only. `06_geneset_tests.csv` contains S5 rows only when a cell type is not estimable, to record why.
