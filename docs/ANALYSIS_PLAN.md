# Pre-specified analysis plan

**Project:** Patient-level single-cell analysis of platinum resistance in HGSOC (GSE165897), with validation in TCGA-OV
**Author:** Syeed Syed Asif
**Version:** 1.0, 23 September 2026
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
| P1 | CAFs | Cheon et al. (2014) collagen-remodelling signature | Up |
| P2 | Malignant | `resistance_up_core` (project gene-set JSON) | Up |

A primary result is positive only if Holm-adjusted p < 0.05 in the expected direction. Effect sizes are reported as the mean log2 fold change across set genes, with the number of patients per group.

## 5. Positive control

`stress_IEG_state_GSE165897` (FOS, JUN, NR4A1) should be higher in PN than TN malignant cells (paired, patient as blocking factor). If this is not reproduced, the pipeline is checked before any other result is interpreted.

## 6. Secondary analyses (exploratory; Benjamini-Hochberg within each family)

1. P1 and P2 repeated with CRS labels.
2. P1 with `matrix_CAF` in place of the Cheon signature; P2 with inferCNV-confirmed malignant cells.
3. Other project gene sets (`sensitivity_up`, `drug_efflux`, `partial_EMT_Xu2022`, `extrinsic_TME`) in malignant cells, CAFs and macrophages.
4. Resistance x timepoint interaction: paired design, `~ patient + group:timepoint`.
5. Gene-level differential expression per cell type. Reported as exploratory only; no gene-level result is claimed as a finding.
6. Cell-type composition: `speckle::propeller` on per-sample proportions (logit transform).
7. Cell-cell communication in TN samples, inferred per patient. Pre-specified axes: malignant to CAF, CAF to malignant, macrophage to malignant. Patient-level interaction scores compared between groups (Wilcoxon).

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
