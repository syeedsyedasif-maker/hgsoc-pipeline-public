# PLAN: Simulated HGSOC scRNA-seq → deconvolution → spatial scaffold

**Purpose.** A *fully simulated*, reproducible pipeline that mirrors the four dissertation aims for
cisplatin resistance in High-Grade Serous Ovarian Cancer (HGSOC). Because the data is simulated, we
**know the right answer** ("ground truth"). Every stage is scored against that planted truth, so the
scaffold proves it can recover known signal *under realistic confounds* before real data arrives.

> **This validates pipeline MACHINERY, not biology.** The genes are synthetic. We never interpret them
> biologically (no GSEA/pathways). We only ask: *does each method recover the signal we planted?*

When real data is ready: delete `00_simulate.R` and the `truth/` checks, point `01` to `04` at the real
central object, and the rest of the pipeline still runs.

---

## 1. One data ecosystem: the central object

**Choice: `SingleCellExperiment` (SCE).** One-line justification: `splatter` emits an SCE natively and
the whole Bioconductor stack (scater/scran, matrix extraction for every downstream tool) speaks it, so
we avoid a lossy Seurat round-trip. All tool-specific formats are produced by small **adapter
functions** (`R/adapters.R`), not by ad-hoc conversions scattered through the scripts.

Central artefact: `data/sce.rds` (counts + logcounts + `colData` with `celltype`, `site` +
`rowData` with chromosome arm and CNV/DE truth columns).

---

## 2. Biological model the simulation encodes

### Cell types (splatter "groups")
| Group | Role | In sc reference? | In bulk/spatial? |
|---|---|---|---|
| `malignant_ascites` | tumour, ascites site | ✅ | ✅ |
| `malignant_omental` | tumour, omental site | ✅ | ✅ |
| `fibroblast` (CAF) | stroma | ✅ | ✅ |
| `macrophage` | immune | ✅ | ✅ |
| `endothelial` | stroma | ✅ | ✅ |
| `Tcell` | immune | ✅ | ✅ |
| `adipocyte` | stroma (fat) | ❌ **omitted** | ✅ |

Malignant is about 45% of cells (split across the two malignant groups). Stroma and immune make up the
rest. Exact proportions live in `config.yml`.

### Two sites (`ascites`, `omental`): a property of **every** cell
- **Malignant** identity is site-specific by construction, so it becomes two splatter groups. This is
  what lets us plant DE *between sites within the tumour* (the Aim 2 resistance signal).
- **Stroma/immune** share one expression programme across sites, but their **abundance differs by
  site** (e.g. omentum is CAF/adipocyte-rich; ascites is spheroid-like, malignant-dominated). Site is
  assigned as metadata to hit the per-site composition in `config.yml`.

### Planted ground truth (written to `truth/`)
1. **Resistance DE programme.** Up in `malignant_omental` vs `malignant_ascites`, created with
   splatter `de.prob`/`de.facLoc`. Truth = the gene set + fold-changes. (Aim 2 target.)
2. **CNV map.** Malignant groups carry contiguous chromosome-arm gains/losses (clonal, shared by both
   malignant groups); stroma stays diploid. Truth = arm→gain/loss map + per-cell aneuploidy label.
   (Aim 1 target.)
3. **Cell-type identity DE.** Every group gets distinct expression so clustering/deconvolution *can*
   separate them (expected, not the headline signal).
4. **Per-patient proportions** (Aim 3) and **spatial layout** (Aim 4), written as tables/matrices.

### Genes → pseudo-genome
splatter gives counts only. We assign genes **in order** to a fixed ordered list of chromosome arms
(`1p,1q,…,22q,Xp,Xq`), giving each a `chr/start/stop`. This produces the **inferCNV gene-ordering
file** and defines the contiguous blocks we scale for CNV.

---

## 3. Realism confounds, so the pipeline can't win trivially

All **toggleable** in `config.yml` (`confounds:` block). Modelled on Hippen 2023 deconvolution artifacts.

| Confound | What it mimics | How we implement it | Where it bites |
|---|---|---|---|
| `batch_shift` | sc reference differs from the bulk assay batch | per-gene log-normal multiplier applied to **bulk/spatial only** | deconvolution (03/04) |
| `mrna_enrichment` | poly-A capture vs rRNA depletion | scale up a defined gene subset in **bulk only** | deconvolution (03) |
| `adipocyte_blindspot` | dissociation loss of fragile cells | adipocyte present in **bulk+spatial**, **absent from sc reference** | scored as a *result* (03) |

The adipocyte blind spot is a **key result, not a bug**: we explicitly score that deconvolution fails
to detect a population it has no reference for.

---

## 4. Parameter table (defaults; full set in `config.yml`)

| Parameter | Default | Tiny (thin-slice/test) |
|---|---|---|
| genes | 2000 | 300 |
| cells | ~3500 | ~400 |
| spatial grid | 30 × 30 | 10 × 10 |
| patients (survival) | 120 | 30 |
| threads | 1 (reproducible) | 1 |
| global seed | 1000 | 1000 |
| CNV gain / loss factor | 2.0 / 0.5 | same |
| resistance fold threshold (truth) | log2(1.5) | same |

Defaults target **a few minutes on a laptop**. Each script also sets its own `set.seed()`.

---

## 5. Package list: install source and fallback

| Package | Source | Used for | Fallback if install fails |
|---|---|---|---|
| splatter | Bioconductor | data generation | none (required; **installed ✓**) |
| SingleCellExperiment, SummarizedExperiment | Bioconductor | central object | none (**installed ✓**) |
| scater | Bioconductor | normalise/PCA | base PCA (**installed ✓**) |
| scran, igraph, bluster | Bioconductor | graph clustering | base `kmeans` on PCs |
| simsurv | CRAN | survival times | hand-rolled exponential hazard (**installed ✓**) |
| survival | CRAN | Cox/KM | none (**installed ✓**) |
| ggplot2, pheatmap | CRAN | figures | base graphics (**installed ✓**) |
| **infercnv** | Bioconductor (+ JAGS) | Aim 1 CNV calling | **`R/cnv_fallback.R`**: windowed expression vs diploid reference (the inferCNV idea, no JAGS) |
| **BayesPrism** | GitHub | Aim 3 deconvolution | **`R/deconvolve.R`**: NNLS on cell-type signatures |
| **spacexr (RCTD)** | GitHub | Aim 4 spatial deconv | **`R/rctd_fallback.R`**: per-spot NNLS |

**Update (tools now installed):** Rtools43 and JAGS 4.3.2 were later installed, and inferCNV,
BayesPrism, and RCTD were validated end to end on the `tiny` profile (7/7 PASS). At full scale only
inferCNV is practical, so the `default` profile keeps inferCNV primary and uses the fallbacks for Aims
3 and 4 (BayesPrism takes about 2.5 h on 120 patients, RCTD runs out of memory on the 900-spot grid).
The named tools stay as the primary path, gated by `tools:` toggles in `config.yml` plus a
`requireNamespace()` check. The fallbacks are clearly labelled `[FALLBACK]` in every log line and
reproduce the logic of each named tool. See the README for the per-tool detail and the `03b`
comparative benchmark.

---

## 6. File tree

```
sim_pipeline/
├─ PLAN.md  README.md  config.yml  install.R  Dockerfile  run_all.R
├─ sessionInfo.txt  renv.lock                      # reproducibility (generated)
├─ R/
│  ├─ utils.R         # load_config, seeding, logging, PASS/FAIL printer
│  ├─ adapters.R      # SCE → inferCNV / BayesPrism / SpatialRNA formats
│  ├─ metrics.R       # ARI, confusion, precision/recall/F1, RMSE, Pearson, spatial corr
│  ├─ deconvolve.R    # BayesPrism wrapper + NNLS fallback
│  ├─ cnv_fallback.R  # inferCNV wrapper + windowed-expression fallback
│  └─ rctd_fallback.R # RCTD wrapper + per-spot NNLS fallback
├─ 00_simulate.R      # generate truth (splatter + CNV overlay + confounds)
├─ 01_qc_cluster.R    # Aim 1: QC/cluster + CNV → malignant vs stroma
├─ 02_dge.R           # Aim 2: DE omental vs ascites malignant
├─ 03_deconvolution.R # Aim 3: pseudobulk + survival + deconvolution
├─ 04_spatial.R       # Aim 4: mock slide + spot deconvolution + colocalisation
├─ data/   truth/   results/   figs/             # all generated
```

---

## 7. Validation plan (each stage prints one PASS/FAIL line)

| Stage | Recovers | Metric vs truth | PASS rule (default) |
|---|---|---|---|
| 00 | n/a | self-checks: truth files written, dims, CNV signal present | all truth files exist |
| 01 | malignant vs stroma | confusion matrix + ARI vs true labels | ARI ≥ 0.7 |
| 02 | resistance genes | precision/recall/F1 vs planted set; volcano | F1 ≥ 0.7 |
| 03 | cell-type fractions | Pearson + RMSE **per type**; adipocyte blind-spot; Cox+KM | mean Pearson ≥ 0.7 & adipocyte flagged undetected |
| 04 | mal-CAF colocalisation | spatial cross-correlation / permutation test vs layout | core>margin CAF enrichment, p<0.05 |

(120 patients in Aim 3 is a **demo, not a powered** survival analysis. This is stated in the output
and the README.)

---

## 8. Build order

1. **Thin vertical slice**: all five scripts run end-to-end on **tiny** data with fallbacks, each
   printing its PASS/FAIL line. Prove the skeleton before deepening.
2. **Deepen phase by phase** (00→04): richer CNV, confounds, proper clustering, simsurv, spatial
   layout, real statistics, running and validating after each.
3. **Reproducibility wrap**: `install.R` (tested), `Dockerfile`, `renv.lock`, `sessionInfo.txt`,
   `README.md`, `run_all.R`.

Single-threaded by default; `threads` exposed in config. `set.seed()` in every script.
