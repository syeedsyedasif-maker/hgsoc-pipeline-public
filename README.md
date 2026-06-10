# Simulated HGSOC pipeline: scRNA-seq, DE, deconvolution, spatial

A fully simulated, reproducible scaffold for a dissertation on cisplatin resistance in High-Grade
Serous Ovarian Cancer (HGSOC). It generates single-cell data with a known ground truth, runs the four
dissertation aims, and scores each one against that truth. When real data is ready, you delete the
simulator and point the same scripts at your real dataset.

### The four aims
1. **Aim 1, find the tumour.** Separate malignant (aneuploid) cells from stroma using a CNV signal.
2. **Aim 2, resistance programme.** Recover the genes differentially expressed in omental versus
   ascites malignant cells (the planted cisplatin-resistance signal).
3. **Aim 3, deconvolution and survival.** Estimate cell-type fractions from bulk pseudobulk under
   realistic confounds, and test whether composition predicts survival.
4. **Aim 4, spatial structure.** Recover malignant and CAF colocalisation in a mock tissue slide.

> **Important: this validates machinery, not biology.** Every gene here is synthetic. The pipeline
> tests whether the analysis machinery recovers a signal that was deliberately planted, under
> realistic noise. It says nothing about real biology. Do not run pathway, GSEA, or literature
> interpretation on these genes. The only question each stage answers is: given data where the answer
> is known, does the method find it?

---

## What you get (at a glance)

Five numbered scripts, one per aim, each printing a single `[PASS]/[FAIL]` line that compares its
output to the planted truth:

| Script | Aim | Recovers | Validation (PASS rule) |
|---|---|---|---|
| `00_simulate.R` | n/a | plants the truth | CNV signal present, resistance set non-empty |
| `01_qc_cluster.R` | 1 | malignant vs stroma | Adjusted Rand Index >= 0.70 vs true labels |
| `02_dge.R` | 2 | resistance DE programme | F1 >= 0.70 vs planted gene set |
| `03_deconvolution.R` | 3 | cell-type fractions, survival | mean Pearson >= 0.70, adipocyte blind spot shown, Cox direction |
| `04_spatial.R` | 4 | malignant and CAF colocalisation | CAF concentrated in core (permutation p < 0.05) |

Latest runs, both 7/7 PASS:
- `default` (inferCNV primary, NNLS for Aims 3 and 4, about 6 min): ARI 0.99, F1 0.85, deconvolution
  Pearson 0.91, adipocyte 8% missed as expected, malignant survival hazard ratio above 1, spatial
  permutation p 0.001.
- `tiny` (all three named tools, inferCNV with BayesPrism and RCTD, about 7 min): 7/7, deconvolution
  Pearson 0.98. See "Primary tools vs fallbacks" below for the comparison and the full-scale caveats.

---

## Quick start

### Option A: local R (R 4.3.x with Bioconductor 3.18)
```r
# from inside sim_pipeline/
Rscript install.R          # core packages (binaries; no compiler needed)
# (optional) INSTALL_OPTIONAL=1 Rscript install.R   # also build inferCNV/BayesPrism/RCTD
Rscript run_all.R tiny     # about 7 min: runs all three named tools end to end
Rscript run_all.R          # full run, about 6 min (inferCNV primary; NNLS for Aims 3 and 4)
```
Outputs land in `truth/`, `results/`, `figs/`, plus `sessionInfo.txt` and `results/validation_log.tsv`.

### Option B: Docker (fully pinned, includes the real primary tools)
```bash
docker build -t hgsoc-sim .          # installs JAGS plus inferCNV/BayesPrism/spacexr
docker run --rm hgsoc-sim            # runs run_all.R and prints the PASS/FAIL summary
```
The Linux image has a compiler, so the optional tools build there even though the Windows host they
were developed on could not build them.

---

## The biological model the simulation encodes

- Cell types (splatter groups): `malignant`, `fibroblast` (CAF), `macrophage`, `endothelial`,
  `Tcell`, `adipocyte`. About 45% malignant, the rest stroma and immune.
- Two sites (`ascites`, `omental`) on every cell. Malignant is split into `malignant_ascites` and
  `malignant_omental` so DE can be planted between them. Stroma and immune share one programme across
  sites but differ in abundance by site (omentum is CAF and adipocyte rich, ascites is spheroid and
  malignant dominated).
- Planted truths (written to `truth/`):
  - Resistance DE programme, up in `malignant_omental` versus `malignant_ascites` (the Aim 2 signal).
  - CNV map: malignant cells carry clonal chromosome-arm gains and losses, stroma stays diploid (Aim 1).
  - Per-patient proportions and spatial layout (Aims 3 and 4).

---

## The realism confounds (so the pipeline cannot win trivially)

Modelled on the deconvolution artifacts described in Hippen 2023. All toggleable in `config.yml`:

| Confound | Mimics | Where it bites |
|---|---|---|
| batch shift | sc reference differs from bulk/spatial assay batch | a per-gene multiplier is applied to bulk/spatial only, so reference and bulk are not biologically identical, which makes deconvolution harder |
| mRNA-enrichment | poly-A capture vs rRNA depletion | a defined gene subset is over-represented in the bulk |
| adipocyte blind spot | fragile cells lost in dissociation | adipocyte is in the bulk and spatial data but omitted from the sc reference |

The adipocyte blind spot is a result, not a bug. Script `03` scores that deconvolution cannot detect a
population it has no reference for, and prints the missed fraction (about 8%).

---

## What each validation proves (and how to read the figure)

- `01`, `figs/01_cnv_score_by_celltype.png`: a per-cell aneuploidy score (smoothed expression along
  the genome versus the diploid stroma reference). Malignant cells sit above the dashed threshold,
  stroma below. Reference-baseline check: `ref_fpr` is how often the reference cells themselves trip
  the threshold. If it were high the reference would be unreliable and every call would be suspect.
  Here it is about 1%.
- `02`, `figs/02_volcano.png`: omental versus ascites DE. Planted resistance genes (red) should sit in
  the significant tails. Both malignant groups carry the same CNV, so copy number cancels out of this
  contrast.
- `03`, `figs/03_fraction_scatter.png` and `figs/03_km_by_malignant_tertile.png`: estimated versus
  true cell-type fraction (points should hug the diagonal), and Kaplan-Meier survival split by
  estimated malignant tertile (higher malignant fraction gives worse survival, matching the planted
  hazard). Note on fractions: bulk deconvolution recovers the fraction of expression each type
  contributes, not the raw fraction of cells, because cell types carry different amounts of mRNA. We
  score against the expression-fraction truth (`results/03_true_fractions.csv`). 120 patients is a
  demonstration of the machinery, not a powered survival study.
- `04`, `figs/04_spatial_maps.png`: planted versus inferred malignant and CAF maps. CAF should light
  up the desmoplastic core. The PASS test is a permutation test asking whether inferred CAF is
  concentrated in the malignant-dense core versus the equally-malignant margin. The raw
  fraction-versus-fraction correlation is negative, which is expected: fractions are compositional and
  malignant spans both core and margin. That is why we use the permutation test, not a correlation.

---

## Primary tools vs fallbacks

Three named tools are the intended primary methods. They are installed (Rtools43 with JAGS 4.3.2) and
were validated end to end on the `tiny` profile: 7/7 PASS with `inferCNV`, `BayesPrism`, and `RCTD` all
running, no fallbacks. Each tool has a clearly labelled `[FALLBACK]` that reproduces the same logic, so
the pipeline runs anywhere.

| Aim | Primary tool | Fallback (always runs) |
|---|---|---|
| 1 | inferCNV | `R/cnv_fallback.R`: smoothed expression along the genome versus a diploid reference (inferCNV's core idea without the HMM) |
| 3 | BayesPrism | `R/deconvolve.R`: marker-based NNLS on cell-type signatures |
| 4 | RCTD (spacexr) | `R/rctd_fallback.R`: per-spot NNLS (a spot is a small bulk sample) |

### Named tool versus fallback: do they agree?
The named tools and the fallbacks recover the same planted signal. The named tools match or slightly
beat the fallbacks, which cross-validates the faster fallbacks.

| Aim | Metric | Fallback | Named tool | Scale tested |
|---|---|---|---|---|
| 1 | ARI vs true labels | 0.979 (windowed expr) | 0.988 (inferCNV) | full, about 5 min |
| 2 | F1 (no external tool) | 0.846 | 0.846 | DE is a plain t-test |
| 3 | mean Pearson | 0.909 (NNLS) | 0.981 (BayesPrism) | full |
| 4 | core-vs-margin perm p | 0.001 (NNLS) | 0.001 (RCTD) | tiny |

### Performance at full scale (honest caveats)
- inferCNV: about 5 min on 3500 cells. Practical, so it is primary in the default profile.
- BayesPrism: correct, but about 2.5 h single-threaded on 120 patients. Too slow for a few-minutes
  default, so the full default uses the NNLS fallback. Enable `use_bayesprism: true` for a real run
  and raise `threads`. The fallback matches it closely (0.91 versus 0.98 Pearson).
- RCTD: runs out of memory on the 30x30 (900-spot) grid in this environment. This is a hard native
  crash, not catchable by the fallback guard. Use it on `tiny`, shrink `spatial_grid`, or use a
  machine with more RAM.

So the `tiny` profile turns all three named tools on (small enough to be fast and fit in memory). The
`default` profile keeps inferCNV primary and uses the fallbacks for Aims 3 and 4.

Switching: edit the `tools:` block in `config.yml` (per profile). Each script checks
`requireNamespace()` and falls back automatically (logging `[FALLBACK]`) if a tool is missing, so
toggling is safe. inferCNV needs JAGS at runtime. If JAGS is a user-local install, `ensure_jags()` in
`R/cnv_fallback.R` auto-detects it, and `.Renviron` sets `JAGS_HOME` for runs started with `cd`.

---

## Comparative deconvolution benchmark (`03b`)

The default Aim 3 builds bulk from the same cells used as the single-cell reference. That is a
pseudobulk benchmark, and Hippen et al. 2023 argue it is too easy: real bulk and single-cell data are
not biologically identical, so naive methods can score well on pseudobulk and then drop on real bulk.
`03b_deconv_benchmark.R` adds the missing test. It builds a second bulk as an independent draw (a
separate splatter batch with its own per-gene platform shift, a different cell-type composition, and
two extra cell types dropped from the reference), then deconvolves both bulks against the same
reference using NNLS, MuSiC, and InstaPrism.

On the `medium` profile the expected pattern appears. Mean Pearson against the true fractions:

| Method | pseudobulk | independent-draw | change |
|---|---|---|---|
| NNLS (naive) | 0.81 | 0.45 | -0.36 |
| MuSiC (naive) | 0.87 | 0.78 | -0.09 |
| InstaPrism (Bayesian) | 0.92 | 0.93 | +0.01 |

The naive methods drop under the reference-to-bulk mismatch. The Bayesian method holds. InstaPrism is
a documented fast approximation of BayesPrism and is the engine that produced the Bayesian numbers.
True BayesPrism was run once for confirmation but did not complete within the 900 second per-condition
budget at this scale, so it is recorded as "did not complete" and the InstaPrism numbers stand. This
is a hypothesis test, not a build gate. The separation depends on the platform mismatch being
substantial; that setting is recorded in `truth/benchmark_notes.md`. See also `figs/03b_benchmark.png`.
Run it with `SIM_PROFILE=medium Rscript 03b_deconv_benchmark.R`.

---

## Configuration

Everything is in `config.yml`, with profiles `default` (full), `tiny` (fast smoke test), and `medium`
(used by the `03b` benchmark). Select a profile with `Rscript run_all.R tiny` or the `SIM_PROFILE`
environment variable. Key settings: data sizes, the global `seed`, `threads` (single-threaded by
default for reproducibility), DE and CNV strength, the three confound toggles, the survival hazard
model, the spatial layout, the tool toggles, the benchmark block, and the PASS/FAIL thresholds.

---

## Outputs

```
truth/    resistance_genes.csv, cnv_map.csv, cell_metadata.csv, gene_order.tsv,
          composition_by_site.csv, benchmark_notes.md
results/  per-stage metrics CSVs, validation_log.tsv (all PASS/FAIL lines), peek_* intermediates,
          03b_benchmark.csv
figs/     one figure per aim, plus 03b_benchmark.png
data/     sce.rds (the central SingleCellExperiment), sce_clustered.rds
```

---

## Swapping in real data

The central object is a single `SingleCellExperiment` (`data/sce.rds`) with `colData` columns
`celltype` and `site`. Tool-specific formats are produced only by the adapters in `R/adapters.R`. To
go live:
1. Delete `00_simulate.R` and replace `data/sce.rds` with your real SCE (same `celltype` and `site`
   columns, gene coordinates in `rowData` for Aim 1).
2. Delete the VALIDATION blocks at the end of `01` to `04` (clearly marked) and the `truth/` reads. The
   analysis above them stands on its own.
3. Turn on the `tools:` toggles once inferCNV, BayesPrism, and RCTD are installed.

---

## Reproducibility notes

- R 4.3.1 with Bioconductor 3.18, single-threaded, `set.seed()` in every script.
- Exact package versions are pinned in `renv.lock` (`renv::restore()` to reproduce; 255 packages
  including inferCNV, BayesPrism, spacexr). A snapshot of the run environment is in `sessionInfo.txt`.
- Built on Windows 11 with Rtools43 and JAGS 4.3.2 installed. The fallbacks were validated first (no
  compiler needed), then all three named tools were installed and validated on `tiny`.
- inferCNV runtime note: it imports `rjags`, which needs JAGS. JAGS here is a user-local install
  (`%LOCALAPPDATA%\Programs\JAGS`) with no registry entry, and `rjags` expects `JAGS_HOME` to be the
  JAGS version directory (it appends `/x64/bin` itself). `ensure_jags()` in `R/cnv_fallback.R`
  auto-detects this, and `.Renviron` sets it for runs started with `cd`. inferCNV and Seurat also need
  Matrix 1.6.4 or newer.
- The Dockerfile reproduces the environment with the primary tools (the Linux image has a compiler).
- See `PLAN.md` for the full design rationale, parameter table, and package and fallback table.

---

## Acknowledgement

This pipeline was designed and directed by Syeed Syed Asif as part of an MSc dissertation. The study
design, biological rationale, and interpretation of results are the author's own. The R implementation
was written with the assistance of Anthropic's Claude (via Claude Code) under the author's direction.

---

## License

MIT License. See the `LICENSE` file.
