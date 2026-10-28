#!/usr/bin/env bash
# Read-only evidence collector. Writes only WRITEUP_STATE.md and writeup_bundle/ in the repo root.
# Run from hgsoc-pipeline-public/ with hgsoc_pipeline active.
set -u
REPO=/data/New_server_layout_example/Papers_and_Projects/syeed_MScProject/hgsoc-pipeline-public
SCR=$(dirname "$(readlink -f "$0")")
cd "$REPO"
OUT=WRITEUP_STATE.md
GIT_STATUS_BEFORE="$(git status --short)"   # captured before we create anything
CCC_PY=/userhome/syeed/envs/hgsoc_ccc/bin/python

{
echo "# WRITEUP_STATE — evidence pack"
echo
echo "Generated $(date '+%Y-%m-%d %H:%M %Z') on $(hostname) by a read-only collector. Nothing was edited, rerun, committed or deleted."
echo "Every number below is copied from a file on disk; the source file is named next to it. Where a planned output does not exist, it says **NOT RUN**."
echo

# ---------------------------------------------------------------- 1
echo "## 1. Repo state"
echo
echo "- Repo: \`hgsoc-pipeline-public/\`"
echo "- Branch: \`$(git branch --show-current)\`  (HEAD \`$(git rev-parse --short HEAD)\`, $(git log -1 --format=%ci))"
echo
echo "\`git log --oneline --since=2026-09-01\`"
echo '```'; git log --oneline --since=2026-09-01; echo '```'
echo
echo "\`git status --short\` (captured before this pack was written)"
echo '```'; echo "$GIT_STATUS_BEFORE"; echo '```'
echo
echo "The one modified file, \`site_load_log.txt\`, is a re-run of \`00_load_real_data_site.R\` (GSE180661) that **failed** at the biomaRt step (\"Could not reach Ensembl via biomaRt\", exit status 1). The committed version of the log records the earlier successful run."
echo

# ---------------------------------------------------------------- 2
echo "## 2. ANALYSIS_PLAN.md (verbatim)"
echo
echo "- Path: \`docs/ANALYSIS_PLAN.md\`"
echo "- Version (from file header): **$(grep -m1 -o 'Version:\*\* [0-9.]*' docs/ANALYSIS_PLAN.md | sed 's/Version:\*\* //')**"
echo "- Last commit: \`$(git log -1 --format='%h %ci %s' -- docs/ANALYSIS_PLAN.md)\`"
echo "- All commits touching it:"
echo '```'; git log --format='%h %ci %s' -- docs/ANALYSIS_PLAN.md; echo '```'
echo
echo '````markdown'; cat docs/ANALYSIS_PLAN.md; echo '````'
echo
echo "\`DEVIATIONS.md\` (referenced by the plan), verbatim:"
echo
echo '````markdown'; cat DEVIATIONS.md; echo '````'
echo

# ---------------------------------------------------------------- 3
echo "## 3. gse165897_labels.csv (verbatim)"
echo
echo "- Path: \`metadata/gse165897_labels.csv\` — $(wc -l < metadata/gse165897_labels.csv) line(s), $(stat -c %s metadata/gse165897_labels.csv) bytes"
echo "- Last commit: \`$(git log -1 --format='%h %ci %s' -- metadata/gse165897_labels.csv)\`"
echo
echo '```'; cat metadata/gse165897_labels.csv; echo '```'
echo
echo "**The file contains only the header row. No patient has a PFI or CRS label.**"
echo
echo "Verification against Zhang et al. 2022 Table 1: **none found.**"
echo "- \`git log --all -i --grep=verif --grep=Zhang --grep='Table 1'\`: $(git log --all -i --grep='verif' --grep='Zhang' --grep='Table 1' --oneline | wc -l) matching commits."
echo "- The only text mentioning Table 1 is the plan itself (\"Values are taken from Zhang et al. (2022), Table 1\"), which describes intent, not a completed check. DEVIATIONS.md does not mention the labels."
echo

# ---------------------------------------------------------------- 4
echo "## 4. Stage 2 code check"
echo
echo "**No script in the repo runs a resistant-vs-sensitive DE.** A search of all \`.R\`/\`.py\`/\`.ipynb\` files for \`PFI\`, \`crs\`, \`resistant\`, \`fry\`, \`voom\`, \`edgeR\`, \`limma\`, \`propeller\`, \`liana\` returns no analysis code (only \`coxph\` in \`03_deconvolution.R\`)."
echo
echo "The only Stage 2 DE script is \`02_dge.R\` (last modified $(stat -c '%y' 02_dge.R | cut -d. -f1)), which tests **post-NACT vs treatment-naive** malignant cells. Its aggregation and test lines:"
echo
echo '```r'
sed -n '47,60p' 02_dge.R | awk '{printf "%3d  %s\n", NR+46, $0}'
echo '   ...'
sed -n '81,94p' 02_dge.R | awk '{printf "%3d  %s\n", NR+80, $0}'
echo '```'
echo
echo "- **Pseudobulk by patient × treatment_phase** (raw counts summed with \`rowSums\`, groups < 20 cells dropped), **not** a per-cell test."
echo "- Test: CPM + log2, then a **paired t-test written in base R** (\`pt()\`, \`matrixStats::rowSds\`), BH FDR via \`p.adjust\`. **No DE package is used** (no edgeR, limma/voom or limma::fry, which the plan specifies)."
echo "- Cell type: author EOC annotation (\`celltype == \"malignant\"\`); CAFs are not tested."
echo

# ---------------------------------------------------------------- 5
echo "## 5. Inventory: results/, figs/, data/"
echo
echo "Labels: **REAL** = produced from GSE165897, GSE180661 or TCGA-OV. **SIM** = synthetic/simulation output. **REAL cells, SIM layout** = real GSE165897 cells placed on a simulated slide with a planted architecture."
echo
echo "| File | Size | Modified | Origin | Produced by |"
echo "|---|---:|---|---|---|"
find results figs data -type f -printf '%p\t%s\t%TY-%Tm-%Td %TH:%TM\n' | sort | while IFS=$'\t' read -r p s t; do
  hs=$(numfmt --to=iec --suffix=B "$s")
  case "$p" in
    data/sce.rds) o="REAL"; b="00_load_real_data.R (GSE165897)";;
    data/sce_clustered.rds) o="REAL"; b="01_qc_cluster.R";;
    data/sce_site_comparison.rds) o="REAL"; b="00_load_real_data_site.R (GSE180661); predates commit 70a0ecf";;
    results/infercnv/*) o="REAL"; b="inferCNV via 01_qc_cluster.R";;
    results/01_*|figs/01_*) o="REAL"; b="01_qc_cluster.R (log says profile: default)";;
    results/02_*|results/peek_02*|figs/02_*) o="REAL"; b="02_dge.R (log says profile: default)";;
    results/03_*|figs/03_km*) o="REAL"; b="03_deconvolution.R (profile: real)";;
    results/04_*|figs/04_*) o="REAL cells, SIM layout"; b="04_spatial.R (mock slide)";;
    results/05_*|results/tcga_ov_*|figs/tcga_ov_*) o="REAL"; b="05_signature_vs_celltype.R";;
    results/tcga_celltype_scores_long.csv|figs/tcga_celltype_score_distributions.png) o="REAL"; b="distribution_plot.R";;
    results/tcga_mean_celltype_proportions.csv|figs/tcga_celltype_pie.png) o="REAL"; b="pie_chart_plot.R";;
    results/tcga_celltype_survival_scatter_data.csv|figs/tcga_celltype_vs_survival_scatter.png|figs/tcga_celltype_kde.png) o="REAL"; b="**no script in repo** (proportion_survival_plot.R is 0 bytes)";;
    figs/03b_benchmark.png|figs/03_fraction_scatter.png) o="SIM"; b="initial commit 9fb71b7 (synthetic pipeline)";;
    results/validation_log.tsv) o="MIXED"; b="appended by every stage; includes 2026-07-21 rows of unclear profile";;
    *) o="?"; b="?";;
  esac
  echo "| \`$p\` | $hs | $t | $o | $b |"
done
echo
echo "\`truth/\` (not requested, listed for completeness) holds the synthetic planted truth: **SIM**."
echo

# ---------------------------------------------------------------- 6
echo "## 6. Real-data results"
echo
echo "Analyses in the plan that have **no output on disk**: every PFI/CRS contrast, P1, P2, the positive control, limma::fry gene-set tests, speckle::propeller, liana-py, and the plan's TCGA-OV validation. All results on disk date from 2026-08-05 to 2026-08-19, **before** the plan was committed (2026-09-23)."
echo
cat "$SCR/section6.md"
echo

# ---------------------------------------------------------------- 7
echo "## 7. PASS / FAIL / FALLBACK lines and run-log tails"
echo
echo "### 7a. results/validation_log.tsv (verbatim)"
echo '```'; cat results/validation_log.tsv; echo '```'
echo
echo "### 7b. Every PASS/FAIL/FALLBACK line in any log"
echo '```'
grep -Hn '\[PASS\]\|\[FAIL\]\|\[FALLBACK\]\|PASSED\|FAILED' results/0*_run_log.txt site_load_log.txt verify_site_log.txt install_log.txt raw_data/GSE180661/download_log.txt 2>/dev/null
echo '```'
echo
for f in results/01_run_log.txt results/02_run_log.txt results/03_run_log.txt results/04_run_log.txt results/05_run_log.txt site_load_log.txt verify_site_log.txt raw_data/GSE180661/download_log.txt install_log.txt; do
  echo "### 7c. Last 30 lines: \`$f\` ($(stat -c '%y' "$f" | cut -d. -f1))"
  echo '```'; tail -n 30 "$f"; echo '```'
  echo
done

# ---------------------------------------------------------------- 8
echo "## 8. Environment"
echo
echo "### sessionInfo files on disk"
echo
echo "| File | Next to | R | Platform | Key packages |"
echo "|---|---|---|---|---|"
for f in sessionInfo.txt results/tcga_ov_sessionInfo_stage05.txt; do
  rv=$(head -1 "$f"); pl=$(sed -n 2p "$f" | sed 's/Platform: //')
  pk=$(grep -oE '\b(GSVA|SingleCellExperiment|SummarizedExperiment|scater|scuttle|splatter|survival|ggplot2|edgeR|limma|Seurat|infercnv|GenomeInfoDb|matrixStats|jsonlite|yaml|data\.table|Matrix)_[0-9.-]+' "$f" | sort -u | tr '\n' ' ')
  echo "| \`$f\` | $(dirname "$f") | $rv | $pl | $pk |"
done
echo
echo "- Root \`sessionInfo.txt\` (modified $(stat -c '%y' sessionInfo.txt | cut -d. -f1)) is from **Windows 11, R 4.3.1 ucrt, time zone Asia/Dubai**, with \`splatter\` attached: it belongs to the synthetic pipeline, not to any real-data run on this server."
echo "- \`results/tcga_ov_sessionInfo_stage05.txt\` is the only sessionInfo for a real-data run (Stage 05, Linux conda, R 4.3.1)."
echo "- **No sessionInfo exists for real-data stages 00, 01, 02, 03, 04.**"
echo "- \`renv.lock\`: R $(grep -A1 '"R": {' renv.lock | grep -o '[0-9.]\+' | head -1), Bioconductor $(sed -n '12p' renv.lock | grep -o '[0-9.]\+')."
echo
echo "### Current hgsoc_pipeline env"
echo '```'
R --version | head -1
Rscript -e 'for (p in c("SingleCellExperiment","infercnv","GSVA","survival","edgeR","limma","scater","scran","matrixStats","BayesPrism","spacexr","speckle")) cat(sprintf("%-22s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e) "NOT INSTALLED")))' 2>/dev/null
echo '```'
echo
echo "### hgsoc_ccc env (\`/userhome/syeed/envs/hgsoc_ccc\`)"
echo '```'
"$CCC_PY" --version 2>&1
"$CCC_PY" -I -c "import importlib.metadata as m; [print(f'{p:10s}', m.version(p)) for p in ['liana','scanpy','anndata','mudata','pandas','numpy']]" 2>&1
echo '```'
echo "Matches \`env_hgsoc_ccc.yml\` / \`requirements_ccc.txt\` (python 3.13.2, liana 1.10.0). The env exists but **no liana output exists anywhere on disk**."
echo

# ---------------------------------------------------------------- 9
cat "$SCR/section9.md"
} > "$OUT"

# ---------------------------------------------------------------- bundle
B=writeup_bundle
mkdir -p "$B/figs_real" "$B/figs_simulation" "$B/tables" "$B/logs" "$B/env" "$B/plan"
for f in figs/*.png; do
  case "$f" in figs/03b_benchmark.png|figs/03_fraction_scatter.png) cp -p "$f" "$B/figs_simulation/";; *) cp -p "$f" "$B/figs_real/";; esac
done
cp -p results/infercnv/infercnv_subclusters.png "$B/figs_real/01_infercnv_subclusters_heatmap.png"
cp -p results/01_metrics.csv results/01_confusion_malignant.csv results/02_de_results.csv results/02_plausibility_check.csv \
      results/02_README.md results/03_metrics.csv results/03_estimated_fractions.csv results/04_metrics.csv \
      results/tcga_ov_signature_celltype_correlations.csv results/tcga_ov_corr_rho_GSVA.csv results/tcga_ov_corr_rho_ssGSEA.csv \
      results/tcga_mean_celltype_proportions.csv results/validation_log.tsv "$B/tables/"
cp -p results/0*_run_log.txt site_load_log.txt verify_site_log.txt "$B/logs/"
cp -p sessionInfo.txt "$B/env/sessionInfo_root_windows_synthetic.txt"
cp -p results/tcga_ov_sessionInfo_stage05.txt env_hgsoc_ccc.yml requirements_ccc.txt env_hgsoc_pipeline.yml "$B/env/"
cp -p docs/ANALYSIS_PLAN.md DEVIATIONS.md metadata/gse165897_labels.csv "$B/plan/"
cp -p "$OUT" "$B/"
rm -f "$B/writeup_bundle.zip"
(cd "$B" && zip -qr writeup_bundle.zip . -x writeup_bundle.zip)
echo "wrote $REPO/$OUT and $REPO/$B/writeup_bundle.zip"
