# GSE165897 diagnostics summary

The diagnostics ran on 28 October 2026 and this summary was written on 29 October 2026. They are read-only checks on the Stage 01 output (`data/sce_clustered.rds`) and the raw GEO files. No pipeline stage was rerun.
The source tables are in `results/diagnostics/`, made by `results/diagnostics/diagnostics_2026-10-28.R`. That folder is gitignored, so every number below is copied from those CSVs.

## 1. inferCNV reference set and settings

- **Where it runs:** Stage 01 (`01_qc_cluster.R`) calls `run_cnv()` → `call_malignant_infercnv()` in `R/cnv_fallback.R`. The run on 5 August 2026 used inferCNV as the primary tool, not the fallback.
- **Reference cells:** every author cell type except malignant. That is endothelial (79), fibroblast (CAF-1/2/3, 5,660), macrophage (6,917) and T cell (20,664), **33,320 cells in all**. They are **pooled across all 11 patients and both timepoints**: the inferCNV annotation file groups cells by cell type, not by patient.
- **`infercnv::run()` settings:** `cutoff = 0.1`, `cluster_by_groups = TRUE`, `denoise = TRUE`, **`HMM = FALSE`**, `num_threads = 1`, `no_plot = TRUE`. All other arguments are package defaults: `analysis_mode = "subclusters"` with Leiden partitioning, and `sd_amplifier = 1.5`.
- **How a cell is called malignant.** The call is not inferCNV's own (there is no HMM). The pipeline scores each cell as the mean |denoised expression − reference mean| and calls it malignant above the reference mean + 3 SD. That threshold is set from the same reference cells it is then judged on. So the reference false-positive rate is low partly by design.

## 2. inferCNV calls vs author malignant labels

"Truth" here is the author annotation (EOC = malignant). Precision = TP / called malignant. Recall = TP / author malignant. Ref. FPR = FP / author non-malignant. ARI is computed with the pipeline's `adjusted_rand_index()` from `R/metrics.R`, the same function Stage 01 uses.

**Overall: precision 0.969, recall 0.639, ARI 0.643** (5,629 TP, 178 FP, 3,177 FN, 33,142 TN). This fails Stage 01's ARI ≥ 0.70 check; DEVIATIONS.md D1 records why that does not change the analysis plan.

| Patient | Timepoint | Cells | Author malignant | Called malignant | TP | FP | FN | TN | Precision | Recall | Ref. FPR | ARI |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ALL | ALL | 42126 | 8806 | 5807 | 5629 | 178 | 3177 | 33142 | 0.969 | 0.639 | 0.0053 | 0.643 |
| EOC1005 | post-NACT | 1026 | 418 | 330 | 329 | 1 | 89 | 607 | 0.997 | 0.787 | 0.0016 | 0.678 |
| EOC1005 | treatment-naive | 2798 | 487 | 432 | 412 | 20 | 75 | 2291 | 0.954 | 0.846 | 0.0087 | 0.835 |
| EOC136 | post-NACT | 3409 | 2427 | 2064 | 2063 | 1 | 364 | 981 | 1.000 | 0.850 | 0.0010 | 0.614 |
| EOC136 | treatment-naive | 2333 | 205 | 175 | 153 | 22 | 52 | 2106 | 0.874 | 0.746 | 0.0103 | 0.759 |
| EOC153 | post-NACT | 1863 | 116 | 57 | 47 | 10 | 69 | 1737 | 0.825 | 0.405 | 0.0057 | 0.499 |
| EOC153 | treatment-naive | 949 | 330 | 119 | 119 | 0 | 211 | 619 | 1.000 | 0.361 | 0.0000 | 0.270 |
| EOC227 | post-NACT | 1550 | 174 | 180 | 158 | 22 | 16 | 1354 | 0.878 | 0.908 | 0.0160 | 0.852 |
| EOC227 | treatment-naive | 1252 | 57 | 64 | 55 | 9 | 2 | 1186 | 0.859 | 0.965 | 0.0075 | 0.896 |
| EOC3 | post-NACT | 1715 | 223 | 61 | 54 | 7 | 169 | 1485 | 0.885 | 0.242 | 0.0047 | 0.302 |
| EOC3 | treatment-naive | 2465 | 203 | 19 | 16 | 3 | 187 | 2259 | 0.842 | 0.079 | 0.0013 | 0.121 |
| EOC349 | post-NACT | 1260 | 66 | 56 | 55 | 1 | 11 | 1193 | 0.982 | 0.833 | 0.0008 | 0.887 |
| EOC349 | treatment-naive | 1076 | 30 | 20 | 19 | 1 | 11 | 1045 | 0.950 | 0.633 | 0.0010 | 0.746 |
| EOC372 | post-NACT | 4321 | 131 | 100 | 80 | 20 | 51 | 4170 | 0.800 | 0.611 | 0.0048 | 0.672 |
| EOC372 | treatment-naive | 577 | 281 | 162 | 162 | 0 | 119 | 296 | 1.000 | 0.577 | 0.0000 | 0.344 |
| EOC443 | post-NACT | 3450 | 103 | 78 | 74 | 4 | 29 | 3343 | 0.949 | 0.718 | 0.0012 | 0.805 |
| EOC443 | treatment-naive | 1903 | 1086 | 771 | 771 | 0 | 315 | 817 | 1.000 | 0.710 | 0.0000 | 0.447 |
| EOC540 | post-NACT | 1665 | 18 | 8 | 8 | 0 | 10 | 1647 | 1.000 | 0.444 | 0.0000 | 0.609 |
| EOC540 | treatment-naive | 1692 | 202 | 169 | 134 | 35 | 68 | 1455 | 0.793 | 0.663 | 0.0235 | 0.636 |
| EOC733 | post-NACT | 2866 | 372 | 133 | 130 | 3 | 242 | 2491 | 0.977 | 0.349 | 0.0012 | 0.430 |
| EOC733 | treatment-naive | 1652 | 1550 | 629 | 629 | 0 | 921 | 102 | 1.000 | 0.406 | 0.0000 | -0.032 |
| EOC87 | post-NACT | 1188 | 191 | 93 | 88 | 5 | 103 | 992 | 0.946 | 0.461 | 0.0050 | 0.508 |
| EOC87 | treatment-naive | 1116 | 136 | 87 | 73 | 14 | 63 | 966 | 0.839 | 0.537 | 0.0143 | 0.566 |

How to read the table:
- **Recall is the weak point.** Precision is 0.79–1.00 in every sample, but recall ranges from 0.079 (EOC3 treatment-naive) to 0.965.
- **ARI misleads when one class dominates.** EOC733 treatment-naive has ARI −0.032 even though precision is 1.000, because 1,550 of its 1,652 cells are author-malignant. Report precision and recall for single samples, not ARI.

## 3. Cell-count reconciliation against Zhang et al. 2022

**Zhang et al. (2022) report 51,786 cells after QC: 8,806 malignant, 8,045 stromal and 34,935 immune.** The raw GEO files contain exactly that:
- `GSE165897_UMIcounts_HGSOC.tsv.gz` has 51,786 cell columns.
- `GSE165897_cellInfo_HGSOC.tsv.gz` has 51,786 rows.
- The cell-info group totals are EOC 8,806, Stromal 8,045 and Immune 34,935.

**The analysed object has 42,126 cells, so 9,660 are missing.** None were lost to QC; Stage 01 QC removed 0 cells. All 9,660 belong to cell subtypes that `00_load_real_data.R` drops on purpose, because the pipeline's fixed six-type vocabulary has no slot for them:
- **Immune, 7,354:** B cells 2,635, NK 1,744, mast cells 804, DC-1 622, plasma cells 586, DC-2 365, pDC 310, ILC 288.
- **Stromal, 2,306:** all mesothelial cells.

**No malignant and no fibroblast cells are missing.** All 8,806 EOC cells and all 5,660 CAF-1/2/3 cells are in the analysed object, as are all endothelial cells, macrophages and T cells. The analysed stromal total (5,739) is lower than Zhang's 8,045 only because mesothelial cells are excluded.

By author cell type:

| Author group | Author subtype | Raw GEO | Analysed | Missing | Pipeline label |
|---|---|---|---|---|---|
| EOC | EOC_C1 | 972 | 972 | 0 | malignant |
| EOC | EOC_C10 | 397 | 397 | 0 | malignant |
| EOC | EOC_C11 | 1530 | 1530 | 0 | malignant |
| EOC | EOC_C12 | 1225 | 1225 | 0 | malignant |
| EOC | EOC_C2 | 737 | 737 | 0 | malignant |
| EOC | EOC_C3 | 196 | 196 | 0 | malignant |
| EOC | EOC_C4 | 449 | 449 | 0 | malignant |
| EOC | EOC_C5 | 813 | 813 | 0 | malignant |
| EOC | EOC_C6 | 628 | 628 | 0 | malignant |
| EOC | EOC_C7 | 460 | 460 | 0 | malignant |
| EOC | EOC_C8 | 1267 | 1267 | 0 | malignant |
| EOC | EOC_C9 | 132 | 132 | 0 | malignant |
| Immune | B-cells | 2635 | 0 | 2635 | DROPPED by 00_load_real_data.R |
| Immune | NK | 1744 | 0 | 1744 | DROPPED by 00_load_real_data.R |
| Immune | Mast-cells | 804 | 0 | 804 | DROPPED by 00_load_real_data.R |
| Immune | DC-1 | 622 | 0 | 622 | DROPPED by 00_load_real_data.R |
| Immune | Plasma-cells | 586 | 0 | 586 | DROPPED by 00_load_real_data.R |
| Immune | DC-2 | 365 | 0 | 365 | DROPPED by 00_load_real_data.R |
| Immune | pDC | 310 | 0 | 310 | DROPPED by 00_load_real_data.R |
| Immune | ILC | 288 | 0 | 288 | DROPPED by 00_load_real_data.R |
| Immune | Macrophages | 6917 | 6917 | 0 | macrophage |
| Immune | T-cells | 20664 | 20664 | 0 | Tcell |
| Stromal | Mesothelial | 2306 | 0 | 2306 | DROPPED by 00_load_real_data.R |
| Stromal | CAF-1 | 2842 | 2842 | 0 | fibroblast |
| Stromal | CAF-2 | 1558 | 1558 | 0 | fibroblast |
| Stromal | CAF-3 | 1260 | 1260 | 0 | fibroblast |
| Stromal | Endothelial | 79 | 79 | 0 | endothelial |
| EOC | ALL | 8806 | 8806 | 0 |  |
| Immune | ALL | 34935 | 27581 | 7354 |  |
| Stromal | ALL | 8045 | 5739 | 2306 |  |
| TOTAL | ALL | 51786 | 42126 | 9660 |  |

By sample. Every sample loses only cells of the dropped subtypes:

| Sample | Patient | Timepoint | Raw GEO | Analysed | Missing | Missing by subtype |
|---|---|---|---|---|---|---|
| EOC1005_interval_Tumor | EOC1005 | post-NACT | 1288 | 1026 | 262 | B-cells=11; DC-1=50; DC-2=8; ILC=8; Mast-cells=16; Mesothelial=152; NK=11; pDC=4; Plasma-cells=2 |
| EOC1005_primary_Peritoneum | EOC1005 | treatment-naive | 3071 | 2798 | 273 | B-cells=76; DC-1=5; DC-2=59; ILC=8; Mast-cells=3; Mesothelial=105; NK=5; pDC=3; Plasma-cells=9 |
| EOC136_interval_Omentum | EOC136 | post-NACT | 3602 | 3409 | 193 | B-cells=7; DC-1=39; DC-2=3; ILC=1; Mast-cells=27; Mesothelial=98; NK=6; pDC=5; Plasma-cells=7 |
| EOC136_primary_Mesentery | EOC136 | treatment-naive | 2501 | 2333 | 168 | B-cells=14; DC-2=26; Mast-cells=28; Mesothelial=58; NK=15; pDC=11; Plasma-cells=16 |
| EOC153_interval_Omentum | EOC153 | post-NACT | 2401 | 1863 | 538 | B-cells=270; DC-1=21; DC-2=49; ILC=47; Mast-cells=35; Mesothelial=23; NK=47; pDC=33; Plasma-cells=13 |
| EOC153_primary_Omentum | EOC153 | treatment-naive | 1380 | 949 | 431 | B-cells=28; DC-1=14; DC-2=16; ILC=6; Mast-cells=7; Mesothelial=309; NK=9; pDC=32; Plasma-cells=10 |
| EOC227_interval_Omentum | EOC227 | post-NACT | 2015 | 1550 | 465 | B-cells=241; DC-2=13; ILC=9; Mast-cells=3; Mesothelial=13; NK=169; pDC=7; Plasma-cells=10 |
| EOC227_primary_Omentum | EOC227 | treatment-naive | 2178 | 1252 | 926 | B-cells=831; DC-2=20; ILC=3; Mast-cells=1; Mesothelial=23; NK=21; pDC=2; Plasma-cells=25 |
| EOC3_interval_Omentum | EOC3 | post-NACT | 2503 | 1715 | 788 | B-cells=52; DC-1=59; DC-2=19; ILC=24; Mast-cells=209; Mesothelial=131; NK=156; pDC=37; Plasma-cells=101 |
| EOC3_primary_Peritoneum | EOC3 | treatment-naive | 3142 | 2465 | 677 | B-cells=22; DC-1=70; DC-2=1; ILC=6; Mast-cells=5; Mesothelial=294; NK=233; pDC=27; Plasma-cells=19 |
| EOC349_interval_Omentum | EOC349 | post-NACT | 1471 | 1260 | 211 | B-cells=73; DC-2=98; ILC=25; Mast-cells=3; Mesothelial=3; NK=6; Plasma-cells=3 |
| EOC349_primary_Peritoneum | EOC349 | treatment-naive | 1483 | 1076 | 407 | B-cells=21; DC-1=16; DC-2=2; ILC=13; Mesothelial=329; NK=16; pDC=6; Plasma-cells=4 |
| EOC372_interval_Peritoneum | EOC372 | post-NACT | 4671 | 4321 | 350 | B-cells=47; DC-2=10; ILC=13; Mast-cells=30; Mesothelial=39; NK=10; pDC=3; Plasma-cells=198 |
| EOC372_primary_Peritoneum | EOC372 | treatment-naive | 711 | 577 | 134 | B-cells=5; DC-1=1; DC-2=2; Mesothelial=64; NK=5; pDC=1; Plasma-cells=56 |
| EOC443_interval_Omentum | EOC443 | post-NACT | 4463 | 3450 | 1013 | B-cells=195; DC-1=106; DC-2=5; ILC=93; Mast-cells=136; Mesothelial=55; NK=345; pDC=47; Plasma-cells=31 |
| EOC443_primary_Omentum | EOC443 | treatment-naive | 2122 | 1903 | 219 | B-cells=4; DC-1=26; DC-2=2; ILC=2; Mast-cells=5; Mesothelial=110; NK=17; pDC=15; Plasma-cells=38 |
| EOC540_interval_Omentum | EOC540 | post-NACT | 2043 | 1665 | 378 | B-cells=18; DC-1=12; DC-2=1; ILC=9; Mast-cells=212; Mesothelial=7; NK=113; Plasma-cells=6 |
| EOC540_primary_Omentum | EOC540 | treatment-naive | 1921 | 1692 | 229 | B-cells=7; DC-1=21; DC-2=5; ILC=12; Mast-cells=34; Mesothelial=42; NK=81; pDC=21; Plasma-cells=6 |
| EOC733_interval_Omentum | EOC733 | post-NACT | 4077 | 2866 | 1211 | B-cells=670; DC-1=130; DC-2=4; ILC=2; Mast-cells=12; Mesothelial=5; NK=361; pDC=20; Plasma-cells=7 |
| EOC733_primary_Peritoneum | EOC733 | treatment-naive | 1742 | 1652 | 90 | DC-1=35; DC-2=1; Mesothelial=41; NK=8; pDC=4; Plasma-cells=1 |
| EOC87_interval_Omentum | EOC87 | post-NACT | 1491 | 1188 | 303 | B-cells=12; DC-1=17; DC-2=2; ILC=5; Mast-cells=33; Mesothelial=150; NK=66; pDC=10; Plasma-cells=8 |
| EOC87_primary_Peritoneum | EOC87 | treatment-naive | 1510 | 1116 | 394 | B-cells=31; DC-2=19; ILC=2; Mast-cells=5; Mesothelial=255; NK=44; pDC=22; Plasma-cells=16 |

## 4. Cells per patient × timepoint (author labels)

These are the analysed cells per sample. Plan §3 sets a minimum of 20 cells per pseudobulk. Bold values with ⚠ are below it.

| Patient | Timepoint | Malignant (EOC) | CAF (CAF-1/2/3) |
|---|---|---|---|
| EOC1005 | post-NACT | 418 | 54 |
| EOC1005 | treatment-naive | 487 | 133 |
| EOC136 | post-NACT | 2427 | 300 |
| EOC136 | treatment-naive | 205 | 904 |
| EOC153 | post-NACT | 116 | 441 |
| EOC153 | treatment-naive | 330 | 91 |
| EOC227 | post-NACT | 174 | **13** ⚠ |
| EOC227 | treatment-naive | 57 | **7** ⚠ |
| EOC3 | post-NACT | 223 | 260 |
| EOC3 | treatment-naive | 203 | 57 |
| EOC349 | post-NACT | 66 | **1** ⚠ |
| EOC349 | treatment-naive | 30 | 25 |
| EOC372 | post-NACT | 131 | 2106 |
| EOC372 | treatment-naive | 281 | 77 |
| EOC443 | post-NACT | 103 | 197 |
| EOC443 | treatment-naive | 1086 | 305 |
| EOC540 | post-NACT | **18** ⚠ | 358 |
| EOC540 | treatment-naive | 202 | 133 |
| EOC733 | post-NACT | 372 | **15** ⚠ |
| EOC733 | treatment-naive | 1550 | **2** ⚠ |
| EOC87 | post-NACT | 191 | 141 |
| EOC87 | treatment-naive | 136 | 40 |

- **Malignant:** only EOC540 post-NACT is below 20 (18 cells). This is why Stage 02's paired analysis has 10 pairs, not 11. All 11 treatment-naive samples pass, so P2 can use every patient.
- **CAF:** five samples are below 20: EOC227 (post-NACT 13, treatment-naive 7), EOC349 (post-NACT 1) and EOC733 (post-NACT 15, treatment-naive 2). So P1 (CAFs, treatment-naive) can use at most 9 patients. Under the PFI labels that is 6 resistant vs 3 sensitive. Under CRS it is 2 vs 7, which is not estimable (DEVIATIONS.md D2).

## 5. R and package versions

- **Diagnostics (28 October 2026):** `hgsoc_pipeline` conda env on rbgo-server2. R 4.3.1, Bioconductor 3.18, SingleCellExperiment 1.24.0, SummarizedExperiment 1.32.0, data.table 1.17.8, Matrix 1.6-5. ARI used the pipeline's base-R `adjusted_rand_index()`.
- **Stage 01 inferCNV run (5 August 2026):** no sessionInfo was saved, so the inferCNV version used then is **not recorded**. The env now has infercnv 1.20.0.
- **Same env, now:** edgeR 4.0.16, limma 3.58.1, GSVA 1.50.0, speckle 1.2.0. speckle was installed on 28 October 2026, after the diagnostics ran; it was the only package added. `env_hgsoc_pipeline.yml` is the current export.
- **Load warnings:** many packages warn that they were built under R 4.3.2 or 4.3.3 while running under R 4.3.1.
