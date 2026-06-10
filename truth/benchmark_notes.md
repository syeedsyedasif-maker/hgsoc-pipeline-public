# Comparative deconvolution benchmark (03b): notes

## Hypothesis (Hippen 2023)
Naive reference-based methods (MuSiC, with NNLS as a simple baseline) score well on pseudobulk, where
the bulk is built from the same cells as the single-cell reference. They are expected to degrade when
the bulk is an independent draw: different patients or platform, a different composition, and cell
types that are present in the bulk but missing from the reference. Bayesian methods of the BayesPrism
family are expected to hold up better.

## What was observed
Profile and scale: medium (1500 genes, 2500 cells in 2 batches, 80 bulk samples per condition).
Mean Pearson per recoverable cell type, pseudobulk then independent-draw, with the change:

```
     method  pseudobulk  independent-draw  degradation
       NNLS       0.805             0.451        0.354
      MuSiC       0.871             0.782        0.089
 InstaPrism       0.917             0.930       -0.013
 BayesPrism          NA                NA           NA   (did not complete, see below)
```

The naive methods dropped under the mismatch: NNLS by 0.354 (Pearson 0.805 to 0.451, RMSE 0.166 to
0.320), MuSiC by 0.089. The Bayesian method did not drop (InstaPrism 0.917 to 0.930). Mean naive drop
0.221 versus Bayesian -0.013. This matches the pattern described in Hippen 2023.

## Platform-mismatch setting (stated for transparency)
The independent-draw bulk is splatter Batch 2 with a per-gene batch effect set to
`independent_batch_facLoc = 0.70` and `independent_batch_facScale = 0.50` (the `benchmark` block in
`config.yml`). This is a substantial shift, chosen to represent a realistic scRNA-versus-bulk platform
gap. The separation between methods depends on this choice: with a weak mismatch the naive methods do
not degrade and all methods look similar. An earlier run at `facLoc = 0.30` showed no degradation
(NNLS even scored slightly higher on the independent draw). The result is therefore conditional on a
substantial mismatch, and that choice is recorded here rather than hidden.

## Engines and scale
- Reference downsampled to 150 cells per type; 200 marker genes selected.
- InstaPrism is the engine behind every Bayesian number above. It is a documented fast approximation
  of BayesPrism.
- True BayesPrism was attempted once for confirmation. It did not complete within the 900 second
  per-condition budget at medium scale (see the runtime table below). Every "BayesPrism-class" number
  here is therefore InstaPrism, not BayesPrism.
- Per-method time budget: 900 seconds. Anything over is logged as "did not complete" rather than left
  to hang, which is itself a practical finding.

```
     method              bulk  mean_pearson  mean_rmse  runtime_sec  status
       NNLS        pseudobulk         0.805      0.166          0.1  ok
       NNLS  independent-draw         0.451      0.320          0.0  ok
      MuSiC        pseudobulk         0.871      0.114          3.1  ok
      MuSiC  independent-draw         0.782      0.132          2.7  ok
 InstaPrism        pseudobulk         0.917      0.100          0.5  ok
 InstaPrism  independent-draw         0.930      0.116          0.3  ok
 BayesPrism        pseudobulk            NA         NA        900.4  did not complete (> 900s)
 BayesPrism  independent-draw            NA         NA        900.3  did not complete (> 900s)
```

## Design choice: batches, not a fresh different-seed simulation
The independent target is splatter Batch 2. It is a separate sample with its own per-gene technical
profile (the batch effect, standing in for a different patient or platform), plus a different
composition and extra reference dropouts. A fresh simulation with a different seed would re-randomise
which gene means which cell type. The reference would then be meaningless and every method would score
near zero, which is not a useful comparison. Batches keep gene identity fixed, so the contrast is
interpretable.

## Caveat
This validates machinery on simulated data, not biology. It reproduces the shape of the Hippen 2023
robustness contrast under controlled conditions. It is not evidence about real tissue.
