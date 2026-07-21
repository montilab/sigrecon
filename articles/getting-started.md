# Getting Started with sigrecon

## The signature recontextualization problem

Given a gene signature $`X`$ describing a genetic or chemical
perturbation in one biological context (a cell line, a tissue), what is
the corresponding gene signature of that *same* perturbation in a
*different* context? Any method that takes a ranked list of genes from a
source context and outputs a ranked list of genes for a target context
is performing **signature recontextualization**.

`sigrecon` provides two baseline recontextualization methods –
projection-based scoring
([`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md))
and network propagation
([`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md) +
[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md))
– plus benchmarking tools
([`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md))
to evaluate how well a recontextualized signature recovers the true
signature in the target context.

This vignette walks through both methods end to end on a small, bundled
real-data example, requiring no downloads or external setup.

## Installation

``` r

BiocManager::install("montilab/sigrecon", dependencies = TRUE)
```

``` r

library(sigrecon)
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
```

## The demo data

`sigrecon` bundles small, real-data recontextualization demos for each
of the datasets used in its benchmarking study (SciPlex, Tahoe,
DrugMatrix, Perturb-seq) – see
[`?demo_datasets`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
for all four. This vignette uses the SciPlex triple:

- **`demo_sciplex_sigs`** – real SciPlex K562 signatures (top
  differentially expressed genes per drug), treated here as **source**
  signatures, as if defined in a different biological context than the
  one we’re recontextualizing into.
- **`demo_sciplex_se`** – real SciPlex A549 expression data,
  representing the **target** context we want to recontextualize the
  source signatures into.
- **`demo_sciplex_true_sigs`** – the real A549 signatures for the same
  drugs, used as **ground truth** to check how well recontextualization
  worked.

``` r

dim(demo_sciplex_se)
#> [1] 1025   48
length(demo_sciplex_sigs)
#> [1] 23
demo_sciplex_sigs[["Panobinostat (LBH589)"]][1:5]
#> [1] "ENSG00000160963" "ENSG00000259124" "ENSG00000214456" "ENSG00000216863"
#> [5] "ENSG00000103196"
```

Both `demo_sciplex_se` and `demo_sciplex_true_sigs` are built from a set
of A549 drugs *disjoint* from the drugs actually being evaluated – so
nothing about the specific perturbations you’re trying to
recontextualize leaks into the expression data or network used to
reconstruct them.

## Method 1: `projectCor()`

[`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md)
reconstructs a signature in a new context by scoring each sample against
the source signature (via GSVA or eigengene scoring) and then ranking
genes by their correlation with that per-sample score. The
highest-correlated genes become the recontextualized signature.

``` r

recon_projectcor <- projectCor(demo_sciplex_se, demo_sciplex_sigs, score = "gsva")
#> ℹ No assay name provided; using default assay 'logcounts'

# One recontextualized signature per source signature
length(recon_projectcor)
#> [1] 23
recon_projectcor[["Panobinostat (LBH589)"]][1:5]
#> [1] "ENSG00000198804" "ENSG00000210082" "ENSG00000211459" "ENSG00000212907"
#> [5] "ENSG00000198899"
```

## Method 2: network propagation

Network propagation works differently:
[`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
first learns a gene co-expression network from the target-context
expression data, then
[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
propagates the source signature across that network (via random walk
with restart) to find the genes most “reachable” from the seed genes –
these become the recontextualized signature.

``` r

target_expr <- t(assay(demo_sciplex_se))
network <- wgcna.adj(target_expr, cor.type = "signed hybrid", diag_zero = TRUE, igraph = TRUE)
#> Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use =
#> "all.obs", : bicor: zero MAD in variable 'x'. Pearson correlation was used for
#> individual columns with zero (or missing) MAD.
#> Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use =
#> "all.obs", : bicor: zero MAD in variable 'y'. Pearson correlation was used for
#> individual columns with zero (or missing) MAD.
#>    Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k.  max.k.
#> 1      1   0.0248  0.649         0.9590 1.62e+02  1.54e+02 301.000
#> 2      2   0.4190 -1.570         0.8980 4.13e+01  3.49e+01 125.000
#> 3      3   0.8380 -1.910         0.9840 1.37e+01  9.72e+00  60.800
#> 4      4   0.8200 -2.030         0.9140 5.37e+00  3.11e+00  34.100
#> 5      5   0.8720 -1.960         0.9670 2.40e+00  1.09e+00  20.900
#> 6      6   0.8640 -1.920         0.9280 1.18e+00  4.20e-01  13.500
#> 7      7   0.8970 -1.860         0.9690 6.27e-01  1.72e-01   9.100
#> 8      8   0.9090 -1.810         0.9660 3.53e-01  7.62e-02   6.370
#> 9      9   0.3490 -2.660         0.3040 2.10e-01  3.50e-02   4.590
#> 10    10   0.3500 -2.510         0.2950 1.30e-01  1.63e-02   3.400
#> 11    12   0.3460 -2.330         0.2670 5.57e-02  3.88e-03   1.990
#> 12    14   0.9550 -1.450         0.9430 2.71e-02  1.00e-03   1.260
#> 13    16   0.2500 -1.790         0.0529 1.46e-02  2.60e-04   0.843
#> 14    18   0.2800 -2.360         0.0889 8.53e-03  6.76e-05   0.609
#> 15    20   0.2540 -2.070         0.0421 5.37e-03  1.88e-05   0.466
#> Optimal power selected: 5
#> Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use =
#> "all.obs", : bicor: zero MAD in variable 'x'. Pearson correlation was used for
#> individual columns with zero (or missing) MAD.

recon_netprop <- network_sig(network, seeds = demo_sciplex_sigs, sig = "rwr")
#> Using weighted graph
#> Reached Convergence. Iteration step: 57
recon_netprop[["Panobinostat (LBH589)"]][1:5]
#> [1] "ENSG00000151746" "ENSG00000170776" "ENSG00000187079" "ENSG00000196935"
#> [5] "ENSG00000154380"
```

## Benchmarking against the ground truth

[`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md)
compares a recontextualized signature against both the unmodified source
signature and the true target-context signature, reporting:

- **`jacc`** – Jaccard overlap between the recontextualized and true
  signature
- **`NES`** – rank-based enrichment (via `fgsea`) of the
  recontextualized signature within the true signature’s full ranked
  gene list
- **`displaced`/`kept`** – how many source signature genes were replaced
  vs. retained during recontextualization

``` r

eval_df <- sig_eval_table(
  source_sigs = demo_sciplex_sigs,
  pred_sigs   = recon_projectcor,
  true_sigs   = demo_sciplex_true_sigs,
  source      = "sciplex_k562",
  target      = "sciplex_a549"
)

head(eval_df[, c("gene", "jacc", "NES", "padj")])
#>                         gene       jacc       NES         padj
#> 1    Abexinostat (PCI-24781) 0.10497238 2.8044745 5.514516e-12
#> 2 Alvespimycin (17-DMAG) HCl 0.08571429 0.8068178 7.919192e-01
#> 3                      AR-42 0.07602339 2.3954446 3.518487e-07
#> 4        Belinostat (PXD101) 0.01769912 0.9169023 6.261664e-01
#> 5                   CUDC-101 0.05714286 2.4557251 1.055620e-06
#> 6                   CUDC-907 0.17058824 3.0029322 6.638770e-14
```

### Is recontextualization actually helping?

The natural baseline to compare against is doing *nothing* – using the
unmodified source signature directly as the “prediction.” Computing the
same evaluation with `pred_sigs = demo_sciplex_sigs` gives that
no-change baseline:

``` r

no_change_df <- sig_eval_table(
  source_sigs = demo_sciplex_sigs,
  pred_sigs   = demo_sciplex_sigs,
  true_sigs   = demo_sciplex_true_sigs,
  source      = "sciplex_k562",
  target      = "sciplex_a549"
)
#> Warning in fgsea_wrapper(ref = ref, data = data, scoreType = scoreType, : No
#> overlap between 'data' and 'ref', returning NULL

data.frame(
  method    = c("no_change", "projectCor"),
  mean_jacc = c(mean(no_change_df$jacc), mean(eval_df$jacc)),
  mean_NES  = c(mean(no_change_df$NES, na.rm = TRUE), mean(eval_df$NES, na.rm = TRUE))
)
#>       method  mean_jacc mean_NES
#> 1  no_change 0.07674496 1.747552
#> 2 projectCor 0.08377736 1.458653
```

On this demo, results are mixed:
[`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md)
improves average Jaccard overlap with the true A549 signature, but
average enrichment (NES) is actually slightly *lower* than just reusing
the K562 signature unmodified. This is a realistic outcome, not a bug –
recontextualization does not always outperform the no-change baseline
for every method/metric/dataset combination, and this demo’s disjoint,
sample-limited design (see
[`?demo_datasets`](https://montilab.github.io/sigrecon/reference/demo_datasets.md))
makes it a harder, noisier test than the full-scale datasets in the
package’s benchmarking study. Try
[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)’s
output (`recon_netprop`) in place of `recon_projectcor` above to compare
methods, or see the full results for all datasets and methods in the
package’s companion analysis repository.
[`sigrecon::paired_eval_table()`](https://montilab.github.io/sigrecon/reference/paired_eval_table.md)
formalizes the no-change-vs-recontextualized comparison shown here,
including paired significance testing, for larger benchmarking
workflows.

## Next steps

- [`?demo_datasets`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  for the other three bundled demo triples (Tahoe, DrugMatrix,
  Perturb-seq)
- [`?get_dataset`](https://montilab.github.io/sigrecon/reference/get_dataset.md)
  and
  [`?list_datasets`](https://montilab.github.io/sigrecon/reference/list_datasets.md)
  for the full-size datasets these demos are subsetted from
- [`?sig_filter_fn`](https://montilab.github.io/sigrecon/reference/sig_filter_fn.md)
  for converting your own differential expression results into the
  signature format used throughout this vignette
- The package README covers bringing your own signature and expression
  data in more detail, including gene-ID-namespace and normalization
  gotchas.
