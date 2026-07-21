# Network Propagation in Depth

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

The [Getting
Started](https://montilab.github.io/sigrecon/articles/getting-started.md)
vignette uses
[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
with its defaults. This vignette goes deeper: how
[`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
builds the network
[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
operates on, the two different algorithms behind `sig = "corr"`
vs. `sig = "rwr"`, and what each of
[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)’s
options (`p`, `avg_p`, `bootstrap`, `n_bootstraps`, `limit`) actually
controls.

We’ll use the bundled SciPlex demo data throughout.

``` r

dim(demo_sciplex_se)
#> [1] 1025   48
length(demo_sciplex_sigs)
#> [1] 23
```

## Building the network with `wgcna.adj()`

[`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
builds a weighted gene co-expression network from an expression matrix
(samples as rows, genes as columns – note this is **transposed**
relative to `SummarizedExperiment`’s gene-rows/sample-columns
convention):

``` r

target_expr <- t(assay(demo_sciplex_se))

network <- wgcna.adj(
  target_expr,
  cor.type  = "signed hybrid", # correlation sign handling -- see below
  diag_zero = TRUE,            # zero out self-loops
  igraph    = TRUE             # return an igraph object, not a raw matrix
)
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

network
#> IGRAPH 636095a UNW- 1025 348313 -- 
#> + attr: name (v/c), weight (e/n)
#> + edges from 636095a (vertex names):
#>  [1] ENSG00000160963--ENSG00000259124 ENSG00000160963--ENSG00000108821
#>  [3] ENSG00000160963--ENSG00000216863 ENSG00000160963--ENSG00000100027
#>  [5] ENSG00000160963--ENSG00000183570 ENSG00000160963--ENSG00000169282
#>  [7] ENSG00000160963--ENSG00000179314 ENSG00000160963--ENSG00000123104
#>  [9] ENSG00000160963--ENSG00000143344 ENSG00000160963--ENSG00000236830
#> [11] ENSG00000160963--ENSG00000158486 ENSG00000160963--ENSG00000197959
#> [13] ENSG00000160963--ENSG00000164330 ENSG00000160963--ENSG00000214595
#> [15] ENSG00000160963--ENSG00000010030 ENSG00000160963--ENSG00000184005
#> + ... omitted several edges
```

A few of
[`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)’s
arguments worth knowing about:

- **`cor.type`**: `"unsigned"` treats positive and negative correlations
  as equally strong connections; `"signed"` only connects genes that
  move *together*; `"signed hybrid"` (used here) is a common middle
  ground – negative correlations are down-weighted rather than either
  kept at full strength or discarded entirely.
- **`cor.fn`**: `"bicor"` (biweight midcorrelation, the default) is more
  robust to outliers than Pearson (`"cor"`) – relevant for pseudobulk
  expression data with a handful of replicates per condition.
- **`min.sft`/`beta`/`powers`**: WGCNA networks raise the correlation
  matrix to a soft-thresholding power to approximate a scale-free
  topology.
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
  automatically searches `powers` for the smallest one reaching
  `min.sft` fit; pass `beta` directly to skip that search if you already
  know the power you want.
- **`diag_zero`**: whether to zero the diagonal (self-correlation)
  before returning – almost always `TRUE`, since a gene being perfectly
  correlated with itself is not informative for propagation.

## `sig = "corr"` vs. `sig = "rwr"`: two different algorithms

[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
supports two mechanically distinct ways of finding a seed signature’s
neighbors in the network:

- **`sig = "corr"`**: directly uses the network’s edge weights. For each
  candidate gene, it averages that gene’s correlation to every seed
  gene, then keeps the top `limit` genes by that average. This is a
  single lookup – no iteration, no walking across the network – so it’s
  fast, but it only “sees” genes with a direct (first-degree) connection
  to the seeds.
- **`sig = "rwr"`** (random walk with restart): simulates a walker that
  starts at the seed genes and, at each step, either moves to a
  network-neighbor or teleports back to a seed gene (with probability
  `p`, see below). Run to convergence, this produces a **stationary
  probability** for every gene in the network – how “reachable” it is
  from the seeds, accounting for indirect (multi-hop) paths. The top
  genes by stationary probability become the recontextualized signature.

``` r

one_seed <- demo_sciplex_sigs["Panobinostat (LBH589)"]

sig_corr <- network_sig(network, seeds = one_seed, sig = "corr", limit = 30)
sig_rwr  <- network_sig(network, seeds = one_seed, sig = "rwr", limit = 30)
#> Using weighted graph
#> Reached Convergence. Iteration step: 29

# How much do the two methods agree on which genes are "close" to the seed?
length(intersect(sig_corr[[1]], sig_rwr[[1]]))
#> [1] 15
jaccard(sig_corr[[1]], sig_rwr[[1]])
#> [1] 0.3333333
```

The two methods substantially disagree here – `sig = "rwr"` can surface
genes several hops away that `sig = "corr"` never considers, since
`"corr"` is restricted to genes directly connected to a seed gene by an
edge. Neither is “more correct” in general; which one recovers the true
target-context signature better is an empirical question this package’s
benchmarking tools
([`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md))
are built to answer per-dataset.

## `sig = "rwr"` options in depth

### `p`: the restart probability

`p` controls how far the walk tends to travel from the seed genes before
teleporting back. Low `p` lets the walk wander further across the
network (more indirect, multi-hop genes get picked up); high `p` keeps
it close to the seeds (results converge toward the seeds’ immediate
neighborhood, similar to `sig = "corr"`).

``` r

sig_p_low  <- network_sig(network, seeds = one_seed, sig = "rwr", p = 0.05, limit = 30)
#> Using weighted graph
#> Reached Convergence. Iteration step: 47
sig_p_high <- network_sig(network, seeds = one_seed, sig = "rwr", p = 0.9, limit = 30)
#> Using weighted graph
#> Reached Convergence. Iteration step: 5

jaccard(sig_p_low[[1]], sig_p_high[[1]])
#> [1] 0.01694915
```

### `avg_p`: ensembling over a range of restart values

Picking a single “right” `p` for a given network is not always obvious.
Setting `avg_p = TRUE` instead runs the random walk at `avg_p_length`
restart values evenly spaced across the `avg_p_vals` range and averages
the resulting stationary probabilities – trading a single,
possibly-arbitrary choice of `p` for a smoother estimate less sensitive
to that choice.

``` r

sig_avg_p <- network_sig(
  network,
  seeds = one_seed,
  sig = "rwr",
  avg_p = TRUE,
  avg_p_vals = c(0.05, 0.5),
  avg_p_length = 5,
  limit = 30
)
#> Using weighted graph
#> Reached Convergence. Iteration step: 47
#> Reached Convergence. Iteration step: 22
#> Reached Convergence. Iteration step: 16
#> Reached Convergence. Iteration step: 13
#> Reached Convergence. Iteration step: 11
length(sig_avg_p[[1]])
#> [1] 30
```

### `bootstrap`/`n_bootstraps`: significance-based selection instead of top-`limit`

By default (`bootstrap = FALSE`),
[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
just keeps the top `limit` genes by stationary probability – an
arbitrary cutoff. Setting `bootstrap = TRUE` instead runs the *same*
random walk from `n_bootstraps` randomly-chosen gene sets (matched in
size to the real seed set) to build an empirical null distribution, then
keeps only genes whose real stationary probability exceeds the 99th
percentile of that null – a significance-based cutoff rather than a
fixed count. This means the resulting signature can be smaller (or
larger) than `limit`, and costs roughly `n_bootstraps` extra random
walks per seed set, so it’s noticeably slower.

``` r

sig_bootstrap <- network_sig(
  network,
  seeds = one_seed,
  sig = "rwr",
  bootstrap = TRUE,
  n_bootstraps = 30, # kept small here for vignette speed; real analyses may want more
  limit = 30 # unused when bootstrap = TRUE and p != 1, but harmless to leave set
)
#> Using weighted graph
#> Reached Convergence. Iteration step: 29
#> Reached Convergence. Iteration step: 56
length(sig_bootstrap[[1]])
#> [1] 54
```

## Under the hood: `random_walk()` and seed matrices

`network_sig(sig = "rwr")` is itself built on the lower-level, exported
[`random_walk()`](https://montilab.github.io/sigrecon/reference/random_walk.md)
function, which takes a **seed matrix** – a (genes x seed-sets) matrix
marking which genes are seeds for which signature – rather than a named
list of gene vectors. Building one by hand shows exactly what
[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
constructs internally:

``` r

gene_names <- igraph::V(network)$name
seed_genes <- intersect(one_seed[[1]], gene_names)

seed_mat <- matrix(
  0, nrow = length(gene_names), ncol = 1,
  dimnames = list(gene_names, names(one_seed))
)
seed_mat[seed_genes, 1] <- 1

raw_probs <- random_walk(network, seed_mat, restart = 0.1)
#> Using weighted graph
#> Reached Convergence. Iteration step: 29
dim(raw_probs)
#> [1] 1025    1

# The top genes by stationary probability match extract_sig_mat()'s logic
# inside network_sig() -- this is the un-thresholded, raw output.
head(sort(raw_probs[, 1], decreasing = TRUE), 5)
#> ENSG00000151746 ENSG00000170776 ENSG00000187079 ENSG00000196935 ENSG00000154380 
#>     0.005362367     0.004703069     0.004574624     0.004490471     0.004369237
```

[`random_walk()`](https://montilab.github.io/sigrecon/reference/random_walk.md)
returns raw stationary probabilities for every gene, with no top-`limit`
cutoff or bootstrap-based thresholding applied – that gene selection
step is what
[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
adds on top.

## Next steps

- The [Getting
  Started](https://montilab.github.io/sigrecon/articles/getting-started.md)
  vignette for the full recontextualization + benchmarking workflow
- [`?network_sig`](https://montilab.github.io/sigrecon/reference/network_sig.md),
  [`?random_walk`](https://montilab.github.io/sigrecon/reference/random_walk.md),
  [`?wgcna.adj`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
  for full parameter references
- [`?demo_datasets`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  for the other three bundled demo triples
