# Getting Started with sigRecon

## The signature recontextualization problem

Given a gene signature $`X`$ describing a genetic or chemical
perturbation in one biological context (a cell line, a tissue), what is
the corresponding gene signature of that *same* perturbation in a
*different* context? Any method that takes a ranked list of genes from a
source context and outputs a ranked list of genes for a target context
is performing **signature recontextualization**.

`sigrecon` provides three baseline recontextualization methods –
projection-based scoring
([`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md)),
network propagation
([`netProp()`](https://montilab.github.io/sigrecon/reference/netProp.md)),
and a DESeq2-based mean method (`recontextualize(method = "mean")`) –
plus benchmarking tools
([`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md))
to evaluate how well a recontextualized signature recovers the true
signature in the target context.

This vignette walks through all three methods end to end on a small,
bundled real-data example, requiring no downloads or external setup.

## Installation

``` r

BiocManager::install("montilab/sigrecon", dependencies = TRUE)
```

``` r

library(sigrecon)
library(SummarizedExperiment)
```

## Demo Data

`sigrecon` bundles small, real-data recontextualization demos for each
of the datasets used in its benchmarking study (SciPlex, Tahoe,
DrugMatrix, Perturb-seq) – see
[`?demo_datasets`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
for all four. This vignette uses the SciPlex demo data:

- **`demo_sciplex_sigs`** – SciPlex K562 signatures (top differentially
  expressed genes per drug), treated here as **source** signatures, as
  if defined in a different biological context than the one we’re
  recontextualizing into.
- **`demo_sciplex_se`** – SciPlex A549 expression data, representing the
  **target** context we want to recontextualize the source signatures
  into.
- **`demo_sciplex_true_sigs`** – SciPlex A549 signatures for the same
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

Note `demo_sciplex_se` contains a *disjoint* set of drug profiles from
those in `demo_sciplex_sigs` and `demo_sciplex_true_sigs` – so nothing
about the specific perturbations you’re trying to recontextualize leaks
into the expression data used to predict new signatures.

## Method 1: Projection-based

[`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md)
predicts a signature in a new context by scoring each sample against the
source signature (via GSVA or eigengene scoring) and then ranking genes
by their correlation with that per-sample score. The highest-correlated
genes become the recontextualized signature.

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

## Method 2: Network Propagation

Network propagation works differently:
[`netProp()`](https://montilab.github.io/sigrecon/reference/netProp.md)
first learns a gene co-expression network from the target-context
expression data, then propagates the source signature across that
network (via random walk with restart) to find the genes most
“reachable” from the seed genes – these become the recontextualized
signature.

``` r

recon_netprop <- netProp(demo_sciplex_se, seeds = demo_sciplex_sigs, sig = "rwr")
#> Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use =
#> "all.obs", : bicor: zero MAD in variable 'x'. Pearson correlation was used for
#> individual columns with zero (or missing) MAD.
#> Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use =
#> "all.obs", : bicor: zero MAD in variable 'y'. Pearson correlation was used for
#> individual columns with zero (or missing) MAD.
#>    Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k. max.k.
#> 1      1    0.119 -0.993        0.33900 156.0000  1.42e+02 265.00
#> 2      2    0.919 -1.640        0.93300  39.3000  3.02e+01 112.00
#> 3      3    0.959 -1.580        0.96400  13.3000  7.90e+00  59.50
#> 4      4    0.962 -1.520        0.96900   5.5500  2.34e+00  37.90
#> 5      5    0.874 -1.590        0.90900   2.7000  7.98e-01  26.50
#> 6      6    0.918 -1.500        0.95300   1.4700  2.87e-01  19.90
#> 7      7    0.918 -1.460        0.94600   0.8820  1.14e-01  15.70
#> 8      8    0.935 -1.390        0.94600   0.5680  4.71e-02  12.80
#> 9      9    0.907 -1.350        0.88500   0.3900  2.01e-02  10.90
#> 10    10    0.933 -1.320        0.91500   0.2820  9.13e-03   9.43
#> 11    12    0.879 -1.260        0.87200   0.1690  2.00e-03   7.51
#> 12    14    0.226 -1.520        0.03280   0.1170  4.82e-04   6.78
#> 13    16    0.213 -1.810       -0.00155   0.0889  1.20e-04   6.28
#> 14    18    0.217 -1.720        0.00168   0.0723  3.10e-05   5.87
#> 15    20    0.175 -1.430        0.03040   0.0614  8.29e-06   5.51
#> Optimal power selected: 2
#> Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use =
#> "all.obs", : bicor: zero MAD in variable 'x'. Pearson correlation was used for
#> individual columns with zero (or missing) MAD.
recon_netprop[["Panobinostat (LBH589)"]][1:5]
#> [1] "ENSG00000166206" "ENSG00000144036" "ENSG00000130449" "ENSG00000248905"
#> [5] "ENSG00000067225"
```

[`netProp()`](https://montilab.github.io/sigrecon/reference/netProp.md)
is a one-call wrapper around two lower-level functions,
[`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
(builds the network) and
[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
(propagates across it). Building the network is the expensive step, so
if you need to propagate many different seed sets across the *same*
network, build it once and pass it back in via `ig` instead of calling
[`netProp()`](https://montilab.github.io/sigrecon/reference/netProp.md)
repeatedly:

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

recon_netprop2 <- netProp(ig = network, seeds = demo_sciplex_sigs, sig = "rwr")
identical(recon_netprop, recon_netprop2)
#> [1] FALSE
```

See the [Network Propagation in
Depth](https://montilab.github.io/sigrecon/articles/2-network-propagation-details.md)
vignette for more on
[`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)/[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)’s
individual parameters.

## Method 3: Mean (DESeq2)

The mean method takes a different approach entirely: rather than
recontextualizing a signature you already have, it *derives* one
directly from target-context expression data, by running DESeq2 between
perturbed and control samples and keeping the top upregulated genes. It
needs raw counts and sample metadata identifying which samples are
perturbed vs. control (and by what), so it doesn’t fit the SciPlex demo
signatures used above – here’s a minimal, fabricated example showing the
shape of the inputs:

``` r

set.seed(42)
n_genes <- 30

# g1/g2 are upregulated 4x in pertA samples; everything else is unperturbed noise.
ctrl_base <- rnbinom(n_genes, mu = 200, size = 10)
pert_base <- ctrl_base
pert_base[1:2] <- pert_base[1:2] * 4

sample_counts <- function(base, n) {
  sapply(seq_len(n), function(i) rnbinom(length(base), mu = pmax(base, 1), size = 10))
}

counts <- cbind(sample_counts(ctrl_base, 4), sample_counts(pert_base, 3))
rownames(counts) <- paste0("g", seq_len(n_genes))
colnames(counts) <- c(paste0("ctrl", 1:4), paste0("pertA", 1:3))

coldata <- S4Vectors::DataFrame(
  perturbation = c(rep("control", 4), rep("pertA", 3)),
  condition    = c(rep("Control", 4), rep("Perturbed", 3)),
  row.names    = colnames(counts)
)

se_mean <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)

recon_mean <- recontextualize(
  method = "mean",
  se = se_mean,
  perturbation_col = "perturbation",
  condition_col = "condition",
  alpha = 1,
  limit = 2
)
recon_mean
#> $pertA
#> [1] "g2" "g1"
```

`g1`/`g2` are the two genes consistently higher in `pertA` samples than
in controls, so they come out as the derived signature for `pertA`.

## One entry point: `recontextualize()`

Each method above
([`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md),
[`netProp()`](https://montilab.github.io/sigrecon/reference/netProp.md),
the mean method) can also be run through a single dispatcher,
[`recontextualize()`](https://montilab.github.io/sigrecon/reference/recontextualize.md),
by setting `method = "projectCor"`, `"networkProp"`, or `"mean"`. This
is convenient when method choice is itself a parameter of a larger
workflow (e.g. looping over multiple methods in a benchmarking script)
rather than something fixed at call-site.

``` r

identical(
  recontextualize(method = "projectCor", se = demo_sciplex_se, sigs = demo_sciplex_sigs, score = "gsva"),
  projectCor(demo_sciplex_se, demo_sciplex_sigs, score = "gsva")
)
#> ℹ No assay name provided; using default assay 'logcounts'
#> ℹ No assay name provided; using default assay 'logcounts'
#> [1] TRUE

identical(
  recontextualize(method = "networkProp", se = demo_sciplex_se, seeds = demo_sciplex_sigs, sig = "rwr"),
  netProp(demo_sciplex_se, seeds = demo_sciplex_sigs, sig = "rwr")
)
#> Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use =
#> "all.obs", : bicor: zero MAD in variable 'x'. Pearson correlation was used for
#> individual columns with zero (or missing) MAD.
#> Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use =
#> "all.obs", : bicor: zero MAD in variable 'y'. Pearson correlation was used for
#> individual columns with zero (or missing) MAD.
#>    Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k. max.k.
#> 1      1    0.119 -0.993        0.33900 156.0000  1.42e+02 265.00
#> 2      2    0.919 -1.640        0.93300  39.3000  3.02e+01 112.00
#> 3      3    0.959 -1.580        0.96400  13.3000  7.90e+00  59.50
#> 4      4    0.962 -1.520        0.96900   5.5500  2.34e+00  37.90
#> 5      5    0.874 -1.590        0.90900   2.7000  7.98e-01  26.50
#> 6      6    0.918 -1.500        0.95300   1.4700  2.87e-01  19.90
#> 7      7    0.918 -1.460        0.94600   0.8820  1.14e-01  15.70
#> 8      8    0.935 -1.390        0.94600   0.5680  4.71e-02  12.80
#> 9      9    0.907 -1.350        0.88500   0.3900  2.01e-02  10.90
#> 10    10    0.933 -1.320        0.91500   0.2820  9.13e-03   9.43
#> 11    12    0.879 -1.260        0.87200   0.1690  2.00e-03   7.51
#> 12    14    0.226 -1.520        0.03280   0.1170  4.82e-04   6.78
#> 13    16    0.213 -1.810       -0.00155   0.0889  1.20e-04   6.28
#> 14    18    0.217 -1.720        0.00168   0.0723  3.10e-05   5.87
#> 15    20    0.175 -1.430        0.03040   0.0614  8.29e-06   5.51
#> Optimal power selected: 2
#> Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use =
#> "all.obs", : bicor: zero MAD in variable 'x'. Pearson correlation was used for
#> individual columns with zero (or missing) MAD.
#> Using weighted graph
#> Reached Convergence. Iteration step: 16
#> Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use =
#> "all.obs", : bicor: zero MAD in variable 'x'. Pearson correlation was used for
#> individual columns with zero (or missing) MAD.
#> Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use =
#> "all.obs", : bicor: zero MAD in variable 'y'. Pearson correlation was used for
#> individual columns with zero (or missing) MAD.
#>    Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k. max.k.
#> 1      1    0.119 -0.993        0.33900 156.0000  1.42e+02 265.00
#> 2      2    0.919 -1.640        0.93300  39.3000  3.02e+01 112.00
#> 3      3    0.959 -1.580        0.96400  13.3000  7.90e+00  59.50
#> 4      4    0.962 -1.520        0.96900   5.5500  2.34e+00  37.90
#> 5      5    0.874 -1.590        0.90900   2.7000  7.98e-01  26.50
#> 6      6    0.918 -1.500        0.95300   1.4700  2.87e-01  19.90
#> 7      7    0.918 -1.460        0.94600   0.8820  1.14e-01  15.70
#> 8      8    0.935 -1.390        0.94600   0.5680  4.71e-02  12.80
#> 9      9    0.907 -1.350        0.88500   0.3900  2.01e-02  10.90
#> 10    10    0.933 -1.320        0.91500   0.2820  9.13e-03   9.43
#> 11    12    0.879 -1.260        0.87200   0.1690  2.00e-03   7.51
#> 12    14    0.226 -1.520        0.03280   0.1170  4.82e-04   6.78
#> 13    16    0.213 -1.810       -0.00155   0.0889  1.20e-04   6.28
#> 14    18    0.217 -1.720        0.00168   0.0723  3.10e-05   5.87
#> 15    20    0.175 -1.430        0.03040   0.0614  8.29e-06   5.51
#> Optimal power selected: 2
#> Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use =
#> "all.obs", : bicor: zero MAD in variable 'x'. Pearson correlation was used for
#> individual columns with zero (or missing) MAD.
#> Using weighted graph
#> Reached Convergence. Iteration step: 16
#> [1] FALSE

identical(
  recontextualize(
    method = "mean", se = se_mean,
    perturbation_col = "perturbation", condition_col = "condition",
    alpha = 1, limit = 2
  ),
  recon_mean
)
#> converting counts to integer mode
#> -- note: fitType='parametric', but the dispersion trend was not well captured by the
#>    function: y = a/x + b, and a local regression fit was automatically substituted.
#>    specify fitType='local' or 'mean' to avoid this message next time.
#> [1] TRUE
```

For `method = "networkProp"`,
[`recontextualize()`](https://montilab.github.io/sigrecon/reference/recontextualize.md)
also accepts a pre-built `ig` (just like
[`netProp()`](https://montilab.github.io/sigrecon/reference/netProp.md))
to skip network construction and reuse a cached network:

``` r

identical(
  recontextualize(method = "networkProp", ig = network, seeds = demo_sciplex_sigs, sig = "rwr"),
  netProp(ig = network, seeds = demo_sciplex_sigs, sig = "rwr")
)
#> [1] FALSE
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
the K562 signature unmodified. This is a realistic outcome,
recontextualization does not always outperform the no-change baseline
for every method/metric/dataset combination (see the paper for more
details on benchmarking results).

[`sigrecon::paired_eval_table()`](https://montilab.github.io/sigrecon/reference/paired_eval_table.md)
formalizes the no-change-vs-recontextualized comparison shown here,
including paired significance testing, for larger benchmarking
workflows.

## Next Steps

- The [Network Propagation in
  Depth](https://montilab.github.io/sigrecon/articles/2-network-propagation-details.md)
  vignette for
  [`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)’s
  options
- The
  [BYOD](https://montilab.github.io/sigrecon/articles/3-bring-your-own-signature.md)
  vignette for bringing your own perturbational data and signatures.
- The [Benchmarking
  Philosophy](https://montilab.github.io/sigrecon/articles/4-benchmarking-and-splits.md)
  vignette for how we benchmark methods by taking into account target
  data availability.
