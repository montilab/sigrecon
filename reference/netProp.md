# Network-propagation based Recontextualization, from expression data to signature.

Learns a gene co-expression network from target-context expression data
with
[`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md),
then propagates seed signatures across it with
[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md).
This is the one-call convenience wrapper for the two-step
[`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md) +
[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
workflow; use
[`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)/[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
directly (or pass a pre-built `ig`) when you need to reuse the same
network across many propagation calls, since building the network is the
expensive step.

## Usage

``` r
netProp(
  se = NULL,
  seeds,
  ig = NULL,
  sig = c("rwr", "corr"),
  avg_p = FALSE,
  avg_p_vals = c(1e-04, 0.1),
  avg_p_length = 5,
  p = 0.1,
  bootstrap = FALSE,
  n_bootstraps = 1000,
  limit = 30,
  assay_name = NULL,
  nfeatures = NULL,
  min.sft = 0.85,
  beta = NULL,
  cores = 1,
  cor.fn = c("bicor", "cor"),
  cor.type = c("unsigned", "signed hybrid", "signed"),
  powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
  diag_zero = TRUE
)
```

## Arguments

- se:

  A SummarizedExperiment object containing target-context expression
  data. Not needed if `ig` is supplied.

- seeds:

  Either a single unnamed gene "TP53", a named list of genes, or a list
  of named lists of genes.

- ig:

  An optional pre-built igraph, e.g. from a previous
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
  call. If supplied, network construction is skipped and `se` is
  ignored.

- sig:

  A string specifying the type of network signature: random walk,
  correlation etc.

- avg_p:

  A boolean specifying whether to ensemble random walk results over a
  range of restart values

- avg_p_vals:

  A numeric vector specifying the start and end of a arithmetic sequence
  to explore restart values.

- avg_p_length:

  A numeric specifying how many values within `avg_p_vals` to include in
  the ensemble

- p:

  A numeric specifying the restart value for random walk, default=0.1

- bootstrap:

  A boolean specifying whether to use empirical distributions of
  stationary values to find significant genes.

- n_bootstraps:

  A numeric specifying the number of bootstraps to perform.

- limit:

  A numeric specifying the number of genes to be included in the network
  signature. Default is 30.

- assay_name:

  Assay name in `se` to use for network construction. Defaults to the
  first assay.

- nfeatures:

  Number of variable genes to use when learning the graph from `se`.
  Defaults to `min(10000, nrow(assay(se)))`.

- min.sft:

  Minimum scale-free topology fitting index used by
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md).

- beta:

  Optional soft-thresholding power passed to
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md).

- cores:

  Number of CPU cores passed to
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md).

- cor.fn:

  Correlation function passed to
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md).

- cor.type:

  Correlation network type passed to
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md).

- powers:

  Candidate power values passed to
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
  when `beta` is `NULL`.

- diag_zero:

  Whether to zero the diagonal of the learned adjacency matrix before
  converting it to igraph.

## Value

A named list of recontextualized gene signatures.
