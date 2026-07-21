# Network-propagation based Recontextualization.

Network-propagation based Recontextualization.

## Usage

``` r
network_sig(
  ig,
  seeds,
  sig = c("corr", "rwr"),
  avg_p = FALSE,
  avg_p_vals = c(1e-04, 0.1),
  avg_p_length = 5,
  p = 0.1,
  bootstrap = FALSE,
  n_bootstraps = 1000,
  limit = 30
)
```

## Arguments

- ig:

  network given as an igraph

- seeds:

  Either a single unnamed gene "TP53", a named list of genes, or a list
  of named lists of genes.

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

## Value

vector of gene strings
