# Create (gene x seed) prior matrix based on seed signatures.

This function creates a binary matrix representing seed genes in the
context of a graph.

## Usage

``` r
seed_matrix(ig, seeds, bootstrap = FALSE, n_bootstraps = 1000)
```

## Arguments

- ig:

  An igraph object representing the network that has the seed genes as
  vertices

- seeds:

  Either a single unnamed gene "TP53", a named list of genes, or a list
  of named lists of genes.

- bootstrap:

  A boolean specifying whether to use empirical distributions of
  stationary values to find significant genes.

- n_bootstraps:

  A numeric specifying the number of bootstraps to perform.

## Value

A binary matrix where rows represent genes in the graph and columns
represent seed sets.
