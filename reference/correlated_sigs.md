# Recontextualize seed signatures with correlation based neighbors

Recontextualize seed signatures with correlation based neighbors

## Usage

``` r
correlated_sigs(corr_mat, seeds, limit = 30)
```

## Arguments

- corr_mat:

  Correlation Matrix

- seeds:

  Either a single unnamed gene "TP53", a named list of genes, or a list
  of named lists of genes.

- limit:

  Number of genes to keep in the output, or a vector of lengths. Default
  is 30.
