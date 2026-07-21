# Extracts a signature from a (gene x seed) matrix of stationary probability values. This is the recontextualized signature. If doing ks.test, you don't need to find the top_n. Just find ks.test(original, recontextualized ranking) before and after.

Extracts a signature from a (gene x seed) matrix of stationary
probability values. This is the recontextualized signature. If doing
ks.test, you don't need to find the top_n. Just find ks.test(original,
recontextualized ranking) before and after.

## Usage

``` r
extract_sig_mat(
  mat,
  bootstraps = NULL,
  sig_bins = NULL,
  percentile = 0.99,
  limit = 30
)
```

## Arguments

- mat:

  (n_gene, n_seed) matrix of stationary probability values from rwr_mat

- bootstraps:

  A (n_gene, n_bins\*n_bootstraps) matrix specifying results from
  bootstrapped random walks.

- sig_bins:

  A named list describing the length of each perturbation.

- percentile:

  A double between 0,1 indicating the proportion cutoff for
  bootstrap-based signature derivation.

- limit:

  Number of genes to keep in the output, or a vector of lengths. Default
  is 30.

## Value

A named list of genesets. Each list element is the recontextualized
signature for that seed.
