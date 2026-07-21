# Construct WGCNA Adjacency Matrix from Correlation Matrix

This function constructs a weighted gene co-expression network adjacency
matrix using WGCNA, starting from a pre-computed correlation matrix.

## Usage

``` r
wgcna.power(cor_mat, cores = 1, diag_zero = TRUE)
```

## Arguments

- cor_mat:

  A correlation matrix of gene expression data.

- cores:

  Number of CPU cores to use for parallel computing. Default is 1.

- diag_zero:

  If TRUE, sets the diagonal of the adjacency matrix to zero. Default is
  TRUE.

## Value

An adjacency matrix representing the gene co-expression network.

## Details

This function performs the following steps:

1.  Selects the optimal soft-thresholding power using
    WGCNA::pickSoftThreshold.fromSimilarity().

2.  Constructs the adjacency matrix using
    WGCNA::adjacency.fromSimilarity().

3.  Optionally sets the diagonal to zero and/or converts the result to
    an igraph object.
