# Construct WGCNA Adjacency Matrix

This function constructs a weighted gene co-expression network adjacency
matrix using WGCNA.

## Usage

``` r
wgcna.adj(
  mat,
  min.sft = 0.85,
  beta = NULL,
  cores = 1,
  cor.fn = c("bicor", "cor"),
  cor.type = c("unsigned", "signed hybrid", "signed"),
  powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
  igraph = FALSE,
  diag_zero = FALSE
)
```

## Arguments

- mat:

  A matrix with rows as samples and columns as genes

- min.sft:

  Minimum scale-free topology fitting index R^2 to pick
  soft-thresholding power. Default is 0.85.

- beta:

  Soft-thresholding power. If NULL, it will be automatically selected.
  Default is NULL.

- cores:

  Number of CPU cores to use for parallel computing. Default is 1.

- cor.fn:

  Correlation function to use. Either "bicor" (biweight midcorrelation)
  or "cor" (Pearson correlation). Default is "bicor".

- cor.type:

  Type of correlation network. Options are "unsigned", "signed hybrid",
  or "signed". Default is "unsigned".

- powers:

  Vector of soft-thresholding powers to try. Default is c(seq(1, 10, by
  = 1), seq(12, 20, by = 2)).

- igraph:

  If TRUE, returns an igraph object instead of a matrix. Default is
  FALSE.

- diag_zero:

  If TRUE, sets the diagonal of the adjacency matrix to zero. Default is
  FALSE.

## Value

An adjacency matrix or an igraph object representing the gene
co-expression network.

## Details

This function performs the following steps:

1.  Prepares the expression data.

2.  Selects the soft-thresholding power (if not provided).

3.  Constructs the adjacency matrix using WGCNA::adjacency().

4.  Optionally converts the result to an igraph object.
