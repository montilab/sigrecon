# Find Common Genes with Highest Median Absolute Deviation (MAD) Across ExpressionSets

This function identifies common genes with the highest MAD across
multiple ExpressionSet objects, up to a specified limit.

## Usage

``` r
common_mad_genes(esets, limit = 2500, parallel = FALSE, filter_zero = FALSE)
```

## Arguments

- esets:

  A list of ExpressionSet objects to compare.

- limit:

  An integer specifying the maximum number of common MAD genes to
  return. Default is 2500.

- parallel:

  A logical indicating whether to use parallel processing. Default is
  FALSE.

- filter_zero:

  A logical indicating whether to filter out genes with zero variance.
  Default is FALSE.

## Value

A character vector of common gene names with highest MAD.

## Details

The function ranks genes in each ExpressionSet by their MAD, then
iteratively selects genes that are present in all ExpressionSets until
reaching the specified limit or exhausting all common genes.
