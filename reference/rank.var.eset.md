# Rank Genes in an ExpressionSet/SummarizedExperiment by Variability

This function ranks genes in an ExpressionSet/SummarizedExperiment
object based on their variability, as measured by a specified function
(default is median absolute deviation).

## Usage

``` r
rank.var.eset(eset, fn = mad, filter_zero = FALSE)
```

## Arguments

- eset:

  An ExpressionSet/SummarizedExperiment object containing gene
  expression data.

- fn:

  A function to measure variability. Default is `mad` (median absolute
  deviation).

- filter_zero:

  A logical indicating whether to filter out genes with zero variance.
  Default is FALSE.

## Value

A pvector object containing ranked gene names.

## Details

The function calculates the variability of each gene using the specified
function (default: mad). Genes are then sorted in descending order of
variability. If `filter_zero` is TRUE, genes with zero variance are
removed before ranking.
