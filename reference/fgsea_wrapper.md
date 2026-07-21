# fgsea wrapper for gene symbol vectors

fgsea wrapper for gene symbol vectors

## Usage

``` r
fgsea_wrapper(ref, data, scoreType = "std", eps = 1e-50)
```

## Arguments

- ref:

  Vector of gene symbols representing the full ranked list (with
  weights)

- data:

  Vector of gene symbols representing the gene set/pathway to test

- scoreType:

  GSEA score type: "std" (default), "pos", or "neg"

- eps:

  Precision for p-value calculation (default: 1e-50)

## Value

A list with fgsea results (ES, NES, pval, padj, leadingEdge, size)

## Examples

``` r
# ref = all genes ranked by importance
ref <- c("TP53", "MYC", "BRCA1", "KRAS", "EGFR", "PTEN", "AKT1", "PIK3CA")
# data = gene set to test for enrichment
data <- c("TP53", "BRCA1", "EGFR")
result <- fgsea_wrapper(ref, data)
#> Error in fgsea_wrapper(ref, data): could not find function "fgsea_wrapper"
```
