# Vectorized fgsea wrapper for multiple gene set comparisons

Vectorized fgsea wrapper for multiple gene set comparisons

## Usage

``` r
v.fgsea(ref_vecs, data_vecs, scoreType = "std", eps = 1e-50, BPPARAM = NULL)
```

## Arguments

- ref_vecs:

  Named list of ranked gene symbol vectors (each is a full ranked list)

- data_vecs:

  Named list of gene sets to test (each is a pathway/gene set)

- scoreType:

  GSEA score type: "std" (default), "pos", or "neg"

- eps:

  Precision for p-value calculation (default: 1e-50)

- BPPARAM:

  A BiocParallelParam object specifying parallel backend. If NULL,
  automatically selects appropriate backend based on platform. See
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  for options.

## Value

A data frame with fgsea results for each ref-data pair

## Examples

``` r
library(BiocParallel)

# Multiple ranked lists (e.g., from different conditions)
ref_vecs <- list(
  condition1 = c("TP53", "MYC", "BRCA1", "KRAS", "EGFR", "PTEN", "AKT1"),
  condition2 = c("MYC", "KRAS", "TP53", "EGFR", "BRCA1", "PTEN", "AKT1")
)

# Gene sets to test in each condition
data_vecs <- list(
  condition1 = c("TP53", "BRCA1", "EGFR"),
  condition2 = c("TP53", "BRCA1", "EGFR")
)

# Sequential execution
results <- v.fgsea(ref_vecs, data_vecs)

# Parallel execution (Unix/Mac)
results <- v.fgsea(
  ref_vecs, data_vecs,
  BPPARAM = MulticoreParam(workers = 2)
)

# Parallel execution (Windows)
results <- v.fgsea(
  ref_vecs, data_vecs,
  BPPARAM = SnowParam(workers = 2)
)

# With reproducibility
results <- v.fgsea(
  ref_vecs, data_vecs,
  BPPARAM = MulticoreParam(workers = 2, RNGseed = 123)
)
```
