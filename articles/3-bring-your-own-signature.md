# Bringing Your Own Signature and Expression Data

``` r

library(sigrecon)
library(SummarizedExperiment)
```

The [Getting
Started](https://montilab.github.io/sigrecon/articles/getting-started.md)
vignette uses `sigrecon`’s bundled demo data. This vignette instead
walks through taking **your own** differential expression (DE) results
and expression data all the way through
[`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md)/[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
and
[`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md).

## Extracting DEGs from a differential expression table

`sigrecon` represents a signature as a named list, one element per
perturbation, each a list with `up` (top significant upregulated genes)
and `up_full` (all genes, ranked).
[`sig_filter_fn()`](https://montilab.github.io/sigrecon/reference/sig_filter_fn.md)
builds this from a DE results table.

### DESeq2 Table

[`DESeq2::results()`](https://rdrr.io/pkg/DESeq2/man/results.html)
output uses `log2FoldChange` and `padj` column names, and typically has
gene IDs as rownames rather than a column – move them into a column
first. Here’s a small simulated example with two perturbations:

``` r

set.seed(1)
genes <- paste0("ENSG", sprintf("%06d", 1:40))

deseq2_table <- data.frame(
  ensembl_id     = rep(genes, times = 2),
  product        = rep(c("drugA", "drugB"), each = 40),
  log2FoldChange = c(rnorm(40, mean = 1), rnorm(40, mean = -1)),
  padj           = runif(80, 0, 0.2)
)

deseq2_sigs <- sig_filter_fn(
  deseq2_table,
  perts      = c("drugA", "drugB"),
  pert_col   = "product",
  log2fc_col = "log2FoldChange",
  pval_col   = "padj",
  geneid_col = "ensembl_id",
  limit      = 20
)
#> [1] "drugA"
#> [1] "drugB"

str(deseq2_sigs$drugA, max.level = 1)
#> List of 2
#>  $ up     : chr [1:5] "ENSG000007" "ENSG000022" "ENSG000033" "ENSG000026" ...
#>  $ up_full: chr [1:40] "ENSG000007" "ENSG000022" "ENSG000039" "ENSG000015" ...
```

### Seurat `FindMarkers()` Table

Seurat’s `FindMarkers()` output uses `avg_log2FC` and `p_val_adj` –
these happen to be
[`sig_filter_fn()`](https://montilab.github.io/sigrecon/reference/sig_filter_fn.md)’s
defaults, so no column-name arguments are needed if your table already
has a perturbation-label column named `"product"` and gene IDs as
rownames-turned-column named `"ensembl_id"` (otherwise, override
`pert_col`/`geneid_col` the same way):

``` r

seurat_table <- data.frame(
  ensembl_id = rep(genes, times = 2),
  product    = rep(c("drugA", "drugB"), each = 40),
  avg_log2FC = c(rnorm(40, mean = 1), rnorm(40, mean = -1)),
  p_val_adj  = runif(80, 0, 0.2)
)

seurat_sigs <- sig_filter_fn(seurat_table, perts = c("drugA", "drugB"), limit = 20)
#> [1] "drugA"
#> [1] "drugB"
str(seurat_sigs$drugA, max.level = 1)
#> List of 2
#>  $ up     : chr [1:6] "ENSG000027" "ENSG000002" "ENSG000032" "ENSG000009" ...
#>  $ up_full: chr [1:40] "ENSG000027" "ENSG000002" "ENSG000022" "ENSG000040" ...
```

[`sig_filter_fn()`](https://montilab.github.io/sigrecon/reference/sig_filter_fn.md)
can be used on any differential expression output (limma, DESeq2, Mast)
as long as you tell it which columns hold the perturbation label, log2
fold-change, p-value, and gene ID via
`pert_col`/`log2fc_col`/`pval_col`/`geneid_col`.

## Next steps

- The [Getting
  Started](https://montilab.github.io/sigrecon/articles/getting-started.md)
  vignette for a walkthrough using real bundled data, including a
  no-change baseline comparison
- [`?demo_datasets`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  for real source/target signature pairs across four datasets, useful as
  a reference for what well-formed input looks like
