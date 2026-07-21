# Bringing Your Own Signature and Expression Data

``` r

library(sigrecon)
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
```

The [Getting
Started](https://montilab.github.io/sigrecon/articles/getting-started.md)
vignette uses `sigrecon`’s bundled demo data. This vignette instead
walks through taking **your own** differential expression (DE) results
and expression data all the way through
[`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md)/[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
and
[`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md).

## 1. From a DE table to a signature

`sigrecon` represents a signature as a named list, one element per
perturbation, each a list with `up` (top significant upregulated genes)
and `up_full` (all genes, ranked).
[`sig_filter_fn()`](https://montilab.github.io/sigrecon/reference/sig_filter_fn.md)
builds this from a DE results table.

### DESeq2-shaped input

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

### Seurat `FindMarkers()`-shaped input

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

**The lesson**:
[`sig_filter_fn()`](https://montilab.github.io/sigrecon/reference/sig_filter_fn.md)
doesn’t care which DE tool produced your table – DESeq2,
`limma::topTable()`, Seurat’s `FindMarkers()` (any test, including
MAST), or something else entirely – as long as you tell it which columns
hold the perturbation label, log2 fold-change, p-value, and gene ID via
`pert_col`/`log2fc_col`/`pval_col`/`geneid_col`. Mixing DE tools across
datasets in the same analysis (e.g. DESeq2 for one cell line, MAST for
another) just means passing different arguments per call – there’s no
need to reconcile the tables into one schema first.

## 2. Gene ID namespace harmonization

A second, easy-to-miss issue: your DE table’s gene IDs need to be in the
**same namespace** as whatever you’re comparing against (a target
signature, or your expression data’s rownames) – both Ensembl IDs, or
both HGNC gene symbols, not a mix.
[`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md)/[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
match genes by exact string, so a namespace mismatch doesn’t error – it
silently produces empty or near-empty results.

If your DE table has gene symbols but you need Ensembl IDs (or vice
versa), join against a symbol/ID mapping table before calling
[`sig_filter_fn()`](https://montilab.github.io/sigrecon/reference/sig_filter_fn.md).
A common source for this mapping is an annotation package (e.g.
`org.Hs.eg.db`) or, if you built your expression data from a `Seurat`
object, its own feature metadata:

``` r

# Seurat objects often carry a symbol <-> Ensembl mapping in feature metadata
id_map <- seurat_obj@assays$RNA@meta.features |>
  tibble::rownames_to_column("ensembl_id") |>
  dplyr::select(ensembl_id, feature_name)

diff_table <- diff_table |>
  dplyr::left_join(id_map, by = c("gene_symbol" = "feature_name"))
```

Some symbols map to more than one Ensembl ID or vice versa; inspect
`id_map` for duplicates before joining, and decide how to resolve them
(e.g. keep the most-variable feature) rather than letting the join
silently fan out rows.

## 3. Recontextualizing your signature

Once you have a signature in the `up`/`up_full` format, the rest of the
workflow is the same as the [Getting
Started](https://montilab.github.io/sigrecon/articles/getting-started.md)
vignette – just with your own expression data in place of the bundled
demo. Building a small simulated `SummarizedExperiment` over the same
gene universe here to keep this vignette self-contained:

``` r

n_samples <- 20
expr <- matrix(
  rnorm(length(genes) * n_samples, mean = 5), nrow = length(genes), ncol = n_samples,
  dimnames = list(genes, paste0("sample", 1:n_samples))
)

my_se <- SummarizedExperiment(
  assays  = list(logcounts = expr),
  colData = DataFrame(sample_id = colnames(expr))
)
```

[`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md)/[`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)’s
`sigs`/`seeds` arguments want just the `up` gene vectors, not the full
`up`/`up_full` list:

``` r

my_sigs_up <- lapply(deseq2_sigs, function(x) x$up)

recon_projectcor <- projectCor(my_se, my_sigs_up, score = "gsva")
#> ℹ No assay name provided; using default assay 'logcounts'
#> ! Some gene sets have size one. Consider setting minSize > 1

my_mat <- t(assay(my_se))
my_network <- wgcna.adj(my_mat, cor.type = "signed hybrid", diag_zero = TRUE, igraph = TRUE)
#>    Power SFT.R.sq slope truncated.R.sq  mean.k. median.k.   max.k.
#> 1      1    0.243  1.61         0.1380 7.39e+00  7.35e+00 8.680000
#> 2      2    0.417 -1.73         0.6080 2.11e+00  2.06e+00 3.040000
#> 3      3    0.197 -7.55         0.0306 7.38e-01  7.05e-01 1.320000
#> 4      4    0.188 -5.36         0.0331 2.96e-01  2.58e-01 0.640000
#> 5      5    0.141 -3.92        -0.1030 1.32e-01  1.07e-01 0.331000
#> 6      6    0.193 -3.91         0.0434 6.34e-02  4.89e-02 0.178000
#> 7      7    0.156 -3.11        -0.0831 3.25e-02  2.14e-02 0.109000
#> 8      8    0.207 -3.20         0.0449 1.76e-02  9.50e-03 0.071500
#> 9      9    0.232 -4.92         0.0136 9.91e-03  4.28e-03 0.048200
#> 10    10    0.316 -5.45         0.1240 5.79e-03  2.00e-03 0.033100
#> 11    12    0.335 -4.56         0.1790 2.14e-03  4.84e-04 0.016000
#> 12    14    0.368 -5.15         0.1880 8.60e-04  1.22e-04 0.007900
#> 13    16    0.297 -4.71         0.1770 3.65e-04  3.00e-05 0.003930
#> 14    18    0.338 -4.44         0.2250 1.61e-04  7.57e-06 0.001960
#> 15    20    0.249 -3.53         0.2380 7.32e-05  1.95e-06 0.000979
#> Using the following power: 6
recon_netprop <- network_sig(my_network, seeds = my_sigs_up, sig = "rwr")
#> Using weighted graph
#> Reached Convergence. Iteration step: 46
```

## 4. Benchmarking

[`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md)’s
`true_sigs` argument expects the full `up`/`up_full` shape (it needs the
full ranked list to compute rank-based enrichment), so pass the
[`sig_filter_fn()`](https://montilab.github.io/sigrecon/reference/sig_filter_fn.md)
output directly rather than the `up`-only list used above:

``` r

eval_df <- sig_eval_table(
  source_sigs = my_sigs_up,
  pred_sigs   = recon_projectcor,
  true_sigs   = seurat_sigs,
  source      = "my_source",
  target      = "my_target"
)

head(eval_df[, c("gene", "jacc", "NES", "padj")])
#>    gene jacc      NES      padj
#> 1 drugA    0 1.158982 0.4194831
#> 2 drugB    0 1.057345 0.4194831
```

(This example uses `seurat_sigs` as a stand-in “ground truth” purely to
demonstrate the call shape – in a real analysis, `true_sigs` would be an
independently-derived signature for the actual target context you’re
recontextualizing into, not another simulated signature over the same
random data.)

## Next steps

- [`?sig_filter_fn`](https://montilab.github.io/sigrecon/reference/sig_filter_fn.md)
  for the full parameter reference
- The [Getting
  Started](https://montilab.github.io/sigrecon/articles/getting-started.md)
  vignette for a walkthrough using real bundled data, including a
  no-change baseline comparison
- [`?demo_datasets`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  for real source/target signature pairs across four datasets, useful as
  a reference for what well-formed input looks like
