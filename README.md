<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# Signature Recontextualization

The signature recontextualization problem describes a simple goal in computational biology:
Given gene signature X of a genetic or chemical perturbation in model organism Y, what is the corresponding gene signature of the same perturbation in model organism Z?

Any method (gene-regulatory network or deep-learning based) that is able to input a ranked list of genes from one biological context and output another ranked list of genes for another context performs this task of 'signature recontextualization'.

This repository contains benchmarking tasks and data for evaluating signature recontextualization, plus implementations of two of our own methods: projection-based scoring (`projectCor()`) and network propagation (`wgcna.adj()` + `network_sig()`).

## Installation

```r
BiocManager::install("montilab/sigrecon", dependencies = TRUE)
```

## Quick start

The examples below run entirely on a small, bundled real-data example. `demo_sciplex_sigs` represent differentially expressed gene signatures (DEGs) for drug perturbations of K562 cells, treated here as the "source" biological context. `demo_sciplex_se` is expression data for A549 cells representing the "target" context. `demo_sciplex_true_sigs` represent DEGs of the same drug perturbations in A549, used as ground truth to benchmark recontextualization quality. We provide demo data for each of the datasets (perturb-seq, sciplex, drugmatrix, tahoe) described in the benchmarking study. Full pseudobulk datasets can be found on zenodo.

### 1. Recontextualize with `projectCor()`

`projectCor()` reconstructs a signature in a new context by scoring samples against the source signature (via GSVA or eigengenes) and ranking genes by their correlation with that score.

```r
library(sigrecon)

recon_projectcor <- projectCor(demo_sciplex_se, demo_sciplex_sigs, score = "gsva")
```

### 2. Recontextualize with network propagation

`wgcna.adj()` learns a gene co-expression network from the target-context expression data; `network_sig()` then propagates the source signature across that network (random walk with restart) to find its target-context analog.

```r
target_expr <- t(SummarizedExperiment::assay(demo_sciplex_se))
network <- wgcna.adj(target_expr, cor.type = "signed hybrid", diag_zero = TRUE, igraph = TRUE)

recon_netprop <- network_sig(network, seeds = demo_sciplex_sigs, sig = "rwr")
```

### 3. Benchmark against ground truth

`sig_eval_table()` compares a reconstructed signature against both the unmodified source signature and the true target-context signature, reporting Jaccard overlap and rank-based enrichment (NES via fgsea).

```r
eval_df <- sig_eval_table(
  source_sigs = demo_sciplex_sigs,
  pred_sigs   = recon_projectcor,
  true_sigs   = demo_sciplex_true_sigs,
  source      = "sciplex_k562",
  target      = "sciplex_a549"
)

head(eval_df[, c("gene", "jacc", "NES", "padj")])
```

## More demo data

`demo_sciplex_*` is one of four real-data demo triples bundled with the package, each following the same `demo_<dataset>_se` / `demo_<dataset>_sigs` / `demo_<dataset>_true_sigs` shape — see `?demo_datasets` for details on all four (SciPlex, Tahoe, DrugMatrix, Perturb-seq).

For the full-size datasets these demos are subsetted from (and others), see `list_datasets()` and `get_dataset()`:

```r
list_datasets()
tahoe_nci_h23 <- get_dataset("tahoe.nci_h23")
```

## Bring your own signature or expression data

`projectCor()` and `network_sig()` work on any expression data and any gene signature, not just the bundled demo data — the demo objects above are just a fast way to see the functions run. This section covers what shape your own data needs to be in.

### Signature format

Both functions expect a **named list of character vectors of gene IDs** — one element per perturbation/condition, e.g.:

```r
my_sigs <- list(
  drugA = c("GENE1", "GENE7", "GENE12", ...),
  drugB = c("GENE3", "GENE9", ...)
)
```

If you're starting from a differential expression results table instead (e.g. output from `DESeq2::results()`, `limma::topTable()`, or Seurat's `FindMarkers()`), `sig_filter_fn()` converts it into this format, selecting the top upregulated genes per perturbation:

```r
# diff_table: one row per gene per perturbation, with a perturbation-label
# column, a log2 fold-change column, a p-value column, and a gene-ID column.
# Column names below match sig_filter_fn()'s defaults ("product",
# "avg_log2FC", "p_val_adj", "ensembl_id"); pass pert_col/log2fc_col/
# pval_col/geneid_col to match your own table's column names instead
# (Seurat::FindMarkers() output already uses "avg_log2FC"/"p_val_adj").
my_sigs <- sig_filter_fn(
  diff_table,
  perts = c("drugA", "drugB"),
  limit = 100
)

# my_sigs$drugA$up      -- top 100 significantly upregulated genes
# my_sigs$drugA$up_full -- all genes, ranked by log2FC * -log10(padj)
```

`sig_eval_table()`'s `true_sigs` argument expects this full `list(up = ..., up_full = ...)` shape; `projectCor()`/`network_sig()`'s `sigs`/`seeds` arguments just want the plain gene vectors (`lapply(my_sigs, function(x) x$up)`).

### Expression data format

| Function | Expected input | Orientation | Notes |
|---|---|---|---|
| `projectCor(se, sigs, score)` | `SummarizedExperiment` | genes as rows, samples as columns (standard Bioconductor convention) | `stopifnot(is(se, "SummarizedExperiment"))` — a matrix alone will error |
| `wgcna.adj(mat, ...)` (feeds `network_sig()`) | plain numeric matrix | **samples as rows, genes as columns** (transposed from `SummarizedExperiment` convention) | |

In both cases, **row/column gene identifiers must be in the same namespace as your signature's gene IDs** (e.g. both Ensembl IDs, or both HGNC symbols) — `projectCor()`/`network_sig()` match genes by exact string, so a namespace mismatch silently produces empty or near-empty results rather than an error.

A couple of practical gotchas learned from building this package's own demo datasets:
- **Normalize raw counts first.** `wgcna.adj()`'s correlation-based network construction behaves poorly on raw counts (skewed distributions, genes with zero variance producing `NA` correlations that make `igraph::graph_from_adjacency_matrix()` error outright). Log2-CPM (or similar) normalization first avoids this. If your data is already normalized/log-transformed (e.g. microarray intensities), skip this step.
- **Drop zero-variance genes** before building a network, for the same reason: `apply(mat, 2, var) > 0` (recall `mat` is samples × genes here).
- **Starting from a `Seurat` object?** Extract expression and wrap it as a `SummarizedExperiment`:
  ```r
  library(Seurat)
  library(SummarizedExperiment)

  expr_mat <- as.matrix(GetAssayData(seurat_obj, layer = "data"))  # normalized, genes x samples
  my_se <- SummarizedExperiment(
    assays = list(logcounts = expr_mat),
    colData = DataFrame(seurat_obj@meta.data)
  )
  ```

### Putting it together

```r
recon_projectcor <- projectCor(my_se, lapply(my_sigs, function(x) x$up), score = "gsva")

my_mat <- t(SummarizedExperiment::assay(my_se))
my_network <- wgcna.adj(my_mat, cor.type = "signed hybrid", diag_zero = TRUE, igraph = TRUE)
recon_netprop <- network_sig(my_network, seeds = lapply(my_sigs, function(x) x$up), sig = "rwr")
```