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

`projectCor()` and `network_sig()` work on any `SummarizedExperiment` and any named list of gene-symbol/ID vectors, not just the bundled demo data — the demo objects above are just a fast way to see the functions run. `sig_filter_fn()` converts a differential expression results table (e.g. from `DESeq2`, `limma`, or Seurat's `FindMarkers()`) into the signature list format used throughout this package.