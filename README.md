<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# sigRecon: Signature Recontextualization

The signature recontextualization problem describes a simple goal in computational biology:
Given a gene signature X of a genetic or chemical perturbation in model organism Y, what is the corresponding gene signature of the same perturbation in model organism Z?

Any method that is able to input a ranked list of genes from one biological context and output another ranked list of genes for another context performs this task of 'signature recontextualization'.

This repository contains benchmarking tasks and data for evaluating signature recontextualization (as reported in [sigrecon_benchmarking](https://github.com/montilab/sigRecon_Benchmarking)), plus implementations of two of our own methods: projection-based scoring (`projectCor()`) and network propagation (`network_sig()`).

## Installation

Note: BiocManager installer is used to handle dependences. This package is currently not hosted on Bioconductor.
```r
BiocManager::install("montilab/sigrecon", dependencies = TRUE)
```

## Quick Start

The example below runs entirely on a small, bundled real-data example: `demo_sciplex_sigs` (source-context signatures), `demo_sciplex_se` (target-context expression), and `demo_sciplex_true_sigs` (target-context ground truth). Demo data is bundled for each dataset in the benchmarking study (Perturb-seq, SciPlex, DrugMatrix, Tahoe); full pseudobulk datasets are on Zenodo.

```r
library(sigrecon)

# Recontextualize with projectCor
recon_projectcor <- projectCor(demo_sciplex_se, demo_sciplex_sigs, score = "gsva")

# Or with network propagation
target_expr <- t(SummarizedExperiment::assay(demo_sciplex_se))
network <- wgcna.adj(target_expr, cor.type = "signed hybrid", diag_zero = TRUE, igraph = TRUE)
recon_netprop <- network_sig(network, seeds = demo_sciplex_sigs, sig = "rwr")

# Benchmark against the true target-context signature
eval_df <- sig_eval_table(
  source_sigs = demo_sciplex_sigs,
  pred_sigs   = recon_projectcor,
  true_sigs   = demo_sciplex_true_sigs,
  source      = "sciplex_k562",
  target      = "sciplex_a549"
)
head(eval_df[, c("gene", "jacc", "NES", "padj")])
```
