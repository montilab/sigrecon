<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# Signature Recontextualization

**Signature recontextualization**: given a gene signature for a perturbation in one biological context, predict the corresponding signature for that same perturbation in a different context. This package provides two baseline methods — projection-based scoring (`projectCor()`) and network propagation (`wgcna.adj()` + `network_sig()`) — plus benchmarking tools and data to evaluate recontextualization methods generally.

## Installation

```r
BiocManager::install("montilab/sigrecon", dependencies = TRUE)
```

## Quick start

Runs entirely on bundled real-data demo objects (`demo_sciplex_sigs`/`_se`/`_true_sigs`; one such triple exists per dataset — Perturb-seq, SciPlex, DrugMatrix, Tahoe):

```r
library(sigrecon)

# Recontextualize with projectCor() (projection-based scoring)
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

See `vignette("getting-started", package = "sigrecon")` for a full walkthrough, `?demo_datasets` for the other three demo triples, and `list_datasets()`/`get_dataset()` for the full-size source datasets (also on Zenodo).

## Bring your own data

`projectCor()`/`network_sig()` work on any data, not just the demos. Signatures are a named list of gene-ID vectors (`sig_filter_fn()` builds these from a DE table). Expression: `projectCor()` wants a `SummarizedExperiment` (genes × samples); `wgcna.adj()` wants a plain matrix, **transposed** (samples × genes) — normalize raw counts and drop zero-variance genes first. Gene IDs must share a namespace with your signature. A full worked example is planned as a vignette; see `?sig_filter_fn`, `?projectCor`, `?wgcna.adj` in the meantime.