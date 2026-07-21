# Package index

## Recontextualization methods

Methods that take a source-context signature and reconstruct it in a
target context.

- [`recontextualize()`](https://montilab.github.io/sigrecon/reference/recontextualize.md)
  : Recontextualize signatures with one of either networkProp,
  projectCor, or mean.
- [`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md)
  : Projection based Recontextualization
- [`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
  : Network-propagation based Recontextualization.
- [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
  : Construct WGCNA Adjacency Matrix

## Benchmarking

Evaluating a recontextualized signature against source/true signatures.

- [`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md)
  : Evaluating recontextualization methods
- [`paired_eval_table()`](https://montilab.github.io/sigrecon/reference/paired_eval_table.md)
  : Evaluate recontextualization methods relative to no-change baseline
- [`jaccard()`](https://montilab.github.io/sigrecon/reference/jaccard.md)
  : Calculate Jaccard Similarity Between Two Sets
- [`v.jaccard()`](https://montilab.github.io/sigrecon/reference/v.jaccard.md)
  : Calculate Jaccard Similarity for Multiple Pairs of Sets
- [`v.fgsea()`](https://montilab.github.io/sigrecon/reference/v.fgsea.md)
  : Vectorized fgsea wrapper for multiple gene set comparisons
- [`fishers_meta_p()`](https://montilab.github.io/sigrecon/reference/fishers_meta_p.md)
  : Combine P-values Using Fisher's Method
- [`count_displaced_genes()`](https://montilab.github.io/sigrecon/reference/count_displaced_genes.md)
  : Count number of displaced seeds

## Datasets

Bundled demo data and on-demand access to the full-size datasets.

- [`demo_sciplex_se`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  [`demo_sciplex_sigs`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  [`demo_sciplex_true_sigs`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  [`demo_tahoe_se`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  [`demo_tahoe_sigs`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  [`demo_tahoe_true_sigs`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  [`demo_drugmatrix_se`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  [`demo_drugmatrix_sigs`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  [`demo_drugmatrix_true_sigs`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  [`demo_perturbseq_se`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  [`demo_perturbseq_sigs`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  [`demo_perturbseq_true_sigs`](https://montilab.github.io/sigrecon/reference/demo_datasets.md)
  : Demo datasets
- [`get_dataset()`](https://montilab.github.io/sigrecon/reference/get_dataset.md)
  : Fetch a bundled perturbational signature dataset
- [`list_datasets()`](https://montilab.github.io/sigrecon/reference/list_datasets.md)
  : List available on-demand datasets

## Utilities

- [`make_bpparam()`](https://montilab.github.io/sigrecon/reference/make_bpparam.md)
  : Create appropriate BiocParallel backend
- [`sig_filter_fn()`](https://montilab.github.io/sigrecon/reference/sig_filter_fn.md)
  : Extract DEGs from a Differential Expression Table
