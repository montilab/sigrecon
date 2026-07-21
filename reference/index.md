# Package index

## Recontextualization methods

Baseline methods that take a source-context signature and reconstruct it
in a target context.

- [`projectCor()`](https://montilab.github.io/sigrecon/reference/projectCor.md)
  : Reconstruct Gene Signatures Using Projection Scores
- [`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)
  : Finds a simulated network signature
- [`recontextualize()`](https://montilab.github.io/sigrecon/reference/recontextualize.md)
  : Recontextualize signatures with a selected baseline method

## Network utilities

Building and propagating across gene co-expression networks.

- [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
  : Construct WGCNA Adjacency Matrix
- [`random_walk()`](https://montilab.github.io/sigrecon/reference/random_walk.md)
  : Perform a random walk with restart (personalized page rank) on an
  igraph given a seed matrix, and return stationary probabilties.
  Stripped down and corrected version of dnet:
  https://rdrr.io/cran/dnet/src/R/dRWR.r

## Benchmarking & evaluation

Evaluating a recontextualized signature against source/true signatures.

- [`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md)
  : Evaluate Signature Prediction
- [`paired_eval_table()`](https://montilab.github.io/sigrecon/reference/paired_eval_table.md)
  : Pair a Recontextualization Evaluation Table With a No-Change
  Baseline
- [`sig_filter_fn()`](https://montilab.github.io/sigrecon/reference/sig_filter_fn.md)
  : Filter Significant Genes by Perturbation
- [`recon_eval_df()`](https://montilab.github.io/sigrecon/reference/recon_eval_df.md)
  : Evaluate recontextualized signatures
- [`jaccard()`](https://montilab.github.io/sigrecon/reference/jaccard.md)
  : Calculate Jaccard Similarity Between Two Sets
- [`jaccard_matrix()`](https://montilab.github.io/sigrecon/reference/jaccard_matrix.md)
  : Create a Jaccard Similarity Matrix for Multiple Sets
- [`v.jaccard()`](https://montilab.github.io/sigrecon/reference/v.jaccard.md)
  : Calculate Jaccard Similarity for Multiple Pairs of Sets
- [`v.fgsea()`](https://montilab.github.io/sigrecon/reference/v.fgsea.md)
  : Vectorized fgsea wrapper for multiple gene set comparisons
- [`fgsea_wrapper()`](https://montilab.github.io/sigrecon/reference/fgsea_wrapper.md)
  : fgsea wrapper for gene symbol vectors
- [`fishers_meta_p()`](https://montilab.github.io/sigrecon/reference/fishers_meta_p.md)
  : Combine P-values Using Fisher's Method
- [`count_displaced_genes()`](https://montilab.github.io/sigrecon/reference/count_displaced_genes.md)
  : Count number of displaced seeds Note: This not commutative.
  Num_displaced is the number of genes that are from the original seed
  that are not in the recon.
- [`bin_presence()`](https://montilab.github.io/sigrecon/reference/bin_presence.md)
  : Check Presence in Bins for Numbers 1-100

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

## Internal

Unexported helpers, documented for maintainers but not part of the
package’s public API.

- [`all_dists_nodesets()`](https://montilab.github.io/sigrecon/reference/all_dists_nodesets.md)
  : Find all pairwise distances between nodes in an igraph
- [`annotate_prob_vec()`](https://montilab.github.io/sigrecon/reference/annotate_prob_vec.md)
  : Returns a dataframe object with the stationary probability value and
  column indicating whether gene was a seed
- [`binarize_matrices()`](https://montilab.github.io/sigrecon/reference/binarize_matrices.md)
  : Binarize Multiple Matrices Based on Cutoffs
- [`col_normalize()`](https://montilab.github.io/sigrecon/reference/col_normalize.md)
  : Normalize Columns of a Matrix
- [`common_mad_genes()`](https://montilab.github.io/sigrecon/reference/common_mad_genes.md)
  : Find Common Genes with Highest Median Absolute Deviation (MAD)
  Across ExpressionSets
- [`common_signature_filter()`](https://montilab.github.io/sigrecon/reference/common_signature_filter.md)
  : Filter Signatures for Common Nodes Across Graphs
- [`consensus_net()`](https://montilab.github.io/sigrecon/reference/consensus_net.md)
  : Create a Consensus Network from Multiple Networks
- [`correlated_sigs()`](https://montilab.github.io/sigrecon/reference/correlated_sigs.md)
  : Recontextualize seed signatures with correlation based neighbors
- [`diag_zero()`](https://montilab.github.io/sigrecon/reference/diag_zero.md)
  : Set Diagonal of a Matrix to Zero
- [`extract_sig_mat()`](https://montilab.github.io/sigrecon/reference/extract_sig_mat.md)
  : Extracts a signature from a (gene x seed) matrix of stationary
  probability values. This is the recontextualized signature. If doing
  ks.test, you don't need to find the top_n. Just find ks.test(original,
  recontextualized ranking) before and after.
- [`ggempty()`](https://montilab.github.io/sigrecon/reference/ggempty.md)
  : An empty ggplot
- [`ggeplot()`](https://montilab.github.io/sigrecon/reference/ggeplot.md)
  : Enrichment plot implemented in ggplot
- [`largest_connected_subgraph()`](https://montilab.github.io/sigrecon/reference/largest_connected_subgraph.md)
  : Extract Largest Connected Subgraph
- [`pvector`](https://montilab.github.io/sigrecon/reference/pvector.md)
  : A push/pop capable vector
- [`rank.var.eset()`](https://montilab.github.io/sigrecon/reference/rank.var.eset.md)
  : Rank Genes in an ExpressionSet/SummarizedExperiment by Variability
- [`seed_matrix()`](https://montilab.github.io/sigrecon/reference/seed_matrix.md)
  : Create (gene x seed) prior matrix based on seed signatures.
- [`seurat_common_var_genes()`](https://montilab.github.io/sigrecon/reference/seurat_common_var_genes.md)
  : Find Common Variable Genes Across Seurat Objects
- [`sft.check()`](https://montilab.github.io/sigrecon/reference/sft.check.md)
  : Check and Return Soft-Thresholding Power from WGCNA
- [`v.correlated_sigs()`](https://montilab.github.io/sigrecon/reference/v.correlated_sigs.md)
  : Recontextualize seed signatures with correlation based neighbors
- [`wgcna.power()`](https://montilab.github.io/sigrecon/reference/wgcna.power.md)
  : Construct WGCNA Adjacency Matrix from Correlation Matrix
