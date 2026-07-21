# Recontextualize signatures with a selected baseline method

Recontextualize signatures with a selected baseline method

## Usage

``` r
recontextualize(
  method = c("networkProp", "projectCor", "mean"),
  se = NULL,
  seeds = NULL,
  sigs = NULL,
  score = c("gsva", "eigen"),
  sig = c("rwr", "corr"),
  avg_p = FALSE,
  avg_p_vals = c(1e-04, 0.1),
  avg_p_length = 5,
  p = 0.1,
  bootstrap = FALSE,
  n_bootstraps = 1000,
  limit = NULL,
  nfeatures = NULL,
  min.sft = 0.85,
  beta = NULL,
  cores = 1,
  cor.fn = c("bicor", "cor"),
  cor.type = c("unsigned", "signed hybrid", "signed"),
  powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
  diag_zero = TRUE,
  perturbation_col = NULL,
  condition_col = NULL,
  perturbed_label = "Perturbed",
  control_label = "Control",
  design_vars = NULL,
  assay_name = NULL,
  alpha = 0.05,
  min_count = NULL,
  min_samples = 1
)
```

## Arguments

- method:

  Baseline recontextualization method. One of `"networkProp"`,
  `"projectCor"`, or `"mean"`.

- se:

  A SummarizedExperiment object used by `"projectCor"`, `"mean"`, and to
  learn the `"networkProp"` graph.

- seeds:

  Seed signatures for the `"networkProp"` method.

- sigs:

  Signature list for the `"projectCor"` method. For `"mean"`, names are
  interpreted as perturbation labels and lengths are used as output
  sizes when `limit` is not supplied.

- score:

  Scoring method used by `"projectCor"`.

- sig:

  Network signature mode passed to
  [`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md).
  Defaults to `"rwr"`.

- avg_p:

  Whether to ensemble network propagation over multiple restart
  probabilities.

- avg_p_vals:

  Range of restart probabilities used when `avg_p = TRUE`.

- avg_p_length:

  Number of restart probabilities to average when `avg_p = TRUE`.

- p:

  Restart probability passed to
  [`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md).

- bootstrap:

  Whether to use bootstrap-based extraction in
  [`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md).

- n_bootstraps:

  Number of bootstrap replicates passed to
  [`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md).

- limit:

  Number of genes to keep for each output signature. If `NULL`,
  `"networkProp"` and `"projectCor"` keep the original signature lengths
  and `"mean"` uses the lengths of `sigs`.

- nfeatures:

  Number of variable genes to use when learning the `"networkProp"`
  graph from `se`. Defaults to `min(10000, nrow(assay(se)))`.

- min.sft:

  Minimum scale-free topology fitting index used by
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
  when learning the `"networkProp"` graph.

- beta:

  Optional soft-thresholding power passed to
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
  for `"networkProp"`.

- cores:

  Number of CPU cores passed to
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
  for `"networkProp"`.

- cor.fn:

  Correlation function passed to
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
  for `"networkProp"`.

- cor.type:

  Correlation network type passed to
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
  for `"networkProp"`.

- powers:

  Candidate power values passed to
  [`wgcna.adj()`](https://montilab.github.io/sigrecon/reference/wgcna.adj.md)
  for `"networkProp"` when `beta` is `NULL`.

- diag_zero:

  Whether to zero the diagonal of the learned `"networkProp"` adjacency
  matrix before converting it to igraph.

- perturbation_col:

  Column in `colData(se)` identifying the perturbation label for each
  profile. Used by `"mean"`.

- condition_col:

  Column in `colData(se)` indicating whether a profile is perturbed or
  control. Used by `"mean"`.

- perturbed_label:

  Value in `condition_col` corresponding to perturbed profiles.

- control_label:

  Value in `condition_col` corresponding to control profiles.

- design_vars:

  Optional character vector of additional covariates to place before
  condition in the DESeq2 design formula for `"mean"`.

- assay_name:

  Assay name in `se` to use for `"mean"`. Defaults to the first assay.

- alpha:

  Adjusted p-value cutoff used by the `"mean"` method when selecting
  upregulated genes.

- min_count:

  Optional count threshold used to filter low-expression genes before
  DESeq2 for `"mean"`. If `NULL`, no pre-filtering is applied.

- min_samples:

  Minimum number of samples that must satisfy `min_count` for a gene to
  be retained when `min_count` is provided.

## Value

A named list of recontextualized gene signatures.
