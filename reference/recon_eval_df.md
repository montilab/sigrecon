# Evaluate recontextualized signatures

Evaluate recontextualized signatures

## Usage

``` r
recon_eval_df(
  ig,
  seed_name,
  source_sigs,
  dest_sigs,
  restart = 0.75,
  avg_p = FALSE,
  avg_p_vals = c(1e-04, 0.1),
  avg_p_length = 5,
  bootstrap = FALSE,
  n_bootstraps = 1000,
  recon = TRUE,
  use_weights = TRUE,
  weights.pwr = 1,
  normalize = c("row", "column", "laplacian"),
  save = FALSE,
  save_path = "",
  limit = 30
)
```

## Arguments

- ig:

  igraph

- seed_name:

  name of list of genesets to be used to label dataframe

- source_sigs:

  list of genesets (to be recontextualized)

- dest_sigs:

  list of genesets (the ground truth)

- restart:

  restart value

- avg_p:

  A boolean specifying whether to ensemble random walk results over a
  range of restart values

- avg_p_vals:

  A numeric vector specifying the start and end of a geometric sequence
  to explore restart values.

- avg_p_length:

  A numeric specifying how many values within `avg_p_vals` to include in
  the ensemble

- bootstrap:

  A boolean specifying whether to use empirical distributions of
  stationary values to find significant genes.

- n_bootstraps:

  A numeric specifying the number of bootstraps to perform.

- recon:

  Boolean indicating whether recontextualization occurs

- use_weights:

  Boolean indicating whether to use weights in the KS Test

- weights.pwr:

  Power to raise weights to

- normalize:

  Normalization strategy to employ

- save:

  Boolean indicating whether to save recontextualized signatures

- save_path:

  file path to save recontextualized signatures

- limit:

  Number of genes to keep in the output, or a vector of lengths
