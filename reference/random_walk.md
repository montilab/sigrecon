# Perform a random walk with restart (personalized page rank) on an igraph given a seed matrix, and return stationary probabilties. Stripped down and corrected version of dnet: https://rdrr.io/cran/dnet/src/R/dRWR.r

Perform a random walk with restart (personalized page rank) on an igraph
given a seed matrix, and return stationary probabilties. Stripped down
and corrected version of dnet: https://rdrr.io/cran/dnet/src/R/dRWR.r

## Usage

``` r
random_walk(
  ig,
  seed_mat,
  restart = 0.1,
  epsilon = NULL,
  normalize = c("row", "column", "laplacian", "none")
)
```

## Arguments

- ig:

  igraph object

- seed_mat:

  (Gene, num_seeds) matrix with prior weights for each gene in a seed
  set. See seed_matrix.

- restart:

  the restart probability for RWR

- epsilon:

  Exploration factor

- normalize:

  Normalization strategy

## Value

It returns a sparse matrix with stationary probabilities.
