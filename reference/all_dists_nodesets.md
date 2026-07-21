# Find all pairwise distances between nodes in an igraph

This function computes the pairwise distances between nodes within each
set of nodes provided in a list, based on the structure of an input
graph.

## Usage

``` r
all_dists_nodesets(g, node_list)
```

## Arguments

- g:

  An igraph object representing the graph.

- node_list:

  A list where each element is a vector of node names or IDs.

## Value

A list of vectors, each containing the pairwise distances for the
corresponding node set.

## Details

The function first calculates the full distance matrix for the graph
using igraph::distances(). Then, for each set of nodes in `node_list`,
it extracts the relevant submatrix and returns the lower triangular
part, which represents all pairwise distances within that set.
