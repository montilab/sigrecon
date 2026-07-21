# Filter Signatures for Common Nodes Across Graphs

This function takes a list of graphs and a list of node signatures, and
filters each signature to include only nodes that are present in the
largest connected component of all provided graphs.

## Usage

``` r
common_signature_filter(graphs, signatures)
```

## Arguments

- graphs:

  A list of igraph objects

- signatures:

  A list of node signatures (vectors of node names)

## Value

A list of filtered signatures containing only nodes common to all graphs

## Details

For each signature:

1.  It checks which nodes are in the largest connected component of each
    graph.

2.  It finds the intersection of these nodes across all graphs.

3.  It returns this intersection as the new filtered signature.
