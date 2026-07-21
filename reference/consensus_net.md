# Create a Consensus Network from Multiple Networks

This function creates a consensus network from a list of adjacency
matrices by averaging and thresholding the connections.

## Usage

``` r
consensus_net(list_of_nets, threshold = 0.9)
```

## Arguments

- list_of_nets:

  A list of adjacency matrices, each representing a network.

- threshold:

  A numeric value between 0 and 1 representing the consensus threshold.
  Default is 0.9.

## Value

An adjacency matrix representing the consensus network.

## Details

The function performs the following steps:

1.  Averages all input adjacency matrices.

2.  Applies a threshold to the averaged matrix.

3.  Returns a binary adjacency matrix where 1 indicates a consensus
    connection and 0 indicates no consensus.
