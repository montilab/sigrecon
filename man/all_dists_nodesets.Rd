% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{all_dists_nodesets}
\alias{all_dists_nodesets}
\title{Find all pairwise distances between nodes in an igraph}
\usage{
all_dists_nodesets(g, node_list)
}
\arguments{
\item{g}{An igraph object representing the graph.}

\item{node_list}{A list where each element is a vector of node names or IDs.}
}
\value{
A list of vectors, each containing the pairwise distances for the corresponding node set.
}
\description{
This function computes the pairwise distances between nodes within each set of nodes
provided in a list, based on the structure of an input graph.
}
\details{
The function first calculates the full distance matrix for the graph using igraph::distances().
Then, for each set of nodes in \code{node_list}, it extracts the relevant submatrix and returns
the lower triangular part, which represents all pairwise distances within that set.
}
