#' Set Diagonal of a Matrix to Zero
#'
#' @description
#' This function sets the diagonal elements of a matrix to zero.
#'
#' @param matrix A square matrix.
#'
#' @return A matrix of the same dimensions as the input, with diagonal elements set to zero.
diag_zero <- function(matrix) {
  diag(matrix) <- 0
  return(matrix)
}

#' Create a Consensus Network from Multiple Networks
#'
#' @description
#' This function creates a consensus network from a list of adjacency matrices by averaging
#' and thresholding the connections.
#'
#' @param list_of_nets A list of adjacency matrices, each representing a network.
#' @param threshold A numeric value between 0 and 1 representing the consensus threshold. Default is 0.9.
#'
#' @return An adjacency matrix representing the consensus network.
#'
#' @details
#' The function performs the following steps:
#' 1. Averages all input adjacency matrices.
#' 2. Applies a threshold to the averaged matrix.
#' 3. Returns a binary adjacency matrix where 1 indicates a consensus connection and 0 indicates no consensus.
#'
#' @importFrom purrr reduce
consensus_net <- function(list_of_nets, threshold=0.9) {
  ### List of nets: A list of adjacency matrices
  avg_adj <- reduce(list_of_nets, `+`) / length(list_of_nets)
  threshed_adj <- (avg_adj > threshold) * 1
  return(threshed_adj)
}

# WGCNA wrappers

#' Check and Return Soft-Thresholding Power from WGCNA
#'
#' @description
#' This function checks the soft-thresholding power estimate from WGCNA's scale-free topology analysis
#' and returns an appropriate power value.
#'
#' @param sft A list object returned by WGCNA's pickSoftThreshold or pickSoftThreshold.fromSimilarity function.
#'
#' @return An integer representing the soft-thresholding power to be used.
#'
#' @details
#' If a valid power estimate is found in the input, it is returned. If the power estimate is NA,
#' a default value of 6 is returned. The function prints a message indicating which power is being used.
sft.check <- function(sft) {
  beta <- sft$powerEstimate
  if (is.na(beta)) {
    beta <- 6 # Default
    cat("Using the following power:", beta, "\n")
  } else {
    cat("Optimal power selected:", beta, "\n")
  }
  return(beta)
}

#' Construct WGCNA Adjacency Matrix from Correlation Matrix
#'
#' @description
#' This function constructs a weighted gene co-expression network adjacency matrix using WGCNA,
#' starting from a pre-computed correlation matrix.
#'
#' @param cor_mat A correlation matrix of gene expression data.
#' @param cores Number of CPU cores to use for parallel computing. Default is 1.
#' @param diag_zero If TRUE, sets the diagonal of the adjacency matrix to zero. Default is TRUE.
#'
#' @return An adjacency matrix representing the gene co-expression network.
#'
#' @details
#' This function performs the following steps:
#' 1. Selects the optimal soft-thresholding power using WGCNA::pickSoftThreshold.fromSimilarity().
#' 2. Constructs the adjacency matrix using WGCNA::adjacency.fromSimilarity().
#' 3. Optionally sets the diagonal to zero and/or converts the result to an igraph object.
#'
#' @importFrom doParallel registerDoParallel
wgcna.power <- function(cor_mat,
                        cores=1,
                        diag_zero=TRUE) {

  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("The 'WGCNA' package is required for wgcna.power(). Install it with BiocManager::install('WGCNA').")
  }

  # Set parallel computing environment
  doParallel::registerDoParallel(cores=cores)

  # Pick soft threshold via scale-free fit
  sft <- WGCNA::pickSoftThreshold.fromSimilarity(similarity=cor_mat)

  # Check selected power
  beta <- sft.check(sft)

  # Construct co-expression similarity
  adj <- WGCNA::adjacency.fromSimilarity(similarity=cor_mat,
                                         power=beta)

  if(diag_zero) {
    adj <- diag_zero(adj)
  }

  return(adj)
}

#' Construct WGCNA Adjacency Matrix
#'
#' @description
#' This function constructs a weighted gene co-expression network adjacency matrix using WGCNA.
#'
#' @param mat A matrix with rows as samples and columns as genes
#' @param min.sft Minimum scale-free topology fitting index R^2 to pick soft-thresholding power. Default is 0.85.
#' @param beta Soft-thresholding power. If NULL, it will be automatically selected. Default is NULL.
#' @param cores Number of CPU cores to use for parallel computing. Default is 1.
#' @param cor.fn Correlation function to use. Either "bicor" (biweight midcorrelation) or "cor" (Pearson correlation). Default is "bicor".
#' @param cor.type Type of correlation network. Options are "unsigned", "signed hybrid", or "signed". Default is "unsigned".
#' @param powers Vector of soft-thresholding powers to try. Default is c(seq(1, 10, by = 1), seq(12, 20, by = 2)).
#' @param igraph If TRUE, returns an igraph object instead of a matrix. Default is FALSE.
#' @param diag_zero If TRUE, sets the diagonal of the adjacency matrix to zero. Default is FALSE.
#'
#' @return An adjacency matrix or an igraph object representing the gene co-expression network.
#'
#' @details
#' This function performs the following steps:
#' 1. Prepares the expression data.
#' 2. Selects the soft-thresholding power (if not provided).
#' 3. Constructs the adjacency matrix using WGCNA::adjacency().
#' 4. Optionally converts the result to an igraph object.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom Biobase exprs
#' @importFrom SummarizedExperiment assays
#' @importFrom igraph graph_from_adjacency_matrix
#'
#' @export
wgcna.adj <- function(mat,
                      min.sft=0.85,
                      beta=NULL,
                      cores=1,
                      cor.fn=c("bicor", "cor"),
                      cor.type=c("unsigned", "signed hybrid", "signed"),
                      powers=c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                      igraph=FALSE,
                      diag_zero=FALSE) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("The 'WGCNA' package is required for wgcna.adj(). Install it with BiocManager::install('WGCNA').")
  }

  # Ad hoc namespace changes
  bicor = WGCNA::bicor
  cor = WGCNA::bicor

  # Handle arguments
  args <- as.list(environment())
  cor.fn <- match.arg(cor.fn)
  cor.type <- match.arg(cor.type)
  # Correlation options
  if (cor.fn == "cor") cor.options = list(use="p")
  if (cor.fn == "bicor") cor.options = list(pearsonFallback="individual")

  # Set parallel computing environment
  doParallel::registerDoParallel(cores=cores)

  # Pick soft threshold via scale-free fit
  if (is.null(beta)) {
    sft <- WGCNA::pickSoftThreshold(data=mat,
                                    corFnc=cor.fn,
                                    RsquaredCut=min.sft,
                                    powerVector=powers)

    # Check selected power
    beta <- sft.check(sft)
  }

  # Construct co-expression similarity
  adj <- WGCNA::adjacency(datExpr=mat,
                          power=beta,
                          corFnc=cor.fn,
                          type=cor.type,
                          corOptions=cor.options)
  if(diag_zero) {
    adj <- diag_zero(adj)
  }
  if(igraph) {
    adj <- igraph::graph_from_adjacency_matrix(adj, weighted=TRUE, mode="undirected")
  }
  return(adj)
}


