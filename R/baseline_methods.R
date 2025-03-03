library(GSVA)
library(igraph)
library(Matrix)
library(tidyverse)

#' Reconstruct Gene Signatures Using GSVA Scores
#'
#' @description
#' This function reconstructs gene signatures based on their correlation with GSVA (Gene Set Variation Analysis) scores.
#'
#' @param sce A SingleCellExperiment object containing gene expression data.
#' @param sigs A list of gene signatures, where each element is a character vector of gene names.
#'
#' @return A list of reconstructed gene signatures, with the same structure as the input `sigs`.
#'
#' @details
#' The function performs the following steps:
#' 1. Calculates GSVA scores for the input signatures.
#' 2. Computes the correlation between gene expression and GSVA scores.
#' 3. Ranks genes based on their correlation with each signature's GSVA score.
#' 4. Selects the top-ranking genes to form new signatures of the same length as the original ones.
#'
#' @note
#' This function requires the GSVA, stats, and dplyr packages.
#'
#' @importFrom GSVA gsva
#' @importFrom stats cor
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange slice pull desc
#'
#' @export
gsva_recon <- function(sce,
                       sigs) {

  # GSVA Scores
  gsva_eset <- GSVA::gsva(sce, sigs, mx.diff = TRUE, verbose = FALSE)

  genes_data <- t(sce@assays@data$counts)
  gsva_scores <- t(gsva_eset@assays@data$es)
  corr_mat <- stats::cor(genes_data, gsva_scores, method = "pearson")

  results <- as_tibble(corr_mat, rownames="gene")
  new_sigs <- list()
  for(sig_name in names(sigs)) {
    # Obtain Rank of Genes most correlated to GSVA scores
    results$rank <- rank(dplyr::desc(results[[sig_name]]))

    # Obtain new signature
    sig_length <- length(sigs[[sig_name]])
    results <- results %>% dplyr::arrange(rank)
    new_sig <- results %>% dplyr::slice(1:sig_length) %>% dplyr::pull(gene)
    new_sigs[[sig_name]] <- new_sig
  }

  return(new_sigs)
}


#' Create (gene x seed) prior matrix based on seed signatures.
#'
#' @description
#' This function creates a binary matrix representing seed genes in the context of a graph.
#'
#' @param ig An igraph object representing the network that has the seed genes as vertices
#' @param seeds Either a single unnamed gene "TP53", a named list of genes, or a list of named lists of genes.
#'
#' @return A binary matrix where rows represent genes in the graph and columns represent seed sets.
#'
#' @importFrom igraph V
#'
#' @export
seed_matrix <- function(ig, seeds) {
  # Seeds can be a single character vector. In that case need to list-ify it.
  if (is(seeds, "character") && length(seeds) == 1) {
    seeds <- list(seeds)
    names(seeds) <- seeds
  } else if (is.null(names(seeds))) {
    stop("Seed signature needs name")
  }

  n_nodes <- length(seeds)

  # Create Seed Matrix
  mat <- matrix(0, nrow = length(V(ig)), ncol = n_nodes)
  rownames(mat) <- V(ig)$name
  colnames(mat) <- names(seeds)

  all_genes <- unname(unlist(seeds))

  seed_genes_filter <- all(all_genes %in% rownames(mat))

  if (!seed_genes_filter) {
    print("Seeds Genes are not all contained in the graph.")
    print("Filtering Seed Genes to those contained in the graph.")
    seeds <- lapply(seeds, function(x) x[x %in% rownames(mat)])
  }

  if (n_nodes > 1) {
    for (seed_name in names(seeds)) {
      genes <- seeds[[seed_name]]
      mat[genes, seed_name] <- 1
    }
  } else {
    genes <- seeds[[1]]
    mat[genes, ] <- 1
  }

  return(mat)
}

#' Perform a random walk on an igraph given seeds, return stationary probabilities
#'
#' @param ig igraph object
#' @param seeds Gene Symbol(s) that are seed nodes for the random walk. If random walk is to be done on multiple genes individually or multiple sets of genes, they should be in seperate elements of the list.
#' @param restart Numeric, Probability of restarting at the seed nodes
#' @param epsilon Exploration factor
#' @param normalize Normalization strategy
#' @return (n_gene, n_seeds) matrix of stationary probability values
#'
#' @export
#'
rwr_mat <- function(ig, seeds, restart = 0.75, epsilon = NULL, normalize = c("row", "column", "laplacian", "none")) {
  # Assumes all seeds are present in ig graph
  # browser()
  stopifnot(is(ig, "igraph"))
  stopifnot("name" %in% igraph::vertex_attr_names(ig))

  normalize <- match.arg(normalize)

  seed_mat <- seed_matrix(ig, seeds)

  rwr <- random_walk(ig = ig, seed_mat = seed_mat, restart = restart, epsilon = epsilon, normalize = normalize)

  return(rwr)
}


#' Find the top n genes per seed in a (gene x seed) matrix.
#' This is the recontextualized signature.
#' If doing ks.test, you don't need to find the top_n. Just find ks.test(original, recontextualized ranking) before and after.
#' @param mat (n_gene, n_seed) matrix of stationary probability values from rwr_mat
#' @param limit Number of genes to keep in the output, or a vector of lengths
#' @return A named list of genesets. Each list element is the recontextualized signature for that seed.
#'
#' @export
top_n_mat <- function(mat, limit = 30) {
  top_n <- list()

  if (length(limit) == 1) {
    limits <- rep(limit, ncol(mat))
  } else {
    limits <- limit
  }

  for (col in 1:ncol(mat)) {
    colname <- colnames(mat)[[col]]
    mat_col <- mat[, col]
    limit <- limits[[col]]
    sig <- sort(mat_col, decreasing = TRUE) %>%
      head(limit) %>%
      names()
    top_n[[colname]] <- sig
  }

  return(top_n)
}


#' Returns a dataframe object with the stationary probability value and column indicating whether gene was a seed
#' @param prob_vec (n_gene, 1) matrix
#' @param seeds character vector
#'
#' @export
annotate_prob_vec <- function(prob_vec, seeds) {
  stopifnot(is(prob_vec, "dgCMatrix") | is(prob_vec, "Matrix"))
  stopifnot("Can only annotate (n x 1) vectors" = dim(prob_vec)[[2]] == 1)

  if (is(prob_vec,"dgCMatrix")) {
    df <- as.data.frame(as.matrix(prob_vec))
  } else {
    df <- as.data.frame(prob_vec)
  }

  colnames(df) <- "prob"
  df$seed <- rownames(df) %in% seeds
  df <- df %>% arrange(desc(prob))

  return(df)
}

#' Perform a random walk on an igraph given seeds, return stationary probabilities
#'
#' @param ig igraph object
#' @param seeds Gene Symbol(s) that are seed nodes for the random walk. If random walk is to be done on multiple genes individually or multiple sets of genes, they should be in seperate elements of the list.
#' @param restart Numeric, Probability of restarting at the seed nodes
#' @param normalize Normalization strategy
#' @return List of Annotated Dataframes. Each Dataframe has columns for gene label, probability value, and seed status.
#'
#' @export
#'
rwr_df <- function(ig, seeds, restart = 1e-2, normalize = c("row", "column", "laplacian", "none")) {

  normalize <- match.arg(normalize)
  mat <- rwr_mat(ig = ig, seeds = seeds, restart = restart, normalize = normalize)

  dfs <- sapply(colnames(mat),
                function(x) annotate_prob_vec(mat[, x, drop = FALSE], seeds = seeds[[x]]),
                USE.NAMES = TRUE,
                simplify = FALSE
  )
  return(dfs)
}


#' Perform a random walk with restart (personalized page rank) on an igraph given a seed matrix, and return stationary probabilties.
#' Stripped down and corrected version of dnet: https://rdrr.io/cran/dnet/src/R/dRWR.r
#'
#' @param ig igraph object
#' @param seed_mat (Gene, num_seeds) matrix with prior weights for each gene in a seed set. See seed_matrix.
#' @param restart the restart probability for RWR
#' @param epsilon Exploration factor
#' @param normalize Normlization strategy
#' @return It returns a sparse matrix with stationary probabilities.
#'
#' @importFrom igraph list.edge.attributes as_adjacency_matrix
#' @importFrom Matrix Diagonal colSums rowSums t
#'
#' @export
random_walk <- function(ig, seed_mat, restart = 0.1, epsilon = NULL, normalize = c("row", "column", "laplacian", "none")) {

  # Type checks
  stopifnot(is(ig) == "igraph")
  stopifnot(is(seed_mat, "matrix"))
  normalize <- match.arg(normalize)

  # Get Adjacency matrix
  if ("weight" %in% igraph::list.edge.attributes(ig)) {
    adj_mat <- igraph::as_adjacency_matrix(ig, attr = "weight")
    adj_mat[is.na(adj_mat)] <- 0
    print("Using weighted graph")
  } else {
    adj_mat <- igraph::as_adjacency_matrix(ig, attr = NULL)
    adj_mat[is.na(adj_mat)] <- 0
    print("Using unweighted graph")
  }

  # Normalize adjacency matrix. DNet::dRWR had the multiplication orders flipped for row and column.
  if (normalize == "row") {
    D <- Matrix::Diagonal(x = (Matrix::rowSums(adj_mat))^(-1))
    nadjM <- D %*% adj_mat
  } else if (normalize == "column") {
    D <- Matrix::Diagonal(x = (Matrix::colSums(adj_mat))^(-1))
    nadjM <- adj_mat %*% D
  } else if (normalize == "laplacian") {
    D <- Matrix::Diagonal(x = (Matrix::colSums(adj_mat))^(-0.5))
    nadjM <- D %*% adj_mat %*% D
  } else if (normalize == "none") {
    nadjM <- adjM
  }

  ## Stopping Criteria
  stop_delta <- 1e-5 # L1 norm of successive iterations of Transition Matrix multiplication
  stop_step <- 100 # maximum steps of iterations

  # Normalize Seed Matrix
  norm_seed_mat <- seed_mat %*% Matrix::Diagonal(x = (Matrix::colSums(seed_mat))^(-1))
  colnames(norm_seed_mat) <- colnames(seed_mat)
  norm_seed_mat <- as(norm_seed_mat, "CsparseMatrix")

  ## Initial Variables
  P0 <- norm_seed_mat
  PT <- P0
  r <- restart
  step <- 0
  delta <- 1

  ## Exploration Parameters
  if (!is.null(epsilon)) {
    stopifnot(is(epsilon,"numeric"))

    n_genes_seed_mean <- as.integer(mean(apply(seed_mat, 2, sum)))
    n_explore <- as.integer(n_genes_seed_mean * epsilon)
    r <- 1
    paste0("Stopping Criterion: ", n_explore, " displaced genes.")
  }

  ## Dnet had the matrix multiplication orders flipped.
  ## This order keeps the distribution in columns after each multiplication.
  while (delta > stop_delta && step <= stop_step) {
    PX <- (1 - r) * Matrix::t(Matrix::t(PT) %*% nadjM) + r * P0
    delta <- sum(abs(PX - PT))
    PT <- PX
    step <- step + 1

    ## Write function that counts number of displaced.
    ## Add stopping condition
    # if (!is.null(epsilon)) {
    # }

    if (step > stop_step) {
      print(paste0("Reached maximum iteration steps. Delta: ", delta))
    } else if (delta <= stop_delta) {
      print(paste0("Reached Convergence. Iteration step: ", step))
    }
  }

  return(PX)
}

