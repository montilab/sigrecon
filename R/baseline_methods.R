library(GSVA)
library(igraph)
library(Matrix)
library(abind)
library(parallel)
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

  if (is(genes_data, "sparseMatrix")) {
    genes_data <- as.matrix(genes_data)
  }

  gsva_scores <- t(gsva_eset@assays@data$es)
  corr_mat <- stats::cor(genes_data, gsva_scores, method = "pearson")

  results <- as_tibble(corr_mat, rownames="gene")
  new_sigs <- list()
  for(sig_name in colnames(gsva_scores)) {
    print(sig_name)
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
#' @param bootstrap A boolean specifying whether to use empirical distributions of stationary values to find significant genes.
#' @param n_bootstraps A numeric specifying the number of bootstraps to perform.
#' @return A binary matrix where rows represent genes in the graph and columns represent seed sets.
#'
#' @importFrom igraph V
seed_matrix <- function(ig,
                        seeds,
                        bootstrap = FALSE,
                        n_bootstraps = 1000) {

  # Seeds can be a single character vector. In that case need to list-ify it.
  if (is.character(seeds) && length(seeds) == 1) {
    seeds <- setNames(list(seeds), seeds)
  }
  if (is.null(names(seeds))) stop("Seed signature needs name")

  gene_names <- igraph::V(ig)$name
  n_pbs <- length(seeds)

  # Create Seed Matrix
  if (bootstrap) {
    bins_filter <- bin_presence(lapply(seeds, function(x) length(x)) %>% unlist)
    n_bins <- sum(bins_filter)
    mat <- matrix(0, nrow = length(gene_names), ncol = n_bins * n_bootstraps)
    rownames(mat) <- gene_names
    colnames(mat) <- bins_filter %>% which %>% names %>% rep(each = n_bootstraps)

    # If bootstrap matrix is (n_genes, n_bootstraps x n_bins)
    for (present_bin in seq_along(which(bins_filter))) {
      bin <- which(bins_filter)[present_bin]
      for (i in 1:n_bootstraps) {
        # Sample a random set of genes with size sampled between bin end-points
        last_digit <- i %% 10
        last_digit <- ifelse(last_digit == 0, 10, last_digit)

        sample_size <- (bin - 1) * 10 + last_digit
        genes <- sample(gene_names, sample_size)

        bin_col <- (present_bin-1) * n_bootstraps + i
        mat[genes, bin_col] <- 1
      }
    }
  } else {
    mat <- matrix(0, nrow = length(gene_names), ncol = n_pbs)
    rownames(mat) <- gene_names
    colnames(mat) <- names(seeds)
    # If not bootstrap matrix is (n_genes, n_seeds)
    for (seed_name in names(seeds)) {
      genes <- seeds[[seed_name]]
      mat[genes, seed_name] <- 1
    }
  }

  return(mat)
}

#' Perform a random walk on an igraph given seeds, return stationary probabilities
#'
#' @param ig igraph object
#' @param seeds Either a single unnamed gene "TP53", a named list of genes, or a list of named lists of genes
#' @param restart A numeric specifying the probability of restarting at the seed nodes
#' @param avg_p A boolean specifying whether to ensemble random walk results over a range of restart values
#' @param avg_p_vals A numeric vector specifying the start and end of a arithmetic sequence to explore restart values
#' @param avg_p_length A numeric specifying how many values within `avg_p_vals` to include in the ensemble
#' @param bootstrap A boolean specifying whether to use empirical distributions of stationary values to find significant genes.
#' @param n_bootstraps A numeric specifying the number of bootstraps to perform.
#' @param epsilon Exploration factor
#' @param normalize Normalization strategy
#' @return (n_gene, n_seeds) matrix of stationary probability values
rwr_mat <- function(ig,
                    seeds,
                    restart = 0.75,
                    avg_p = FALSE,
                    avg_p_vals = c(1e-4, 1e-1),
                    avg_p_length = 5,
                    bootstrap = FALSE,
                    n_bootstraps = 1000,
                    epsilon = NULL,
                    normalize = c("row", "column", "laplacian", "none")) {
  # Assumes all seeds are present in ig graph
  # browser()
  stopifnot(is(ig, "igraph"))
  stopifnot("name" %in% igraph::vertex_attr_names(ig))

  normalize <- match.arg(normalize)

  seed_mat <- seed_matrix(ig, seeds, bootstrap = bootstrap, n_bootstraps = n_bootstraps)

  if (avg_p) {
    end <- avg_p_vals[2]
    start <- avg_p_vals[1]
    avg_p_seq <- seq(start, end, length.out = avg_p_length)

    mat <- matrix(0, nrow = nrow(seed_mat), ncol = ncol(seed_mat))
    rownames(mat) <- rownames(seed_mat)
    colnames(mat) <- colnames(seed_mat)

    for (restart_val in avg_p_seq) {
      rwr <- random_walk(ig = ig, seed_mat = seed_mat, restart = restart_val, epsilon = epsilon, normalize = normalize)
      mat <- mat + rwr
    }
    mat <- mat/avg_p_length
  } else {
    mat <- random_walk(ig = ig, seed_mat = seed_mat, restart = restart, epsilon = epsilon, normalize = normalize)
  }

  return(mat)
}


#' Extracts a signature from a (gene x seed) matrix of stationary probability values.
#' This is the recontextualized signature.
#' If doing ks.test, you don't need to find the top_n. Just find ks.test(original, recontextualized ranking) before and after.
#' @param mat (n_gene, n_seed) matrix of stationary probability values from rwr_mat
#' @param bootstraps A (n_gene, n_bins*n_bootstraps) matrix specifying results from bootstrapped random walks.
#' @param sig_bins A named list describing the length of each perturbation.
#' @param percentile A double between 0,1 indicating the proportion cutoff for bootstrap-based signature derivation.
#' @param limit Number of genes to keep in the output, or a vector of lengths. Default is 30.
#' @return A named list of genesets. Each list element is the recontextualized signature for that seed.
extract_sig_mat <- function(mat,
                            bootstraps = NULL,
                            sig_bins = NULL,
                            percentile = 0.99,
                            limit = 30) {

  sigs <- list()

  if (length(limit) == 1) {
    limits <- rep(limit, ncol(mat))
  } else {
    limits <- limit
  }

  if(!is.null(bootstraps)) {
    # Extract top n signatures based on gene-level comparison with bootstraps
    percentiles <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
    colnames(percentiles) <- colnames(mat)
    rownames(percentiles) <- rownames(mat)
    intervals <- unique(colnames(bootstraps))
    for(i in 1:nrow(mat)) {
      for(j in 1:ncol(mat)) {
        # Finding the matching bootstraps depending on the length of pb
        pb_name <- colnames(percentiles)[j]
        sig_bin <- sig_bins[[pb_name]]
        bounds <- t(sapply(strsplit(intervals, "-"), as.numeric))
        interval <- intervals[which(sig_bin >= bounds[,1] & sig_bin <= bounds[,2])]
        bootstrap_idx <- which(colnames(bootstraps) == interval)

        percentiles[i,j] <- mean(bootstraps[i,bootstrap_idx] <= mat[i,j])
      }
    }
    percentiles[percentiles < percentile] <- 0

    for (col in 1:ncol(mat)) {
      colname <- colnames(percentiles)[[col]]
      mat_col <- percentiles[, col]
      mat_col <- mat_col[mat_col != 0]
      sig <- sort(mat_col, decreasing = TRUE) %>% names()
      sigs[[colname]] <- sig
    }

  } else {
    # Manually extract the 'top n' signatures set by limit parameter.
    for (col in 1:ncol(mat)) {
      colname <- colnames(mat)[[col]]
      mat_col <- mat[, col]
      limit <- limits[[col]]
      sig <- sort(mat_col, decreasing = TRUE) %>%
        head(limit) %>%
        names()
      sigs[[colname]] <- sig
    }
  }
  return(sigs)
}


#' Returns a dataframe object with the stationary probability value and column indicating whether gene was a seed
#' @param prob_vec (n_gene, 1) matrix
#' @param seeds character vector
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


#' Finds a simulated network signature
#'
#' @param ig network given as an igraph
#' @param seeds Either a single unnamed gene "TP53", a named list of genes, or a list of named lists of genes.
#' @param sig A string specifying the type of network signature: random walk, correlation etc.
#' @param avg_p A boolean specifying whether to ensemble random walk results over a range of restart values
#' @param avg_p_vals A numeric vector specifying the start and end of a arithmetic sequence to explore restart values.
#' @param avg_p_length A numeric specifying how many values within `avg_p_vals` to include in the ensemble
#' @param p A numeric specifying the restart value for random walk, default=0.1
#' @param bootstrap A boolean specifying whether to use empirical distributions of stationary values to find significant genes.
#' @param n_bootstraps A numeric specifying the number of bootstraps to perform.
#' @param limit A numeric specifying the number of genes to be included in the network signature. Default is 30.
#'
#' @return vector of gene strings
#'
#' @importFrom igraph as_adjacency_matrix
#' @importFrom abind abind
#' @importFrom parallel mclapply
#' @export
network_sig <- function(ig,
                        seeds,
                        sig = c("corr", "rwr"),
                        avg_p = FALSE,
                        avg_p_vals = c(1e-4, 1e-1),
                        avg_p_length = 5,
                        p = 0.1,
                        bootstrap = FALSE,
                        n_bootstraps = 1000,
                        limit = 30) {
  stopifnot(is(seeds, "character") | is(seeds, "list"))
  stopifnot(is(ig, "igraph"))

  gene_names <- igraph::V(ig)$name
  all_genes <- unname(unlist(seeds))
  seed_genes_filter <- all(all_genes %in% gene_names)
  if (!seed_genes_filter) {
    print("Filtering Seed Genes to those contained in the graph.")
    seeds <- lapply(seeds, function(x) x[x %in% gene_names])
  }

  if (sig == "corr") {
    if ("weight" %in% list.edge.attributes(ig)) {
      cor_mat <- igraph::as_adjacency_matrix(ig, attr = "weight")
    } else {
      cor_mat <- igraph::as_adjacency_matrix(ig, attr = NULL)
    }
    net_sig <- correlated_sigs(corr_mat = cor_mat, seeds = seeds, limit = limit)
  } else if (sig == "rwr") {
    obs_mat <- rwr_mat(ig,
                       seeds,
                       restart = p,
                       avg_p = avg_p,
                       avg_p_vals = avg_p_vals,
                       avg_p_length = avg_p_length)
    if(bootstrap & (p != 1)) {
      # If bootstrap matrix is (n_genes, n_bootstraps x n_bins)
      mat_bootstraps <- rwr_mat(ig,
                                seeds,
                                restart = p,
                                avg_p = avg_p,
                                avg_p_vals = avg_p_vals,
                                avg_p_length = avg_p_length,
                                bootstrap = bootstrap,
                                n_bootstraps = n_bootstraps) %>% as.matrix()

      sig_bins <- lapply(seeds, length)
      net_sig <- extract_sig_mat(obs_mat, bootstraps = mat_bootstraps, sig_bins = sig_bins)
    } else {
      net_sig <- extract_sig_mat(obs_mat, bootstraps = NULL, limit = limit)
    }
  }
  return(net_sig)
}


#' Recontextualize seed signatures with correlation based neighbors
#'
#' @param corr_mat Correlation Matrix
#' @param seeds Either a single unnamed gene "TP53", a named list of genes, or a list of named lists of genes.
#' @param limit Number of genes to keep in the output, or a vector of lengths. Default is 30.
correlated_sigs <- function(corr_mat, seeds, limit = 30) {

  # Seeds can be a single character vector. In that case need to list-ify it.
  if (is(seeds, "character") && length(seeds) == 1) {
    seeds <- list(seeds)
    names(seeds) <- seeds
  } else if (is.null(names(seeds))) {
    stop("Seed signature needs name")
  }

  all_genes <- unname(unlist(seeds))
  seed_filter <- all_genes %in% rownames(corr_mat)

  stopifnot("Seed Genes are not all contained in the Correlation matrix." = all(seed_filter))

  top_corrs <- list()
  for (gs_name in names(seeds)) {
    geneset <- seeds[[gs_name]]

    if (length(geneset) > 1) {
      mean_corrs <- apply(corr_mat[, geneset], 1, mean)
      top_corr <- sort(mean_corrs, decreasing = TRUE)[1:limit]
      top_corr <- names(top_corr)
      top_corrs[[gs_name]] <- top_corr
    } else {
      top_corr <- sort(corr_mat[, geneset], decreasing = TRUE)[1:limit]
      top_corrs[[gs_name]] <- names(top_corr)
    }
  }

  return(top_corrs)
}

#' Recontextualize seed signatures with correlation based neighbors
#'
#' @param corr_mats List of correlation matrices
#' @param seeds List of seeds, length has to match corr_mats
#' @param limit Number of genes to keep in the output, or a vector of lengths. Default is 30.
v.correlated_sigs <- function(corr_mats, seeds, limit = 30) {
  # browser()
  sigs <- list()
  for (i in seq_along(corr_mats)) {
    corr_mat <- corr_mats[[i]]
    name <- names(corr_mats)[[i]]
    seed <- seeds[i]

    top_corrs <- correlated_sigs(corr_mat, seed, limit)

    if (is.null(name)) {
      sigs[[i]] <- top_corrs
    } else {
      sigs[[name]] <- top_corrs
    }
  }
  return(sigs)
}

#' Perform a random walk with restart (personalized page rank) on an igraph given a seed matrix, and return stationary probabilties.
#' Stripped down and corrected version of dnet: https://rdrr.io/cran/dnet/src/R/dRWR.r
#'
#' @param ig igraph object
#' @param seed_mat (Gene, num_seeds) matrix with prior weights for each gene in a seed set. See seed_matrix.
#' @param restart the restart probability for RWR
#' @param epsilon Exploration factor
#' @param normalize Normalization strategy
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

