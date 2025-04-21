library(WGCNA)
library(Biobase)
library(SummarizedExperiment)

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

#' Perform Stratified Sampling on an ExpressionSet
#'
#' @description
#' This function performs stratified sampling on an ExpressionSet object, splitting it into training and test sets.
#'
#' @param eset An ExpressionSet object containing gene expression data and phenotype data.
#' @param stratify_column The name of the column in pData(eset) to use for stratification.
#' @param sample_id_column The name of the column in pData(eset) containing sample IDs.
#' @param p The proportion of samples to include in the training set. Default is 0.5.
#'
#' @return A list containing two ExpressionSet objects:
#'   \item{train}{The ExpressionSet containing the training samples}
#'   \item{test}{The ExpressionSet containing the test samples}
#'
#' @details
#' The function stratifies the samples based on the specified column and then randomly selects
#' a proportion p of samples from each stratum for the training set. The remaining samples
#' form the test set.
#'
#' @importFrom Biobase pData
#' @importFrom dplyr group_by sample_frac ungroup pull
eset_stratified_sample <- function(eset, stratify_column, sample_id_column, p=0.5) {
  print(stratify_column)
  pdata <- pData(eset)
  pdata_strat <- pdata %>% group_by({{stratify_column}}) %>% sample_frac(p) %>% ungroup()
  sample_ids <- pdata_strat %>% pull({{sample_id_column}})
  sample_filter <- eset[[sample_id_column]] %in% sample_ids
  train_eset <- eset[,sample_filter]
  test_eset <- eset[,!sample_filter]
  return(list(train = train_eset, test = test_eset))
}

#' Split Dataset, Learn Network, and Save Results
#'
#' @description
#' This function splits an expression set, learns a gene network using a specified method,
#' and saves the results for both training and test sets.
#'
#' @param eset An ExpressionSet object containing gene expression data.
#' @param stratify_col Column name in eset's pData for stratification.
#' @param sample_id_col Column name in eset's pData for sample IDs.
#' @param path Directory path to save the results.
#' @param name Base name for the output files.
#' @param igraph If TRUE, returns igraph objects instead of matrices. Default is FALSE.
#' @param diag_zero If TRUE, sets the diagonal of adjacency matrices to zero. Default is FALSE.
#' @param p Proportion of samples to include in the training set. Default is 0.5.
#' @param n_split Number of times to repeat the split-learn-save process. Default is 10.
#' @param cor.type Type of correlation for WGCNA. Options are "unsigned", "signed hybrid", or "signed". Default is "unsigned".
#' @param learn Learning method to use. Currently only "wgcna" is implemented. Default is "wgcna".
#'
#' @details
#' For each split:
#' 1. The eset is split into training and test sets using stratified sampling.
#' 2. A network is learned on both the training and test sets using the specified method.
#' 3. The resulting networks are saved as RDS files.
#'
#' @note
#' This function requires the eset_stratified_sample and wgcna.adj functions to be available.
#' Currently, only the WGCNA method is fully implemented.
split_learn_save <- function(eset,
                             stratify_col,
                             sample_id_col,
                             path,
                             name,
                             igraph=FALSE,
                             diag_zero=FALSE,
                             p=0.5,
                             n_split=10,
                             cor.type=c("unsigned", "signed hybrid", "signed"),
                             learn=c("wgcna","megena","aracne","shine")) {

  learn <- match.arg(learn)
  cor.type <- match.arg(cor.type)

  for (i in 1:n_split) {
    # Split eset
    esets <- eset_stratified_sample(eset, stratify_col, sample_id_col, p)
    train_eset <- esets$train
    test_eset <- esets$test

    if(learn=="wgcna") {
      # Learn wgcna adjacency matrix on both splits
      train_ig <- wgcna.adj(train_eset, cor.type = cor.type, igraph = igraph, diag_zero = diag_zero)
      test_ig <- wgcna.adj(test_eset, cor.type = cor.type, igraph = igraph, diag_zero = diag_zero)
      igs <- list(train = train_ig, test=test_ig)

      saveRDS(igs, file=file.path(path, paste0(name, "_", i,".rds")))
    }
  }
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
#' @importFrom WGCNA pickSoftThreshold.fromSimilarity adjacency.fromSimilarity
#' @importFrom doParallel registerDoParallel
wgcna.power <- function(cor_mat,
                        cores=1,
                        diag_zero=TRUE) {

  # Set parallel computing environment
  doParallel::registerDoParallel(cores=cores)

  # Pick soft threshold via scale-free fit
  sft <- WGCNA::pickSoftThreshold.fromSimilarity(similarity=cor_mat)

  # Check selected power
  beta <- sft.check(sft)
  print(beta)

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
#' @param eset An ExpressionSet or SummarizedExperiment object containing gene expression data.
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
#' @importFrom WGCNA bicor adjacency pickSoftThreshold
#' @importFrom doParallel registerDoParallel
#' @importFrom Biobase exprs
#' @importFrom SummarizedExperiment assays
#' @importFrom igraph graph_from_adjacency_matrix
#'
#' @export
wgcna.adj <- function(eset,
                      min.sft=0.85,
                      beta=NULL,
                      cores=1,
                      cor.fn=c("bicor", "cor"),
                      cor.type=c("unsigned", "signed hybrid", "signed"),
                      powers=c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                      igraph=FALSE,
                      diag_zero=FALSE) {
  # Ad hoc namespace changes
  bicor = WGCNA::bicor
  cor = WGCNA::bicor

  # Handle arguments
  args <- as.list(environment())
  print(args)
  cor.fn <- match.arg(cor.fn)
  cor.type <- match.arg(cor.type)
  print(cor.fn)
  # Correlation options
  if (cor.fn == "cor") cor.options = list(use="p")
  if (cor.fn == "bicor") cor.options = list(pearsonFallback="individual")
  print(cor.options)

  # Set parallel computing environment
  doParallel::registerDoParallel(cores=cores)

  # Format expression set into sample x gene matrix
  if(is(eset,"SummarizedExperiment")) {
    dat <- t(SummarizedExperiment::assays(eset)[[1]])
  } else {
    dat <- t(Biobase::exprs(eset))
  }

  # Pick soft threshold via scale-free fit
  if (is.null(beta)) {
    sft <- WGCNA::pickSoftThreshold(data=dat,
                                    corFnc=cor.fn,
                                    RsquaredCut=min.sft,
                                    powerVector=powers)

    # Check selected power
    beta <- sft.check(sft)
    print(beta)
  }

  # Construct co-expression similarity
  adj <- WGCNA::adjacency(datExpr=dat,
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


#' Compute Sparse Partial Correlation Network with SILGGM
#'
#' This function performs sparse partial correlation estimation using the SILGGM package,
#' with optional p-value filtering, WGCNA power transformation, and igraph conversion.
#'
#' @param mat A matrix with rows as features and columns as genes
#' @param method Correlation estimation method for SILGGM. Default "B_NW_SL".
#'               See [SILGGM::SILGGM()] for options.
#' @param wgcna_power Logical indicating whether to apply WGCNA soft power thresholding.
#'                    Default TRUE.
#' @param pval_filter Logical indicating whether to filter edges by p-value. Default TRUE.
#' @param pos_filter Logical indicating whether to filter out negative edges. Default TRUE.
#' @param pval Significance threshold for edge filtering. Default 0.05.
#' @param igraph Logical indicating whether to return an igraph object. Default TRUE.
#'
#' @return Either:
#' - Weighted adjacency matrix (if `igraph = FALSE`)
#' - igraph graph object (if `igraph = TRUE`)
#'
#' @details The function performs these steps:
#' 1. Estimates sparse partial correlations using SILGGM
#' 2. Optionally filters edges by non-negativity
#' 3. Optionally filters edges by p-value significance
#' 4. Optionally applies WGCNA soft thresholding
#' 5. Converts to igraph object if requested
#'
#' @importFrom SILGGM SILGGM
#' @importFrom igraph graph_from_adjacency_matrix
#' @export
silggm.adj <- function(mat,
                       method = "B_NW_SL",
                       wgcna_power = TRUE,
                       pval_filter = TRUE,
                       pos_filter = TRUE,
                       pval = 0.05,
                       igraph = TRUE) {
  stopifnot(is(mat, "matrix"))

  silggm_res <- SILGGM::SILGGM(mat, method = method)
  silggm_mat <- as.matrix(silggm_res$partialCor)
  rownames(silggm_mat) <- colnames(silggm_mat) <- colnames(mat)

  if(pos_filter) {
    silggm_mat[silggm_mat <= 0] <- 0
  }

  if (pval_filter) {
    silggm_pval <- as.matrix(silggm_res$p_partialCor)
    silggm_pval_bool <- silggm_pval < p_val
    silggm_mat[!silggm_pval_bool] <- 0
  }

  if (wgcna_power) {
    # WGCNA scale free estimation requires positive edges only.
    silggm_mat[silggm_mat <= 0] <- 0
    silggm_mat <- wgcna.power(silggm_mat)
  }

  if(igraph) {
    silggm_mat <- igraph::graph_from_adjacency_matrix(silggm_mat, weighted=TRUE, mode="undirected")
  }
  return(silggm_mat)
}

