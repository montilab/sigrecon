library(tidyverse)

#' Count number of displaced seeds
#' Note: This not commutative. Num_displaced is the number of genes that are from the original seed that are not in the recon.
#' @param recon_sig Named list of genesets
#' @param seed_sig Named list of genesets
#' @return Named list for which each entry is a list with two elements: displaced seed genes, and non-displaced seed genes
#'
#' @export
count_displaced_genes <- function(recon_sig, seed_sig) {
  # Two lists of genesets should have the same names
  stopifnot(all(names(recon_sig) == names(seed_sig)))

  displaced <- list()
  for (name in names(seed_sig)) {
    seeds_in_recon <- seed_sig[[name]] %in% recon_sig[[name]]
    n_displaced <- sum(!seeds_in_recon)
    n_not_displaced <- sum(seeds_in_recon)
    displaced[[name]] <- list(displaced = n_displaced, not_displaced = n_not_displaced)
  }

  return(displaced)
}

#' Evaluate recontextualized signatures
#'
#' @param ig igraph
#' @param seed_name name of list of genesets to be used to label dataframe
#' @param source_sigs list of genesets (to be recontextualized)
#' @param dest_sigs list of genesets (the ground truth)
#' @param split split number for CCLE TCGA experiment
#' @param experiment benchmarking task type
#' @param restart restart value
#' @param recon Boolean indicating whether recontextualization occurs
#' @param use_weights Boolean indicating whether to use weights in the KS Test
#' @param weights.pwr Power to raise weights to
#' @param normalize Normalization strategy to employ
#' @param limit number of genes to be included in the network signature, default=30
#'
#' @export
recon_eval_df <- function(ig,
                          seed_name,
                          source_sigs,
                          dest_sigs,
                          split = 0,
                          experiment = c("CCLE", "perturb-seq", "drugmatrix", "sciplex"),
                          restart = 0.75,
                          recon = TRUE,
                          use_weights = TRUE,
                          weights.pwr = 1,
                          normalize = c("row", "column", "laplacian"),
                          limit = 30) {

  normalize <- match.arg(normalize)

  # NA Checks
  source_na_filter <- lapply(source_sigs, function(x) sum(is.na(x))) == 0
  source_sigs <- source_sigs[source_na_filter]

  # Read Igraph
  dest_net <- ig

  # Read Signatures
  experiment <- match.arg(experiment)

  if (experiment == "CCLE") {
    dest_train_net <- dest_net$train
    dest_test_net <- dest_net$test
    dest_kd_sigs <- dest_sigs$test
    dest_kd_sigs <- dest_kd_sigs[names(source_sigs)]

    n_genes <- length(V(dest_train_net))
  } else if (experiment == "drugmatrix" | experiment == "perturb-seq" | experiment == "sciplex") {
    dest_full_sigs <- lapply(dest_sigs, function(x) x[["up_full"]])
    dest_kd_sigs <- lapply(dest_sigs, function(x) x[["up"]])
    source_sigs <- lapply(source_sigs, function(x) x[["up"]])
  }

  # Check if source signatures are obtained on the same set of genes as the dest signatures
  stopifnot(all.equal(names(source_sigs), names(dest_kd_sigs)))

  # Find full destination correlation sig for ks tests
  if (experiment == "CCLE") {
    dest_seeds <- names(source_sigs)
    dest_seeds <- as.list(dest_seeds)
    names(dest_seeds) <- dest_seeds
    dest_full_sigs <- network_sig(net = dest_test_net, seeds = dest_seeds, "corr", limit = n_genes)

    # Sanity check that re-obtaining signatures from test network are the same as those generated apriori
    stopifnot(all(unlist(lapply(names(dest_full_sigs), function(x) all.equal(dest_full_sigs[[x]][1:limit], dest_kd_sigs[[x]])))))
  }

  # Baseline: Without Recontextualization
  if (!recon) {
    displaced <- count_displaced_genes(recon_sig = dest_kd_sigs, seed_sig = source_sigs)
    n_displaced <- unname(unlist(lapply(displaced, function(x) x$displaced)))
    n_not_displaced <- unname(unlist(lapply(displaced, function(x) x$not_displaced)))
    jacc_source_dest <- vjaccard(source_sigs, dest_kd_sigs)
    ks_obj <- v.ks.test(data_vecs = source_sigs,
                        ref_vecs = dest_full_sigs,
                        use_weights = use_weights,
                        weights.pwr = weights.pwr)
    Ds <- ks_obj$D_stats
    p_vals <- ks_obj$p_vals

  } else {

    # Stationary Probability matrix (genes, n_knockdown)
    if (experiment == "CCLE") {
      source_kd_recon_probs <- rwr_mat(ig = dest_train_net, seeds = source_sigs, restart = restart, normalize = normalize)
    } else if (experiment == "perturb-seq" | experiment == "drugmatrix" | experiment == "sciplex") {
      source_kd_recon_probs <- rwr_mat(ig = dest_net, seeds = source_sigs, restart = restart, normalize = normalize)
    }

    n_ranks <- nrow(source_kd_recon_probs)

    # Recontextualized Signature
    if (experiment == "CCLE") {
      n_test_sig <- length(dest_kd_sigs[[1]])
    } else {
      n_test_sig <- lapply(source_sigs, length)
    }
    source_kd_recon_sigs_jacc <- top_n_mat(source_kd_recon_probs, n_test_sig)
    saveRDS(source_kd_recon_sigs_jacc, file.path(PATH, paste0("data/sigs/",experiment,"/recon/",seed_name,".rds")))

    # Jaccard Similarity (Top n vs Top n)
    jacc_source_dest <- vjaccard(source_kd_recon_sigs_jacc, dest_kd_sigs)

    # Displacement
    displaced <- count_displaced_genes(recon_sig = source_kd_recon_sigs_jacc, seed_sig = source_sigs)
    n_displaced <- unname(unlist(lapply(displaced, function(x) x$displaced)))
    n_not_displaced <- unname(unlist(lapply(displaced, function(x) x$not_displaced)))

    # KS.Test (Ranks of Top n recontextualized compared to Full list of dest)
    stopifnot(all.equal(names(source_kd_recon_sigs_jacc), names(dest_full_sigs)))
    ks_obj <- v.ks.test(data_vecs = source_kd_recon_sigs_jacc,
                        ref_vecs = dest_full_sigs,
                        use_weights = use_weights,
                        weights.pwr = weights.pwr)
    Ds <- ks_obj$D_stats
    p_vals <- ks_obj$p_vals
  }

  sig_jacc_df <- data.frame(
    cancer = seed_name,
    split = split,
    displaced = n_displaced,
    kept = n_not_displaced,
    rw_p = restart,
    gene = names(dest_kd_sigs),
    jacc = jacc_source_dest,
    k_d = Ds,
    k_p = p_vals
  )

  return(sig_jacc_df)
}


#' Generate Network Signatures from a Path
#'
#' This function reads a network object from a file path and generates network signatures based on provided seeds.
#'
#' @param path A character string specifying the file path to the network object (RDS file).
#' @param seeds A character vector or list of seed genes/nodes.
#' @param sig The signature generation method. Either "corr" for correlation-based or "rwr" for random walk with restart. Default is c("corr", "rwr").
#' @param p A numeric value specifying the restart probability for random walk. Default is 0.1.
#' @param limit An integer specifying the maximum number of nodes to include in the signature. Default is 30.
#'
#' @return A list of network signatures. If the input is a single network, returns a single signature list. If the input is a list of networks, returns a list of signature lists.
#'
#' @export
network_sig_path <- function(path,
                             seeds,
                             sig = c("corr", "rwr"),
                             p = 0.1,
                             limit = 30) {
  net_object <- readRDS(path)

  stopifnot(is(seeds, "character")| is(seeds, "list"))

  if (length(seeds) > 1 & is.null(names(seeds))) {
    stop("Seed signature needs name")
  }

  if (is(net_object, "igraph")) {
    sigs <- network_sig(net_object, seeds, sig, p, limit)
  } else if (is(net_object, "list") || is(net_object, "character")) {
    sigs <- lapply(net_object, function(x) network_sig(x, seeds, sig, p, limit))
  }
  return(sigs)
}


#' Finds a simulated network signature
#'
#' @param net network given as an igraph
#' @param seeds single gene string or vector of gene strings or named list of sets of gene symbols
#' @param sig String, network signature type: random walk, correlation, nearest neighbor
#' @param p Numeric, restart value for random walk, default=0.1
#' @param limit number of genes to be included in the network signature, default=30
#'
#' @return vector of gene strings
#'
#' @importFrom igraph as_adj
#' @export
network_sig <- function(net,
                        seeds,
                        sig = c("corr", "rwr"),
                        p = 0.1,
                        limit = 30) {
  stopifnot(is(seeds, "character")| is(seeds, "list"))

  if (sig == "corr") {
    if ("weight" %in% list.edge.attributes(net)) {
      cor_mat <- igraph::as_adj(net, attr = "weight")
    } else {
      cor_mat <- igraph::as_adj(net, attr = NULL)
    }
    net_sig <- correlated_sigs(corr_mat = cor_mat, seeds = seeds, limit = limit)
  } else if (sig == "rwr") {
    mat <- rwr_mat(net, seeds, restart = p)
    net_sig <- top_n_mat(mat, limit = limit)
  }
  return(net_sig)
}


#' Recontextualize seed signatures with correlation based neighbors
#'
#' @param corr_mat Correlation Matrix
#' @param seeds Seed or List of seeds
#' @param limit Numeric indicating the number of genes to be returned
#'
#' @export
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
#' @param limit Numeric indicating the number of genes to be returned
#'
#' @export
vcorrelated_sigs <- function(corr_mats, seeds, limit = 30) {
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
