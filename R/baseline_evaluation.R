library(tidyverse)

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
    jacc_source_dest <- v.jaccard(source_sigs, dest_kd_sigs)
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
    jacc_source_dest <- v.jaccard(source_kd_recon_sigs_jacc, dest_kd_sigs)

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
