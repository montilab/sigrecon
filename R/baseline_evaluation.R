library(tidyverse)

#' Evaluate recontextualized signatures
#'
#' @param ig igraph
#' @param seed_name name of list of genesets to be used to label dataframe
#' @param source_sigs list of genesets (to be recontextualized)
#' @param dest_sigs list of genesets (the ground truth)
#' @param restart restart value
#' @param avg_p A boolean specifying whether to ensemble random walk results over a range of restart values
#' @param avg_p_vals A numeric vector specifying the start and end of a geometric sequence to explore restart values.
#' @param avg_p_length A numeric specifying how many values within `avg_p_vals` to include in the ensemble
#' @param recon Boolean indicating whether recontextualization occurs
#' @param use_weights Boolean indicating whether to use weights in the KS Test
#' @param weights.pwr Power to raise weights to
#' @param normalize Normalization strategy to employ
#' @param save_path file path to save recontextualized signatures
#' @param limit number of genes to be included in the network signature, default=30
#'
#' @export
recon_eval_df <- function(ig,
                          seed_name,
                          source_sigs,
                          dest_sigs,
                          restart = 0.75,
                          avg_p = FALSE,
                          avg_p_vals = c(1e-4, 1e-1),
                          avg_p_length = 5,
                          recon = TRUE,
                          use_weights = TRUE,
                          weights.pwr = 1,
                          normalize = c("row", "column", "laplacian"),
                          save_path = "",
                          limit = 30) {

  normalize <- match.arg(normalize)
  dest_net <- ig

  dest_full_sigs <- lapply(dest_sigs, function(x) x[["up_full"]])
  dest_short_sigs <- lapply(dest_sigs, function(x) x[["up"]])
  source_sigs <- lapply(source_sigs, function(x) x[["up"]])

  # Check if source signatures are obtained on the same set of genes as the dest signatures
  stopifnot(all.equal(names(source_sigs), names(dest_short_sigs)))

  # Baseline: Without Recontextualization
  if (!recon) {
    displaced <- count_displaced_genes(recon_sig = dest_short_sigs, seed_sig = source_sigs)
    n_displaced <- unname(unlist(lapply(displaced, function(x) x$displaced)))
    n_not_displaced <- unname(unlist(lapply(displaced, function(x) x$not_displaced)))
    jacc_source_dest <- v.jaccard(source_sigs, dest_short_sigs)
    ks_obj <- v.ks.test(data_vecs = source_sigs,
                        ref_vecs = dest_full_sigs,
                        use_weights = use_weights,
                        weights.pwr = weights.pwr)
    Ds <- ks_obj$D_stats
    p_vals <- ks_obj$p_vals

  } else {
    # Recontextualized Signature
    recon_sigs <- network_sig(ig = ig,
                              seeds = source_sigs,
                              sig = "rwr",
                              avg_p = avg_p,
                              avg_p_vals = avg_p_vals,
                              avg_p_length = avg_p_length,
                              p = restart,
                              limit = limit)
    saveRDS(recon_sigs, file.path(save_path, paste0(seed_name,".rds")))

    # Jaccard Similarity
    jacc_source_dest <- v.jaccard(recon_sigs, dest_short_sigs)

    # Displacement
    stopifnot(all.equal(names(recon_sigs), names(source_sigs)))
    displaced <- count_displaced_genes(recon_sig = recon_sigs, seed_sig = source_sigs)
    n_displaced <- unname(unlist(lapply(displaced, function(x) x$displaced)))
    n_not_displaced <- unname(unlist(lapply(displaced, function(x) x$not_displaced)))

    # KS.Test (Ranks of Top n recontextualized compared to Full list of dest)
    stopifnot(all.equal(names(recon_sigs), names(dest_full_sigs)))
    ks_obj <- v.ks.test(data_vecs = recon_sigs,
                        ref_vecs = dest_full_sigs,
                        use_weights = use_weights,
                        weights.pwr = weights.pwr)
    Ds <- ks_obj$D_stats
    p_vals <- ks_obj$p_vals
  }

  eval_df <- data.frame(
    source = seed_name,
    displaced = n_displaced,
    kept = n_not_displaced,
    rw_p = restart,
    gene = names(dest_short_sigs),
    jacc = jacc_source_dest,
    k_d = Ds,
    k_p = p_vals
  )

  return(eval_df)
}
