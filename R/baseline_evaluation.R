library(tidyverse)

#' @title Evaluate Signature Prediction
#'
#' @description This function calculates various evaluation metrics for predicted gene signatures
#' against true (ground truth) gene signatures, given a source signature as input.
#' Metrics include Jaccard index, gene displacement, and Kolmogorov-Smirnov (KS) test statistics
#' for rank agreement.
#'
#' @param source_sigs A named list of genesets.
#' @param pred_sigs A named list of genesets.
#' @param true_sigs A named list of genesets containing both the cutoff signature `up` and the full signature `up_full`.
#' @param source A character string describing the starting biological context.
#'
#' @return A data frame with evaluation results. Each row corresponds to a perturbation and context.
#'
#' @export
sig_eval_table <- function(source_sigs,
                           pred_sigs,
                           true_sigs,
                           source = "source_context") {
  dest_short_sigs <- lapply(true_sigs, function(x) x$up)
  dest_full_sigs <- lapply(true_sigs, function(x) x$up_full)
  jacc_source_dest <- v.jaccard(pred_sigs, dest_short_sigs)

  # Displacement
  displaced <- count_displaced_genes(recon_sig = pred_sigs, seed_sig = source_sigs)
  n_displaced <- unname(unlist(lapply(displaced, function(x) x$displaced)))
  n_not_displaced <- unname(unlist(lapply(displaced, function(x) x$not_displaced)))

  # KS.Test (Ranks of Top n recontextualized compared to Full list of dest)
  ks_obj <- v.ks.test(data_vecs = pred_sigs,
                      ref_vecs = dest_full_sigs)
  Ds <- ks_obj$D_stats
  p_vals <- ks_obj$p_vals

  eval_df <- data.frame(
    source = source,
    displaced = n_displaced,
    kept = n_not_displaced,
    gene = names(pred_sigs),
    jacc = jacc_source_dest,
    k_d = Ds,
    k_p = p_vals
  )
  return(eval_df)
}

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
#' @param bootstrap A boolean specifying whether to use empirical distributions of stationary values to find significant genes.
#' @param n_bootstraps A numeric specifying the number of bootstraps to perform.
#' @param recon Boolean indicating whether recontextualization occurs
#' @param use_weights Boolean indicating whether to use weights in the KS Test
#' @param weights.pwr Power to raise weights to
#' @param normalize Normalization strategy to employ
#' @param save Boolean indicating whether to save recontextualized signatures
#' @param save_path file path to save recontextualized signatures
#' @param limit Number of genes to keep in the output, or a vector of lengths
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
                          bootstrap = FALSE,
                          n_bootstraps = 1000,
                          recon = TRUE,
                          use_weights = TRUE,
                          weights.pwr = 1,
                          normalize = c("row", "column", "laplacian"),
                          save = FALSE,
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
                              bootstrap = bootstrap,
                              n_bootstraps = n_bootstraps,
                              limit = limit)
    if(save) {
      saveRDS(recon_sigs, file.path(save_path, paste0(seed_name,".rds")))
    }

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

#' Combine P-values Using Fisher's Method
#'
#' @description
#' Combines multiple p-values from independent tests into a single meta p-value
#' using Fisher's method (also known as Fisher's combined probability test).
#'
#' @details
#' Fisher's method combines p-values by computing the test statistic:
#' \deqn{X = -2 \sum_{i=1}^{k} \ln(p_i)}
#'
#' Under the null hypothesis (all individual null hypotheses are true), X follows
#' a chi-squared distribution with 2k degrees of freedom, where k is the number
#' of p-values. This method assumes that the tests are independent.
#'
#' The combined p-value represents the probability of observing the given set of
#' p-values (or more extreme) if all null hypotheses are true.
#'
#' @param pvals Numeric vector of p-values to combine. All values must be between
#'   0 and 1 (exclusive of 0). NA values are not allowed.
#'
#' @return A single numeric value representing the combined p-value.
#'
#' @examples
#' # Combine three p-values
#' fishers_meta_p(c(0.01, 0.03, 0.25))
#'
#' # Two significant p-values
#' fishers_meta_p(c(0.001, 0.005))
#'
#' # Mix of significant and non-significant
#' fishers_meta_p(c(0.045, 0.23, 0.67))
#'
#' @export
fishers_meta_p <- function(pvals) {
  # pvals: a vector of p-values, e.g. c(0.01, 0.03, 0.25)
  X <- -2 * sum(log(pvals))
  df <- 2 * length(pvals)
  return(pchisq(X, df=df, lower.tail = FALSE))
}
