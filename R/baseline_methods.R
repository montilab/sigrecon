#' Reconstruct Gene Signatures Using Projection Scores
#'
#' @description
#' This function reconstructs gene signatures based on their correlation with
#' per-sample projection scores computed from the input signatures.
#'
#' @param se A SummarizedExperiment object containing gene expression data.
#' @param sigs A list of gene signatures, where each element is a character vector of gene names.
#' @param score Scoring method used to score samples against input signatures.
#'   Either `"gsva"`, `"AUCell"`, or `"eigen"`. Default is `"gsva"`.
#'
#' @return A list of reconstructed gene signatures, with the same structure as the input `sigs`.
#'
#' @details
#' The function performs the following steps:
#' 1. Calculates projection scores for the input signatures with GSVA, AUCell, or eigengenes.
#' 2. Computes the correlation between gene expression and projection scores.
#' 3. Ranks genes based on their correlation with each signature's projection scores.
#' 4. Selects the top-ranking genes to form new signatures of the same length as the original ones.
#'
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC getAUC
#' @importFrom GSVA gsva gsvaParam
#' @importFrom stats cor
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange slice pull desc
#'
#' @export
projectCor <- function(se,
                       sigs,
                       score = c("gsva", "AUCell", "eigen")) {
  score <- match.arg(score)
  stopifnot(is(se, "SummarizedExperiment"))

  expr_mat <- SummarizedExperiment::assay(se)

  if (score == "gsva") {
    score_param <- GSVA::gsvaParam(se, sigs, maxDiff = TRUE)
    score_res <- GSVA::gsva(score_param, verbose = FALSE)
    if (is(score_res, "SummarizedExperiment")) {
      score_mat <- SummarizedExperiment::assay(score_res)
    } else {
      score_mat <- score_res
    }
  } else if (score == "AUCell") {
    if (is(expr_mat, "sparseMatrix")) {
      expr_mat <- as.matrix(expr_mat)
    }

    rankings <- AUCell::AUCell_buildRankings(exprMat = expr_mat,
                                             plotStats = FALSE,
                                             verbose = FALSE)
    auc <- AUCell::AUCell_calcAUC(geneSets = sigs,
                                  rankings = rankings,
                                  verbose = FALSE)
    score_mat <- AUCell::getAUC(auc)
  } else {
    if (is(expr_mat, "sparseMatrix")) {
      expr_mat <- as.matrix(expr_mat)
    }

    gene_names <- rownames(expr_mat)
    if (is.null(gene_names)) {
      stop("Expression assay must have rownames to compute eigengene scores.")
    }

    sig_names <- names(sigs)
    eig_scores <- setNames(lapply(sig_names, function(sig_name) {
      genes <- intersect(sigs[[sig_name]], gene_names)
      if (length(genes) == 0) {
        message(sprintf(
          "Skipping signature '%s' because no genes were found in the expression assay.",
          sig_name
        ))
        return(NULL)
      }

      sig_expr <- expr_mat[genes, , drop = FALSE]
      gene_sds <- apply(sig_expr, 1, stats::sd)
      keep <- !is.na(gene_sds) & gene_sds > 0
      sig_expr <- sig_expr[keep, , drop = FALSE]

      if (nrow(sig_expr) == 0) {
        stop(sprintf("All genes in signature '%s' have zero variance.", sig_name))
      }

      sig_expr <- t(scale(t(sig_expr), center = TRUE, scale = TRUE))

      eigengene <- if (nrow(sig_expr) == 1) {
        as.numeric(sig_expr[1, ])
      } else {
        pca <- stats::prcomp(t(sig_expr), center = FALSE, scale. = FALSE)
        as.numeric(pca$x[, 1])
      }

      # Align the PC direction with the average standardized module signal.
      avg_signal <- colMeans(sig_expr)
      align_cor <- stats::cor(eigengene, avg_signal)
      if (!is.na(align_cor) && align_cor < 0) {
        eigengene <- -eigengene
      }

      eigengene
    }), sig_names)

    keep_scores <- !vapply(eig_scores, is.null, logical(1))
    eig_scores <- eig_scores[keep_scores]

    if (length(eig_scores) == 0) {
      stop("No signatures had genes present in the expression assay for eigengene scoring.")
    }

    score_mat <- do.call(rbind, eig_scores)
    rownames(score_mat) <- names(eig_scores)
    colnames(score_mat) <- colnames(expr_mat)
  }

  genes_data <- t(expr_mat)
  if (is(genes_data, "sparseMatrix")) {
    genes_data <- as.matrix(genes_data)
  }

  if (is(score_mat, "sparseMatrix")) {
    score_mat <- as.matrix(score_mat)
  }

  proj_scores <- t(score_mat)
  corr_mat <- stats::cor(genes_data, proj_scores, method = "pearson")
  if (is.null(dim(corr_mat))) {
    corr_mat <- matrix(
      corr_mat,
      ncol = 1,
      dimnames = list(colnames(genes_data), colnames(proj_scores))
    )
  }

  results <- tibble::as_tibble(corr_mat, rownames = "gene")
  new_sigs <- list()
  for(sig_name in colnames(proj_scores)) {
    # Obtain rank of genes most correlated to the projection scores
    results$rank <- rank(dplyr::desc(results[[sig_name]]))

    # Obtain new signature
    sig_length <- length(sigs[[sig_name]])
    results <- results %>% dplyr::arrange(rank)
    new_sig <- results %>% dplyr::slice(1:sig_length) %>% dplyr::pull(gene)
    new_sigs[[sig_name]] <- new_sig
  }

  return(new_sigs)
}

gsva_recon <- function(se,
                       sigs) {
  projectCor(se = se, sigs = sigs, score = "gsva")
}

.normalize_recontext_input <- function(x, arg_name) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.character(x) && length(x) == 1) {
    x <- stats::setNames(list(x), x)
  }

  if (!is.list(x)) {
    stop(sprintf("'%s' must be a named list or a single character string.", arg_name))
  }

  if (is.null(names(x)) || any(names(x) == "")) {
    stop(sprintf("'%s' must be a named list.", arg_name))
  }

  x
}

.recontext_input_or <- function(primary, fallback) {
  if (!is.null(primary)) {
    return(primary)
  }

  fallback
}

.resolve_recontext_limits <- function(items,
                                      limit = NULL,
                                      default = 30) {
  item_names <- names(items)
  n_items <- length(items)

  if (is.null(limit)) {
    limits <- lengths(items)
    limits[limits == 0] <- default
    names(limits) <- item_names
    return(limits)
  }

  if (length(limit) == 1) {
    limits <- rep(limit, n_items)
    names(limits) <- item_names
    return(limits)
  }

  if (!is.null(names(limit))) {
    if (!all(item_names %in% names(limit))) {
      stop("Named 'limit' must include every requested signature.")
    }

    limits <- unname(limit[item_names])
    names(limits) <- item_names
    return(limits)
  }

  if (length(limit) != n_items) {
    stop("'limit' must have length 1, match the number of signatures, or be a named vector.")
  }

  names(limit) <- item_names
  limit
}

.mean_recontextualize <- function(se,
                                  sigs = NULL,
                                  limit = NULL,
                                  perturbation_col,
                                  condition_col,
                                  perturbed_label = "Perturbed",
                                  control_label = "Control",
                                  design_vars = NULL,
                                  assay_name = NULL,
                                  alpha = 0.05,
                                  min_count = NULL,
                                  min_samples = 1) {
  stopifnot(is(se, "SummarizedExperiment"))

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("The 'DESeq2' package is required for method = 'mean'.")
  }

  if (is.null(assay_name)) {
    assay_names <- names(SummarizedExperiment::assays(se))
    assay_name <- assay_names[[1]]
  }

  if (is.null(assay_name) || !assay_name %in% names(SummarizedExperiment::assays(se))) {
    stop("'assay_name' must identify an assay present in 'se'.")
  }

  counts <- SummarizedExperiment::assays(se)[[assay_name]]
  if (is(counts, "sparseMatrix")) {
    counts <- as.matrix(counts)
  }

  if (!all(counts == round(counts))) {
    stop("The 'mean' method requires integer count data in the selected assay.")
  }

  coldata <- as.data.frame(SummarizedExperiment::colData(se))
  if (is.null(rownames(coldata))) {
    rownames(coldata) <- colnames(counts)
  }

  if (!all(colnames(counts) %in% rownames(coldata))) {
    stop("colData rownames must match assay column names for method = 'mean'.")
  }
  coldata <- coldata[colnames(counts), , drop = FALSE]

  needed_cols <- c(perturbation_col, condition_col, design_vars)
  missing_cols <- setdiff(needed_cols, colnames(coldata))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required colData columns for method = 'mean': %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  controls <- rownames(coldata)[as.character(coldata[[condition_col]]) == control_label]
  if (length(controls) == 0) {
    stop(sprintf("No control samples found where '%s' == '%s'.", condition_col, control_label))
  }

  sigs <- .normalize_recontext_input(sigs, "sigs")
  if (is.null(sigs)) {
    pert_names <- unique(as.character(
      coldata[[perturbation_col]][as.character(coldata[[condition_col]]) == perturbed_label]
    ))
    pert_names <- pert_names[!is.na(pert_names)]
    sigs <- stats::setNames(vector("list", length(pert_names)), pert_names)
  }

  limits <- .resolve_recontext_limits(sigs, limit = limit)
  mean_sigs <- list()

  for (pert_name in names(sigs)) {
    perturbed_samples <- rownames(coldata)[
      as.character(coldata[[condition_col]]) == perturbed_label &
        as.character(coldata[[perturbation_col]]) == pert_name
    ]

    if (length(perturbed_samples) == 0) {
      message(sprintf(
        "Skipping signature '%s' because no perturbed samples were found in 'se'.",
        pert_name
      ))
      next
    }

    sample_ids <- c(perturbed_samples, controls)
    subset_counts <- counts[, sample_ids, drop = FALSE]
    subset_meta <- coldata[sample_ids, , drop = FALSE]

    if (!is.null(min_count)) {
      keep <- rowSums(subset_counts >= min_count) >= min_samples
      subset_counts <- subset_counts[keep, , drop = FALSE]
    }

    if (nrow(subset_counts) == 0) {
      message(sprintf(
        "Skipping signature '%s' because all genes were filtered before DESeq2.",
        pert_name
      ))
      next
    }

    subset_meta$.recontext_condition <- factor(
      ifelse(rownames(subset_meta) %in% perturbed_samples, perturbed_label, control_label),
      levels = c(control_label, perturbed_label)
    )

    design_formula <- stats::reformulate(c(design_vars, ".recontext_condition"))

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = subset_counts,
      colData = subset_meta,
      design = design_formula
    )
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res <- DESeq2::results(
      dds,
      contrast = c(".recontext_condition", perturbed_label, control_label)
    )

    res_df <- tibble::as_tibble(as.data.frame(res), rownames = "gene")
    res_df <- res_df[!is.na(res_df$log2FoldChange), , drop = FALSE]
    res_df$padj_for_score <- res_df$padj
    res_df$padj_for_score[is.na(res_df$padj_for_score)] <- 1
    res_df$padj_for_score[res_df$padj_for_score <= 0] <- .Machine$double.xmin
    res_df$score <- res_df$log2FoldChange * (-log10(res_df$padj_for_score))

    sig_df <- res_df %>%
      dplyr::filter(log2FoldChange > 0) %>%
      dplyr::filter(padj <= alpha | is.na(alpha)) %>%
      dplyr::arrange(dplyr::desc(score), dplyr::desc(log2FoldChange), gene)

    mean_sigs[[pert_name]] <- sig_df %>%
      dplyr::slice(seq_len(min(limits[[pert_name]], nrow(sig_df)))) %>%
      dplyr::pull(gene)
  }

  mean_sigs
}

.build_networkprop_graph <- function(se,
                                     assay_name = NULL,
                                     nfeatures = NULL,
                                     min.sft = 0.85,
                                     beta = NULL,
                                     cores = 1,
                                     cor.fn = c("bicor", "cor"),
                                     cor.type = c("unsigned", "signed hybrid", "signed"),
                                     powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                                     diag_zero = TRUE) {
  stopifnot(is(se, "SummarizedExperiment"))

  if (is.null(assay_name)) {
    assay_names <- names(SummarizedExperiment::assays(se))
    assay_name <- assay_names[[1]]
  }

  if (is.null(assay_name) || !assay_name %in% names(SummarizedExperiment::assays(se))) {
    stop("'assay_name' must identify an assay present in 'se'.")
  }

  counts <- SummarizedExperiment::assays(se)[[assay_name]]
  if (is(counts, "sparseMatrix")) {
    counts <- as.matrix(counts)
  }

  if (is.null(rownames(counts))) {
    stop("The selected assay must have rownames for network construction.")
  }

  lib_sizes <- colSums(counts)
  if (any(lib_sizes <= 0)) {
    stop("All samples must have positive library sizes for method = 'networkProp'.")
  }

  norm_counts <- sweep(counts, 2, lib_sizes, "/")
  norm_counts <- log1p(norm_counts)

  network_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(log_norm = norm_counts),
    colData = SummarizedExperiment::colData(se)
  )

  if (is.null(nfeatures)) {
    nfeatures <- min(10000L, nrow(network_se))
  } else {
    nfeatures <- min(as.integer(nfeatures), nrow(network_se))
  }

  ranked_genes <- rank.var.eset(network_se)$values
  keep_genes <- ranked_genes[seq_len(nfeatures)]
  network_mat <- t(norm_counts[keep_genes, , drop = FALSE])

  wgcna.adj(
    mat = network_mat,
    min.sft = min.sft,
    beta = beta,
    cores = cores,
    cor.fn = match.arg(cor.fn),
    cor.type = match.arg(cor.type),
    powers = powers,
    igraph = TRUE,
    diag_zero = diag_zero
  )
}

#' Recontextualize signatures with a selected baseline method
#'
#' @param method Baseline recontextualization method. One of `"networkProp"`,
#'   `"projectCor"`, or `"mean"`.
#' @param se A SummarizedExperiment object used by `"projectCor"`, `"mean"`,
#'   and to learn the `"networkProp"` graph.
#' @param seeds Seed signatures for the `"networkProp"` method.
#' @param sigs Signature list for the `"projectCor"` method. For `"mean"`, names
#'   are interpreted as perturbation labels and lengths are used as output sizes
#'   when `limit` is not supplied.
#' @param score Scoring method used by `"projectCor"`.
#' @param sig Network signature mode passed to [network_sig()]. Defaults to `"rwr"`.
#' @param avg_p Whether to ensemble network propagation over multiple restart probabilities.
#' @param avg_p_vals Range of restart probabilities used when `avg_p = TRUE`.
#' @param avg_p_length Number of restart probabilities to average when `avg_p = TRUE`.
#' @param p Restart probability passed to [network_sig()].
#' @param bootstrap Whether to use bootstrap-based extraction in [network_sig()].
#' @param n_bootstraps Number of bootstrap replicates passed to [network_sig()].
#' @param limit Number of genes to keep for each output signature. If `NULL`,
#'   `"networkProp"` and `"projectCor"` keep the original signature lengths and
#'   `"mean"` uses the lengths of `sigs`.
#' @param nfeatures Number of variable genes to use when learning the
#'   `"networkProp"` graph from `se`. Defaults to `min(10000, nrow(assay(se)))`.
#' @param min.sft Minimum scale-free topology fitting index used by
#'   [wgcna.adj()] when learning the `"networkProp"` graph.
#' @param beta Optional soft-thresholding power passed to [wgcna.adj()] for
#'   `"networkProp"`.
#' @param cores Number of CPU cores passed to [wgcna.adj()] for `"networkProp"`.
#' @param cor.fn Correlation function passed to [wgcna.adj()] for `"networkProp"`.
#' @param cor.type Correlation network type passed to [wgcna.adj()] for `"networkProp"`.
#' @param powers Candidate power values passed to [wgcna.adj()] for `"networkProp"`
#'   when `beta` is `NULL`.
#' @param diag_zero Whether to zero the diagonal of the learned `"networkProp"`
#'   adjacency matrix before converting it to igraph.
#' @param perturbation_col Column in `colData(se)` identifying the perturbation
#'   label for each profile. Used by `"mean"`.
#' @param condition_col Column in `colData(se)` indicating whether a profile is
#'   perturbed or control. Used by `"mean"`.
#' @param perturbed_label Value in `condition_col` corresponding to perturbed profiles.
#' @param control_label Value in `condition_col` corresponding to control profiles.
#' @param design_vars Optional character vector of additional covariates to place
#'   before condition in the DESeq2 design formula for `"mean"`.
#' @param assay_name Assay name in `se` to use for `"mean"`. Defaults to the first assay.
#' @param alpha Adjusted p-value cutoff used by the `"mean"` method when selecting
#'   upregulated genes.
#' @param min_count Optional count threshold used to filter low-expression genes
#'   before DESeq2 for `"mean"`. If `NULL`, no pre-filtering is applied.
#' @param min_samples Minimum number of samples that must satisfy `min_count`
#'   for a gene to be retained when `min_count` is provided.
#'
#' @return A named list of recontextualized gene signatures.
#' @export
recontextualize <- function(method = c("networkProp", "projectCor", "mean"),
                            se = NULL,
                            seeds = NULL,
                            sigs = NULL,
                            score = c("gsva", "AUCell", "eigen"),
                            sig = c("rwr", "corr"),
                            avg_p = FALSE,
                            avg_p_vals = c(1e-4, 1e-1),
                            avg_p_length = 5,
                            p = 0.1,
                            bootstrap = FALSE,
                            n_bootstraps = 1000,
                            limit = NULL,
                            nfeatures = NULL,
                            min.sft = 0.85,
                            beta = NULL,
                            cores = 1,
                            cor.fn = c("bicor", "cor"),
                            cor.type = c("unsigned", "signed hybrid", "signed"),
                            powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                            diag_zero = TRUE,
                            perturbation_col = NULL,
                            condition_col = NULL,
                            perturbed_label = "Perturbed",
                            control_label = "Control",
                            design_vars = NULL,
                            assay_name = NULL,
                            alpha = 0.05,
                            min_count = NULL,
                            min_samples = 1) {
  method <- match.arg(method)

  if (is.null(se)) {
    stop("'se' must be supplied to recontextualize().")
  }

  if (method == "networkProp") {
    seeds <- .normalize_recontext_input(.recontext_input_or(seeds, sigs), "seeds")
    if (is.null(seeds)) {
      stop("'seeds' or 'sigs' must be supplied for method = 'networkProp'.")
    }

    if (is.null(limit)) {
      limit <- lengths(seeds)
    }

    ig <- .build_networkprop_graph(
      se = se,
      assay_name = assay_name,
      nfeatures = nfeatures,
      min.sft = min.sft,
      beta = beta,
      cores = cores,
      cor.fn = cor.fn,
      cor.type = cor.type,
      powers = powers,
      diag_zero = diag_zero
    )

    return(network_sig(
      ig = ig,
      seeds = seeds,
      sig = match.arg(sig),
      avg_p = avg_p,
      avg_p_vals = avg_p_vals,
      avg_p_length = avg_p_length,
      p = p,
      bootstrap = bootstrap,
      n_bootstraps = n_bootstraps,
      limit = limit
    ))
  }

  if (method == "projectCor") {
    sigs <- .normalize_recontext_input(.recontext_input_or(sigs, seeds), "sigs")
    if (is.null(sigs)) {
      stop("'sigs' or 'seeds' must be supplied for method = 'projectCor'.")
    }

    return(projectCor(
      se = se,
      sigs = sigs,
      score = match.arg(score)
    ))
  }

  sigs <- .normalize_recontext_input(.recontext_input_or(sigs, seeds), "sigs")
  if (is.null(perturbation_col) || is.null(condition_col)) {
    stop("'perturbation_col' and 'condition_col' are required for method = 'mean'.")
  }

  .mean_recontextualize(
    se = se,
    sigs = sigs,
    limit = limit,
    perturbation_col = perturbation_col,
    condition_col = condition_col,
    perturbed_label = perturbed_label,
    control_label = control_label,
    design_vars = design_vars,
    assay_name = assay_name,
    alpha = alpha,
    min_count = min_count,
    min_samples = min_samples
  )
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
#' @importFrom Matrix sparseMatrix
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
  gene_idx <- setNames(seq_along(gene_names), gene_names)
  n_pbs <- length(seeds)

  # Create Seed Matrix
  if (bootstrap) {
    bins_filter <- bin_presence(lapply(seeds, function(x) length(x)) %>% unlist)
    present_bins <- which(bins_filter)
    n_bins <- length(present_bins)
    n_cols <- n_bins * n_bootstraps
    row_indices <- vector("list", n_cols)
    colnames_mat <- names(present_bins) %>% rep(each = n_bootstraps)

    # If bootstrap matrix is (n_genes, n_bootstraps x n_bins)
    for (present_bin in seq_along(present_bins)) {
      bin <- present_bins[present_bin]
      for (i in seq_len(n_bootstraps)) {
        # Sample a random set of genes with size sampled between bin end-points
        last_digit <- i %% 10
        last_digit <- ifelse(last_digit == 0, 10, last_digit)

        sample_size <- (bin - 1) * 10 + last_digit
        genes <- sample.int(length(gene_names), sample_size)

        bin_col <- (present_bin-1) * n_bootstraps + i
        row_indices[[bin_col]] <- genes
      }
    }

    mat <- Matrix::sparseMatrix(i = unlist(row_indices, use.names = FALSE),
                                j = rep(seq_len(n_cols), lengths(row_indices)),
                                x = 1,
                                dims = c(length(gene_names), n_cols),
                                dimnames = list(gene_names, colnames_mat))
  } else {
    row_indices <- unlist(lapply(seeds, function(genes) unname(gene_idx[genes])), use.names = FALSE)
    col_indices <- rep(seq_len(n_pbs), lengths(seeds))

    mat <- Matrix::sparseMatrix(i = row_indices,
                                j = col_indices,
                                x = 1,
                                dims = c(length(gene_names), n_pbs),
                                dimnames = list(gene_names, names(seeds)))
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
  transition_mat <- prepare_rwr_transition(ig = ig, normalize = normalize)
  mat <- rwr_from_seed_matrix(transition_mat = transition_mat,
                              seed_mat = seed_mat,
                              restart = restart,
                              avg_p = avg_p,
                              avg_p_vals = avg_p_vals,
                              avg_p_length = avg_p_length,
                              epsilon = epsilon)

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
    intervals <- unique(colnames(bootstraps))
    bounds <- t(vapply(strsplit(intervals, "-", fixed = TRUE), as.numeric, numeric(2)))
    bootstrap_idx_by_interval <- split(seq_len(ncol(bootstraps)), colnames(bootstraps))
    percentiles <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))

    for (j in seq_len(ncol(mat))) {
      pb_name <- colnames(mat)[j]
      sig_bin <- sig_bins[[pb_name]]
      interval <- intervals[which(sig_bin >= bounds[, 1] & sig_bin <= bounds[, 2])]
      bootstrap_idx <- bootstrap_idx_by_interval[[interval]]
      obs_col <- as.numeric(mat[, j])
      bootstrap_subset <- bootstraps[, bootstrap_idx, drop = FALSE]

      if (inherits(bootstrap_subset, "Matrix")) {
        counts <- numeric(length(obs_col))
        for (k in seq_len(ncol(bootstrap_subset))) {
          counts <- counts + (as.numeric(bootstrap_subset[, k]) <= obs_col)
        }
        percentiles[, j] <- counts / ncol(bootstrap_subset)
      } else {
        percentiles[, j] <- rowMeans(bootstrap_subset <= obs_col)
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
      mat_col <- as.numeric(mat[, col])
      names(mat_col) <- rownames(mat)
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
  sig <- match.arg(sig)

  if (is.character(seeds) && length(seeds) == 1) {
    seeds <- setNames(list(seeds), seeds)
  }

  gene_names <- igraph::V(ig)$name
  all_genes <- unname(unlist(seeds))
  seed_genes_filter <- all(all_genes %in% gene_names)
  if (!seed_genes_filter) {
    message("Filtering Seed Genes to those contained in the graph.")
    seeds <- lapply(seeds, function(x) x[x %in% gene_names])
  }

  if (is.list(seeds)) {
    empty_sigs <- names(seeds)[lengths(seeds) == 0]
    if (length(empty_sigs) > 0) {
      for (sig_name in empty_sigs) {
        message(sprintf(
          "Skipping signature '%s' because its genes were not found in the igraph network.",
          sig_name
        ))
      }
      seeds <- seeds[lengths(seeds) > 0]
    }
  }

  if (length(seeds) == 0) {
    message("No signatures remained after filtering to genes present in the igraph network.")
    return(list())
  }

  if (sig == "corr") {
    if ("weight" %in% list.edge.attributes(ig)) {
      cor_mat <- igraph::as_adjacency_matrix(ig, attr = "weight")
    } else {
      cor_mat <- igraph::as_adjacency_matrix(ig, attr = NULL)
    }
    net_sig <- correlated_sigs(corr_mat = cor_mat, seeds = seeds, limit = limit)
  } else if (sig == "rwr") {
    transition_mat <- prepare_rwr_transition(ig = ig, normalize = "row")
    obs_seed_mat <- seed_matrix(ig, seeds)
    obs_mat <- rwr_from_seed_matrix(transition_mat = transition_mat,
                                    seed_mat = obs_seed_mat,
                                    restart = p,
                                    avg_p = avg_p,
                                    avg_p_vals = avg_p_vals,
                                    avg_p_length = avg_p_length)
    if(bootstrap & (p != 1)) {
      # If bootstrap matrix is (n_genes, n_bootstraps x n_bins)
      bootstrap_seed_mat <- seed_matrix(ig,
                                        seeds,
                                        bootstrap = bootstrap,
                                        n_bootstraps = n_bootstraps)
      mat_bootstraps <- rwr_from_seed_matrix(transition_mat = transition_mat,
                                             seed_mat = bootstrap_seed_mat,
                                             restart = p,
                                             avg_p = avg_p,
                                             avg_p_vals = avg_p_vals,
                                             avg_p_length = avg_p_length)

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
  stopifnot(is(seed_mat, "matrix") || is(seed_mat, "Matrix"))
  normalize <- match.arg(normalize)

  transition_mat <- prepare_rwr_transition(ig = ig, normalize = normalize)
  return(rwr_from_seed_matrix(transition_mat = transition_mat,
                              seed_mat = seed_mat,
                              restart = restart,
                              epsilon = epsilon))
}

prepare_rwr_transition <- function(ig, normalize = c("row", "column", "laplacian", "none")) {
  normalize <- match.arg(normalize)

  # Get Adjacency matrix
  if ("weight" %in% igraph::edge_attr_names(ig)) {
    adj_mat <- igraph::as_adjacency_matrix(ig, attr = "weight")
    adj_mat[is.na(adj_mat)] <- 0
    message("Using weighted graph")
  } else {
    adj_mat <- igraph::as_adjacency_matrix(ig, attr = NULL)
    adj_mat[is.na(adj_mat)] <- 0
    message("Using unweighted graph")
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
    nadjM <- adj_mat
  }

  nadjM
}

normalize_seed_matrix <- function(seed_mat) {
  norm_seed_mat <- seed_mat %*% Matrix::Diagonal(x = (Matrix::colSums(seed_mat))^(-1))
  colnames(norm_seed_mat) <- colnames(seed_mat)
  as(norm_seed_mat, "CsparseMatrix")
}

rwr_from_seed_matrix <- function(transition_mat,
                                 seed_mat,
                                 restart = 0.1,
                                 avg_p = FALSE,
                                 avg_p_vals = c(1e-4, 1e-1),
                                 avg_p_length = 5,
                                 epsilon = NULL) {
  norm_seed_mat <- normalize_seed_matrix(seed_mat)

  if (avg_p) {
    avg_p_seq <- seq(avg_p_vals[1], avg_p_vals[2], length.out = avg_p_length)
    mat <- Matrix::Matrix(0, nrow = nrow(norm_seed_mat), ncol = ncol(norm_seed_mat),
                          dimnames = dimnames(norm_seed_mat), sparse = TRUE)

    for (restart_val in avg_p_seq) {
      mat <- mat + random_walk_from_transition(transition_mat = transition_mat,
                                               norm_seed_mat = norm_seed_mat,
                                               restart = restart_val,
                                               epsilon = epsilon)
    }
    return(mat / avg_p_length)
  }

  random_walk_from_transition(transition_mat = transition_mat,
                              norm_seed_mat = norm_seed_mat,
                              restart = restart,
                              epsilon = epsilon)
}

random_walk_from_transition <- function(transition_mat,
                                        norm_seed_mat,
                                        restart = 0.1,
                                        epsilon = NULL) {
  stopifnot(is(norm_seed_mat, "Matrix"))

  ## Stopping Criteria
  stop_delta <- 1e-5 # L1 norm of successive iterations of Transition Matrix multiplication
  stop_step <- 100 # maximum steps of iterations

  ## Initial Variables
  P0 <- norm_seed_mat
  PT <- P0
  r <- restart
  step <- 0
  delta <- 1

  ## Exploration Parameters
  if (!is.null(epsilon)) {
    stopifnot(is(epsilon,"numeric"))

    n_genes_seed_mean <- as.integer(mean(Matrix::colSums(norm_seed_mat > 0)))
    n_explore <- as.integer(n_genes_seed_mean * epsilon)
    r <- 1
    paste0("Stopping Criterion: ", n_explore, " displaced genes.")
  }

  ## Dnet had the matrix multiplication orders flipped.
  ## This order keeps the distribution in columns after each multiplication.
  while (delta > stop_delta && step <= stop_step) {
    PX <- (1 - r) * Matrix::t(Matrix::t(PT) %*% transition_mat) + r * P0
    delta <- sum(abs(PX - PT))
    PT <- PX
    step <- step + 1

    ## Write function that counts number of displaced.
    ## Add stopping condition
    # if (!is.null(epsilon)) {
    # }

    if (step > stop_step) {
      message(paste0("Reached maximum iteration steps. Delta: ", delta))
    } else if (delta <= stop_delta) {
      message(paste0("Reached Convergence. Iteration step: ", step))
    }
  }

  return(PX)
}
