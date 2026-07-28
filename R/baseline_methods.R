.normalize_recontext_input <- function(x, arg_name) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.character(x) && length(x) == 1) {
    x <- stats::setNames(list(x), x)
  }

  if (!is.list(x)) {
    stop(sprintf(
      "'%s' must be a named list or a single character string.",
      arg_name
    ))
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

.resolve_recontext_limits <- function(items, limit = NULL, default = 30) {
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
    stop(
      "'limit' must have length 1, match the number of signatures, or be a named vector."
    )
  }

  names(limit) <- item_names
  limit
}

.mean_recontextualize <- function(
  se,
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
  min_samples = 1
) {
  stopifnot(is(se, "SummarizedExperiment"))

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("The 'DESeq2' package is required for method = 'mean'.")
  }

  if (is.null(assay_name)) {
    assay_names <- names(SummarizedExperiment::assays(se))
    assay_name <- assay_names[[1]]
  }

  if (
    is.null(assay_name) ||
      !assay_name %in% names(SummarizedExperiment::assays(se))
  ) {
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

  controls <- rownames(coldata)[
    as.character(coldata[[condition_col]]) == control_label
  ]
  if (length(controls) == 0) {
    stop(sprintf(
      "No control samples found where '%s' == '%s'.",
      condition_col,
      control_label
    ))
  }

  sigs <- .normalize_recontext_input(sigs, "sigs")
  if (is.null(sigs)) {
    pert_names <- unique(as.character(
      coldata[[perturbation_col]][
        as.character(coldata[[condition_col]]) == perturbed_label
      ]
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
      ifelse(
        rownames(subset_meta) %in% perturbed_samples,
        perturbed_label,
        control_label
      ),
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

#' @title Recontextualize signatures with one of either networkProp, projectCor, or mean.
#'
#' @param method Baseline recontextualization method. One of `"networkProp"`,
#'   `"projectCor"`, or `"mean"`.
#' @param se A SummarizedExperiment object used by `"projectCor"`, `"mean"`,
#'   and to learn the `"networkProp"` graph.
#' @param seeds Seed signatures for the `"networkProp"` method.
#' @param sigs Signature list for the `"projectCor"` method. For `"mean"`, names
#'   are interpreted as perturbation labels and lengths are used as output sizes
#'   when `limit` is not supplied.
#' @param ig An optional pre-built igraph for method = `"networkProp"`, e.g.
#'   from a previous [wgcna.adj()] call. If supplied, network construction is
#'   skipped and `se` is not required. See [netProp()].
#' @param score Scoring method used by `"projectCor"`.
#' @param sig Network signature mode passed to [netProp()]. Defaults to `"rwr"`.
#' @param avg_p Whether to ensemble network propagation over multiple restart probabilities.
#' @param avg_p_vals Range of restart probabilities used when `avg_p = TRUE`.
#' @param avg_p_length Number of restart probabilities to average when `avg_p = TRUE`.
#' @param p Restart probability passed to [netProp()].
#' @param bootstrap Whether to use bootstrap-based extraction in [netProp()].
#' @param n_bootstraps Number of bootstrap replicates passed to [netProp()].
#' @param limit Number of genes to keep for each output signature. If `NULL`,
#'   `"networkProp"` and `"projectCor"` keep the original signature lengths and
#'   `"mean"` uses the lengths of `sigs`.
#' @param nfeatures Number of variable genes to use when learning the
#'   `"networkProp"` graph from `se`. Defaults to `min(10000, nrow(assay(se)))`.
#' @param min.sft Minimum scale-free topology fitting index used by
#'   [netProp()] when learning the `"networkProp"` graph.
#' @param beta Optional soft-thresholding power passed to [netProp()] for
#'   `"networkProp"`.
#' @param cores Number of CPU cores passed to [netProp()] for `"networkProp"`.
#' @param cor.fn Correlation function passed to [netProp()] for `"networkProp"`.
#' @param cor.type Correlation network type passed to [netProp()] for `"networkProp"`.
#' @param powers Candidate power values passed to [netProp()] for `"networkProp"`
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
recontextualize <- function(
  method = c("networkProp", "projectCor", "mean"),
  se = NULL,
  seeds = NULL,
  sigs = NULL,
  ig = NULL,
  score = c("gsva", "eigen"),
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
  min_samples = 1
) {
  method <- match.arg(method)

  if (is.null(se) && is.null(ig)) {
    stop("'se' (or 'ig' for method = 'networkProp') must be supplied to recontextualize().")
  }

  if (method == "networkProp") {
    seeds <- .normalize_recontext_input(
      .recontext_input_or(seeds, sigs),
      "seeds"
    )
    if (is.null(seeds)) {
      stop("'seeds' or 'sigs' must be supplied for method = 'networkProp'.")
    }

    if (is.null(limit)) {
      limit <- lengths(seeds)
    }

    return(netProp(
      se = se,
      seeds = seeds,
      ig = ig,
      sig = match.arg(sig),
      avg_p = avg_p,
      avg_p_vals = avg_p_vals,
      avg_p_length = avg_p_length,
      p = p,
      bootstrap = bootstrap,
      n_bootstraps = n_bootstraps,
      limit = limit,
      assay_name = assay_name,
      nfeatures = nfeatures,
      min.sft = min.sft,
      beta = beta,
      cores = cores,
      cor.fn = cor.fn,
      cor.type = cor.type,
      powers = powers,
      diag_zero = diag_zero
    ))
  }

  if (is.null(se)) {
    stop("'se' must be supplied to recontextualize().")
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
    stop(
      "'perturbation_col' and 'condition_col' are required for method = 'mean'."
    )
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
