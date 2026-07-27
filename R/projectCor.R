#' @title Projection based Recontextualization
#'
#' @description
#' This function reconstructs gene signatures based on their correlation with
#' per-sample projection scores computed from the input signatures.
#'
#' @param se A SummarizedExperiment object containing gene expression data.
#' @param sigs A list of gene signatures, where each element is a character vector of gene names.
#' @param score Scoring method used to score samples against input signatures.
#'   Either `"gsva"` or `"eigen"`. Default is `"gsva"`.
#'
#' @return A list of reconstructed gene signatures, with the same structure as the input `sigs`.
#'
#' @details
#' The function performs the following steps:
#' 1. Calculates projection scores for the input signatures with GSVA or eigengenes.
#' 2. Computes the correlation between gene expression and projection scores.
#' 3. Ranks genes based on their correlation with each signature's projection scores.
#' 4. Selects the top-ranking genes to form new signatures of the same length as the original ones.
#'
#' @importFrom stats cor
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange slice pull desc
#'
#' @export
projectCor <- function(se, sigs, score = c("gsva", "eigen")) {
  score <- match.arg(score)
  stopifnot(is(se, "SummarizedExperiment"))

  expr_mat <- SummarizedExperiment::assay(se)

  if (score == "gsva") {
    if (!requireNamespace("GSVA", quietly = TRUE)) {
      stop(
        "The 'GSVA' package is required for score = 'gsva'. Install it with BiocManager::install('GSVA'), or use score = 'eigen' instead."
      )
    }
    score_param <- GSVA::gsvaParam(se, sigs, maxDiff = TRUE)
    score_res <- GSVA::gsva(score_param, verbose = FALSE)
    if (is(score_res, "SummarizedExperiment")) {
      score_mat <- SummarizedExperiment::assay(score_res)
    } else {
      score_mat <- score_res
    }
  } else {
    if (is(expr_mat, "sparseMatrix")) {
      expr_mat <- as.matrix(expr_mat)
    }

    gene_names <- rownames(expr_mat)
    if (is.null(gene_names)) {
      stop("Expression assay must have rownames to compute eigengene scores.")
    }

    sig_names <- names(sigs)
    eig_scores <- setNames(
      lapply(sig_names, function(sig_name) {
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
          stop(sprintf(
            "All genes in signature '%s' have zero variance.",
            sig_name
          ))
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
      }),
      sig_names
    )

    keep_scores <- !vapply(eig_scores, is.null, logical(1))
    eig_scores <- eig_scores[keep_scores]

    if (length(eig_scores) == 0) {
      stop(
        "No signatures had genes present in the expression assay for eigengene scoring."
      )
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
  for (sig_name in colnames(proj_scores)) {
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

gsva_recon <- function(se, sigs) {
  projectCor(se = se, sigs = sigs, score = "gsva")
}
