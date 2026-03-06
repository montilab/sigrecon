#' A push/pop capable vector
#'
#' @description
#' An R6 class that implements a persistent vector with push and pop operations.
#'
#' @importFrom R6 R6Class
#'
#' @keywords internal
pvector <- R6Class("pvector", list(
  #' @field values A vector of values
  values = NULL,

  #' @description
  #' Create a pvector
  #' @param values A vector of values
  #' @return A new pvector
  initialize = function(values=c()) {
    self$values <- values
  },
  #' @description
  #' Print pvector
  #' @return NULL
  print = function() {
    base::print(self$values)
    invisible(self)
  },
  #' @description
  #' Get length of pvector
  #' @return An integer
  length = function() {
    base::length(self$values)
  },
  #' @description
  #' Pop vector
  #' @return Popped value
  pop = function() {
    if (length(self$values) > 0) {
      popped.value <- self$values[1]
      self$values <- self$values[-1]
      return(popped.value)
    }
  },
  #' @description
  #' Push values
  #' @param pushed.values A vector of values
  #' @return NULL
  push = function(pushed.values) {
    self$values <- c(self$values, pushed.values)
  }
))

#' Rank Genes in an ExpressionSet/SummarizedExperiment by Variability
#'
#' This function ranks genes in an ExpressionSet/SummarizedExperiment object based on their variability,
#' as measured by a specified function (default is median absolute deviation).
#'
#' @param eset An ExpressionSet/SummarizedExperiment object containing gene expression data.
#' @param fn A function to measure variability. Default is `mad` (median absolute deviation).
#' @param filter_zero A logical indicating whether to filter out genes with zero variance. Default is FALSE.
#'
#' @return A pvector object containing ranked gene names.
#'
#' @details
#' The function calculates the variability of each gene using the specified function (default: mad).
#' Genes are then sorted in descending order of variability. If `filter_zero` is TRUE,
#' genes with zero variance are removed before ranking.
#'
#' @importFrom Biobase exprs
#' @importFrom SummarizedExperiment assay
#' @importFrom stats mad
rank.var.eset <- function(eset, fn=mad, filter_zero=FALSE) {
  stopifnot(
    "Input object must be an ExpressionSet or a SummarizedExperiment" =
      is(eset, "ExpressionSet") || is(eset, "SummarizedExperiment")
  )

  # Extract expression matrix based on object type
  if (is(eset, "ExpressionSet")) {
    eset.mat <- Biobase::exprs(eset)
  } else if (is(eset, "SummarizedExperiment")) {
    # Default to the first assay if no specific assay name is provided or if there's only one
    # For more robust usage, one might add an 'assay_name' argument.
    eset.mat <- SummarizedExperiment::assay(eset)
  }
  gene.var <- apply(eset.mat, 1, fn)

  if(filter_zero) {
    genes.keep <- gene.var != 0
    non_zero_genes <- rownames(eset.mat)[genes.keep]
    eset.mat <- eset.mat[non_zero_genes,]
  }

  ranked.genes <- pvector$new(names(sort(gene.var, decreasing=TRUE)))
  return(ranked.genes)
}


#' Find Common Genes with Highest Median Absolute Deviation (MAD) Across ExpressionSets
#'
#' This function identifies common genes with the highest MAD across multiple ExpressionSet objects,
#' up to a specified limit.
#'
#' @param esets A list of ExpressionSet objects to compare.
#' @param limit An integer specifying the maximum number of common MAD genes to return. Default is 2500.
#' @param parallel A logical indicating whether to use parallel processing. Default is FALSE.
#' @param filter_zero A logical indicating whether to filter out genes with zero variance. Default is FALSE.
#'
#' @return A character vector of common gene names with highest MAD.
#'
#' @details
#' The function ranks genes in each ExpressionSet by their MAD, then iteratively selects
#' genes that are present in all ExpressionSets until reaching the specified limit or
#' exhausting all common genes.
#'
#' @importFrom parallel detectCores makeCluster parLapply
#' @importFrom doParallel registerDoParallel
#' @importFrom Biobase featureNames
common_mad_genes <- function(esets, limit=2500, parallel=FALSE, filter_zero=FALSE) {
  if(parallel) {
    no_cores <- detectCores() - 1
    registerDoParallel(cores=no_cores)
    cl <- makeCluster(no_cores, type="FORK")
    ranked.genes <- parallel::parLapply(esets, rank.var.eset)
  } else {
    ranked.genes <- lapply(esets, rank.var.eset, filter_zero=filter_zero)
    print("Finished ranking esets")
  }

  genes.selected <- c()
  i <- 1
  while (TRUE) {
    popped <- ranked.genes[[i]]$pop()
    in_all_esets <- all(sapply(esets, function (x) popped %in% featureNames(x)))
    if (popped %in% genes.selected | !in_all_esets) next

    genes.selected <- c(genes.selected, popped)
    i <- i+1

    if (i == length(ranked.genes)+1) i <- 1
    if (length(genes.selected) >= limit) break
  }
  return(genes.selected)
}

#' Find Common Variable Genes Across Seurat Objects
#'
#' This function identifies common variable genes across multiple Seurat objects,
#' up to a specified limit.
#'
#' @param seurat_objs A list of Seurat objects to compare.
#' @param limit An integer specifying the maximum number of common variable genes to return.
#'
#' @return A character vector of common variable gene names.
#'
#' @details
#' The function iterates through the variable features of each Seurat object,
#' selecting genes that are present in all objects. It continues until it reaches
#' the specified limit or exhausts all common variable genes.
#'
#' @importFrom Seurat VariableFeatures
seurat_common_var_genes <- function(seurat_objs, limit) {

  if(seurat_objs[[1]]@version == "5.0.1") {
    pvectors <- lapply(seurat_objs, function(x) pvector$new(Seurat::VariableFeatures(x)))
  } else {
    pvectors <- lapply(seurat_objs, function(x) pvector$new(x@assays$RNA@var.features))
  }
  selected <- c()
  i <- 1
  while (TRUE) {
    popped <- pvectors[[i]]$pop()
    in_all_esets <- all(sapply(seurat_objs, function (x) popped %in% rownames(x)))
    if (popped %in% selected | !in_all_esets) next

    selected <- c(selected, popped)
    i <- i+1

    if (i == length(pvectors)+1) i <- 1
    if (length(selected) >= limit) break
  }
  return(selected)
}

#' Check Presence in Bins for Numbers 1-100
#'
#' Given a numeric vector with values between 1 and 100, returns a logical vector
#' indicating whether at least one value falls into each bin of size 10
#' (i.e., 1-10, 11-20, ..., 91-100).
#'
#' @param x Numeric vector. Values should be between 1 and 100.
#' @return A named logical vector of length 10. Each element is \code{TRUE} if at least one
#'   value in \code{x} falls into the corresponding bin, otherwise \code{FALSE}.
#' @examples
#' vec <- c(3, 15, 27, 45, 58, 99)
#' bin_presence(vec)
#' @export
bin_presence <- function(x) {
  breaks <- seq(1, 100, by = 10)
  labels <- paste(breaks, breaks + 9, sep = "-")
  labels[length(labels)] <- "91-100"
  bins <- cut(x, breaks = c(breaks, 101), right = FALSE, labels = labels)
  present <- labels %in% bins
  names(present) <- labels
  return(present)
}

#' Filter Significant Genes by Perturbation
#'
#' @description
#' Filters and ranks genes by significance for each perturbation in a differential
#' expression table. Returns both all genes ranked by a combined score and
#' significantly up-regulated genes that pass thresholds.
#'
#' @param diff_table A data frame containing differential expression results with
#'   perturbation identifiers, log2 fold changes, and adjusted p-values.
#' @param perts A character vector of perturbation names to filter for.
#' @param alpha Numeric. Adjusted p-value threshold for significance. Default is 0.05.
#' @param limit Integer. Maximum number of top significant genes to return per
#'   perturbation. Default is 100.
#' @param pert_col Character. Name of the column containing perturbation identifiers.
#'   Default is "product".
#' @param log2fc_col Character. Name of the column containing log2 fold change values.
#'   Default is "avg_log2FC".
#' @param pval_col Character. Name of the column containing adjusted p-values.
#'   Default is "p_val_adj".
#' @param geneid_col Character. Name of the column containing gene identifiers.
#'
#' @return A nested list where each perturbation contains:
#'   \itemize{
#'     \item{\code{up}}: Character vector of row names for top significantly
#'       upregulated genes (filtered by alpha and limited by limit parameter)
#'     \item{\code{up_full}}: Character vector of row names for all genes ranked
#'       by combined score (log2FC * -log10(p_val_adj))
#'   }
#'
#' @details
#' The function computes a combined significance score as:
#' \deqn{score = log2FC \times -log10(p_{adj})}
#'
#' For the "up" results, genes must meet three criteria:
#' \enumerate{
#'   \item log2FC > 0 (upregulated)
#'   \item adjusted p-value <= alpha
#'   \item Ranked in top N genes (specified by limit)
#' }
#' @importFrom dplyr filter arrange slice
#' @export
sig_filter_fn <- function(diff_table,
                          perts,
                          alpha = 0.05,
                          limit = 100,
                          pert_col = "product",
                          log2fc_col = "avg_log2FC",
                          pval_col = "p_val_adj",
                          geneid_col = "ensembl_id") {
  results <- list()

  for(pb in perts) {
    print(pb)

    # Create combined score column
    diff_table$logFC_adjpval <- diff_table[[log2fc_col]] * (-log10(diff_table[[pval_col]]))

    # Get all genes for this perturbation, ordered by score
    full_sig <- diff_table %>%
      dplyr::filter(!!sym(pert_col) == !!pb) %>%
      dplyr::arrange(desc(logFC_adjpval)) %>%
      dplyr::pull(geneid_col)

    # Get significant upregulated genes
    sig_sig <- diff_table %>%
      dplyr::filter(!!sym(pert_col) == !!pb) %>%
      dplyr::filter(!!sym(log2fc_col) > 0) %>%
      dplyr::filter(!!sym(pval_col) <= alpha) %>%
      dplyr::arrange(desc(logFC_adjpval)) %>%
      dplyr::slice(1:limit) %>%
      dplyr::pull(geneid_col)

    results[[pb]][["up"]] <- sig_sig
    results[[pb]][["up_full"]] <- full_sig
  }

  return(results)
}

#' @title Run MDMR regression with robust error handling
#' @description A helper function to perform MDMR regression given a gene signature,
#'   a SummarizedExperiment object, and a phenotype label from colData.
#' @param signature A character vector of gene IDs.
#' @param se A SummarizedExperiment object containing expression data and phenotype information.
#' @param phenotype_label A character string specifying the column name in colData(se)
#'   to be used as the phenotype for MDMR (e.g., "AGE").
#' @param assay_label A character string specifying the name of the assay in the se object to be used for MDMR.
#' @param method_label An optional character string for logging/warning purposes,
#'   indicating which method (e.g., "non_recon") is calling this function.
#' @return A list with elements 'stat' (the MDMR statistic) and 'r2' (the R-squared/pr.sq),
#'   or a list with NA_real_ values if an error or invalid input occurs.
#'
#' @importFrom MDMR mdmr
#' @importFrom SummarizedExperiment colData
#'
#' @export
mdmr_eval <- function(signature,
                      se,
                      phenotype_label,
                      assay_label = "DESeq2_log",
                      method_label = "MDMR_run") {

  # Initialize return values
  res <- list(stat = NA_real_, r2 = NA_real_)

  # --- 1. Validate Phenotype Data ---
  if (!phenotype_label %in% colnames(colData(se))) {
    warning(paste0("Skipping ", method_label, ": Phenotype '", phenotype_label, "' not found in colData(se)."))
    return(res)
  }

  pheno_data <- colData(se)[[phenotype_label]]
  if (is.null(pheno_data) || length(unique(pheno_data)) < 2) {
    warning(paste0("Skipping ", method_label, ": Phenotype '", phenotype_label, "' data is insufficient or constant."))
    return(res)
  }

  # Calculate distance matrix for phenotype
  y_d <- dist(as.matrix(pheno_data), method = "manhattan")

  # --- 2. Validate Gene Signature and Expression Data ---
  if (is.null(signature) || length(signature) < 2) {
    warning(paste0("Skipping ", method_label, ": Gene signature insufficient (<2 genes)."))
    return(res)
  }

  all_features_in_se <- rownames(se)
  valid_genes <- intersect(signature, all_features_in_se)

  if (length(valid_genes) < 2) {
    warning(paste0("Skipping ", method_label, ": Too few valid genes (", length(valid_genes), ") from signature found in se object."))
    return(res)
  }

  expression_mat <- se[valid_genes,]@assays@data[[assay_label]]

  if (nrow(expression_mat) < 2 || ncol(expression_mat) < 2) {
    warning(paste0("Skipping ", method_label, ": Expression matrix for signature is insufficient (",
                   nrow(expression_mat), " genes, ", ncol(expression_mat), " samples)."))
    return(res)
  }

  # --- 3. Run MDMR with error handling ---
  mdmr_output <- tryCatch({
    # Transpose for MDMR: samples as rows, genes as columns
    mdmr_res <- mdmr(t(expression_mat), y_d)
    list(stat = mdmr_res$stat[1,], r2 = mdmr_res$pr.sq[1,])
  }, error = function(e) {
    warning(paste0("MDMR failed for ", method_label, ": ", e$message))
    res # Return initialized NAs on error
  })

  return(mdmr_output)
}

#' Create appropriate BiocParallel backend
#'
#' @param workers Number of parallel workers (default: 1 for sequential)
#' @param RNGseed Random seed for reproducibility (default: NULL)
#' @param progress Show progress bar (default: FALSE)
#' @param type Force backend type: "multicore", "snow", or "serial" (default: auto-detect)
#' @return A BiocParallelParam object
#'
#' @examples
#' # Sequential (default)
#' bp <- make_bpparam()
#'
#' # Parallel on Unix/Mac (forking)
#' bp <- make_bpparam(workers = 2, RNGseed = 123)
#'
#' # Parallel on Windows (socket)
#' bp <- make_bpparam(workers = 2, RNGseed = 123)
#'
#' # With progress bar
#' bp <- make_bpparam(workers = 2, progress = TRUE)
#'
#' @importFrom BiocParallel MulticoreParam SnowParam SerialParam
#' @export
make_bpparam <- function(workers = 1,
                         RNGseed = NULL,
                         progress = FALSE,
                         type = NULL) {

  # Sequential execution
  if (workers == 1) {
    return(BiocParallel::SerialParam())
  }

  # Determine backend type
  if (is.null(type)) {
    type <- if (.Platform$OS.type == "unix") "multicore" else "snow"
  }

  # Create appropriate backend
  if (type == "multicore") {
    if (.Platform$OS.type != "unix") {
      warning("MulticoreParam only works on Unix/Mac, falling back to SnowParam")
      type <- "snow"
    }
  }

  BPPARAM <- switch(type,
                    "multicore" = BiocParallel::MulticoreParam(
                      workers = workers,
                      RNGseed = RNGseed,
                      progressbar = progress
                    ),
                    "snow" = BiocParallel::SnowParam(
                      workers = workers,
                      RNGseed = RNGseed,
                      progressbar = progress
                    ),
                    "serial" = BiocParallel::SerialParam(),
                    stop("Unknown backend type: ", type)
  )

  return(BPPARAM)
}

#' fgsea wrapper for gene symbol vectors
#'
#' @param ref Vector of gene symbols representing the full ranked list (with weights)
#' @param data Vector of gene symbols representing the gene set/pathway to test
#' @param scoreType GSEA score type: "std" (default), "pos", or "neg"
#' @param eps Precision for p-value calculation (default: 1e-50)
#' @return A list with fgsea results (ES, NES, pval, padj, leadingEdge, size)
#'
#' @examples
#' # ref = all genes ranked by importance
#' ref <- c("TP53", "MYC", "BRCA1", "KRAS", "EGFR", "PTEN", "AKT1", "PIK3CA")
#' # data = gene set to test for enrichment
#' data <- c("TP53", "BRCA1", "EGFR")
#' result <- fgsea_wrapper(ref, data)
#'
#' @importFrom fgsea fgseaMultilevel
#' @export
fgsea_wrapper <- function(ref,
                          data,
                          scoreType = "std",
                          eps = 1e-50) {

  # Input validation
  if (length(ref) == 0) {
    stop("'ref' must be a non-empty vector of ranked gene symbols")
  }
  if (length(data) == 0) {
    warning("'data' is empty, returning NULL")
    return(NULL)
  }

  # Check for duplicates in ref
  if (any(duplicated(ref))) {
    warning("'ref' contains duplicates, removing them")
    ref <- unique(ref)
  }

  n <- length(ref)

  # Create stats vector: descending rank from 1 to -1
  # Position 1 (most important gene in ref) gets rank 1
  # Last position gets rank -1
  stats <- seq(from = 1, to = -1, length.out = n)
  names(stats) <- ref

  # data is the gene set/pathway to test
  # Only test genes that appear in ref
  pathway_genes <- intersect(data, ref)

  if (length(pathway_genes) == 0) {
    warning("No overlap between 'data' and 'ref', returning NULL")
    return(NULL)
  }

  # Prepare pathways list
  pathways <- list(geneset = pathway_genes)

  # Run fgsea
  result <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    scoreType = scoreType,
    eps = eps
  )

  # Extract and format results
  return(list(
    ES = result$ES[1],
    NES = result$NES[1],
    pval = result$pval[1],
    padj = result$padj[1],
    log2err = result$log2err[1],
    size = result$size[1],
    leadingEdge = result$leadingEdge[[1]]
  ))
}

#' Vectorized fgsea wrapper for multiple gene set comparisons
#'
#' @param ref_vecs Named list of ranked gene symbol vectors (each is a full ranked list)
#' @param data_vecs Named list of gene sets to test (each is a pathway/gene set)
#' @param scoreType GSEA score type: "std" (default), "pos", or "neg"
#' @param eps Precision for p-value calculation (default: 1e-50)
#' @param BPPARAM A BiocParallelParam object specifying parallel backend.
#'   If NULL, automatically selects appropriate backend based on platform.
#'   See \code{BiocParallel::BiocParallelParam} for options.
#' @return A data frame with fgsea results for each ref-data pair
#'
#' @examples
#' library(BiocParallel)
#'
#' # Multiple ranked lists (e.g., from different conditions)
#' ref_vecs <- list(
#'   condition1 = c("TP53", "MYC", "BRCA1", "KRAS", "EGFR", "PTEN", "AKT1"),
#'   condition2 = c("MYC", "KRAS", "TP53", "EGFR", "BRCA1", "PTEN", "AKT1")
#' )
#'
#' # Gene sets to test in each condition
#' data_vecs <- list(
#'   condition1 = c("TP53", "BRCA1", "EGFR"),
#'   condition2 = c("TP53", "BRCA1", "EGFR")
#' )
#'
#' # Sequential execution
#' results <- v.fgsea(ref_vecs, data_vecs)
#'
#' # Parallel execution (Unix/Mac)
#' results <- v.fgsea(
#'   ref_vecs, data_vecs,
#'   BPPARAM = MulticoreParam(workers = 2)
#' )
#'
#' # Parallel execution (Windows)
#' results <- v.fgsea(
#'   ref_vecs, data_vecs,
#'   BPPARAM = SnowParam(workers = 2)
#' )
#'
#' # With reproducibility
#' results <- v.fgsea(
#'   ref_vecs, data_vecs,
#'   BPPARAM = MulticoreParam(workers = 2, RNGseed = 123)
#' )
#'
#' @importFrom fgsea fgseaMultilevel
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam SnowParam
#' @export
v.fgsea <- function(ref_vecs,
                    data_vecs,
                    scoreType = "std",
                    eps = 1e-50,
                    BPPARAM = NULL) {

  # Input validation
  if (!is.list(ref_vecs) || !is.list(data_vecs)) {
    stop("'ref_vecs' and 'data_vecs' must be named lists")
  }

  if (is.null(names(ref_vecs)) || is.null(names(data_vecs))) {
    stop("'ref_vecs' and 'data_vecs' must have names")
  }

  # Find matching names
  common_names <- intersect(names(ref_vecs), names(data_vecs))

  if (length(common_names) == 0) {
    stop("No matching names between 'ref_vecs' and 'data_vecs'")
  }

  # Filter to only matching pairs
  ref_vecs <- ref_vecs[common_names]
  data_vecs <- data_vecs[common_names]

  # Set up parallel backend if not provided
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::SerialParam()
  }

  # Define function to run on each pair
  run_pair <- function(name) {
    ref <- ref_vecs[[name]]
    data <- data_vecs[[name]]

    tryCatch({
      result <- fgsea_wrapper(
        ref = ref,
        data = data,
        scoreType = scoreType,
        eps = eps
      )

      if (is.null(result)) {
        return(data.frame(
          name = name,
          ES = NA,
          NES = NA,
          pval = NA,
          padj = NA,
          log2err = NA,
          size = 0,
          leadingEdge = I(list(character(0))),
          stringsAsFactors = FALSE
        ))
      }

      return(data.frame(
        name = name,
        ES = result$ES,
        NES = result$NES,
        pval = result$pval,
        padj = result$padj,
        log2err = result$log2err,
        size = result$size,
        leadingEdge = I(list(result$leadingEdge)),
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      warning(sprintf("Error processing '%s': %s", name, e$message))
      return(data.frame(
        name = name,
        ES = NA,
        NES = NA,
        pval = NA,
        padj = NA,
        log2err = NA,
        size = 0,
        leadingEdge = I(list(character(0))),
        stringsAsFactors = FALSE
      ))
    })
  }

  # Run analysis with BiocParallel
  results_list <- BiocParallel::bplapply(
    X = common_names,
    FUN = run_pair,
    BPPARAM = BPPARAM
  )

  # Combine results
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL

  # Recompute FDR across all tests
  valid_pvals <- !is.na(results_df$pval)
  if (any(valid_pvals)) {
    results_df$padj[valid_pvals] <- p.adjust(
      results_df$pval[valid_pvals],
      method = "BH"
    )
  }

  return(results_df)
}

#' Calculate Jaccard Similarity for Multiple Pairs of Sets
#'
#' @param sigs1 A named list of genesets (vectors)
#' @param sigs2 A named list of genesets (vectors)
#'
#' @return A vector of Jaccard similarity scores
#'
#' @export
v.jaccard <- function(sigs1, sigs2) {

  if(!setequal(names(sigs1), names(sigs2))) {
    print("Filtering sets to matching pairs.")
    shared_vecs <- intersect(names(sigs1), names(sigs2))
    removed_vecs <- setdiff(names(sigs1), names(sigs2))
    print(paste0("Removed ", length(removed_vecs), " perturbations."))
    sigs1 <- sigs1[shared_vecs]
    sigs2 <- sigs2[shared_vecs]
  }

  stopifnot(all.equal(names(sigs1), names(sigs2)))

  sims <- c()
  for(i in seq_along(sigs1)) {
    sim <- jaccard(sigs1[[i]], sigs2[[i]])
    sims <- c(sims, sim)
  }
  return(sims)
}

#' Calculate Jaccard Similarity Between Two Sets
#'
#' @param a A vector representing a set
#' @param b A vector representing a set
#'
#' @return The Jaccard similarity score between sets a and b
#'
#' @export
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

#' Create a Jaccard Similarity Matrix for Multiple Sets
#'
#' @param sets A list of sets (vectors)
#'
#' @return A matrix of Jaccard similarity scores between all pairs of sets
#'
#' @export
jaccard_matrix <- function(sets) {
  n_sets <- length(sets)
  jaccard_mat <- matrix(data = NA,nrow=n_sets,ncol=n_sets)
  rownames(jaccard_mat) <- names(sets)
  colnames(jaccard_mat) <- names(sets)

  for (i in seq_along(sets)) {
    set_i <- sets[[i]]
    for (j in seq_along(sets)) {
      set_j <- sets[[j]]
      jaccard_mat[i, j] <- jaccard(set_i, set_j)
    }
  }

  return(jaccard_mat)
}

#' Count number of displaced seeds
#' Note: This not commutative. Num_displaced is the number of genes that are from the original seed that are not in the recon.
#' @param recon_sig Named list of genesets
#' @param seed_sig Named list of genesets
#' @return Named list for which each entry is a list with two elements: displaced seed genes, and non-displaced seed genes
#'
#' @export
count_displaced_genes <- function(recon_sig, seed_sig) {

  if(!setequal(names(recon_sig), names(seed_sig))) {
    print("Filtering sets to matching pairs.")
    shared_vecs <- intersect(names(recon_sig), names(seed_sig))
    removed_vecs <- setdiff(names(recon_sig), names(seed_sig))
    print(paste0("Removed ", length(removed_vecs), " perturbations."))
    recon_sig <- recon_sig[shared_vecs]
    seed_sig <- seed_sig[shared_vecs]
  }

  stopifnot(all.equal(names(recon_sig), names(seed_sig)))

  displaced <- list()
  for (name in names(seed_sig)) {
    seeds_in_recon <- seed_sig[[name]] %in% recon_sig[[name]]
    n_displaced <- sum(!seeds_in_recon)
    n_not_displaced <- sum(seeds_in_recon)
    displaced[[name]] <- list(displaced = n_displaced, not_displaced = n_not_displaced)
  }

  return(displaced)
}

#' Normalize Columns of a Matrix
#'
#' @param matrix A numeric matrix to be normalized
#'
#' @return A matrix with columns normalized by their sums
col_normalize <- function(matrix) {
  return(scale(matrix, center=FALSE, scale=colSums(matrix)))
}

#' Binarize Multiple Matrices Based on Cutoffs
#'
#' @param matrix_list A list of matrices to be binarized
#' @param cutoff_list A list of cutoff values, one for each matrix
#'
#' @return A list of binarized matrices
#'
#' @details
#' This function takes a list of matrices and a corresponding list of cutoff values.
#' For each matrix, values less than or equal to the cutoff are set to 0, and values
#' greater than the cutoff are set to 1.
binarize_matrices <- function(matrix_list, cutoff_list) {

  #Two lists should be of equal length, with the same names
  stopifnot(names(matrix_list)==names(cutoff_list))
  binarized_matrices <- list()
  for (i in seq_along(cutoff_list)) {
    cutoff <- cutoff_list[[i]]
    matrix <- matrix_list[[i]]
    name <- names(matrix_list)[[i]]
    matrix_filter <- matrix <= cutoff
    matrix[matrix_filter] <- 0
    matrix[!matrix_filter] <- 1
    binarized_matrices[[name]] <- matrix
  }
  return(binarized_matrices)
}

#' Extract Largest Connected Subgraph
#'
#' @param igraph An igraph object
#'
#' @return An igraph object representing the largest connected component of the input graph
#'
#' @importFrom igraph components V subgraph
largest_connected_subgraph <- function(igraph) {
  component_filter <- components(igraph)$membership == 1
  component_nodes <- V(igraph)$name[component_filter]
  largest_subgraph <- subgraph(igraph, component_nodes)
  return(largest_subgraph)
}

# Prepare signature so that only the intersection across all graphs is used
drop_nas_list <- function(l) {
  #l is a list
  l <- l[!is.na(l)]
  return(l)
}

filter_list <- function(l, filter=1) {
  #l is a list
  l <- l[l==filter]
  return(l)
}

#' Filter Signatures for Common Nodes Across Graphs
#'
#' @param graphs A list of igraph objects
#' @param signatures A list of node signatures (vectors of node names)
#'
#' @return A list of filtered signatures containing only nodes common to all graphs
#'
#' @description
#' This function takes a list of graphs and a list of node signatures, and filters
#' each signature to include only nodes that are present in the largest connected
#' component of all provided graphs.
#'
#' @details
#' For each signature:
#' 1. It checks which nodes are in the largest connected component of each graph.
#' 2. It finds the intersection of these nodes across all graphs.
#' 3. It returns this intersection as the new filtered signature.
#'
#' @importFrom igraph components
#' @importFrom purrr reduce
#' @importFrom magrittr %>%
common_signature_filter <- function(graphs, signatures) {
  new_sigs <- list()
  for (sig_name in names(signatures)) {
    sig <- signatures[[sig_name]]
    node_common <- lapply(graphs, function(x) igraph::components(x)$membership[sig] %>% drop_nas_list %>% filter_list %>% names) %>% purrr::reduce(intersect)
    new_sigs[[sig_name]] <- node_common
  }
  return(new_sigs)
}

#' Extract Signature Nodes from Random Walk Results
#'
#' @param rwr_df A dataframe containing the results of RWR, including columns for
#'               cancer type, signature, node type, and node names.
#' @param sig_label The label of the signature to filter for.
#' @param node_label The label of the node type to filter for.
#' @param cancer_label The label of the cancer type to filter for.
#'
#' @return A vector of node names that are significant for the specified signature,
#'         node type, and cancer type.
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom dplyr pull
get_signature <- function(rwr_df, sig_label, node_label, cancer_label) {
  # rwr_df: Dataframe of stationary probabilities with cancer, signature, and node annotations
  # Return: List of node labels that are significant after doing signature propagation on cancer_label graph with sig_label sig

  sig <- rwr_df %>% dplyr::filter((Node == node_label) & (cancer == cancer_label) & (signature == sig_label)) %>% pull(name)
  return(sig)
}

get_signature_list <- function(rwr_df, sig_label, node_label="Significant") {
  # Returns a named list of significant nodes for each cancer
  cancer_labels <- unique(rwr_df[['cancer']])
  names(cancer_labels) <- cancer_labels
  return(lapply(cancer_labels, get_signature, rwr_df=rwr_df, node_label=node_label, sig_label=sig_label))
}

#' Find all pairwise distances between nodes in an igraph
#'
#' This function computes the pairwise distances between nodes within each set of nodes
#' provided in a list, based on the structure of an input graph.
#'
#' @param g An igraph object representing the graph.
#' @param node_list A list where each element is a vector of node names or IDs.
#'
#' @return A list of vectors, each containing the pairwise distances for the corresponding node set.
#'
#' @details
#' The function first calculates the full distance matrix for the graph using igraph::distances().
#' Then, for each set of nodes in `node_list`, it extracts the relevant submatrix and returns
#' the lower triangular part, which represents all pairwise distances within that set.
all_dists_nodesets <- function(g, node_list) {
  stopifnot(class(g) == "igraph")
  stopifnot(class(node_list) == "list")

  g_dist <- igraph::distances(g)

  all_dists_nodeset <- function(dist_matrix, nodes) {
    node_filter <- nodes %in% rownames(dist_matrix)
    nodes <- nodes[node_filter]
    g_dist_sig <- g_dist[nodes, nodes]
    g_dist_sig <- g_dist_sig[lower.tri(g_dist_sig)]
    return(g_dist_sig)
  }

  dists_nodesets <- lapply(node_list, all_dists_nodeset, dist_matrix=g_dist)
  return(dists_nodesets)
}

#' An empty ggplot
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 ggplot theme_void
ggempty <- function() {
  ggplot() +
    theme_void()
}

#' Enrichment plot implemented in ggplot
#'
#' @param n The length of a ranked list
#' @param positions A vector of positions in the ranked list
#' @param x_axis The x-axis of a running enrichment score
#' @param y_axis The y-axis of a running enrichment score
#' @param title Plot title
#' @return A ggplot object
#'
#' @importFrom ggplot2 qplot aes geom_rug geom_hline geom_vline annotate theme element_text element_blank element_line element_rect
ggeplot <- function(n, positions, x_axis, y_axis, title="") {
  score <- which.max(abs(y_axis))
  qplot(x_axis,
        y_axis,
        main=title,
        ylab="Running Enrichment Score",
        xlab="Position in Ranked List of Genes",
        geom="line")+
    geom_rug(data=data.frame(positions), aes(x=positions), inherit.aes=FALSE)+
    geom_hline(yintercept=0) +
    geom_vline(xintercept=n/2, linetype="dotted") +
    annotate("point", x=x_axis[score], y=y_axis[score], color="red") +
    annotate("text", x=x_axis[score]+n/20, y=y_axis[score], label=round(y_axis[score],2)) +
    annotate("point", x=x_axis[score], y=y_axis[score], color="red") +
    theme(plot.title=element_text(hjust=0.5),
          panel.background=element_blank(),
          axis.line=element_line(color="black"),
          panel.border=element_rect(color="black", fill=NA, size=1))
}
