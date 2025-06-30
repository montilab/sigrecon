library(tidyverse)
library(Seurat)
library(anndata)
library(reticulate)
reticulate::use_condaenv("r-sceasy")
library(SeuratDisk)
library(Matrix)
library(doParallel)
library(DESeq2)
registerDoParallel(15)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATAPATH <- file.path(Sys.getenv("CBM"), "otherStudies/perturbational_data/tahoe")
do_save <- FALSE

rm_low_rnaseq_counts <- function(
    eset,
    min_samples=NULL,
    class_id=NULL,
    assay_id=NULL,
    reads_per = 1000000,
    min_thresh=1
)
{
  ## BEGIN input checks
  stopifnot( methods::is(eset,"ExpressionSet") || methods::is(eset,"SummarizedExperiment") || methods::is(eset, "Seurat") )
  stopifnot( xor(is.null(min_samples),is.null(class_id)) )
  ## END input checks

  if ( methods::is(eset,"ExpressionSet") )
  {
    if ( !is.null(class_id)) {
      stopifnot( class_id %in% colnames(Biobase::pData(eset)) )
      groups <- Biobase::pData(eset)[,class_id]
      min_samples <- max(min_thresh,table(groups))
    }
    rpm <- colSums(Biobase::exprs(eset))/reads_per
    filter_ind <- t(apply(Biobase::exprs(eset), 1,function(x) {x>rpm}))
    filter_ind_rowsums <- apply(filter_ind, 1, sum)
    return(eset[filter_ind_rowsums >= min_samples,])
  }
  else if ( methods::is(eset,"SummarizedExperiment") )
  {
    if ( is.null(assay_id) )
      assay_id <- names(SummarizedExperiment::assays(eset))[1]
    stopifnot( assay_id %in% names(SummarizedExperiment::assays(eset)) )

    if ( !is.null(class_id)) {
      stopifnot( class_id %in% colnames(SummarizedExperiment::colData(eset)) )
      groups <- SummarizedExperiment::colData(eset)[,class_id]
      min_samples <- max(min_thresh,table(groups))
    }
    counts <- SummarizedExperiment::assays(eset)[[assay_id]]
    rpm <- colSums(counts)/reads_per
    filter_ind <- t(apply(counts, 1, function(x) {x>rpm}))
    filter_ind_rowsums <- apply(filter_ind, 1, sum)
    return(eset[filter_ind_rowsums >= min_samples,])
  }
  else if ( methods::is(eset,"Seurat") )
  {
    # Default to RNA assay if not specified
    if (is.null(assay_id)) assay_id <- "RNA"
    stopifnot(assay_id %in% names(eset@assays))
    counts <- Seurat::GetAssayData(eset, assay = assay_id, slot = "counts")
    meta <- eset@meta.data

    # If class_id is specified, use it to determine min_samples
    if (!is.null(class_id)) {
      stopifnot(class_id %in% colnames(meta))
      groups <- meta[, class_id]
      min_samples <- max(min_thresh, table(groups))
    }
    rpm <- colSums(counts) / reads_per
    filter_ind <- t(apply(counts, 1, function(x) {x > rpm}))
    filter_ind_rowsums <- rowSums(filter_ind)
    keep_genes <- which(filter_ind_rowsums >= min_samples)
    # Subset Seurat object to keep only the filtered genes
    return(subset(eset, features = rownames(counts)[keep_genes]))
  }
  else
    stop( "unrecognized oject type: ", class(eset) )
}

if (do_save) {
  # Loading data
  adata <- anndata::read_h5ad(file.path(DATAPATH, "pseudobulk/merged_pseudobulk.h5ad"))
  # Convert dgRMatrix to dgCMatrix and transpose
  mat_t <- as(adata$X, "CsparseMatrix")      # Convert to dgCMatrix
  mat_c_t <- t(mat_t)
  seurat_obj <- CreateSeuratObject(counts = mat_c_t,
                                   meta.data = adata$obs)

  saveRDS(seurat_obj, file.path(DATAPATH, "pseudobulk/merged_pseudobulk.rds"))
} else {
  seurat_obj <- readRDS(file.path(DATAPATH, "pseudobulk/merged_pseudobulk.rds"))
}

highest_dose_filter_metadata <- seurat_obj@meta.data %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  dplyr::mutate(drug_conc = as.numeric(drug_conc)) %>%
  dplyr::group_by(cell_name, drug_name) %>%
  dplyr::filter(drug_conc == max(drug_conc, na.rm = TRUE)) %>%
  dplyr::ungroup()

seurat_obj_f <- seurat_obj[,highest_dose_filter_metadata$cell_id]
seurat_obj_f <- rm_low_rnaseq_counts(seurat_obj_f, min_samples = 3)

celllines <- seurat_obj_f$cell_name %>% unique
drugs <- seurat_obj_f$drug_name %>% unique

sig_dfs <- foreach(cell = celllines, .combine = dplyr::bind_rows) %:%
  foreach (drug = drugs) %dopar% {
    print(paste(drug, cell))

    meta <- seurat_obj_f@meta.data
    meta$cell_id <- rownames(meta)
    drug_cells <- meta %>% filter(cell_name == cell & drug_name == drug)
    control_cells <- meta %>% filter(cell_name == cell & drug_name == "DMSO_TF")

    # Combine for DESeq2
    sub_cells <- c(drug_cells$cell_id, control_cells$cell_id)
    sub_meta <- meta %>% dplyr::filter(cell_id %in% sub_cells)
    sub_meta$condition <- ifelse(sub_meta$drug_name == drug, "treatment", "control")
    if(sum(sub_meta$condition == "treatment") < 2) {
      print(paste0("Skipping ", drug, " for ", cell, ". Fewer than 2 replicates."))
      return(NULL)
    } else if (sum(sub_meta$condition == "control") < 2) {
      print(paste0("Skipping ", drug, " for ", cell, ". Fewer than 2 controls."))
      return(NULL)
    }
    sub_count <- seurat_obj_f[,sub_cells]@assays$RNA$count

    # Run DESeq2
    dds <- DESeqDataSetFromMatrix(countData = sub_count, colData = sub_meta, design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", "treatment", "control"))
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$cell_line <- cell
    res_df$drug <- drug
    res_df <- res_df %>%
      tidyr::drop_na() %>%
      dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
      dplyr::arrange(desc(logFC_pval))
    res_df
  }

# saveRDS(sig_dfs, file.path(PATH, "data/tahoe/tahoe_deseq_dfs.rds"))

tahoe_sigs <- list()
for(sig_df in sig_dfs) {
  cell_line <- unique(sig_df$cell_line) %>% as.character
  drug <- unique(sig_df$drug) %>% as.character

  up_sig <- sig_df %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::slice(1:100) %>% dplyr::pull(gene)
  full_sig <- sig_df %>% pull(gene)


  tahoe_sigs[[cell_line]][[drug]][["up"]] <- up_sig
  tahoe_sigs[[cell_line]][[drug]][["up_full"]] <- full_sig
}

# saveRDS(tahoe_sigs, file.path(PATH, "data/tahoe/tahoe_sigs.rds"))

# Subsetting each cell line to shared drugs (110 drugs for 50 cell lines)
pbs <- lapply(tahoe_sigs, function(x) names(x)) %>% purrr::reduce(., intersect)
tahoe_sigs_filtered <- lapply(tahoe_sigs, function(x) x[pbs])
saveRDS(tahoe_sigs_filtered, file.path(PATH, "data/tahoe/tahoe_sigs_filtered.rds"))

usethis::use_data(tahoe_sigs_filtered)
