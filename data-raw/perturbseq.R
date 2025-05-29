## code to prepare `perturbseq` dataset goes here
library(tidyverse)
library(anndata)
library(reticulate)
library(Seurat)

PATH <- ""

# Load AnnData
# Download from here https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387
rpe1_essential_raw_sc_path <- file.path(PATH, "rpe1_raw_singlecell_01.h5ad")
rpe1_essential_raw_sc_data <- anndata::read_h5ad(rpe1_essential_raw_sc_path)

k562_essential_raw_sc_path <- file.path(PATH, "K562_essential_raw_singlecell_01.h5ad")
k562_essential_raw_sc_data <- read_h5ad(k562_essential_raw_sc_path)

# Changing to Gene Symbols
colnames(rpe1_essential_raw_sc_data) <- rpe1_essential_raw_sc_data$var$gene_name
colnames(k562_essential_raw_sc_data) <- k562_essential_raw_sc_data$var$gene_name
rpe1_essential_raw_sc_data <- rpe1_essential_raw_sc_data[,!duplicated(colnames(rpe1_essential_raw_sc_data))]
k562_essential_raw_sc_data <- k562_essential_raw_sc_data[,!duplicated(colnames(k562_essential_raw_sc_data))]

perturb_seq_sig <- function(ann_data_sc, gene) {
  tryCatch({
    gene_filter <- ann_data_sc$obs$gene == gene
    if(sum(gene_filter) <= 30) {
      print(paste0("Skipping ", gene, ". Fewer than 30 cells."))
      return(NULL)
    }
    control_filter <- ann_data_sc$obs$gene == "non-targeting"
    gene_or_control_filter <- gene_filter | control_filter
    gene_exp <- ann_data_sc[gene_or_control_filter,]

    # Making cell ids unique
    cell_ids <- rownames(gene_exp$X)
    if (any(duplicated(cell_ids))) {
      print("Making cell ids non-duplicated.")
      non_dup_ids <- make.unique(cell_ids, sep = "_")
      rownames(gene_exp$X) <- non_dup_ids
      rownames(gene_exp$obs) <- non_dup_ids
    }
    seurat_obj <- CreateSeuratObject(counts = t(gene_exp$X),
                                     meta.data = gene_exp$obs)
    seurat_obj <- NormalizeData(seurat_obj)

    Idents(seurat_obj) <- seurat_obj$gene
    markers <- FindMarkers(seurat_obj,
                           ident.1 = gene,
                           ident.2 = "non-targeting",
                           assay = "RNA",
                           slot = "data",
                           test.use = "MAST",
                           latent.vars = "gem_group",
                           logfc.threshold = 0,
                           only.pos = FALSE,
                           min.cells.group = 3,
                           verbose = TRUE)
    return(markers)
  }, error = function(e) {
    message(paste0("Error in gene ", gene, ": ", e$message))  # Use message() for cleaner output
    return(NULL)
  })
}

shared_essential_genes <- intersect(k562_essential_raw_sc_data$obs$gene, rpe1_essential_raw_sc_data$obs$gene)

k562_sig <- foreach (gene=shared_essential_genes, .combine=c) %dopar% {
  k562_sigs <- list()
  k562_sigs[[gene]] <- perturb_seq_sig(k562_essential_raw_sc_data,
                                       gene)
  k562_sigs
}

rpe1_sig <- foreach (gene=shared_essential_genes, .combine=c) %dopar% {
  rpe1_sigs <- list()
  rpe1_sigs[[gene]] <- perturb_seq_sig(rpe1_essential_raw_sc_data,
                                       gene)
  rpe1_sigs
}

pull_significant_signatures <- function(dge_tbl, pb_gene, alpha=0.05, limit = 100) {
  stopifnot(all(c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj") %in% colnames(dge_tbl)))

  # First confirm that CRISPR target is significantly downregulated
  lfc_pbgene <- dge_tbl[pb_gene,2]
  pval_pbgene <- dge_tbl[pb_gene,5]
  if(is.na(lfc_pbgene) | is.na(pval_pbgene)) {
    print(paste("Perturbation did not effectively target gene of interest. Skipping", pb_gene))
    return(NULL)
  } else if (!(lfc_pbgene < 0 )) {
    print(paste("Perturbation of target gene does not decrease expression. Skipping", pb_gene))
    return(NULL)
  } else if (!(pval_pbgene < alpha)) {
    print(paste("Perturbation of target gene does not decrease expression significantly. Skipping", pb_gene))
    return(NULL)
  }
  pb_data <- dge_tbl %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::filter(gene == pb_gene)
  pb_logfc <- pb_data[,"avg_log2FC"]
  pb_pval <- pb_data[,"p_val_adj"]
  print(paste(pb_pval, pb_logfc))
  pb_filter <- (pb_logfc < 0) & (pb_pval < alpha)

  if (!pb_filter) {
    print("Perturbation did not effectively target gene of interest. Skipping")
    return(NULL)
  }

  dge_tbl$logFC_adjpval <- dge_tbl$avg_log2FC * (-log10(dge_tbl$p_val_adj))
  sig_sig <- dge_tbl %>%
    dplyr::arrange(desc(logFC_adjpval)) %>%
    dplyr::filter(avg_log2FC > 0) %>%
    dplyr::filter(p_val_adj < alpha) %>%
    dplyr::slice(1:limit) %>%
    rownames
  full_sig <- dge_tbl %>%
    dplyr::arrange(desc(logFC_adjpval)) %>%
    rownames
  return(list(up_full = full_sig, up=sig_sig))
}

k562_sigs <- list()
for(i in seq_along(k562_sig)) {
  sig_name <- names(k562_sig)[[i]]
  print(sig_name)
  sig_tbl <- k562_sig[[i]]
  k562_sigs[[sig_name]] <- pull_significant_signatures(sig_tbl, pb_gene = sig_name)
}

rpe1_sigs <- list()
for(i in seq_along(rpe1_sig)) {
  sig_name <- names(rpe1_sig)[[i]]
  print(sig_name)
  sig_tbl <- rpe1_sig[[i]]
  rpe1_sigs[[sig_name]] <- pull_significant_signatures(sig_tbl, pb_gene = sig_name)
}

k562_strong_weak_filter <- lapply(k562_sigs, function(x) length(x$up) >=5) %>% unlist
rpe1_strong_weak_filter <- lapply(rpe1_sigs, function(x) length(x$up) >=5) %>% unlist
k562_sigs <- k562_sigs[k562_strong_weak_filter]
rpe1_sigs <- rpe1_sigs[rpe1_strong_weak_filter]
shared_pb_genes <- intersect(names(k562_sigs), names(rpe1_sigs))
perturbseq.k562 <- k562_sigs[shared_pb_genes]
perturbseq.rpe1 <- rpe1_sigs[shared_pb_genes]

usethis::use_data(perturbseq.k562, overwrite = TRUE)
usethis::use_data(perturbseq.rpe1, overwrite = TRUE)
