## code to prepare `neurips2023` dataset goes here

library(reticulate)
# reticulate::use_condaenv("r-sceasy")
library(anndata)
library(tidyverse)
library(Matrix)
library(limma)
library(edgeR)
library(Seurat)
library(furrr)
library(future)

PATH <- ""

# Read AnnData object
# Download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945
adata <- anndata::read_h5ad(file.path(PATH, "GSE279945_sc_counts_processed.h5ad"))
seurat_obj <- CreateSeuratObject(counts = t(adata$raw$X),
                                 meta.data = adata$obs)
rownames(seurat_obj) <- colnames(adata$X)

# 1. Filter Obs (https://github.com/openproblems-bio/task_perturbation_prediction/blob/main/src/process_dataset/filter_obs/script.R)

# set up obs_filt
obs_filt <- rep(TRUE, nrow(seurat_obj@meta.data))

# Alvocidib only T cells in only 2 donors, remove
v <- obs_filt & seurat_obj@meta.data$sm_name != "Alvocidib"
# BMS-387032 - one donor with only T cells, two other consistent, but only 2 cell types - leave the 2 cell types in, remove donor 2 with only T cells
obs_filt <- obs_filt & !(seurat_obj@meta.data$sm_name == "BMS-387032" & seurat_obj@meta.data$donor_id == "Donor 2")
# BMS-387032 remove myeloid cells and B cells
obs_filt <- obs_filt & !(seurat_obj@meta.data$sm_name == "BMS-387032" & seurat_obj@meta.data$cell_type %in% c("B cells", "Myeloid cells"))
# CGP 60474 has only T cells left, remove
obs_filt <- obs_filt & seurat_obj@meta.data$sm_name != "CGP 60474"
# Canertinib - the variation of Myeloid cell proportions is very large, skip Myeloid
obs_filt <- obs_filt & !(seurat_obj@meta.data$sm_name == "Canertinib" & seurat_obj@meta.data$cell_type == "Myeloid cells")
# Foretinib - large variation in Myeloid cell proportions (some in T cells), skip Myeloid.
obs_filt <- obs_filt & !(seurat_obj@meta.data$sm_name == "Foretinib" & seurat_obj@meta.data$cell_type == "Myeloid cells")
# Ganetespib (STA-9090) - donor 2 has no Myeloid and small NK cells proportions. Skip Myeloid, remove donor 2
obs_filt <- obs_filt & !(seurat_obj@meta.data$sm_name == "Ganetespib (STA-9090)" & seurat_obj@meta.data$donor_id == "Donor 2")
# IN1451 - donor 2 has no NK or B, remove Donor 2
obs_filt <- obs_filt & !(seurat_obj@meta.data$sm_name == "IN1451" & seurat_obj@meta.data$donor_id == "Donor 2")
# Navitoclax - donor 3 doesn't have B cells and has different T and Myeloid proportions, remove donor 3
obs_filt <- obs_filt & !(seurat_obj@meta.data$sm_name == "Navitoclax" & seurat_obj@meta.data$donor_id == "Donor 3")
# PF-04691502 remove Myeloid (only present in donor 3)
obs_filt <- obs_filt & !(seurat_obj@meta.data$sm_name == "PF-04691502" & seurat_obj@meta.data$cell_type == "Myeloid cells")
# Proscillaridin A;Proscillaridin-A remove Myeloid, since the variation is very high (4x)
obs_filt <- obs_filt & !(seurat_obj@meta.data$sm_name == "Proscillaridin A;Proscillaridin-A" & seurat_obj@meta.data$cell_type == "Myeloid cells")
# R428 - skip NK due to high variation (close to 3x)
obs_filt <- obs_filt & !(seurat_obj@meta.data$sm_name == "R428" & seurat_obj@meta.data$cell_type == "NK cells")
# UNII-BXU45ZH6LI - remove due to large variation across all cell types and missing cell types
obs_filt <- obs_filt & seurat_obj@meta.data$sm_name != "UNII-BXU45ZH6LI"

# filter obs
seurat_obj <- seurat_obj[,obs_filt]

# 2. Compute Pseudobulk
seurat_pb <- Seurat::AggregateExpression(seurat_obj,
                                         group.by = "plate_well_celltype_reannotated",
                                         return.seurat = TRUE)
## Adding Pseudobulk metadata

seurat_pb$plate_name <- lapply(strsplit(seurat_pb$plate_well_celltype_reannotated, "-"),
                               function(x) paste(x[[1]], x[[2]], sep = "-")) %>% unlist %>% unname
seurat_pb$well <- lapply(strsplit(seurat_pb$plate_well_celltype_reannotated, "-"),
                         function(x) x[[3]]) %>% unlist %>% unname
seurat_pb$cell_type <- lapply(strsplit(seurat_pb$plate_well_celltype_reannotated, "-"),
                              function(x) x[[4]]) %>% unlist %>% unname
seurat_pb$plate_well_celltype_reannotated <- paste(seurat_pb$plate_name,seurat_pb$well,seurat_pb$cell_type, sep = "_")
pb_metadata <- dplyr::left_join(seurat_pb@meta.data,
                                seurat_obj@meta.data %>% dplyr::select(-c(orig.ident, nCount_RNA, nFeature_RNA, well,
                                                                          cell_count_by_plate_well, cell_count_by_well_celltype, plate_name,
                                                                          cell_type)),
                                by = "plate_well_celltype_reannotated",
                                multiple = "any")
rownames(pb_metadata) <- pb_metadata$orig.ident
seurat_pb@meta.data <- pb_metadata

# Remove samples with no counts
sum(rowSums(seurat_pb@assays$RNA$counts) == 0)
## There are no samples with no counts

# 3. Filter Var
cat("Filtering variables\n")
# Transform function similar to limma's requirements in the Python script
limma_trafo <- function(value) {
  gsub("[^[:alnum:]_]", "_", value)
}

# Unique cell types for processing
genes_to_keep <- lapply(
  unique(seurat_pb@meta.data$cell_type),
  function(cell_type) {
    input_ct <- seurat_pb[,seurat_pb@meta.data$cell_type == cell_type]

    # Prepare count data for edgeR analysis
    counts <- input_ct@assays$RNA$counts

    d <- DGEList(counts)
    design <- model.matrix(~ 0 + sm_name + plate_name, data = input_ct@meta.data %>% mutate_all(limma_trafo))
    keep <- filterByExpr(d, design)

    rownames(d)[keep]
  }
)


# Calculate the intersection of genes across all cell types
gene_filt <- Reduce(intersect, genes_to_keep)

# filter genes
seurat_pb <- seurat_pb[gene_filt, ]

# 4. Limma
seurat_pb$sm_cell_type <- paste(seurat_pb$sm_name, seurat_pb$cell_type, sep = "_")

# select [cell_type, sm_name] pairs which will be used for DE analysis
non_control_obs <- seurat_pb@meta.data %>%
  select(sm_cell_type, cell_type, sm_name, sm_lincs_id, SMILES, split, control) %>%
  distinct() %>%
  filter(sm_name != "Dimethyl Sulfoxide")

# check which cell_types to run limma for
cell_types <- as.character(unique(non_control_obs$cell_type))

d0 <- seurat_pb@assays$RNA$counts %>%
  edgeR::DGEList() %>%
  edgeR::calcNormFactors()

design_matrix <- model.matrix(~ 0 + sm_cell_type + plate_name, seurat_pb@meta.data %>% mutate_all(limma_trafo))

# Voom transformation and lmFit
v <- limma::voom(d0, design = design_matrix, plot = TRUE)
fit <- limma::lmFit(v, design_matrix)

# run limma DE for each cell type and compound
de_df <- furrr::future_map_dfr(
  seq_len(nrow(non_control_obs)),
  .options = furrr_options(seed = TRUE),
  function(row_i) {
    cat("Computing DE contrasts (", row_i, "/", nrow(non_control_obs), ")\n", sep = "")
    sm_cell_type <- as.character(non_control_obs$sm_cell_type[[row_i]])
    cell_type <- as.character(non_control_obs$cell_type[[row_i]])

    control_name <- paste("Dimethyl Sulfoxide", cell_type, sep = "_")
    # run contrast fit
    contrast_formula <- paste0(
      "sm_cell_type", limma_trafo(sm_cell_type),
      " - ",
      "sm_cell_type", limma_trafo(control_name)
    )
    contr <- limma::makeContrasts(
      contrasts = contrast_formula,
      levels = colnames(coef(fit))
    )

    limma::contrasts.fit(fit, contr) %>%
      limma::eBayes(robust = TRUE) %>%
      limma::topTable(n = Inf, sort = "none") %>%
      rownames_to_column("gene") %>%
      mutate(row_i = row_i)
  }
)

# transform data
de_df2 <- de_df %>%
  mutate(
    # convert gene names to factor
    gene = factor(gene),
    # readjust p-values for multiple testing
    adj.P.Value = p.adjust(P.Value, method = "BH"),
    # compute sign fc × log10 p-values
    sign_log10_pval = sign(logFC) * -log10(ifelse(P.Value == 0, .Machine$double.eps, P.Value)),
    sign_log10_adj_pval = sign(logFC) * -log10(ifelse(adj.P.Value == 0, .Machine$double.eps, adj.P.Value)),
    # determine if gene is DE
    is_de = P.Value < 0.05,
    is_de_adj = adj.P.Value < 0.05,
    # compute clipped sign fc × log10 p-values
    clipped_sign_log10_pval = sign(logFC) * -log10(pmax(1e-16, P.Value)),
  ) %>%
  as_tibble()


cat("DE df:\n")
print(head(de_df2))

# Convert to nested list format

de_df2 <- de_df2 %>% dplyr::select(-c("adj.P.Val"))
de_df2$group <- rownames(non_control_obs)[de_df2$row_i]
de_df2$cell_type <- non_control_obs$cell_type[de_df2$row_i]
de_df2$sm_name <- non_control_obs$sm_name[de_df2$row_i]
de_df2$sm_cell_type <- paste(de_df2$sm_name, de_df2$cell_type, sep = "_")

## Row 558 is a repeated example of B-Cell treated with Tivozanib
de_df2 <- de_df2 %>% dplyr::filter(row_i != 558)
de_df2_list <- de_df2 %>%
  dplyr::group_by(sm_cell_type) %>%
  dplyr::mutate(logFC_adjpval = logFC * abs(sign_log10_adj_pval)) %>%
  dplyr::arrange(desc(logFC_adjpval)) %>%
  dplyr::group_split()

names(de_df2_list) <- lapply(de_df2_list, function(x) unique(x$sm_cell_type)) %>% unlist

b_sigs <- list()
myeloid_sigs <- list()
t_sigs <- list()
nk_sigs <- list()

for(pb_df in de_df2_list) {
  celltype <- unique(pb_df$cell_type)
  sm_name <- unique(pb_df$sm_name) %>% as.character
  sm_celltype <- unique(pb_df$sm_cell_type)

  up_sig <- pb_df %>% dplyr::filter(logFC > 0) %>% dplyr::slice(1:100) %>% dplyr::pull(gene)
  full_sig <- pb_df %>% pull(gene)

  if(celltype == "B cells") {
    b_sigs[[sm_name]][["up"]] <- up_sig %>% as.character
    b_sigs[[sm_name]][["up_full"]] <- full_sig %>% as.character
  } else if (celltype == "Myeloid cells") {
    myeloid_sigs[[sm_name]][["up"]] <- up_sig %>% as.character
    myeloid_sigs[[sm_name]][["up_full"]] <- full_sig %>% as.character
  } else if (celltype == "T cells") {
    t_sigs[[sm_name]][["up"]] <- up_sig %>% as.character
    t_sigs[[sm_name]][["up_full"]] <- full_sig %>% as.character
  } else if (celltype == "NK cells") {
    nk_sigs[[sm_name]][["up"]] <- up_sig %>% as.character
    nk_sigs[[sm_name]][["up_full"]] <- full_sig %>% as.character
  }
}

shared_drugs <- purrr::reduce(list(names(b_sigs),
                                   names(t_sigs),
                                   names(nk_sigs),
                                   names(myeloid_sigs)), intersect)

neurips2023.b <- b_sigs[shared_drugs]
neurips2023.nk <- nk_sigs[shared_drugs]
neurips2023.t <- t_sigs[shared_drugs]
neurips2023.myeloid <- myeloid_sigs[shared_drugs]

usethis::use_data(neurips2023.b, overwrite = TRUE)
usethis::use_data(neurips2023.nk, overwrite = TRUE)
usethis::use_data(neurips2023.t, overwrite = TRUE)
usethis::use_data(neurips2023.myeloid, overwrite = TRUE)
