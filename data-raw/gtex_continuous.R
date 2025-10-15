## code to prepare `gtex_continuous` dataset goes here
library(data.table)
library(dplyr)
library(DESeq2)
library(rlang)
library(Biobase)
options(box.path=file.path(Sys.getenv("CBMGIT"), "MLscripts"))
box::use(R/rm_low_rnaseq_counts)

# On the SCC DATA_PATH can be set to "/restricted/projectnb/gtex/montilab/CBMrepositoryData/GTEX_dbGap/processed_data"
DATA_PATH <- ""
PATH <- ""

# Protected data requires dbGAP application (https://www.gtexportal.org/home/protectedDataAccess)
filepaths <- Sys.glob(paste0(DATA_PATH, "/*.rds"))
deseq_dfs <- list()

for (filepath in filepaths) {
  filename <- basename(filepath)
  tissue_name <- sub("^RNAseq_GTEx_v10_(.*)\\.rds$", "\\1", filename)
  print(paste0("Processing ", tissue_name))

  # 1. Load Data
  expr <- readRDS(filepath)
  # DESeq expects a dense matrix
  expr_matrix <- as.matrix(expr@assays@data$counts)

  # Check if there's enough samples for deseq.
  if(ncol(expr) < 10) {
    print(paste0("Fewer than 10 samples, skipping ", tissue_name))
    deseq_name <- tissue_name
    deseq_dfs[[deseq_name]] <- NULL
    next
  }

  # 4. Convert categorical covariates to factors
  expr$SEX <- factor(expr$SEX)
  expr$SMCENTER <- factor(expr$SMCENTER)

  # 5. Center and scale numerical covariates
  expr$SMRIN <- scale(expr$SMRIN)
  expr$SMTSISCH <- scale(expr$SMTSISCH)

  # 6. Remove samples with missing data in any covariate
  # List of candidate covariates
  covariates <- c("AGE", "SEX", "SMRIN", "SMCENTER", "SMTSISCH")
  expr_metadata <- expr@colData %>% data.frame()
  # Check which covariates have more than one unique value
  valid_covariates <- covariates[sapply(expr_metadata[, covariates], function(x) length(unique(x)) > 1)]
  print(valid_covariates)
  design_formula <- as.formula(
    paste("~", paste(valid_covariates, collapse = " + "))
  )
  complete_idx <- complete.cases(expr_metadata[, valid_covariates])
  expr_metadata <- expr_metadata[complete_idx, ]
  expr_sub <- expr_matrix[, complete_idx]

  print(table(expr_metadata$AGE))

  # DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(
    countData = round(expr_sub),
    colData = expr_metadata,
    design = design_formula
  ) %>% DESeq()

  deseq_df <- results(dds, name = "AGE")
  deseq_df <- deseq_df %>%
    as.data.frame %>%
    tidyr::drop_na() %>%
    dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
    dplyr::arrange(desc(logFC_pval))

  deseq_name <- tissue_name
  deseq_dfs[[deseq_name]] <- deseq_df
}

saveRDS(deseq_dfs, file.path(PATH, "/data/gtex/deseq_cont_dfs_all.rds"))

deseq_cont_dfs_all <- readRDS(file.path(PATH, "/data/gtex/deseq_cont_dfs_all.rds"))

gtex_cont_sigs <- list()
for(tissue in names(deseq_cont_dfs_all)) {
  deseq_df <- deseq_cont_dfs_all[[tissue]]
  tissue_up <- deseq_df %>%
    dplyr::filter(log2FoldChange > 0) %>%
    dplyr::slice(1:100) %>%
    rownames
  tissue_up_full <- deseq_df %>%
    rownames
  gtex_cont_sigs[[tissue]][["up"]] <- tissue_up
  gtex_cont_sigs[[tissue]][["up_full"]] <- tissue_up_full
}
saveRDS(gtex_cont_sigs, file.path(PATH, "all_tissue_sigs.rds"))
gtex.aging <- gtex_cont_sigs
usethis::use_data(gtex.aging)
