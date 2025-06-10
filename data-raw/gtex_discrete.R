## code to prepare `gtex_discrete` dataset goes here
library(data.table)
library(dplyr)
library(DESeq2)

DATA_PATH <- ""
PATH <- ""

# 1. Load data
# Download from https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
expr <- fread(file.path(DATA_PATH, "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz"))
# DESeq expects a dense matrix
expr_matrix <- as.matrix(expr[, -(1:2), with=FALSE])
rownames(expr_matrix) <- expr$Name

## Load sample attributes
## Download from https://gtexportal.org/home/downloads/adult-gtex/metadata
sample_attr <- fread(file.path(DATA_PATH, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"))
## Load subject phenotypes
subject_pheno <- fread(file.path(DATA_PATH, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"))
sample_attr$SUBJID <- lapply(strsplit(sample_attr$SAMPID, "-"),
                             function(x) paste(x[1:2], collapse = "-")) %>% unlist
## Add subject phenotypes to sample metadata
gtex_metadata <- dplyr::left_join(sample_attr, subject_pheno, by = "SUBJID")

# 2. GTEx Analysis Parameters
## Which tissues to filter for.
## Some are in the SMTSD field, others are in the SMTS field.
tissue_filters <- list(blood = c("SMTSD", "Whole Blood"),
                       brain = c("SMTS", "Brain"))
## Which phenotype contrasts to evaluate for each tissue.
## First entry is the phenotype column name, second entry is the 'control', third entry is the 'treatment'.
contrast_filters <- list(
  blood = list(c("AGE", "20-29", "60-69"),
               c("AGE", "20-29", "70-79")),
  brain = list(c("AGE", "20-29", "60-69"),
               c("AGE", "20-29", "70-79"))
)

gtex_sample_ids <- function(gtex_metadata,
                            tissue_column,
                            tissue_name,
                            contrast_phenotype,
                            contrast_value) {

  sample_ids <- gtex_metadata %>%
    dplyr::filter(.data[[tissue_column]] == tissue_name) %>%
    dplyr::filter(.data[[contrast_phenotype]] == contrast_value) %>%
    dplyr::pull(SAMPID)
  return(sample_ids)
}

gtex_deseq_helper <- function(gtex_expr,
                              gtex_metadata,
                              tissue_filters,
                              contrast_filters) {

  stopifnot(all.equal(names(tissue_filters), names(contrast_filters)))

  # Get gtex subsets
  deseq_dfs <- list()
  for(tissue_exp in names(tissue_filters)) {
    tissue_filter <- tissue_filters[[tissue_exp]]
    tissue_contrast_filters <- contrast_filters[[tissue_exp]]

    for (contrast_filter in tissue_contrast_filters) {
      tissue_col_name <- tissue_filter[1]
      tissue_name <- tissue_filter[2]
      contrast_phenotype_name <- contrast_filter[1]
      contrast_control_value <- contrast_filter[2]
      contrast_treatment_value <- contrast_filter[3]

      ctrl_samp_ids <- gtex_sample_ids(gtex_metadata,
                                       tissue_col_name,
                                       tissue_name,
                                       contrast_phenotype_name,
                                       contrast_control_value)
      treat_samp_ids <- gtex_sample_ids(gtex_metadata,
                                        tissue_col_name,
                                        tissue_name,
                                        contrast_phenotype_name,
                                        contrast_treatment_value)

      # 1. Get the intersection of available sample IDs
      all_samp_ids <- intersect(colnames(gtex_expr), c(ctrl_samp_ids, treat_samp_ids))

      # 2. Subset expression matrix to matched and ordered sample IDs
      expr_sub <- gtex_expr[, all_samp_ids]

      # 3. Subset and order metadata to match expression matrix columns
      expr_metadata <- gtex_metadata %>%
        dplyr::filter(SAMPID %in% all_samp_ids) %>%
        dplyr::arrange(match(SAMPID, all_samp_ids))

      # 4. Convert categorical covariates to factors
      expr_metadata$AGE <- factor(expr_metadata$AGE)
      expr_metadata$SEX <- factor(expr_metadata$SEX)
      expr_metadata$SMCENTER <- factor(expr_metadata$SMCENTER)

      # 5. Center and scale numerical covariates
      expr_metadata$SMRIN <- scale(expr_metadata$SMRIN)
      expr_metadata$SMTSISCH <- scale(expr_metadata$SMTSISCH)

      # 6. Remove samples with missing data in any covariate
      complete_idx <- complete.cases(expr_metadata[, c("AGE", "SEX", "SMRIN", "SMCENTER", "SMTSISCH")])
      expr_metadata <- expr_metadata[complete_idx, ]
      expr_sub <- expr_sub[, complete_idx]

      print(tissue_name)
      print(table(expr_metadata$AGE))

      # DESeq2 analysis
      dds <- DESeqDataSetFromMatrix(
        countData = round(expr_sub),
        colData = expr_metadata,
        design = ~ AGE + SEX + SMRIN + SMCENTER + SMTSISCH
      ) %>% DESeq()

      deseq_df <- results(dds, contrast = c("AGE", contrast_treatment_value, contrast_control_value))
      deseq_df <- deseq_df %>%
        as.data.frame %>%
        tidyr::drop_na() %>%
        dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
        dplyr::arrange(desc(logFC_pval))

      deseq_name <- paste(tissue_name, paste(contrast_filter, collapse="_"), sep = "_")
      deseq_dfs[[deseq_name]] <- deseq_df
    }
  }
  return(deseq_dfs)
}

deseq_dfs <- gtex_deseq_helper(gtex_expr = expr_matrix,
                               gtex_metadata = gtex_metadata,
                               tissue_filters = tissue_filters,
                               contrast_filters = contrast_filters)
saveRDS(deseq_dfs, file.path(PATH, "deseq_dfs.rds"))
gtex.blood.aging.20.60 <- list(up = deseq_dfs$`Whole Blood_AGE_20-29_60-69` %>%
                                 dplyr::filter(log2FoldChange > 0) %>%
                                 dplyr::slice(1:100) %>%
                                 rownames,
                               up_full = deseq_dfs$`Whole Blood_AGE_20-29_60-69` %>%
                                 rownames)
gtex.blood.aging.20.70 <- list(up = deseq_dfs$`Whole Blood_AGE_20-29_70-79` %>%
                                 dplyr::filter(log2FoldChange > 0) %>%
                                 dplyr::slice(1:100) %>%
                                 rownames,
                               up_full = deseq_dfs$`Whole Blood_AGE_20-29_70-79` %>%
                                 rownames)
gtex.brain.aging.20.60 <- list(up = deseq_dfs$`Brain_AGE_20-29_60-69` %>%
                                 dplyr::filter(log2FoldChange > 0) %>%
                                 dplyr::slice(1:100) %>%
                                 rownames,
                               up_full = deseq_dfs$`Brain_AGE_20-29_60-69` %>%
                                 rownames)
gtex.brain.aging.20.70 <- list(up = deseq_dfs$`Brain_AGE_20-29_70-79` %>%
                                 dplyr::filter(log2FoldChange > 0) %>%
                                 dplyr::slice(1:100) %>%
                                 rownames,
                               up_full = deseq_dfs$`Brain_AGE_20-29_70-79` %>%
                                 rownames)
vennr::vennr(list(blood_20_60 = gtex.blood.aging.20.60$up, blood_20_70 = gtex.blood.aging.20.70$up,
                  brain_20_60 = gtex.brain.aging.20.60$up, brain_20_70 = gtex.brain.aging.20.70$up))
usethis::use_data(gtex.blood.aging.20.70, overwrite=TRUE)
usethis::use_data(gtex.brain.aging.20.70, overwrite=TRUE)
# usethis::use_data(gtex.blood.aging.20.60, overwrite=TRUE)
# usethis::use_data(gtex.brain.aging.20.60, overwrite=TRUE)

