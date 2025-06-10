## code to prepare `gtex_continuous` dataset goes here
library(data.table)
library(dplyr)
library(DESeq2)
library(rlang)

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
gtex_metadata$AGE_ordinal <- with(gtex_metadata, dplyr::case_when(is.na(AGE) ~ NA,
                                                                  AGE == "20-29" ~ 1,
                                                                  AGE == "30-39" ~ 2,
                                                                  AGE == "40-49" ~ 3,
                                                                  AGE == "50-59" ~ 4,
                                                                  AGE == "60-69" ~ 5,
                                                                  AGE == "70-79" ~ 6,
))

# 2. GTEx Analysis Parameters
## Which tissues to filter for.
## Some are in the SMTSD field, others are in the SMTS field.
tissue_filters <- list(blood = c("SMTSD", "Whole Blood"),
                       brain = c("SMTS", "Brain"))

gtex_sample_ids <- function(gtex_metadata,
                            tissue_column,
                            tissue_name) {

  sample_ids <- gtex_metadata %>%
    dplyr::filter(.data[[tissue_column]] == tissue_name) %>%
    dplyr::pull(SAMPID)
  return(sample_ids)
}

gtex_deseq_helper <- function(gtex_expr,
                              gtex_metadata,
                              tissue_filters) {

  # Get gtex subsets
  deseq_dfs <- list()
  for(tissue_exp in names(tissue_filters)) {
    tissue_filter <- tissue_filters[[tissue_exp]]
    tissue_col_name <- tissue_filter[1]
    tissue_name <- tissue_filter[2]

    samp_ids <- gtex_sample_ids(gtex_metadata,
                                tissue_col_name,
                                tissue_name)

    # 1. Get the intersection of available sample IDs
    all_samp_ids <- intersect(colnames(gtex_expr), samp_ids)

    # 2. Subset expression matrix to matched and ordered sample IDs
    expr_sub <- gtex_expr[, all_samp_ids]

    # 3. Subset and order metadata to match expression matrix columns
    expr_metadata <- gtex_metadata %>%
      dplyr::filter(SAMPID %in% all_samp_ids) %>%
      dplyr::arrange(match(SAMPID, all_samp_ids))

    # 4. Convert categorical covariates to factors
    expr_metadata$SEX <- factor(expr_metadata$SEX)
    expr_metadata$SMCENTER <- factor(expr_metadata$SMCENTER)

    # 5. Center and scale numerical covariates
    expr_metadata$SMRIN <- scale(expr_metadata$SMRIN)
    expr_metadata$SMTSISCH <- scale(expr_metadata$SMTSISCH)

    # 6. Remove samples with missing data in any covariate
    complete_idx <- complete.cases(expr_metadata[, c("AGE_ordinal", "SEX", "SMRIN", "SMCENTER", "SMTSISCH")])
    expr_metadata <- expr_metadata[complete_idx, ]
    expr_sub <- expr_sub[, complete_idx]

    print(tissue_name)
    print(table(expr_metadata$AGE))

    # DESeq2 analysis
    dds <- DESeqDataSetFromMatrix(
      countData = round(expr_sub),
      colData = expr_metadata,
      design = ~ AGE_ordinal + SEX + SMRIN + SMCENTER + SMTSISCH
    ) %>% DESeq()

    deseq_df <- results(dds, name = "AGE_ordinal")
    deseq_df <- deseq_df %>%
      as.data.frame %>%
      tidyr::drop_na() %>%
      dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
      dplyr::arrange(desc(logFC_pval))

    deseq_name <- tissue_name
    deseq_dfs[[deseq_name]] <- deseq_df
  }
  return(deseq_dfs)
}

deseq_dfs <- gtex_deseq_helper(gtex_expr = expr_matrix,
                               gtex_metadata = gtex_metadata,
                               tissue_filters = tissue_filters)
saveRDS(deseq_dfs, file.path(PATH, "deseq_cont_dfs.rds"))
gtex.blood.aging <- list(up = deseq_dfs$`Whole Blood` %>%
                           dplyr::filter(log2FoldChange > 0) %>%
                           dplyr::slice(1:100) %>%
                           rownames,
                         up_full = deseq_dfs$`Whole Blood` %>%
                           rownames)
gtex.brain.aging <- list(up = deseq_dfs$Brain %>%
                           dplyr::filter(log2FoldChange > 0) %>%
                           dplyr::slice(1:100) %>%
                           rownames,
                         up_full = deseq_dfs$Brain %>%
                           rownames)

vennr::vennr(list(blood = gtex.blood.aging$up,brain = gtex.brain.aging$up))
usethis::use_data(gtex.blood.aging, overwrite = TRUE)
usethis::use_data(gtex.brain.aging, overwrite = TRUE)

