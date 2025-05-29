## code to prepare `drugmatrix` dataset goes here
library(tidyverse)
library(edgeR)

PATH <- ""
DATA_PATH <- ""
SAVE_PATH <- ""

# Download data from https://ntp.niehs.nih.gov/data/drugmatrix
liver_elist <- readRDS(file=file.path(DATA_PATH, "liver.rds"))
kidney_elist <- readRDS(file=file.path(DATA_PATH, "kidney.rds"))

# Filtering Metadata
liver_samples <- liver_elist %>% pData %>%
  filter(`dose:ch1` != "0 mg/kg") %>%
  select(title, source_name_ch1, `compound:ch1`,`dose:ch1`, `tissue:ch1`, `time:ch1`, `vehicle:ch1`) %>%
  mutate(new_dose = as.numeric(str_extract(`dose:ch1`, "\\d+")), new_time = as.numeric(str_extract(`time:ch1`, "\\d+"))) %>%
  rownames_to_column("sample_id") %>%
  group_by(`compound:ch1`) %>%
  filter((new_dose == max(new_dose))) %>%
  filter(new_time == max(new_time))

kidney_samples <- kidney_elist  %>% pData %>%
  filter(`dose:ch1` != "0 mg/kg") %>%
  select(title, source_name_ch1, `compound:ch1`,`dose:ch1`, `tissue:ch1`, `time:ch1`, `vehicle:ch1`) %>%
  mutate(new_dose = as.numeric(str_extract(`dose:ch1`, "\\d+")), new_time = as.numeric(str_extract(`time:ch1`, "\\d+"))) %>%
  rownames_to_column("sample_id") %>%
  group_by(`compound:ch1`) %>%
  filter((new_dose == max(new_dose))) %>%
  filter(new_time == max(new_time))

# Shared Drugs
liver_kidney_drugs <- intersect(liver_samples$`compound:ch1`, kidney_samples$`compound:ch1`)
#all_drugs <- intersect(liver_kidney_drugs, cultured_liver_samples$`compound:ch1`)

# Limma - Lung
for(drug in liver_kidney_drugs) {
  drug_samples <- liver_samples %>% filter(`compound:ch1`== drug)
  drug_vehicle <- drug_samples %>% pull(`vehicle:ch1`)
  drug_time <- drug_samples %>% pull(`time:ch1`)
  stopifnot(length(unique(drug_vehicle)) == 1)
  stopifnot(length(unique(drug_time)) == 1)
  drug_vehicle <- unique(drug_vehicle)
  drug_time <- unique(drug_time)
  ctrl_samples <- liver_elist %>% pData %>% filter(`dose:ch1` == "0 mg/kg",
                                                   `vehicle:ch1` == drug_vehicle,
                                                   `time:ch1` == drug_time)
  ctrl_ids <- ctrl_samples %>% rownames
  drug_ids <- drug_samples %>% pull(sample_id)

  eset_sub <- liver_elist[,c(ctrl_ids, drug_ids)]
  sampleNames(eset_sub) <- c(paste("Ctrl", seq(1:length(ctrl_ids))), paste(drug, seq(1:length(drug_ids))))
  design_df <- data.frame(sample = sampleNames(eset_sub),
                          group = c(rep("Ctrl", length(ctrl_ids)),
                                    rep("Drug", length(drug_ids))))
  group <- as.factor(design_df$group)
  boxplot(exprs(eset_sub), col=as.numeric(group), las=2, cex.axis = 0.7, boxwex=0.6, outline=FALSE)

  # Limma
  exprs(eset_sub) <- log2(exprs(eset_sub))
  design <- model.matrix(~ 0 + group, eset_sub)
  colnames(design) <- levels(group)

  # Fit linear model for each gene given a series of arrays
  fit <- limma::lmFit(eset_sub, design)
  # Build comparison and compute the satistics
  cont.matrix <- limma::makeContrasts(Drug-Ctrl, levels=design)
  # Apply the contrast matrix to the linear model
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  # Generate the statistics
  fit2 <- limma::eBayes(fit2, 0.05)

  de_genes_drug_ctrl <- limma::topTable(fit2, coef=1, adjust="fdr", number=Inf)
  de_table <- cbind(de_genes_drug_ctrl,
                    fData(eset_sub)[rownames(de_genes_drug_ctrl),c("Gene Symbol", "ENTREZ_GENE_ID")])
  write.csv(de_table, file.path(PATH, paste0("data/sigs/drugmatrix/tables/", drug, "_liver_sigs.csv")), row.names = TRUE)
}

sig_table_paths <- Sys.glob(file.path(PATH, "data/sigs/drugmatrix/tables/*_liver_sigs.csv"))
drugs <- str_extract(sig_table_paths, "(?<=tables/)(.*?)(?=_liver_sigs\\.csv)")
sig_tables <- lapply(sig_table_paths, read.csv)
names(sig_tables) <- drugs

# Extracting Signatures
fc_threshold <- 0
p_val_threshold <- 0.05
max_signature <- 100
sigs <- list()
for(drug in drugs) {
  table <- sig_tables[[drug]]
  print(drug)
  if(nrow(table) == 0) {
    next
    paste("skipping", drug, "no DEGS.")
  }
  print(table %>% dim)

  table$logFC_adjpval <- table$logFC * (-log10(table$adj.P.Val))
  up_sig <- table %>%
    dplyr::filter(logFC > fc_threshold, adj.P.Val < p_val_threshold) %>%
    dplyr::arrange(desc(logFC_adjpval)) %>%
    pull(Gene.Symbol.1)
  up_full_sig <- table %>%
    dplyr::arrange(desc(logFC_adjpval)) %>%
    pull(Gene.Symbol.1)
  dn_sig <- table %>%
    dplyr::filter(logFC < -fc_threshold, adj.P.Val < p_val_threshold) %>%
    dplyr::arrange(logFC_adjpval) %>%
    pull(Gene.Symbol.1)
  dn_full_sig <- table %>%
    dplyr::arrange(logFC_adjpval) %>%
    pull(Gene.Symbol.1)

  up_sig_f <- up_sig[up_sig != ""][1:max_signature]
  dn_sig_f <- dn_sig[dn_sig != ""][1:max_signature]
  up_full_sig <- up_full_sig[up_full_sig != ""]
  dn_full_sig <- dn_full_sig[dn_full_sig != ""]

  if (identical(up_sig, character(0)) | identical(dn_sig, character(0))) {
    paste("skipping", drug, "no significant genes.")
    next
  } else if (length(up_sig) < 10 | length(dn_sig) < 10) {
    paste("skipping", drug, "not enough genes.")
    next
  }
  print(length(up_sig))
  print(length(dn_sig))
  sigs[[drug]] <- list(up = up_sig_f, dn = dn_sig_f,
                       up_full = up_full_sig, dn_full = dn_full_sig)
}

liver_sigs <- sigs

# Limma - Kidney
if(do_save) {
  for(drug in liver_kidney_drugs) {
    drug_samples <- kidney_samples %>% filter(`compound:ch1`== drug)
    drug_vehicle <- drug_samples %>% pull(`vehicle:ch1`)
    drug_time <- drug_samples %>% pull(`time:ch1`)
    stopifnot(length(unique(drug_vehicle)) == 1)
    stopifnot(length(unique(drug_time)) == 1)
    drug_vehicle <- unique(drug_vehicle)
    drug_time <- unique(drug_time)
    ctrl_samples <- kidney_elist %>% pData %>% filter(`dose:ch1` == "0 mg/kg",
                                                      `vehicle:ch1` == drug_vehicle,
                                                      `time:ch1` == drug_time)
    ctrl_ids <- ctrl_samples %>% rownames
    drug_ids <- drug_samples %>% pull(sample_id)

    eset_sub <- kidney_elist[,c(ctrl_ids, drug_ids)]
    sampleNames(eset_sub) <- c(paste("Ctrl", seq(1:length(ctrl_ids))), paste(drug, seq(1:length(drug_ids))))
    design_df <- data.frame(sample = sampleNames(eset_sub),
                            group = c(rep("Ctrl", length(ctrl_ids)),
                                      rep("Drug", length(drug_ids))))
    group <- as.factor(design_df$group)
    boxplot(exprs(eset_sub), col=as.numeric(group), las=2, cex.axis = 0.7, boxwex=0.6, outline=FALSE)

    # Limma
    exprs(eset_sub) <- log2(exprs(eset_sub))
    design <- model.matrix(~ 0 + group, eset_sub)
    colnames(design) <- levels(group)

    # Fit linear model for each gene given a series of arrays
    fit <- limma::lmFit(eset_sub, design)
    # Build comparison and compute the satistics
    cont.matrix <- limma::makeContrasts(Drug-Ctrl, levels=design)
    # Apply the contrast matrix to the linear model
    fit2 <- limma::contrasts.fit(fit, cont.matrix)
    # Generate the statistics
    fit2 <- limma::eBayes(fit2, 0.05)

    de_genes_drug_ctrl <- limma::topTable(fit2, coef=1, adjust="fdr", number=Inf)
    de_table <- cbind(de_genes_drug_ctrl,
                      fData(eset_sub)[rownames(de_genes_drug_ctrl),c("Gene Symbol", "ENTREZ_GENE_ID")])
    write.csv(de_table, file.path(PATH, paste0("data/sigs/drugmatrix/tables/", drug, "_kidney_sigs.csv")), row.names = TRUE)
  }
}

sig_table_paths <- Sys.glob(file.path(PATH, "data/sigs/drugmatrix/tables/*_kidney_sigs.csv"))
drugs <- str_extract(sig_table_paths, "(?<=tables/)(.*?)(?=_kidney_sigs\\.csv)")
sig_tables <- lapply(sig_table_paths, read.csv)
names(sig_tables) <- drugs

# Extracting Signatures
fc_threshold <- 0
p_val_threshold <- 0.05
max_signature <- 100
sigs <- list()
for(drug in drugs) {
  table <- sig_tables[[drug]]
  print(drug)
  if(nrow(table) == 0) {
    next
    paste("skipping", drug, "no DEGS.")
  }
  print(table %>% dim)
  table$logFC_adjpval <- table$logFC * (-log10(table$adj.P.Val))
  up_sig <- table %>%
    dplyr::filter(logFC > fc_threshold, adj.P.Val < p_val_threshold) %>%
    dplyr::arrange(desc(logFC_adjpval)) %>%
    pull(Gene.Symbol.1)
  up_full_sig <- table %>%
    dplyr::arrange(desc(logFC_adjpval)) %>%
    pull(Gene.Symbol.1)
  dn_sig <- table %>%
    dplyr::filter(logFC < -fc_threshold, adj.P.Val < p_val_threshold) %>%
    dplyr::arrange(logFC_adjpval) %>%
    pull(Gene.Symbol.1)
  dn_full_sig <- table %>%
    dplyr::arrange(logFC_adjpval) %>%
    pull(Gene.Symbol.1)

  up_sig_f <- up_sig[up_sig != ""][1:max_signature]
  dn_sig_f <- dn_sig[dn_sig != ""][1:max_signature]
  up_full_sig <- up_full_sig[up_full_sig != ""]
  dn_full_sig <- dn_full_sig[dn_full_sig != ""]

  if (identical(up_sig, character(0)) | identical(dn_sig, character(0))) {
    paste("skipping", drug, "no significant genes.")
    next
  } else if (length(up_sig) < 10 | length(dn_sig) < 10) {
    paste("skipping", drug, "not enough genes.")
    next
  }
  print(length(up_sig))
  print(length(dn_sig))
  sigs[[drug]] <- list(up = up_sig_f, dn = dn_sig_f,
                       up_full = up_full_sig, dn_full = dn_full_sig)
}
kidney_sigs <- sigs

# For upload to github
kidney_sigs <- lapply(kidney_sigs, function(x) x[c("up", "up_full")])
liver_sigs <- lapply(liver_sigs, function(x) x[c("up", "up_full")])
shared_drugs <- intersect(names(kidney_sigs), names(liver_sigs))
kidney_sigs <- kidney_sigs[shared_drugs]
liver_sigs <- liver_sigs[shared_drugs]

drugmatrix.liver <- liver_sigs
drugmatrix.kidney <- kidney_sigs
usethis::use_data(drugmatrix.liver, overwrite= TRUE)
usethis::use_data(drugmatrix.kidney, overwrite = TRUE)
