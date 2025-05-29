## code to prepare `sciplex` dataset goes here
library(tidyverse)
library(data.table)
library(Seurat)
library(doParallel)
options(Seurat.object.assay.version = 'v5')

DATA_PATH <- ""

# Loading Sci-Plex Data
# Download from here https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398
a549 <- readRDS(file.path(DATA_PATH, "a549.rds"))
k562 <- readRDS(file.path(DATA_PATH, "k562.rds"))
mcf7 <- readRDS(file.path(DATA_PATH, "mcf7.rds"))

# Getting Groups
stopifnot(setdiff(unique(a549$product_name), unique(mcf7$product_name)) == 0)
stopifnot(setdiff(unique(a549$product_name), unique(k562$product_name)) == 0)
products <- a549$product_name %>% unique %>% as.vector
products <- products[products != "Vehicle"]
a549_ids <- list()
k562_ids <- list()
mcf7_ids <- list()

for (product in products) {
  a549_ids[[product]] <- a549[[]] %>% rownames_to_column(var= "id") %>% filter(dose == 10000, product_name == product) %>% pull(id)
  k562_ids[[product]] <- k562[[]] %>% rownames_to_column(var= "id") %>% filter(dose == 10000, product_name == product) %>% pull(id)
  mcf7_ids[[product]] <- mcf7[[]] %>% rownames_to_column(var= "id") %>% filter(dose == 10000, product_name == product) %>% pull(id)
}

print("Starting a549")
a549_sigs <- foreach(product = products, .combine=dplyr::bind_rows) %dopar% {
  drug_cell_ids <- a549_ids[[product]]
  if(length(drug_cell_ids) <= 30) {
    print(paste0("Skipping ", product, ". Fewer than 30 cells."))
    NULL
  } else {
    control_cell_ids <- a549[[]] %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name == "Vehicle") %>% pull(id)
    stopifnot(drug_cell_ids %in% colnames(a549))
    stopifnot(control_cell_ids %in% colnames(a549))
    sig_df <- Seurat::FindMarkers(a549,
                                  ident.1 = drug_cell_ids,
                                  ident.2 = control_cell_ids,
                                  test.use = "MAST",
                                  logfc.threshold = 0,
                                  only.pos = FALSE)
    sig_df$product <- product
    sig_df
  }
}

print("Starting K562")
k562_sigs <- foreach(product = products, .combine=dplyr::bind_rows) %dopar% {
  drug_cell_ids <- k562_ids[[product]]
  if(length(drug_cell_ids) <= 30) {
    print(paste0("Skipping ", product, ". Fewer than 30 cells."))
    NULL
  } else {
    control_cell_ids <- k562[[]] %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name == "Vehicle") %>% pull(id)
    stopifnot(drug_cell_ids %in% colnames(k562))
    stopifnot(control_cell_ids %in% colnames(k562))
    sig_df <- Seurat::FindMarkers(k562,
                                  ident.1 = drug_cell_ids,
                                  ident.2 = control_cell_ids,
                                  test.use = "MAST",
                                  logfc.threshold = 0,
                                  only.pos = FALSE)
    sig_df$product <- product
    sig_df
  }
}

print("Starting MCF7")
mcf7_sigs <- foreach(product = products, .combine=dplyr::bind_rows) %dopar% {
  drug_cell_ids <- mcf7_ids[[product]]
  if(length(drug_cell_ids) <= 30) {
    print(paste0("Skipping ", product, ". Fewer than 30 cells."))
    NULL
  } else {
    control_cell_ids <- mcf7[[]] %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name == "Vehicle") %>% pull(id)
    stopifnot(drug_cell_ids %in% colnames(mcf7))
    stopifnot(control_cell_ids %in% colnames(mcf7))
    sig_df <- Seurat::FindMarkers(mcf7,
                                  ident.1 = drug_cell_ids,
                                  ident.2 = control_cell_ids,
                                  test.use = "MAST",
                                  logfc.threshold = 0,
                                  only.pos = FALSE)
    sig_df$product <- product
    sig_df
  }
}

# Filtering signatures
# stopifnot(all.equal(unique(mcf7_sigs$product), unique(k562_sigs$product)))
# stopifnot(all.equal(unique(a549_sigs$product), unique(k562_sigs$product)))

sig_filter_fn <- function(diff_table, products, alpha=0.05, limit=100) {
  results <- list()
  for(drug in products) {
    print(drug)
    diff_table$logFC_adjpval <- diff_table$avg_log2FC * (-log10(diff_table$p_val_adj))
    full_sig <- diff_table %>%
      dplyr::filter(product==drug) %>%
      dplyr::arrange(desc(logFC_adjpval)) %>%
      rownames
    sig_sig <- diff_table %>%
      dplyr::filter(product==drug) %>%
      dplyr::filter(avg_log2FC > 0) %>%
      dplyr::filter(p_val_adj <= alpha) %>%
      dplyr::arrange(desc(logFC_adjpval)) %>%
      dplyr::slice(1:limit) %>%
      rownames
    results[[drug]][["up"]] <- sig_sig
    results[[drug]][["up_full"]] <- full_sig
  }
  return(results)
}

# Filtering for significant signatures
mcf7_sig_list <- sig_filter_fn(mcf7_sigs, unique(mcf7_sigs$product))
k562_sig_list <- sig_filter_fn(k562_sigs, unique(k562_sigs$product))
a549_sig_list <- sig_filter_fn(a549_sigs, unique(a549_sigs$product))
# Removing extra string padding
mcf7_sig_list <- lapply(mcf7_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
k562_sig_list <- lapply(k562_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
a549_sig_list <- lapply(a549_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
# Removing NAs
mcf7_sig_list <- lapply(mcf7_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))
k562_sig_list <- lapply(k562_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))
a549_sig_list <- lapply(a549_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))

lapply(mcf7_sig_list, function(x) length(x$up)) %>% unlist %>% unname %>% fivenum
lapply(k562_sig_list, function(x) length(x$up)) %>% unlist %>% unname %>% fivenum
lapply(a549_sig_list, function(x) length(x$up)) %>% unlist %>% unname %>% fivenum
lapply(mcf7_sig_list, function(x) length(x$up_full)) %>% unlist %>% unname %>% fivenum
lapply(k562_sig_list, function(x) length(x$up_full)) %>% unlist %>% unname %>% fivenum
lapply(a549_sig_list, function(x) length(x$up_full)) %>% unlist %>% unname %>% fivenum
# Filtering to perturbations with a reasonable number of DEGS (arbitrarily set as 5)
mcf7_sig_list <- mcf7_sig_list[lapply(mcf7_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
k562_sig_list <- k562_sig_list[lapply(k562_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
a549_sig_list <- a549_sig_list[lapply(a549_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]

# Filtering to shared products
shared_products <- purrr::reduce(list(names(mcf7_sig_list), names(k562_sig_list), names(a549_sig_list)), intersect)
sciplex.mcf7 <- mcf7_sig_list[shared_products]
sciplex.k562 <- k562_sig_list[shared_products]
sciplex.a549 <- a549_sig_list[shared_products]

usethis::use_data(sciplex.mcf7, overwrite = TRUE)
usethis::use_data(sciplex.k562, overwrite = TRUE)
usethis::use_data(sciplex.a549, overwrite = TRUE)
