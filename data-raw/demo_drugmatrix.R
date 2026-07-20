# Builds demo_drugmatrix_se / demo_drugmatrix_sigs /
# demo_drugmatrix_true_sigs by subsetting real DrugMatrix liver
# expression and the real, already-published
# drugmatrix.kidney/drugmatrix.liver signatures.
#
# Narrative: demo_drugmatrix_sigs are "source" signatures from kidney,
# demo_drugmatrix_se is real target-context (liver) expression data, and
# demo_drugmatrix_true_sigs is the real target-context ground truth
# (liver signatures). scripts/03_DrugMatrix/projectCor/01_project_cor.R
# computes both directions (kidney->liver and liver->kidney) -- this
# picks one arbitrarily.
#
# No data leakage: demo_drugmatrix_se's perturbation samples are a set of
# liver compounds DISJOINT from the ~23 evaluated in
# demo_drugmatrix_sigs/demo_drugmatrix_true_sigs (see
# data-raw/demo_sciplex.R for the full rationale). The disjoint set is
# drawn from the full liver compound pool (not just the 39 with
# kidney+liver signatures), since background samples don't need
# precomputed signatures.
#
# Prerequisites (not part of the repo):
#   * Raw DrugMatrix liver ExpressionSet, available locally at
#     <repo_root>/pb_data/liver.rds. Source:
#     https://ntp.niehs.nih.gov/data/drugmatrix
#   * drugmatrix.kidney / drugmatrix.liver signature objects, fetched via
#     sigrecon::get_dataset() (see R/data.R and #77).

library(Biobase)
library(SummarizedExperiment)

devtools::load_all(".")

pb_path <- file.path("..", "pb_data", "liver.rds")
stopifnot(
  "Raw DrugMatrix liver ExpressionSet not found -- see prerequisites above" =
    file.exists(pb_path)
)
liver_eset <- readRDS(pb_path)
featureNames(liver_eset) <- make.unique(fData(liver_eset)$`Gene Symbol`)

drugmatrix.kidney <- get_dataset("drugmatrix.kidney")
drugmatrix.liver <- get_dataset("drugmatrix.liver")

# Drugs with both a source-context (kidney) and true target-context
# (liver) signature available, and present in the liver ExpressionSet.
common_drugs <- sort(intersect(names(drugmatrix.kidney), names(drugmatrix.liver)))
common_drugs <- intersect(common_drugs, unique(na.omit(pData(liver_eset)$`compound:ch1`)))
stopifnot(length(common_drugs) > 0)

# Sample down to ~23 evaluated drugs (matching the SciPlex/Tahoe demo
# scale).
set.seed(42)
evaluated_drugs <- sort(sample(common_drugs, min(23, length(common_drugs))))

# Disjoint background: a *different* set of ~23 liver compounds
# (excluding the evaluated drugs), drawn from the full liver compound
# pool, to give demo_drugmatrix_se real perturbation-induced expression
# variance without containing any sample for the drugs being evaluated.
set.seed(43)
background_pool <- setdiff(unique(na.omit(pData(liver_eset)$`compound:ch1`)), evaluated_drugs)
background_drugs <- sort(sample(background_pool, min(23, length(background_pool))))

# Gene universe: union of the up-genes needed for every source/true
# signature we're keeping (see data-raw/demo_sciplex.R for why this is
# used instead of variance-based gene selection).
genes_needed <- unique(unlist(c(
  lapply(evaluated_drugs, function(d) drugmatrix.kidney[[d]]$up),
  lapply(evaluated_drugs, function(d) drugmatrix.liver[[d]]$up)
)))
genes_needed <- intersect(genes_needed, featureNames(liver_eset))

# Control samples (vehicle, 0 mg/kg) + max-dose/max-time samples for the
# disjoint background drugs, matching the filter_sample_metadata()
# convention in scripts/03_DrugMatrix/01_DrugMatrix_sigs.R (the
# signatures were derived from max-dose/max-time comparisons). Cap
# controls (279 available) to keep the demo small, matching the other
# demo datasets' scale.
pd <- pData(liver_eset)
set.seed(42)
ctrl_ids <- rownames(pd)[pd$`dose:ch1` == "0 mg/kg"]
ctrl_ids <- sort(sample(ctrl_ids, min(15, length(ctrl_ids))))

drug_pd <- pd[pd$`compound:ch1` %in% background_drugs & !is.na(pd$`compound:ch1`), ]
drug_pd$new_dose <- as.numeric(gsub("[^0-9.]", "", drug_pd$`dose:ch1`))
drug_pd$new_time <- as.numeric(gsub("[^0-9.]", "", drug_pd$`time:ch1`))
drug_ids <- unlist(lapply(background_drugs, function(d) {
  sub_pd <- drug_pd[drug_pd$`compound:ch1` == d, ]
  sub_pd <- sub_pd[sub_pd$new_dose == max(sub_pd$new_dose, na.rm = TRUE), ]
  sub_pd <- sub_pd[sub_pd$new_time == max(sub_pd$new_time, na.rm = TRUE), ]
  rownames(sub_pd)
}))

liver_sub <- liver_eset[genes_needed, c(ctrl_ids, drug_ids)]

# Assay data is already log-transformed (microarray), unlike the raw
# counts in the SciPlex/Tahoe pseudobulk objects -- no additional
# normalization needed.
logcounts <- exprs(liver_sub)
pdata_sub <- pData(liver_sub)[, c("compound:ch1", "dose:ch1", "tissue:ch1", "time:ch1", "vehicle:ch1")]
colnames(pdata_sub) <- c("compound", "dose", "tissue", "time", "vehicle")
pdata_sub$compound <- ifelse(is.na(pdata_sub$compound), "Control", pdata_sub$compound)

demo_drugmatrix_se <- SummarizedExperiment(
  assays = list(logcounts = logcounts),
  colData = DataFrame(pdata_sub)
)

final_genes <- rownames(demo_drugmatrix_se)
stopifnot(
  "demo_drugmatrix_se must not contain any evaluated drug's samples" =
    length(intersect(unique(demo_drugmatrix_se$compound), evaluated_drugs)) == 0
)

demo_drugmatrix_sigs <- setNames(
  lapply(evaluated_drugs, function(d) intersect(drugmatrix.kidney[[d]]$up, final_genes)),
  evaluated_drugs
)

demo_drugmatrix_true_sigs <- setNames(
  lapply(evaluated_drugs, function(d) {
    list(
      up = intersect(drugmatrix.liver[[d]]$up, final_genes),
      up_full = intersect(drugmatrix.liver[[d]]$up_full, final_genes)
    )
  }),
  evaluated_drugs
)

usethis::use_data(demo_drugmatrix_se, overwrite = TRUE)
usethis::use_data(demo_drugmatrix_sigs, overwrite = TRUE)
usethis::use_data(demo_drugmatrix_true_sigs, overwrite = TRUE)
