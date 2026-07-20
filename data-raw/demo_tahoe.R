# Builds demo_tahoe_se / demo_tahoe_sigs / demo_tahoe_true_sigs by
# subsetting real Tahoe A498 pseudobulk expression and the real,
# already-published tahoe.nci_h23/tahoe.a498 signatures.
#
# Narrative: demo_tahoe_sigs are "source" signatures from NCI-H23 (as
# scripts/06_Tahoe/projectCor/01_project_cor.R uses NCI-H23 as its
# source: `source_sig <- lapply(tahoe.nci_h23, function(x) x$up)`),
# demo_tahoe_se is real target-context (A498) expression data, and
# demo_tahoe_true_sigs is the real target-context ground truth (A498
# signatures).
#
# No data leakage: demo_tahoe_se's perturbation samples are a set of
# A498 drugs DISJOINT from the ~23 drugs actually evaluated in
# demo_tahoe_sigs/demo_tahoe_true_sigs (see data-raw/demo_sciplex.R for
# the full rationale).
#
# Prerequisites (not part of the repo):
#   * Raw merged Tahoe pseudobulk Seurat object across all cell lines and
#     drugs, as produced by the Tahoe processing pipeline (see
#     scripts/06_Tahoe/, data-raw/tahoe_pseudobulk.py), available locally
#     at <repo_root>/pb_data/merged_pseudobulk_filtered.rds. Source:
#     https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M
#   * tahoe.nci_h23 / tahoe.a498 signature objects, fetched via
#     sigrecon::get_dataset() (see R/data.R and #77).

library(Seurat)
library(SummarizedExperiment)

devtools::load_all(".")

pb_path <- file.path("..", "pb_data", "merged_pseudobulk_filtered.rds")
stopifnot(
  "Raw Tahoe merged pseudobulk not found -- see prerequisites above" =
    file.exists(pb_path)
)
merged <- readRDS(pb_path)

tahoe.nci_h23 <- get_dataset("tahoe.nci_h23")
tahoe.a498 <- get_dataset("tahoe.a498")

# Drugs with both a source-context (NCI-H23) and true target-context
# (A498) signature available, and present in the A498 pseudobulk data.
common_drugs <- sort(intersect(names(tahoe.nci_h23), names(tahoe.a498)))
a498_drugs <- unique(merged$drug_name[merged$cell_name == "A498"])
common_drugs <- intersect(common_drugs, a498_drugs)
stopifnot(length(common_drugs) > 0)

# Sample down to ~23 evaluated drugs (matching the SciPlex demo's scale)
# out of 107 available, for a reproducible small demo.
set.seed(42)
evaluated_drugs <- sort(sample(common_drugs, min(23, length(common_drugs))))

# Disjoint background: a *different* set of ~23 A498 drugs (excluding
# DMSO_TF and the evaluated drugs) to give demo_tahoe_se real
# perturbation-induced expression variance without containing any sample
# for the drugs being evaluated.
set.seed(43)
a498_drugs <- unique(merged$drug_name[merged$cell_name == "A498"])
background_pool <- setdiff(a498_drugs, c("DMSO_TF", evaluated_drugs))
background_drugs <- sort(sample(background_pool, min(23, length(background_pool))))

# Gene universe: union of the up-genes needed for every source/true
# signature we're keeping (see data-raw/demo_sciplex.R for why this is
# used instead of variance-based gene selection).
genes_needed <- unique(unlist(c(
  lapply(evaluated_drugs, function(d) tahoe.nci_h23[[d]]$up),
  lapply(evaluated_drugs, function(d) tahoe.a498[[d]]$up)
)))
genes_needed <- intersect(genes_needed, rownames(merged))

# Subset to A498, DMSO_TF (control) + the disjoint background drugs'
# samples only -- none of evaluated_drugs' samples are included.
keep_samples <- merged$cell_name == "A498" &
  merged$drug_name %in% c("DMSO_TF", background_drugs)
a498_sub <- subset(merged, cells = colnames(merged)[keep_samples])

counts <- as.matrix(GetAssayData(a498_sub, layer = "counts"))[genes_needed, ]

# Drop zero-variance genes and log2-CPM normalize (see demo_sciplex.R).
keep <- apply(counts, 1, var) > 0
counts <- counts[keep, ]
cpm <- t(t(counts) / colSums(counts)) * 1e6
logcounts <- log2(cpm + 1)

demo_tahoe_se <- SummarizedExperiment(
  assays = list(logcounts = logcounts, counts = counts),
  colData = DataFrame(a498_sub@meta.data)
)

final_genes <- rownames(demo_tahoe_se)
stopifnot(
  "demo_tahoe_se must not contain any evaluated drug's samples" =
    length(intersect(unique(demo_tahoe_se$drug_name), evaluated_drugs)) == 0
)

demo_tahoe_sigs <- setNames(
  lapply(evaluated_drugs, function(d) intersect(tahoe.nci_h23[[d]]$up, final_genes)),
  evaluated_drugs
)

demo_tahoe_true_sigs <- setNames(
  lapply(evaluated_drugs, function(d) {
    list(
      up = intersect(tahoe.a498[[d]]$up, final_genes),
      up_full = intersect(tahoe.a498[[d]]$up_full, final_genes)
    )
  }),
  evaluated_drugs
)

usethis::use_data(demo_tahoe_se, overwrite = TRUE)
usethis::use_data(demo_tahoe_sigs, overwrite = TRUE)
usethis::use_data(demo_tahoe_true_sigs, overwrite = TRUE)
