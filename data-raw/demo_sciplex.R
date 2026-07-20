# Builds the small, real-data quick-start demo objects (demo_sciplex_se,
# demo_sciplex_sigs, demo_sciplex_true_sigs) by subsetting actual SciPlex
# A549 pseudobulk expression and the real, already-published
# sciplex.a549/sciplex.k562 signatures.
#
# Narrative: demo_sciplex_sigs are "source" signatures as if defined in
# another biological context (SciPlex K562), demo_sciplex_se is real
# target-context (SciPlex A549) expression data, and
# demo_sciplex_true_sigs is the real target-context ground truth (SciPlex
# A549 signatures) -- exactly the recontextualization setup used
# throughout scripts/04_SciPlex/.
#
# No data leakage: demo_sciplex_se's perturbation samples are a set of
# A549 drugs DISJOINT from the 23 drugs actually evaluated in
# demo_sciplex_sigs/demo_sciplex_true_sigs. This matches the real
# pipeline's ctrl_se convention (scripts/04_SciPlex/projectCor/01_project_cor.R
# builds ctrl_se from control-only samples so the recontextualization
# input never contains the perturbation being predicted); here we also
# add a disjoint *other*-drug background (rather than control-only)
# because control-only leaves too few samples (2 Vehicle replicates) to
# build a meaningful WGCNA network for the demo.
#
# Prerequisites (not part of the repo -- these are the raw/derived
# datasets underlying the sciplex.* entries in R/data.R):
#   * Raw SciPlex A549 pseudobulk Seurat object, as produced by the
#     SciPlex processing pipeline (see scripts/04_SciPlex/), available
#     locally at <repo_root>/pb_data/a549_filtered_pb.rds. Source:
#     https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398
#   * sciplex.a549 / sciplex.k562 signature objects, fetched via
#     sigrecon::get_dataset() (see R/data.R and #77).

library(Seurat)
library(SummarizedExperiment)

devtools::load_all(".")

pb_path <- file.path("..", "pb_data", "a549_filtered_pb.rds")
stopifnot(
  "Raw SciPlex A549 pseudobulk not found -- see prerequisites above" =
    file.exists(pb_path)
)
a549 <- readRDS(pb_path)

sciplex.a549 <- get_dataset("sciplex.a549")
sciplex.k562 <- get_dataset("sciplex.k562")

# Drugs with both a source-context (K562) and true target-context (A549)
# signature available, and present in the A549 pseudobulk object. These
# are the drugs actually evaluated -- demo_sciplex_se must not contain
# their expression samples (see no-leakage note above).
common_drugs <- sort(intersect(names(sciplex.a549), names(sciplex.k562)))
common_drugs <- intersect(common_drugs, unique(a549$product_name))
stopifnot(length(common_drugs) > 0)

# Sample down to ~23 evaluated drugs for a reproducible small demo.
set.seed(42)
evaluated_drugs <- sort(sample(common_drugs, min(23, length(common_drugs))))

# Disjoint background: a *different* set of ~23 A549 drugs (excluding
# Vehicle and the evaluated drugs) to give demo_sciplex_se real
# perturbation-induced expression variance without containing any sample
# for the drugs being evaluated.
set.seed(43)
background_pool <- setdiff(unique(a549$product_name), c("Vehicle", evaluated_drugs))
background_drugs <- sort(sample(background_pool, min(23, length(background_pool))))

# Gene universe: union of the up-genes needed for every source/true
# signature we're keeping, so projectCor()/network_sig() have every
# signature gene available in demo_sciplex_se. This deliberately diverges
# from the real pipeline's variance-based gene selection
# (rank.var.eset() top-10000) -- picking top-variable genes here would
# risk missing genes the demo signatures actually need, silently
# weakening the demo. The tradeoff: demo_sciplex_se's correlation/network
# structure reflects genes selected for DE-relevance to the evaluated
# drugs, not a generic most-variable-genes panel.
genes_needed <- unique(unlist(c(
  lapply(evaluated_drugs, function(d) sciplex.k562[[d]]$up),
  lapply(evaluated_drugs, function(d) sciplex.a549[[d]]$up)
)))
genes_needed <- intersect(genes_needed, rownames(a549))

# Subset to Vehicle (control) + the disjoint background drugs' samples
# only -- none of evaluated_drugs' samples are included.
keep_samples <- a549$product_name %in% c("Vehicle", background_drugs)
a549_sub <- subset(a549, cells = colnames(a549)[keep_samples])

counts <- as.matrix(GetAssayData(a549_sub, layer = "counts"))[genes_needed, ]

# Drop zero-variance genes (breaks correlation-based network construction)
# and log2-CPM normalize, matching the log-expression convention used for
# network learning elsewhere in the package (e.g. the GTEx DESeq2_log
# assay in scripts/07_GTEX/network_prop/01_gtex_net.R).
keep <- apply(counts, 1, var) > 0
counts <- counts[keep, ]
cpm <- t(t(counts) / colSums(counts)) * 1e6
logcounts <- log2(cpm + 1)

demo_sciplex_se <- SummarizedExperiment(
  assays = list(logcounts = logcounts, counts = counts),
  colData = DataFrame(a549_sub@meta.data)
)

final_genes <- rownames(demo_sciplex_se)
stopifnot(
  "demo_sciplex_se must not contain any evaluated drug's samples" =
    length(intersect(unique(demo_sciplex_se$product_name), evaluated_drugs)) == 0
)

demo_sciplex_sigs <- setNames(
  lapply(evaluated_drugs, function(d) intersect(sciplex.k562[[d]]$up, final_genes)),
  evaluated_drugs
)

demo_sciplex_true_sigs <- setNames(
  lapply(evaluated_drugs, function(d) {
    list(
      up = intersect(sciplex.a549[[d]]$up, final_genes),
      up_full = intersect(sciplex.a549[[d]]$up_full, final_genes)
    )
  }),
  evaluated_drugs
)

usethis::use_data(demo_sciplex_se, overwrite = TRUE)
usethis::use_data(demo_sciplex_sigs, overwrite = TRUE)
usethis::use_data(demo_sciplex_true_sigs, overwrite = TRUE)
