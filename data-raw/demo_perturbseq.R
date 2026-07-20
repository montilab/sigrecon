# Builds demo_perturbseq_se / demo_perturbseq_sigs /
# demo_perturbseq_true_sigs by subsetting real Perturb-seq RPE1
# pseudobulk expression and the real, already-published
# perturbseq.k562/perturbseq.rpe1 signatures.
#
# Narrative: demo_perturbseq_sigs are "source" signatures from K562,
# demo_perturbseq_se is real target-context (RPE1) expression data, and
# demo_perturbseq_true_sigs is the real target-context ground truth
# (RPE1 signatures) -- matching the K562/RPE1 recontextualization setup
# in scripts/02_PerturbSeq/.
#
# No data leakage: demo_perturbseq_se's perturbation samples are a set
# of RPE1 knockdowns DISJOINT from the ~23 evaluated in
# demo_perturbseq_sigs/demo_perturbseq_true_sigs (see
# data-raw/demo_sciplex.R for the full rationale).
#
# The raw counts CSVs are ~2.9GB each, so this reads only the columns
# (samples) actually needed via data.table::fread(select = ...) rather
# than loading the full file, unlike scripts/02_PerturbSeq/*/01_*.R which
# load the whole thing (fine for an HPC job, not for a demo-data build).
#
# Prerequisites (not part of the repo):
#   * Raw Perturb-seq RPE1 pseudobulk counts + metadata CSVs, as produced
#     by the Perturb-seq processing pipeline (see scripts/02_PerturbSeq/),
#     available locally at <repo_root>/pb_data/rpe1_processed_pb.csv and
#     rpe1_processed_pb_metadata.csv. Source:
#     https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387
#   * perturbseq.k562 / perturbseq.rpe1 signature objects, fetched via
#     sigrecon::get_dataset() (see R/data.R and #77).

library(data.table)
library(SummarizedExperiment)

devtools::load_all(".")

counts_path <- file.path("..", "pb_data", "rpe1_processed_pb.csv")
meta_path <- file.path("..", "pb_data", "rpe1_processed_pb_metadata.csv")
stopifnot(
  "Raw Perturb-seq RPE1 pseudobulk CSVs not found -- see prerequisites above" =
    file.exists(counts_path) && file.exists(meta_path)
)

rpe1_meta <- read.csv(meta_path, row.names = 1)

perturbseq.k562 <- get_dataset("perturbseq.k562")
perturbseq.rpe1 <- get_dataset("perturbseq.rpe1")

# Knockdowns with both a source-context (K562) and true target-context
# (RPE1) signature available, and present in the RPE1 metadata.
common_genes <- sort(intersect(names(perturbseq.k562), names(perturbseq.rpe1)))
common_genes <- intersect(common_genes, unique(rpe1_meta$gene))
stopifnot(length(common_genes) > 0)

# Sample down to ~23 evaluated knockdowns (matching the
# SciPlex/Tahoe/DrugMatrix demo scale) out of 1297 available.
set.seed(42)
evaluated_genes <- sort(sample(common_genes, min(23, length(common_genes))))

# Disjoint background: a *different* set of ~23 RPE1 knockdowns
# (excluding non-targeting and the evaluated genes) to give
# demo_perturbseq_se real perturbation-induced expression variance
# without containing any sample for the knockdowns being evaluated.
set.seed(43)
background_pool <- setdiff(unique(rpe1_meta$gene), c("non-targeting", evaluated_genes))
background_genes <- sort(sample(background_pool, min(23, length(background_pool))))

# Cap replicates per knockdown (some genes have 50+) and controls, so the
# demo stays small.
sample_ids <- unlist(lapply(background_genes, function(g) {
  ids <- rownames(rpe1_meta)[rpe1_meta$gene == g]
  ids[seq_len(min(3, length(ids)))]
}))
ctrl_ids <- rownames(rpe1_meta)[rpe1_meta$gene == "non-targeting"]
ctrl_ids <- ctrl_ids[seq_len(min(15, length(ctrl_ids)))]

all_ids <- c(ctrl_ids, sample_ids)

rpe1_dt <- fread(counts_path, select = c("gene_name", all_ids))
counts <- as.matrix(rpe1_dt[, -1, with = FALSE])
rownames(counts) <- rpe1_dt$gene_name

# Gene universe: union of the up-genes needed for every source/true
# signature we're keeping (see data-raw/demo_sciplex.R for why this is
# used instead of variance-based gene selection).
genes_needed <- unique(unlist(c(
  lapply(evaluated_genes, function(g) perturbseq.k562[[g]]$up),
  lapply(evaluated_genes, function(g) perturbseq.rpe1[[g]]$up)
)))
genes_needed <- intersect(genes_needed, rownames(counts))
counts <- counts[genes_needed, ]

# Drop zero-variance genes and log2-CPM normalize (see demo_sciplex.R).
keep <- apply(counts, 1, var) > 0
counts <- counts[keep, ]
cpm <- t(t(counts) / colSums(counts)) * 1e6
logcounts <- log2(cpm + 1)

col_data <- rpe1_meta[all_ids, , drop = FALSE]
col_data$is_control <- ifelse(col_data$gene == "non-targeting", "NTC", "Perturbed")

demo_perturbseq_se <- SummarizedExperiment(
  assays = list(logcounts = logcounts, counts = counts),
  colData = DataFrame(col_data)
)

final_genes <- rownames(demo_perturbseq_se)
stopifnot(
  "demo_perturbseq_se must not contain any evaluated knockdown's samples" =
    length(intersect(unique(demo_perturbseq_se$gene), evaluated_genes)) == 0
)

demo_perturbseq_sigs <- setNames(
  lapply(evaluated_genes, function(g) intersect(perturbseq.k562[[g]]$up, final_genes)),
  evaluated_genes
)

demo_perturbseq_true_sigs <- setNames(
  lapply(evaluated_genes, function(g) {
    list(
      up = intersect(perturbseq.rpe1[[g]]$up, final_genes),
      up_full = intersect(perturbseq.rpe1[[g]]$up_full, final_genes)
    )
  }),
  evaluated_genes
)

usethis::use_data(demo_perturbseq_se, overwrite = TRUE)
usethis::use_data(demo_perturbseq_sigs, overwrite = TRUE)
usethis::use_data(demo_perturbseq_true_sigs, overwrite = TRUE)
