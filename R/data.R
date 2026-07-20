#' Registry of on-demand perturbational signature datasets
#'
#' @description
#' Internal registry describing the datasets available via get_dataset().
#' Each dataset is a nested list of genesets (`up` = top 100 DEGs,
#' `up_full` = full ranked gene list), hosted as a release asset on
#' `montilab/sigrecon` and fetched on first use.
#' @noRd
.sigrecon_dataset_registry <- data.frame(
  name = c(
    "drugmatrix.kidney", "drugmatrix.liver",
    "sciplex.a549", "sciplex.k562", "sciplex.mcf7",
    "perturbseq.k562", "perturbseq.rpe1",
    "neurips2023.b", "neurips2023.nk", "neurips2023.t", "neurips2023.myeloid",
    "gtex.blood", "gtex.brain.hippo",
    "tahoe.a498", "tahoe.hct15", "tahoe.hec_1_a", "tahoe.lovo",
    "tahoe.miapaca_2", "tahoe.nci_h23", "tahoe.panc03.27",
    "tahoe.snu_1", "tahoe.snu_423", "tahoe.sw48"
  ),
  description = c(
    "DEGs between drug (max dose) and control samples in DrugMatrix kidney data. 39 drugs, mouse genes.",
    "DEGs between drug (max dose) and control samples in DrugMatrix liver data. 39 drugs, mouse genes.",
    "DEGs between drug (max dose) and control samples in SciPlex A549 cell line data. 23 drugs, EnsemblIDs.",
    "DEGs between drug (max dose) and control samples in SciPlex K562 cell line data. 23 drugs, EnsemblIDs.",
    "DEGs between drug (max dose) and control samples in SciPlex MCF7 cell line data. 23 drugs, EnsemblIDs.",
    "DEGs between drug (max dose) and control samples in Perturb-seq K562 cell line data. 1297 CRISPRi knockdowns, HGNC symbols.",
    "DEGs between drug (max dose) and control samples in Perturb-seq RPE1 cell line data. 1297 CRISPRi knockdowns, HGNC symbols.",
    "DEGs between drug and control (DMSO) samples in NeurIPS 2023 human PBMC B cells. 135 shared drugs, HGNC symbols.",
    "DEGs between drug and control (DMSO) samples in NeurIPS 2023 human PBMC NK cells. 135 shared drugs, HGNC symbols.",
    "DEGs between drug and control (DMSO) samples in NeurIPS 2023 human PBMC T cells. 135 shared drugs, HGNC symbols.",
    "DEGs between drug and control (DMSO) samples in NeurIPS 2023 human PBMC myeloid cells. 135 shared drugs, HGNC symbols.",
    "DEGs between old and young whole blood bulk samples (age as numerical covariate) from GTEx v10.",
    "DEGs between old and young brain hippocampus bulk samples (age as numerical covariate) from GTEx v10.",
    "DEGs between drug and control samples in Tahoe A498 cell line data. 108 drugs, HGNC/Ensembl mix.",
    "DEGs between drug and control samples in Tahoe HCT15 cell line data. 105 drugs, HGNC/Ensembl mix.",
    "DEGs between drug and control samples in Tahoe HEC-1-A cell line data. 109 drugs, HGNC/Ensembl mix.",
    "DEGs between drug and control samples in Tahoe LOVO cell line data. 109 drugs, HGNC/Ensembl mix.",
    "DEGs between drug and control samples in Tahoe MIAPACA-2 cell line data. 109 drugs, HGNC/Ensembl mix.",
    "DEGs between drug and control samples in Tahoe NCI-H23 cell line data. 108 drugs, HGNC/Ensembl mix.",
    "DEGs between drug and control samples in Tahoe PANC03.27 cell line data. 106 drugs, HGNC/Ensembl mix.",
    "DEGs between drug and control samples in Tahoe SNU-1 cell line data. 105 drugs, HGNC/Ensembl mix.",
    "DEGs between drug and control samples in Tahoe SNU-423 cell line data. 108 drugs, HGNC/Ensembl mix.",
    "DEGs between drug and control samples in Tahoe SW48 cell line data. 104 drugs, HGNC/Ensembl mix."
  ),
  source = c(
    "https://ntp.niehs.nih.gov/data/drugmatrix",
    "https://ntp.niehs.nih.gov/data/drugmatrix",
    "https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398",
    "https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398",
    "https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398",
    "https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387",
    "https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945",
    "https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression",
    "https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression",
    "https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M",
    "https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M",
    "https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M",
    "https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M",
    "https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M",
    "https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M",
    "https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M",
    "https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M",
    "https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M",
    "https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M"
  ),
  stringsAsFactors = FALSE
)

.sigrecon_release_base_url <- "https://github.com/montilab/sigrecon/releases/download/data-v1/"

#' List available on-demand datasets
#'
#' @description
#' Returns a table of the perturbational signature datasets available via
#' [get_dataset()], with a short description and original data source for each.
#'
#' @return A data.frame with columns `name`, `description`, `source`.
#'
#' @examples
#' list_datasets()
#'
#' @export
list_datasets <- function() {
  .sigrecon_dataset_registry
}

#' Fetch a bundled perturbational signature dataset
#'
#' @description
#' Downloads (and locally caches, via `BiocFileCache`) one of the datasets
#' listed in [list_datasets()]. Datasets are hosted as release assets on
#' `montilab/sigrecon` rather than bundled with the package, so the first
#' call to `get_dataset()` for a given name requires network access;
#' subsequent calls reuse the local cache.
#'
#' @param name Dataset name. See [list_datasets()] for available names.
#' @param force Logical. If `TRUE`, re-download even if a cached copy exists.
#'   Default is `FALSE`.
#'
#' @return A named list of genesets. Each element has an `up` component
#'   (top 100 DEGs) and an `up_full` component (the full ranked gene list).
#'
#' @examples
#' \dontrun{
#' tahoe_nci_h23 <- get_dataset("tahoe.nci_h23")
#' }
#'
#' @export
get_dataset <- function(name, force = FALSE) {
  if (!is.character(name) || length(name) != 1) {
    stop("'name' must be a single character string. See list_datasets() for available names.")
  }

  if (!name %in% .sigrecon_dataset_registry$name) {
    stop(sprintf(
      "Unknown dataset '%s'. See list_datasets() for available names.",
      name
    ))
  }

  if (!requireNamespace("BiocFileCache", quietly = TRUE)) {
    stop("The 'BiocFileCache' package is required for get_dataset(). Install it with BiocManager::install('BiocFileCache').")
  }

  url <- paste0(.sigrecon_release_base_url, name, ".rda")
  bfc <- BiocFileCache::BiocFileCache(ask = FALSE)

  if (force) {
    hits <- BiocFileCache::bfcquery(bfc, url, field = "rname", exact = TRUE)
    if (nrow(hits) > 0) {
      BiocFileCache::bfcremove(bfc, hits$rid)
    }
  }

  path <- BiocFileCache::bfcrpath(bfc, url)

  e <- new.env()
  loaded <- load(path, envir = e)
  if (length(loaded) != 1) {
    stop(sprintf(
      "Expected a single object in the downloaded file for '%s', found: %s",
      name, paste(loaded, collapse = ", ")
    ))
  }

  get(loaded, envir = e)
}

#' Demo datasets
#'
#' @description
#' `sigrecon` bundles four small, real-data recontextualization demos
#' (one per major dataset covered by `scripts/`), each a
#' `demo_<dataset>_se` / `demo_<dataset>_sigs` / `demo_<dataset>_true_sigs`
#' triple with the same shape: `demo_<dataset>_se` is target-context
#' expression, `demo_<dataset>_sigs` are real source-context signatures
#' (as if defined in a different biological context), and
#' `demo_<dataset>_true_sigs` is the real target-context ground truth,
#' for benchmarking recontextualization quality with `sig_eval_table()`.
#' Gene sets are restricted to the union of genes needed by the kept
#' signatures (not variance-based selection -- see
#' `data-raw/demo_sciplex.R` for why), and perturbations are sampled down
#' to ~23 for speed. All are fast, self-contained examples for
#' `projectCor()`, `wgcna.adj()`/`network_sig()`, and `sig_eval_table()`
#' -- see the package README.
#'
#' No data leakage: each `demo_<dataset>_se`'s perturbation samples are a
#' set of perturbations DISJOINT from the ~23 evaluated in the
#' corresponding `_sigs`/`_true_sigs` (control samples plus a *different*
#' set of ~23 "background" perturbations, not the evaluated ones) -- so
#' the network/projection input never contains a sample for the
#' perturbation being predicted. See the "No data leakage" note in each
#' `data-raw/demo_*.R` script.
#'
#' | Dataset | Source context | Target context (`_se`) | Built by |
#' |---|---|---|---|
#' | SciPlex | K562 | A549 | `data-raw/demo_sciplex.R` |
#' | Tahoe | NCI-H23 | A498 | `data-raw/demo_tahoe.R` |
#' | DrugMatrix | kidney | liver | `data-raw/demo_drugmatrix.R` |
#' | Perturb-seq | K562 | RPE1 | `data-raw/demo_perturbseq.R` |
#'
#' @name demo_datasets
NULL

#' @rdname demo_datasets
#' @format `demo_sciplex_se`: a `SummarizedExperiment` with 1025 genes
#'   (Ensembl IDs) and 48 samples (Vehicle control + 23 background drugs,
#'   disjoint from the 23 evaluated in `demo_sciplex_sigs`/
#'   `demo_sciplex_true_sigs`), built from real SciPlex A549 pseudobulk
#'   expression
#'   (\url{https://figshare.com/articles/dataset/sciPlex_dataset/24681285?file=43381398}).
#'   Has `logcounts` (log2 CPM, default assay) and `counts` (raw) assays.
#'   `colData` has `orig.ident`, `product_name` (drug, or "Vehicle"), `replicate`.
"demo_sciplex_se"

#' @rdname demo_datasets
#' @format `demo_sciplex_sigs`: named list (by drug) of character vectors
#'   of Ensembl gene IDs -- subset of `sciplex.k562`.
"demo_sciplex_sigs"

#' @rdname demo_datasets
#' @format `demo_sciplex_true_sigs`: named list (by drug), each element a
#'   list with `up`/`up_full` -- subset of `sciplex.a549`.
"demo_sciplex_true_sigs"

#' @rdname demo_datasets
#' @format `demo_tahoe_se`: a `SummarizedExperiment` with 1808 genes (HGNC
#'   symbols) and 55 samples (DMSO_TF control + 23 background drugs,
#'   disjoint from the 23 evaluated in `demo_tahoe_sigs`/
#'   `demo_tahoe_true_sigs`), built from real Tahoe A498 pseudobulk
#'   expression
#'   (\url{https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M}).
#'   Has `logcounts` (log2 CPM, default assay) and `counts` (raw) assays.
"demo_tahoe_se"

#' @rdname demo_datasets
#' @format `demo_tahoe_sigs`: named list (by drug) of character vectors --
#'   subset of `tahoe.nci_h23` (source context).
"demo_tahoe_sigs"

#' @rdname demo_datasets
#' @format `demo_tahoe_true_sigs`: named list (by drug), each element a
#'   list with `up`/`up_full` -- subset of `tahoe.a498` (target context).
"demo_tahoe_true_sigs"

#' @rdname demo_datasets
#' @format `demo_drugmatrix_se`: a `SummarizedExperiment` with 1846 genes
#'   (gene symbols) and 84 samples (15 control + max-dose/time samples
#'   for 23 background drugs, disjoint from the 23 evaluated in
#'   `demo_drugmatrix_sigs`/`demo_drugmatrix_true_sigs`), built from real
#'   DrugMatrix liver microarray expression
#'   (\url{https://ntp.niehs.nih.gov/data/drugmatrix}). Single `logcounts`
#'   assay (already log-transformed microarray intensities, no further
#'   normalization applied). `colData` has `compound`, `dose`, `tissue`,
#'   `time`, `vehicle`.
"demo_drugmatrix_se"

#' @rdname demo_datasets
#' @format `demo_drugmatrix_sigs`: named list (by drug) of character
#'   vectors -- subset of `drugmatrix.kidney` (source context).
"demo_drugmatrix_sigs"

#' @rdname demo_datasets
#' @format `demo_drugmatrix_true_sigs`: named list (by drug), each element
#'   a list with `up`/`up_full` -- subset of `drugmatrix.liver` (target
#'   context).
"demo_drugmatrix_true_sigs"

#' @rdname demo_datasets
#' @format `demo_perturbseq_se`: a `SummarizedExperiment` with 1382 genes
#'   (HGNC symbols) and 84 samples (15 non-targeting control + up to 3
#'   replicates each for 23 background CRISPRi knockdowns, disjoint from
#'   the 23 evaluated in `demo_perturbseq_sigs`/
#'   `demo_perturbseq_true_sigs`), built from real Perturb-seq RPE1
#'   pseudobulk expression
#'   (\url{https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387}).
#'   Has `logcounts` (log2 CPM, default assay) and `counts` (raw) assays.
#'   `colData` has `gene` (knockdown target, or "non-targeting"),
#'   `gem_group`, `is_control`.
"demo_perturbseq_se"

#' @rdname demo_datasets
#' @format `demo_perturbseq_sigs`: named list (by knockdown target) of
#'   character vectors -- subset of `perturbseq.k562` (source context).
"demo_perturbseq_sigs"

#' @rdname demo_datasets
#' @format `demo_perturbseq_true_sigs`: named list (by knockdown target),
#'   each element a list with `up`/`up_full` -- subset of
#'   `perturbseq.rpe1` (target context).
"demo_perturbseq_true_sigs"
