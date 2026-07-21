# Fetch a bundled perturbational signature dataset

Downloads (and locally caches, via `BiocFileCache`) one of the datasets
listed in
[`list_datasets()`](https://montilab.github.io/sigrecon/reference/list_datasets.md).
Datasets are hosted as release assets on `montilab/sigrecon` rather than
bundled with the package, so the first call to `get_dataset()` for a
given name requires network access; subsequent calls reuse the local
cache.

## Usage

``` r
get_dataset(name, force = FALSE)
```

## Arguments

- name:

  Dataset name. See
  [`list_datasets()`](https://montilab.github.io/sigrecon/reference/list_datasets.md)
  for available names.

- force:

  Logical. If `TRUE`, re-download even if a cached copy exists. Default
  is `FALSE`.

## Value

A named list of genesets. Each element has an `up` component (top 100
DEGs) and an `up_full` component (the full ranked gene list).

## Examples

``` r
if (FALSE) { # \dontrun{
tahoe_nci_h23 <- get_dataset("tahoe.nci_h23")
} # }
```
