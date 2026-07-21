# Find Common Variable Genes Across Seurat Objects

This function identifies common variable genes across multiple Seurat
objects, up to a specified limit.

## Usage

``` r
seurat_common_var_genes(seurat_objs, limit)
```

## Arguments

- seurat_objs:

  A list of Seurat objects to compare.

- limit:

  An integer specifying the maximum number of common variable genes to
  return.

## Value

A character vector of common variable gene names.

## Details

The function iterates through the variable features of each Seurat
object, selecting genes that are present in all objects. It continues
until it reaches the specified limit or exhausts all common variable
genes.
