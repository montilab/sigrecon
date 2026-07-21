# Count number of displaced seeds

Quantifies the number of genes from the original source signature that
are still in the recontextualized signature `not_displaced`, and the
number of genes that are `displaced`. Note: This operation is not
commutative.

## Usage

``` r
count_displaced_genes(recon_sig, seed_sig)
```

## Arguments

- recon_sig:

  Named list of genesets

- seed_sig:

  Named list of genesets

## Value

Named list for which each entry is a list with two elements: displaced
seed genes, and non-displaced seed genes
