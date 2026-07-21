# Count number of displaced seeds Note: This not commutative. Num_displaced is the number of genes that are from the original seed that are not in the recon.

Count number of displaced seeds Note: This not commutative.
Num_displaced is the number of genes that are from the original seed
that are not in the recon.

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
