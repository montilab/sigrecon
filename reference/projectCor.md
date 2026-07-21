# Projection based Recontextualization

This function reconstructs gene signatures based on their correlation
with per-sample projection scores computed from the input signatures.

## Usage

``` r
projectCor(se, sigs, score = c("gsva", "eigen"))
```

## Arguments

- se:

  A SummarizedExperiment object containing gene expression data.

- sigs:

  A list of gene signatures, where each element is a character vector of
  gene names.

- score:

  Scoring method used to score samples against input signatures. Either
  `"gsva"` or `"eigen"`. Default is `"gsva"`.

## Value

A list of reconstructed gene signatures, with the same structure as the
input `sigs`.

## Details

The function performs the following steps:

1.  Calculates projection scores for the input signatures with GSVA or
    eigengenes.

2.  Computes the correlation between gene expression and projection
    scores.

3.  Ranks genes based on their correlation with each signature's
    projection scores.

4.  Selects the top-ranking genes to form new signatures of the same
    length as the original ones.
