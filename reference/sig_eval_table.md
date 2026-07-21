# Evaluate Signature Prediction

This function calculates various evaluation metrics for predicted gene
signatures against true (ground truth) gene signatures, given a source
signature as input. Metrics include Jaccard index, gene displacement,
and Kolmogorov-Smirnov (KS) test statistics for rank agreement.

## Usage

``` r
sig_eval_table(
  source_sigs,
  pred_sigs,
  true_sigs,
  splits = FALSE,
  split_file = NULL,
  split_pb_col = "drug",
  split_type = NULL,
  source = "source_context",
  target = "target_context",
  BPPARAM = NULL
)
```

## Arguments

- source_sigs:

  A named list of genesets.

- pred_sigs:

  A named list of genesets.

- true_sigs:

  A named list of genesets containing both the cutoff signature `up` and
  the full signature `up_full`.

- splits:

  Boolean, indicating whether `pred_sigs` is a named list of splits,
  each split containing a separate list of geneset predictions. When
  `split_file` is supplied, split names in `pred_sigs` must be present
  as columns in the split table.

- split_file:

  Optional path to a `.csv` file containing a split table with a
  perturbation column and split logical columns.

- split_pb_col:

  Name of the perturbation column in `split_file`. Default is `"drug"`.

- split_type:

  Optional split subset to evaluate. If `"10th"`, evaluates
  perturbations with `FALSE` in the relevant split column. If `"90th"`,
  evaluates perturbations with `TRUE` in the relevant split column.

- source:

  A character string describing the starting biological context.

- target:

  A character string describing the target biological context.

- BPPARAM:

  A BiocParallelParam object for parallel processing. If NULL, uses
  SerialParam.

## Value

A data frame with evaluation results. Each row corresponds to a
perturbation and context.
