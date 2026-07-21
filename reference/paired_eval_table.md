# Pair a Recontextualization Evaluation Table With a No-Change Baseline

Merges a
[`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md)
result with a saved no-change evaluation table, adds paired delta
metrics, prints median performance and paired Wilcoxon meta-test
summaries, and returns the paired table.

## Usage

``` r
paired_eval_table(combined_df, no_change_eval_path)
```

## Arguments

- combined_df:

  A data frame returned by
  [`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md)
  for a recontextualized/predicted signature set.

- no_change_eval_path:

  File path to an `.rds` containing the no-change
  [`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md)
  output.

## Value

A tibble containing paired no-change (`*_FALSE`) and predicted
(`*_TRUE`) metrics, plus `kept_alpha`, `NES_delta`, and `Jacc_delta`.
