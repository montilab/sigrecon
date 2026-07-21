# Filter Significant Genes by Perturbation

Filters and ranks genes by significance for each perturbation in a
differential expression table. Returns both all genes ranked by a
combined score and significantly up-regulated genes that pass
thresholds.

## Usage

``` r
sig_filter_fn(
  diff_table,
  perts,
  alpha = 0.05,
  limit = 100,
  pert_col = "product",
  log2fc_col = "avg_log2FC",
  pval_col = "p_val_adj",
  geneid_col = "ensembl_id"
)
```

## Arguments

- diff_table:

  A data frame containing differential expression results with
  perturbation identifiers, log2 fold changes, and adjusted p-values.

- perts:

  A character vector of perturbation names to filter for.

- alpha:

  Numeric. Adjusted p-value threshold for significance. Default is 0.05.

- limit:

  Integer. Maximum number of top significant genes to return per
  perturbation. Default is 100.

- pert_col:

  Character. Name of the column containing perturbation identifiers.
  Default is "product".

- log2fc_col:

  Character. Name of the column containing log2 fold change values.
  Default is "avg_log2FC".

- pval_col:

  Character. Name of the column containing adjusted p-values. Default is
  "p_val_adj".

- geneid_col:

  Character. Name of the column containing gene identifiers.

## Value

A nested list where each perturbation contains:

- `up`: Character vector of row names for top significantly upregulated
  genes (filtered by alpha and limited by limit parameter)

- `up_full`: Character vector of row names for all genes ranked by
  combined score (log2FC \* -log10(p_val_adj))

## Details

The function computes a combined significance score as: \$\$score =
log2FC \times -log10(p\_{adj})\$\$

For the "up" results, genes must meet three criteria:

1.  log2FC \> 0 (upregulated)

2.  adjusted p-value \<= alpha

3.  Ranked in top N genes (specified by limit)
