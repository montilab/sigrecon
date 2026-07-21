# Benchmarking Design: Splits and Evaluation Metrics

``` r

library(sigrecon)
```

The [Getting
Started](https://montilab.github.io/sigrecon/articles/getting-started.md)
vignette runs
[`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md)
once and glances at a few of its output columns. This vignette covers
the evaluation design in full: what every output column means, the
ctrl/10th/90th split design used to benchmark methods at varying amounts
of target-context information, and how to compare against a no-change
baseline with
[`paired_eval_table()`](https://montilab.github.io/sigrecon/reference/paired_eval_table.md).

## `sig_eval_table()` output columns

``` r

recon <- projectCor(demo_sciplex_se, demo_sciplex_sigs, score = "gsva")
#> ℹ No assay name provided; using default assay 'logcounts'

eval_df <- sig_eval_table(
  source_sigs = demo_sciplex_sigs,
  pred_sigs   = recon,
  true_sigs   = demo_sciplex_true_sigs,
  source      = "sciplex_k562",
  target      = "sciplex_a549"
)

str(eval_df, max.level = 1)
#> 'data.frame':    23 obs. of  13 variables:
#>  $ source     : chr  "sciplex_k562" "sciplex_k562" "sciplex_k562" "sciplex_k562" ...
#>  $ target     : chr  "sciplex_a549" "sciplex_a549" "sciplex_a549" "sciplex_a549" ...
#>  $ displaced  : int  82 40 68 13 12 70 74 46 63 17 ...
#>  $ kept       : int  18 8 16 2 2 29 26 10 6 8 ...
#>  $ gene       : chr  "Abexinostat (PCI-24781)" "Alvespimycin (17-DMAG) HCl" "AR-42" "Belinostat (PXD101)" ...
#>  $ jacc       : num  0.105 0.0857 0.076 0.0177 0.0571 ...
#>  $ ES         : num  0.554 0.195 0.499 0.314 0.83 ...
#>  $ NES        : num  2.804 0.807 2.395 0.917 2.456 ...
#>  $ pval       : num  1.20e-12 7.58e-01 1.22e-07 5.72e-01 4.13e-07 ...
#>  $ padj       : num  5.51e-12 7.92e-01 3.52e-07 6.26e-01 1.06e-06 ...
#>  $ log2err    : num  0.9101 0.0576 0.6901 0.0727 0.675 ...
#>  $ size       : int  94 41 69 12 12 85 87 24 61 14 ...
#>  $ leadingEdge:List of 23
#>   ..- attr(*, "class")= chr "AsIs"
```

Each row is one perturbation. Columns:

| Column | Meaning |
|----|----|
| `source`, `target` | Labels you passed in, identifying the source/target contexts |
| `gene` | The perturbation name (row identifier) |
| `displaced`, `kept` | How many source-signature genes were replaced vs. retained by recontextualization. `kept / (kept + displaced)` is the fraction of the source signature recontextualization left unchanged. |
| `jacc` | Jaccard overlap between the recontextualized signature and the true target-context signature’s `up` genes |
| `ES`, `NES` | fgsea enrichment score / normalized enrichment score of the recontextualized signature against the true signature’s full ranked gene list (`up_full`) – how concentrated the recontextualized genes are near the top of the true ranking |
| `pval`, `padj` | Significance of that enrichment (raw and BH-adjusted) |
| `log2err`, `size` | fgsea’s estimated log2 error on the p-value, and the number of recontextualized genes that overlapped the ranked list at all |
| `leadingEdge` | The specific genes driving the enrichment signal (a list-column – one character vector per row) |

Higher `jacc` and `NES` mean the recontextualized signature more closely
resembles the true target-context signature – i.e. better
recontextualization.

``` r

# leadingEdge is a list-column: one entry per row
eval_df$leadingEdge[[1]]
#>  [1] "ENSG00000160963" "ENSG00000103740" "ENSG00000259124" "ENSG00000173727"
#>  [5] "ENSG00000184611" "ENSG00000244468" "ENSG00000186115" "ENSG00000172264"
#>  [9] "ENSG00000197959" "ENSG00000144868" "ENSG00000069667" "ENSG00000169918"
#> [13] "ENSG00000110076" "ENSG00000104967" "ENSG00000132718" "ENSG00000150347"
#> [17] "ENSG00000257261" "ENSG00000078018" "ENSG00000184524" "ENSG00000171431"
#> [21] "ENSG00000214456" "ENSG00000248810" "ENSG00000183023" "ENSG00000105409"
#> [25] "ENSG00000141750" "ENSG00000228509" "ENSG00000238212" "ENSG00000158813"
#> [29] "ENSG00000163531" "ENSG00000240521" "ENSG00000106772" "ENSG00000136275"
#> [33] "ENSG00000129244" "ENSG00000150625" "ENSG00000107295" "ENSG00000138028"
#> [37] "ENSG00000134986" "ENSG00000161270" "ENSG00000007516" "ENSG00000198910"
#> [41] "ENSG00000144285" "ENSG00000197872" "ENSG00000124302" "ENSG00000178217"
#> [45] "ENSG00000187672" "ENSG00000123104" "ENSG00000184005" "ENSG00000115419"
#> [49] "ENSG00000109654" "ENSG00000140450" "ENSG00000102003" "ENSG00000021645"
#> [53] "ENSG00000221946" "ENSG00000082014" "ENSG00000141526" "ENSG00000151572"
#> [57] "ENSG00000143603" "ENSG00000182771" "ENSG00000267761" "ENSG00000240875"
#> [61] "ENSG00000006747" "ENSG00000155657" "ENSG00000036530"
```

## The ctrl / 10th / 90th split design

A recontextualization method that’s simply handed the true
target-context signature as extra input isn’t really being tested – it’s
cheating. The package’s benchmarking study evaluates methods under three
conditions that control how much target-context information is
available:

- **`ctrl`**: none of the target perturbations’ information is available
  to the method beyond the source signature itself (what every example
  so far in this vignette series has done).
- **`ctrl_10th`**: the method is allowed 1/10th of the target
  perturbations as reference (e.g. to tune a network or model), then
  evaluated on the *other* 9/10th – the harder, held-out majority.
- **`ctrl_90th`**: the reverse – allowed 9/10th, evaluated on the
  remaining 1/10th.

In the real analysis scripts, which perturbations fall in which split is
decided once and saved to a `drug_splits.csv`/`pb_splits.csv` file: one
column per split (`split_1`, `split_2`, …), each a logical vector
marking a fixed ~10%/90% partition of perturbations, aligned to a
perturbation-label column. The same partition serves both regimes by
flipping which side is the reference and which is evaluated:
`split_type = "10th"` uses the `TRUE` (small) side as reference and
evaluates the `FALSE` (majority) side; `split_type = "90th"` uses the
`FALSE` side as reference and evaluates the `TRUE` side.
`sig_eval_table(splits = TRUE, split_file = ..., split_pb_col = ..., split_type = c("10th", "90th"))`
reads that file directly and filters which perturbations get evaluated
for you.

### A worked example

`splits = TRUE` expects `pred_sigs` to be a *named list of splits*, each
element itself a signature list, with names matching the split file’s
column names. Building a small one here from the demo signatures
(reusing the same predictions per split for simplicity – in a real
study, each split’s predictions come from a method run only on that
split’s reference perturbations):

``` r

set.seed(1)
drugs <- names(demo_sciplex_sigs)

split_tbl <- data.frame(
  drug     = drugs,
  split_1  = sample(c(TRUE, FALSE), length(drugs), replace = TRUE, prob = c(0.1, 0.9)),
  split_2  = sample(c(TRUE, FALSE), length(drugs), replace = TRUE, prob = c(0.1, 0.9))
)
table(split_tbl$split_1)
#> 
#> FALSE  TRUE 
#>    19     4

split_file <- tempfile(fileext = ".csv")
write.csv(split_tbl, split_file, row.names = FALSE)

pred_sigs_splits <- list(split_1 = recon, split_2 = recon)
```

``` r

eval_10th <- sig_eval_table(
  source_sigs  = demo_sciplex_sigs,
  pred_sigs    = pred_sigs_splits,
  true_sigs    = demo_sciplex_true_sigs,
  source       = "sciplex_k562",
  target       = "sciplex_a549",
  splits       = TRUE,
  split_file   = split_file,
  split_pb_col = "drug",
  split_type   = "10th"
)
#> Processing Split split_1
#> Processing Split split_2

# One extra "split" column identifying which split each row came from
table(eval_10th$split)
#> 
#>  1  2 
#> 19 23
nrow(eval_10th)
#> [1] 42
```

`split_type = "10th"` evaluated only the perturbations where `split_1`/
`split_2` was `FALSE` (the withheld 9/10ths) – fewer rows than the 23
perturbations in the full `demo_sciplex_sigs`, and split across the two
`split_1`/`split_2` conditions.

## Comparing against a no-change baseline

The natural baseline for “did recontextualization help” is doing nothing
– using the unmodified source signature as the prediction.
[`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md)
computes this the same way as any other prediction, just passing the
source signatures as `pred_sigs`:

``` r

no_change_df <- sig_eval_table(
  source_sigs = demo_sciplex_sigs,
  pred_sigs   = demo_sciplex_sigs,
  true_sigs   = demo_sciplex_true_sigs,
  source      = "sciplex_k562",
  target      = "sciplex_a549"
)
#> Warning in fgsea_wrapper(ref = ref, data = data, scoreType = scoreType, : No
#> overlap between 'data' and 'ref', returning NULL
```

[`paired_eval_table()`](https://montilab.github.io/sigrecon/reference/paired_eval_table.md)
formalizes the no-change-vs-recontextualized comparison: it joins the
two evaluation tables on `source`/`gene`, computes
`kept_alpha`/`NES_delta`/`Jacc_delta`, and prints a paired significance
summary (median performance, then a paired Wilcoxon test per source
combined across perturbations via Fisher’s method). It takes the
no-change table as a **file path**, not a data frame in memory, so save
it first:

``` r

no_change_path <- tempfile(fileext = ".rds")
saveRDS(no_change_df, no_change_path)

paired_df <- paired_eval_table(eval_df, no_change_path)
#> # A tibble: 1 × 2
#>   NES_median jacc_median
#>        <dbl>       <dbl>
#> 1       2.32      0.0857
#> # A tibble: 1 × 2
#>   KS_meta_p Jacc_meta_p
#>       <dbl>       <dbl>
#> 1     0.957       0.618
head(paired_df[, c("gene", "jacc_FALSE", "jacc_TRUE", "NES_FALSE", "NES_TRUE", "Jacc_delta", "NES_delta")])
#> # A tibble: 6 × 7
#>   gene              jacc_FALSE jacc_TRUE NES_FALSE NES_TRUE Jacc_delta NES_delta
#>   <chr>                  <dbl>     <dbl>     <dbl>    <dbl>      <dbl>     <dbl>
#> 1 Abexinostat (PCI…     0.117     0.105      2.40     2.80     -0.0123    0.404 
#> 2 Alvespimycin (17…     0.0179    0.0857     0.674    0.807     0.0679    0.133 
#> 3 AR-42                 0.157     0.0760     2.47     2.40     -0.0812   -0.0791
#> 4 Belinostat (PXD1…     0.0550    0.0177     2.28     0.917    -0.0373   -1.37  
#> 5 CUDC-101              0.121     0.0571     2.52     2.46     -0.0641   -0.0624
#> 6 CUDC-907              0.0874    0.171      2.85     3.00      0.0832    0.158
```

`_FALSE` columns are the no-change baseline, `_TRUE` columns are the
recontextualized prediction (matching
[`sig_eval_table()`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md)’s
internal `splits`/`filter_eval_inputs` naming convention, not a literal
boolean value) – `Jacc_delta`/`NES_delta` being positive means
recontextualization improved on the baseline for that perturbation.

## Next steps

- The [Getting
  Started](https://montilab.github.io/sigrecon/articles/getting-started.md)
  vignette for the full recontextualization workflow
- The [Network Propagation in
  Depth](https://montilab.github.io/sigrecon/articles/network-propagation-details.md)
  vignette for
  [`network_sig()`](https://montilab.github.io/sigrecon/reference/network_sig.md)’s
  options
- [`?sig_eval_table`](https://montilab.github.io/sigrecon/reference/sig_eval_table.md),
  [`?paired_eval_table`](https://montilab.github.io/sigrecon/reference/paired_eval_table.md)
  for full parameter references
