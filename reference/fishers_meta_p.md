# Combine P-values Using Fisher's Method

Combines multiple p-values from independent tests into a single meta
p-value using Fisher's method (also known as Fisher's combined
probability test).

## Usage

``` r
fishers_meta_p(pvals)
```

## Arguments

- pvals:

  Numeric vector of p-values to combine. All values must be between 0
  and 1 (exclusive of 0). NA values are not allowed.

## Value

A single numeric value representing the combined p-value.

## Details

Fisher's method combines p-values by computing the test statistic: \$\$X
= -2 \sum\_{i=1}^{k} \ln(p_i)\$\$

Under the null hypothesis (all individual null hypotheses are true), X
follows a chi-squared distribution with 2k degrees of freedom, where k
is the number of p-values. This method assumes that the tests are
independent.

The combined p-value represents the probability of observing the given
set of p-values (or more extreme) if all null hypotheses are true.

## Examples

``` r
# Combine three p-values
fishers_meta_p(c(0.01, 0.03, 0.25))
#> [1] 0.004170318

# Two significant p-values
fishers_meta_p(c(0.001, 0.005))
#> [1] 6.603036e-05

# Mix of significant and non-significant
fishers_meta_p(c(0.045, 0.23, 0.67))
#> [1] 0.1270948
```
