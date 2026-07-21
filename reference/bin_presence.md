# Check Presence in Bins for Numbers 1-100

Given a numeric vector with values between 1 and 100, returns a logical
vector indicating whether at least one value falls into each bin of size
10 (i.e., 1-10, 11-20, ..., 91-100).

## Usage

``` r
bin_presence(x)
```

## Arguments

- x:

  Numeric vector. Values should be between 1 and 100.

## Value

A named logical vector of length 10. Each element is `TRUE` if at least
one value in `x` falls into the corresponding bin, otherwise `FALSE`.

## Examples

``` r
vec <- c(3, 15, 27, 45, 58, 99)
bin_presence(vec)
#> Error in bin_presence(vec): could not find function "bin_presence"
```
