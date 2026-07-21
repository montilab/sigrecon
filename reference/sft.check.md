# Check and Return Soft-Thresholding Power from WGCNA

This function checks the soft-thresholding power estimate from WGCNA's
scale-free topology analysis and returns an appropriate power value.

## Usage

``` r
sft.check(sft)
```

## Arguments

- sft:

  A list object returned by WGCNA's pickSoftThreshold or
  pickSoftThreshold.fromSimilarity function.

## Value

An integer representing the soft-thresholding power to be used.

## Details

If a valid power estimate is found in the input, it is returned. If the
power estimate is NA, a default value of 6 is returned. The function
prints a message indicating which power is being used.
