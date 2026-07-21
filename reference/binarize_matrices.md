# Binarize Multiple Matrices Based on Cutoffs

Binarize Multiple Matrices Based on Cutoffs

## Usage

``` r
binarize_matrices(matrix_list, cutoff_list)
```

## Arguments

- matrix_list:

  A list of matrices to be binarized

- cutoff_list:

  A list of cutoff values, one for each matrix

## Value

A list of binarized matrices

## Details

This function takes a list of matrices and a corresponding list of
cutoff values. For each matrix, values less than or equal to the cutoff
are set to 0, and values greater than the cutoff are set to 1.
