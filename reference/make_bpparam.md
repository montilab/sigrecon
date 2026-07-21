# Create appropriate BiocParallel backend

Create appropriate BiocParallel backend

## Usage

``` r
make_bpparam(workers = 1, RNGseed = NULL, progress = FALSE, type = NULL)
```

## Arguments

- workers:

  Number of parallel workers (default: 1 for sequential)

- RNGseed:

  Random seed for reproducibility (default: NULL)

- progress:

  Show progress bar (default: FALSE)

- type:

  Force backend type: "multicore", "snow", or "serial" (default:
  auto-detect)

## Value

A BiocParallelParam object

## Examples

``` r
# Sequential (default)
bp <- make_bpparam()

# Parallel on Unix/Mac (forking)
bp <- make_bpparam(workers = 2, RNGseed = 123)

# Parallel on Windows (socket)
bp <- make_bpparam(workers = 2, RNGseed = 123)

# With progress bar
bp <- make_bpparam(workers = 2, progress = TRUE)
```
