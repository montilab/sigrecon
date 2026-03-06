# tests/testthat/test-fgsea-biocparallel.R

# Helper function to create test data
create_test_data <- function(n_tests = 10, seed = 999) {
  set.seed(seed)
  
  ref_vecs <- lapply(1:n_tests, function(i) {
    paste0("GENE", sample(1:500, 200))
  })
  names(ref_vecs) <- paste0("test", 1:n_tests)
  
  data_vecs <- lapply(1:n_tests, function(i) {
    sample(ref_vecs[[i]], 30)
  })
  names(data_vecs) <- paste0("test", 1:n_tests)
  
  list(ref_vecs = ref_vecs, data_vecs = data_vecs)
}

# Helper to get appropriate number of workers for testing
get_test_workers <- function(requested = 2) {
  check_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(check_limit) && check_limit == "TRUE") {
    return(2L)
  }
  return(requested)
}

# Test: BiocParallel reproducibility ----
test_that("v.fgsea produces reproducible results with BiocParallel", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("BiocParallel")
  
  # Create test data
  test_data <- create_test_data(n_tests = 10)
  
  # Determine workers
  n_workers <- get_test_workers(2)
  
  # Run twice with same seed and backend
  bp <- BiocParallel::MulticoreParam(workers = n_workers, RNGseed = 123)
  
  results_1 <- v.fgsea(
    ref_vecs = test_data$ref_vecs,
    data_vecs = test_data$data_vecs,
    BPPARAM = bp
  )
  
  results_2 <- v.fgsea(
    ref_vecs = test_data$ref_vecs,
    data_vecs = test_data$data_vecs,
    BPPARAM = BiocParallel::MulticoreParam(workers = n_workers, RNGseed = 123)
  )
  
  # Sort for comparison
  results_1 <- results_1[order(results_1$name), ]
  results_2 <- results_2[order(results_2$name), ]
  
  # Test reproducibility
  expect_equal(results_1$NES, results_2$NES, 
               info = "NES should be identical with same seed")
  expect_equal(results_1$pval, results_2$pval,
               info = "P-values should be identical with same seed")
  expect_equal(results_1$ES, results_2$ES,
               info = "ES should be identical with same seed")
})

# Test: Results consistency between sequential and parallel ----
test_that("v.fgsea gives identical results between sequential and parallel", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("BiocParallel")
  
  # Create test data
  test_data <- create_test_data(n_tests = 10)
  n_workers <- get_test_workers(4)
  
  # Run sequential with seed
  results_seq <- v.fgsea(
    ref_vecs = test_data$ref_vecs,
    data_vecs = test_data$data_vecs,
    BPPARAM = BiocParallel::SerialParam(RNGseed = 123)
  )
  
  # Run parallel with same seed
  results_par <- v.fgsea(
    ref_vecs = test_data$ref_vecs,
    data_vecs = test_data$data_vecs,
    BPPARAM = BiocParallel::MulticoreParam(workers = n_workers, RNGseed = 123)
  )
  
  # Sort for comparison
  results_seq <- results_seq[order(results_seq$name), ]
  results_par <- results_par[order(results_par$name), ]
  
  # Test exact equality
  expect_equal(results_seq$NES, results_par$NES, 
               info = "NES should be identical with same seed")
  expect_equal(results_seq$pval, results_par$pval,
               info = "P-values should be identical with same seed")
  expect_equal(results_seq$ES, results_par$ES,
               info = "ES should be identical with same seed")
  expect_equal(results_seq$padj, results_par$padj,
               info = "Adjusted p-values should be identical")
})

# Test: Parallel speedup ----
test_that("v.fgsea parallel execution is faster than sequential", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("BiocParallel")
  skip_on_cran()  # Skip timing tests on CRAN
  skip_on_os("windows")  # MulticoreParam doesn't parallelize on Windows
  
  # Create larger test data for meaningful timing
  test_data <- create_test_data(n_tests = 500)
  n_workers <- get_test_workers(4)
  
  # Time sequential
  time_seq <- system.time({
    results_seq <- v.fgsea(
      ref_vecs = test_data$ref_vecs,
      data_vecs = test_data$data_vecs,
      BPPARAM = BiocParallel::SerialParam()
    )
  })
  
  # Time parallel
  time_par <- system.time({
    results_par <- v.fgsea(
      ref_vecs = test_data$ref_vecs,
      data_vecs = test_data$data_vecs,
      BPPARAM = BiocParallel::MulticoreParam(workers = n_workers)
    )
  })
  
  # Parallel should be faster (with some tolerance for overhead)
  speedup <- time_seq[3] / time_par[3]
  expect_gt(speedup, 1.2,
            label = sprintf("Parallel should be faster (speedup: %.2f)", speedup))
})

# Test: SnowParam works on all platforms ----
test_that("v.fgsea works with SnowParam", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("BiocParallel")
  skip_on_cran()  # Skip due to potential resource constraints
  
  # Create test data
  test_data <- create_test_data(n_tests = 5)
  n_workers <- get_test_workers(2)  # Use fewer workers for SnowParam
  
  # Run with SnowParam (works on all platforms including Windows)
  results <- v.fgsea(
    ref_vecs = test_data$ref_vecs,
    data_vecs = test_data$data_vecs,
    BPPARAM = BiocParallel::SnowParam(workers = n_workers, RNGseed = 123)
  )
  
  # Should complete without error and return expected structure
  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 5)
  expect_true(all(c("name", "ES", "NES", "pval", "padj") %in% names(results)))
})

# Test: SerialParam works ----
test_that("v.fgsea works with SerialParam", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("BiocParallel")
  
  # Create test data
  test_data <- create_test_data(n_tests = 5)
  
  # Run with SerialParam
  results <- v.fgsea(
    ref_vecs = test_data$ref_vecs,
    data_vecs = test_data$data_vecs,
    BPPARAM = BiocParallel::SerialParam()
  )
  
  # Should complete without error and return expected structure
  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 5)
  expect_true(all(c("name", "ES", "NES", "pval", "padj") %in% names(results)))
  expect_true(all(!is.na(results$NES)))
})

# Test: NULL BPPARAM defaults to SerialParam ----
test_that("v.fgsea defaults to SerialParam when BPPARAM is NULL", {
  skip_if_not_installed("fgsea")
  
  # Create test data
  test_data <- create_test_data(n_tests = 3)
  
  # Run with NULL BPPARAM
  results <- v.fgsea(
    ref_vecs = test_data$ref_vecs,
    data_vecs = test_data$data_vecs,
    BPPARAM = NULL
  )
  
  # Should complete without error
  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 3)
})

# Test: Input validation ----
test_that("v.fgsea validates inputs correctly", {
  skip_if_not_installed("fgsea")
  
  ref_vecs <- list(test1 = c("GENE1", "GENE2", "GENE3"))
  data_vecs <- list(test1 = c("GENE1", "GENE2"))
  
  # Should error with non-list inputs
  expect_error(
    v.fgsea(ref_vecs = "not a list", data_vecs = data_vecs),
    "must be named lists"
  )
  
  expect_error(
    v.fgsea(ref_vecs = ref_vecs, data_vecs = "not a list"),
    "must be named lists"
  )
  
  # Should error with unnamed lists
  expect_error(
    v.fgsea(
      ref_vecs = list(c("GENE1", "GENE2")),
      data_vecs = list(c("GENE1"))
    ),
    "must have names"
  )
  
  # Should error with no matching names
  expect_error(
    v.fgsea(
      ref_vecs = list(test1 = c("GENE1", "GENE2")),
      data_vecs = list(test2 = c("GENE1"))
    ),
    "No matching names"
  )
})

# Test: FDR correction ----
test_that("v.fgsea correctly adjusts p-values", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("BiocParallel")
  
  # Create test data
  test_data <- create_test_data(n_tests = 10)
  
  results <- v.fgsea(
    ref_vecs = test_data$ref_vecs,
    data_vecs = test_data$data_vecs,
    BPPARAM = BiocParallel::SerialParam()
  )
  
  # padj should be >= pval (FDR correction is conservative)
  expect_true(all(results$padj >= results$pval, na.rm = TRUE),
              info = "Adjusted p-values should be >= raw p-values")
  
  # Manually verify FDR correction
  valid_pvals <- !is.na(results$pval)
  expected_padj <- p.adjust(results$pval[valid_pvals], method = "BH")
  expect_equal(results$padj[valid_pvals], expected_padj,
               info = "FDR correction should match p.adjust with BH method")
})

# Test: Output structure ----
test_that("v.fgsea returns correct output structure", {
  skip_if_not_installed("fgsea")
  
  # Create test data
  test_data <- create_test_data(n_tests = 3)
  
  results <- v.fgsea(
    ref_vecs = test_data$ref_vecs,
    data_vecs = test_data$data_vecs,
    BPPARAM = BiocParallel::SerialParam()
  )
  
  # Check data frame structure
  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 3)
  
  # Check required columns
  required_cols <- c("name", "ES", "NES", "pval", "padj", "log2err", "size", "leadingEdge")
  expect_true(all(required_cols %in% names(results)),
              info = "All required columns should be present")
  
  # Check column types
  expect_type(results$name, "character")
  expect_type(results$ES, "double")
  expect_type(results$NES, "double")
  expect_type(results$pval, "double")
  expect_type(results$padj, "double")
  expect_type(results$log2err, "double")
  expect_type(results$size, "integer")
  expect_type(results$leadingEdge, "list")
})

# Test: Different scoreTypes ----
test_that("v.fgsea works with different scoreTypes", {
  skip_if_not_installed("fgsea")
  
  # Create test data
  test_data <- create_test_data(n_tests = 5)
  
  # Test "std" (default)
  results_std <- v.fgsea(
    ref_vecs = test_data$ref_vecs,
    data_vecs = test_data$data_vecs,
    scoreType = "std",
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_true(all(!is.na(results_std$NES)))
  
  # Test "pos"
  results_pos <- v.fgsea(
    ref_vecs = test_data$ref_vecs,
    data_vecs = test_data$data_vecs,
    scoreType = "pos",
    BPPARAM = BiocParallel::SerialParam()
  )
  # All ES should be >= 0 for "pos"
  expect_true(all(results_pos$ES >= 0 | is.na(results_pos$ES)),
              info = "ES should be non-negative for scoreType='pos'")
  
  # Test "neg"
  results_neg <- v.fgsea(
    ref_vecs = test_data$ref_vecs,
    data_vecs = test_data$data_vecs,
    scoreType = "neg",
    BPPARAM = BiocParallel::SerialParam()
  )
  # All ES should be <= 0 for "neg"
  expect_true(all(results_neg$ES <= 0 | is.na(results_neg$ES)),
              info = "ES should be non-positive for scoreType='neg'")
})

# Test: Error handling ----
test_that("v.fgsea handles errors gracefully", {
  skip_if_not_installed("fgsea")
  
  # Create test data with one problematic entry
  ref_vecs <- list(
    good1 = paste0("GENE", 1:100),
    bad = character(0),  # Empty reference
    good2 = paste0("GENE", 1:100)
  )
  
  data_vecs <- list(
    good1 = paste0("GENE", 1:20),
    bad = paste0("GENE", 1:20),
    good2 = paste0("GENE", 1:20)
  )
  
  # Should warn but not fail completely
  expect_warning(
    results <- v.fgsea(
      ref_vecs = ref_vecs,
      data_vecs = data_vecs,
      BPPARAM = BiocParallel::SerialParam()
    ),
    "Error processing"
  )
  
  # Should still return results for good entries
  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 3)
  
  # Bad entry should have NA values
  bad_row <- results[results$name == "bad", ]
  expect_true(all(is.na(bad_row[c("ES", "NES", "pval")])))
})