test_that("projectCor reconstructs signatures from a SummarizedExperiment", {
  mat <- matrix(c(5, 3, 4, 2, 1,
                  2, 6, 5, 3, 2,
                  4, 2, 6, 1, 3,
                  3, 5, 2, 4, 6),
                nrow = 5, byrow = FALSE,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  sigs <- list(sigA = c("g1", "g2"), sigB = c("g3", "g4", "g5"))

  result <- projectCor(se, sigs)

  expect_identical(names(result), names(sigs))
  expect_identical(lengths(result), lengths(sigs))
  expect_identical(result, gsva_recon(se, sigs))
})

test_that("projectCor supports eigengene scoring", {
  mat <- matrix(c(1, 2, 3, 4,
                  2, 4, 6, 8,
                  4, 3, 2, 1,
                  1, 1, 2, 2,
                  2, 1, 2, 1),
                nrow = 5, byrow = TRUE,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  sigs <- list(sigA = c("g1", "g2"))

  result <- projectCor(se, sigs, score = "eigen")

  expect_identical(names(result), "sigA")
  expect_identical(length(result$sigA), 2L)
  expect_true(setequal(result$sigA, c("g1", "g2")))
})

test_that("projectCor eigengene scoring skips zero-overlap signatures", {
  mat <- matrix(c(1, 2, 3, 4,
                  2, 4, 6, 8,
                  4, 3, 2, 1,
                  1, 1, 2, 2,
                  2, 1, 2, 1),
                nrow = 5, byrow = TRUE,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  sigs <- list(sigA = c("g1", "g2"), sigB = c("missing1", "missing2"))

  expect_message(
    result <- projectCor(se, sigs, score = "eigen"),
    "Skipping signature 'sigB' because no genes were found in the expression assay."
  )

  expect_identical(names(result), "sigA")
  expect_identical(length(result$sigA), 2L)
})
