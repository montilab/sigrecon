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
