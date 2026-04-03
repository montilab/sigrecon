test_that("ridge_benchmark_r2 fits a ridge model on a geneset", {
  skip_if_not_installed("glmnet")

  mat <- matrix(c(
    0, 0, 5, 5,
    0, 0, 4, 4,
    5, 5, 0, 0,
    1, 1, 1, 1
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(c("g1", "g2", "g3", "g4"), paste0("s", 1:4)))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = mat),
    colData = data.frame(perturbed = c(FALSE, FALSE, TRUE, TRUE))
  )

  r2 <- ridge_benchmark_r2(se = se, geneset = c("g1", "g2"), pb_col = "perturbed")

  expect_true(is.finite(r2))
  expect_gte(r2, 0)
})

test_that("sig_eval_table includes ridge benchmark columns when requested", {
  skip_if_not_installed("glmnet")

  mat <- matrix(c(
    0, 0, 5, 5,
    0, 0, 4, 4,
    5, 5, 0, 0,
    1, 1, 1, 1
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(c("g1", "g2", "g3", "g4"), paste0("s", 1:4)))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = mat),
    colData = data.frame(perturbed = c(FALSE, FALSE, TRUE, TRUE))
  )

  source_sigs <- list(drugA = c("g1", "g2"))
  pred_sigs <- list(drugA = c("g1", "g2"))
  true_sigs <- list(drugA = list(up = c("g1", "g2"), up_full = c("g1", "g2", "g3", "g4")))

  eval_df <- sig_eval_table(
    source_sigs = source_sigs,
    pred_sigs = pred_sigs,
    true_sigs = true_sigs,
    se = se,
    ridge_benchmark = TRUE,
    pb_col = "perturbed",
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_true(all(c("source_r2", "pred_r2", "true_r2") %in% colnames(eval_df)))
  expect_true(all(is.finite(as.matrix(eval_df[, c("source_r2", "pred_r2", "true_r2")]))))
})
