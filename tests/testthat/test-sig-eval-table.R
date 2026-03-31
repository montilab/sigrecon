test_that("sig_eval_table filters perturbations by split file for 10th and 90th splits", {
  source_sigs <- list(
    drugA = c("g1", "g2"),
    drugB = c("g2", "g3"),
    drugC = c("g3", "g4")
  )

  pred_split_1 <- list(
    drugA = c("g1", "g2"),
    drugB = c("g2", "g4"),
    drugC = c("g3", "g4")
  )
  pred_split_2 <- list(
    drugA = c("g1", "g3"),
    drugB = c("g2", "g3"),
    drugC = c("g1", "g4")
  )

  pred_sigs <- list(split_1 = pred_split_1, split_2 = pred_split_2)

  true_sigs <- list(
    drugA = list(up = c("g1", "g2"), up_full = c("g1", "g2", "g3", "g4")),
    drugB = list(up = c("g2", "g3"), up_full = c("g2", "g3", "g1", "g4")),
    drugC = list(up = c("g3", "g4"), up_full = c("g3", "g4", "g1", "g2"))
  )

  split_tbl <- data.frame(
    drug = c("drugA", "drugB", "drugC"),
    split_1 = c(FALSE, TRUE, FALSE),
    split_2 = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  split_file <- tempfile(fileext = ".csv")
  utils::write.csv(split_tbl, split_file, row.names = FALSE)

  eval_10th <- sig_eval_table(
    source_sigs = source_sigs,
    pred_sigs = pred_sigs,
    true_sigs = true_sigs,
    splits = TRUE,
    split_file = split_file,
    split_type = "10th",
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_equal(eval_10th$gene[eval_10th$split == 1], c("drugA", "drugC"))
  expect_equal(eval_10th$gene[eval_10th$split == 2], "drugB")

  eval_90th <- sig_eval_table(
    source_sigs = source_sigs,
    pred_sigs = pred_sigs,
    true_sigs = true_sigs,
    splits = TRUE,
    split_file = split_file,
    split_type = "90th",
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_equal(eval_90th$gene[eval_90th$split == 1], "drugB")
  expect_equal(eval_90th$gene[eval_90th$split == 2], c("drugA", "drugC"))
})

test_that("sig_eval_table requires split columns to match names(pred_sigs)", {
  source_sigs <- list(drugA = c("g1", "g2"))
  pred_sigs <- list(split_1 = list(drugA = c("g1", "g2")))
  true_sigs <- list(drugA = list(up = c("g1", "g2"), up_full = c("g1", "g2", "g3")))

  split_tbl <- data.frame(
    drug = "drugA",
    wrong_split = FALSE,
    stringsAsFactors = FALSE
  )
  split_file <- tempfile(fileext = ".csv")
  utils::write.csv(split_tbl, split_file, row.names = FALSE)

  expect_error(
    sig_eval_table(
      source_sigs = source_sigs,
      pred_sigs = pred_sigs,
      true_sigs = true_sigs,
      splits = TRUE,
      split_file = split_file,
      split_type = "10th",
      BPPARAM = BiocParallel::SerialParam()
    ),
    "Split columns in 'split_file' must exactly match names\\(pred_sigs\\) when 'splits = TRUE'\\."
  )
})

test_that("sig_eval_table respects a custom split perturbation column", {
  source_sigs <- list(
    geneA = c("g1", "g2"),
    geneB = c("g2", "g3")
  )
  pred_sigs <- list(
    split_1 = list(
      geneA = c("g1", "g2"),
      geneB = c("g2", "g3")
    )
  )
  true_sigs <- list(
    geneA = list(up = c("g1", "g2"), up_full = c("g1", "g2", "g3")),
    geneB = list(up = c("g2", "g3"), up_full = c("g2", "g3", "g1"))
  )

  split_tbl <- data.frame(
    gene = c("geneA", "geneB"),
    split_1 = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  split_file <- tempfile(fileext = ".csv")
  utils::write.csv(split_tbl, split_file, row.names = FALSE)

  eval_df <- sig_eval_table(
    source_sigs = source_sigs,
    pred_sigs = pred_sigs,
    true_sigs = true_sigs,
    splits = TRUE,
    split_file = split_file,
    split_pb_col = "gene",
    split_type = "10th",
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_equal(eval_df$gene, "geneA")
})
