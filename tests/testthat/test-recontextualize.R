test_that("recontextualize dispatches to netProp and projectCor", {
  counts <- matrix(
    c(50, 52, 10, 11, 12, 13,
      48, 51, 9, 10, 11, 12,
      8, 9, 45, 47, 10, 11,
      7, 8, 44, 46, 9, 10,
      15, 16, 15, 16, 30, 31),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(paste0("g", 1:5), paste0("s", 1:6))
  )
  se_net <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts))
  seeds <- list(sig1 = c("g1", "g2"))

  norm_counts <- log1p(sweep(counts, 2, colSums(counts), "/"))
  se_norm <- SummarizedExperiment::SummarizedExperiment(assays = list(log_norm = norm_counts))
  var_genes <- sigrecon:::rank.var.eset(se_norm)$values[1:4]
  test_ig <- wgcna.adj(
    mat = t(norm_counts[var_genes, , drop = FALSE]),
    beta = 1,
    cor.type = "signed hybrid",
    igraph = TRUE,
    diag_zero = TRUE
  )

  expect_equal(
    recontextualize(
      method = "networkProp",
      se = se_net,
      seeds = seeds,
      sig = "corr",
      limit = 2,
      nfeatures = 4,
      beta = 1,
      cor.type = "signed hybrid"
    ),
    network_sig(
      ig = test_ig,
      seeds = seeds,
      sig = "corr",
      limit = 2
    )
  )

  expr <- matrix(c(10, 12, 9, 11,
                   8, 9, 7, 8,
                   2, 3, 9, 8,
                   3, 2, 8, 9),
                 nrow = 4, byrow = TRUE,
                 dimnames = list(paste0("g", 1:4), paste0("s", 1:4)))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = expr))
  sigs <- list(sigA = c("g1", "g2"))

  expect_equal(
    recontextualize(
      method = "projectCor",
      se = se,
      sigs = sigs,
      score = "eigen"
    ),
    projectCor(
      se = se,
      sigs = sigs,
      score = "eigen"
    )
  )
})

test_that("recontextualize accepts a pre-built ig for method = 'networkProp'", {
  test_mat <- matrix(c(0,1,0,1,0,1,0,1,0), nrow=3, ncol=3, byrow=TRUE)
  test_ig <- igraph::graph_from_adjacency_matrix(test_mat, mode="undirected")
  igraph::V(test_ig)$name <- c("g1", "g2", "g3")

  seeds <- list(keep = c("g1", "g2"))

  expect_equal(
    recontextualize(method = "networkProp", ig = test_ig, seeds = seeds, sig = "corr", limit = 2),
    netProp(ig = test_ig, seeds = seeds, sig = "corr", limit = 2)
  )
})

test_that("recontextualize supports mean signatures with DESeq2 and batch variables", {
  skip_if_not_installed("DESeq2")

  counts <- matrix(
    c(
      20, 22, 19, 21, 240, 250, 235, 25, 24, 26,
      18, 20, 17, 19, 180, 175, 185, 22, 21, 23,
      15, 16, 14, 15, 18, 19, 17, 210, 220, 205,
      12, 11, 13, 12, 15, 14, 16, 170, 165, 175
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c("g1", "g2", "g3", "g4"),
      c("ctrl1", "ctrl2", "ctrl3", "ctrl4", "pertA1", "pertA2", "pertA3", "pertB1", "pertB2", "pertB3")
    )
  )

  coldata <- S4Vectors::DataFrame(
    perturbation = c("control", "control", "control", "control", "pertA", "pertA", "pertA", "pertB", "pertB", "pertB"),
    condition = c("Control", "Control", "Control", "Control", "Perturbed", "Perturbed", "Perturbed", "Perturbed", "Perturbed", "Perturbed"),
    batch = factor(c("b1", "b2", "b1", "b2", "b1", "b2", "b1", "b2", "b1", "b2")),
    row.names = colnames(counts)
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = coldata
  )

  result <- recontextualize(
    method = "mean",
    se = se,
    sigs = list(pertA = c("seed1", "seed2"), pertB = c("seed3", "seed4")),
    perturbation_col = "perturbation",
    condition_col = "condition",
    design_vars = "batch",
    alpha = 1
  )

  expect_identical(names(result), c("pertA", "pertB"))
  expect_equal(length(result$pertA), 2L)
  expect_equal(length(result$pertB), 2L)
  expect_true(setequal(result$pertA, c("g1", "g2")))
  expect_true(setequal(result$pertB, c("g3", "g4")))
})
