# test seeds
hnsc_seed <- c("MTOR","MUC4","TP63","TSC2","MUC16","PTPRD","PTPRT","NOTCH1","TGFBR2")
brca_seed <- c("RB1","AKT1","CDH1","FLNA","IRS4","TP53","CCND1","EP300","FBLN2","GATA3")

# Create toy igraph with vertex names and edge weights
all_genes <- union(hnsc_seed, brca_seed)
tcga_brca_wgcna <- igraph::make_full_graph(length(all_genes))
igraph::V(tcga_brca_wgcna)$name <- all_genes
igraph::E(tcga_brca_wgcna)$weight <- runif(igraph::ecount(tcga_brca_wgcna), min = 0.1, max = 1)

combined_seed <- list(hnsc=hnsc_seed, brca=brca_seed)

# Seed Matrix Function Tests

test_that("Seed Matrix Function Works", {
  # Single Gene
  result <- seed_matrix(tcga_brca_wgcna, "TP53")
  expect_equal(result["TP53",], 1)
  expect_true(all(result[rownames(result) != "TP53",] == 0))

  # Named List of Gene (One geneset)
  result <- seed_matrix(tcga_brca_wgcna, list(hnsc=hnsc_seed))
  expect_true(all(result[hnsc_seed,] == 1))
  expect_true(all(result[!(rownames(result) %in% hnsc_seed),] == 0))

  # List of Named List of genes (List of genesets)
  hnsc_atomized <- as.list(hnsc_seed)
  names(hnsc_atomized) <- hnsc_atomized
  result <- seed_matrix(tcga_brca_wgcna, hnsc_atomized)
  expect_true(all(result[!(rownames(result) %in% hnsc_seed),] == 0))
  expect_true(all(result[hnsc_seed,] == diag(x=1, nrow=length(hnsc_seed))))

  # List of Named List of of more diverse genes (List of genesets)
  result <- seed_matrix(tcga_brca_wgcna, combined_seed)
  expect_true(all(result[hnsc_seed,"hnsc"] == 1))
  expect_false(all(result[hnsc_seed,"brca"] == 1))
  expect_true(all(result[brca_seed,"brca"] == 1))
  expect_false(all(result[brca_seed,"hnsc"] == 1))
})

# Random Walk Tests

test_that("Random Walk Components Work", {
  # Test mat describes a three node graph with edges (1,2), (2,3).
  test_mat <- matrix(c(0,1,0,1,0,1,0,1,0), nrow=3, ncol=3, byrow=TRUE)
  e1 <- matrix(c(1,0,0), nrow=3, ncol=1)
  e2 <- matrix(c(0,1,0), nrow=3, ncol=1)
  e3 <- matrix(c(0,0,1), nrow=3, ncol=1)
  norm_test_mat <- Matrix::Diagonal(x=Matrix::rowSums(test_mat)^(-1)) %*% test_mat

  # Checking Normalization
  expect_true(all(norm_test_mat == matrix(c(0,1,0,0.5,0,0.5,0,1,0), nrow=3, ncol=3, byrow=TRUE)))

  # Checking Matrix Multiplication
  expect_true(all(t(t(e1) %*% norm_test_mat) == e2))
  expect_true(all(t(t(e2) %*% norm_test_mat) == c(0.5,0,0.5)))
  expect_true(all(t(t(e3) %*% norm_test_mat) == e2))
  expect_true(all(t(cbind(e1,e2,e3)) %*% norm_test_mat == t(cbind(e2,c(0.5,0,0.5),e2))))

  # Checking Random Walk Function

  ## If Restart value is 1, you stay at seeds, should converge.
  expect_true(all(random_walk(graph_from_adjacency_matrix(test_mat, mode="undirected"), cbind(e1, e2, e3), normalize="row", restart=1) == cbind(e1, e2, e3)))

  ## For seed_gene matrix, if you have restart value of 1, you stay at seeds, should converge.
  result <- random_walk(tcga_brca_wgcna, seed_mat = seed_matrix(tcga_brca_wgcna, combined_seed), normalize="row", restart=1)
  expect_true(all(result[hnsc_seed,"hnsc"] == 1/length(hnsc_seed)))
  expect_true(all(result[brca_seed,"brca"] == 1/length(brca_seed)))

})

