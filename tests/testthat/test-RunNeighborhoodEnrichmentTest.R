library(testthat)
library(semla)

# Define test data
se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))

se_mcolon <- se_mcolon |>
  ScaleData(verbose = FALSE) |>
  RunPCA(verbose = FALSE) |>
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) |>
  FindClusters(verbose = FALSE)

# Test
test_that("RadialDistance function returns valid results", {
  
  # Run Label Assortativity Analysis
  res <- RunNeighborhoodEnrichmentTest(object = se_mcolon,
                                   column_name = "seurat_clusters",
                                   n_permutations = 2,
                                   nCores = 1)
  
  expect_s3_class(res, "tbl_df")
})
