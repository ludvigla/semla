library(testthat)
library(dplyr)

# Test if the function returns a matrix
test_that("LoadAndMergeMatrices returns a matrix with correct dimensions", {
  paths <- c(system.file("extdata/mousebrain/filtered_feature_bc_matrix.h5", package = "semla"),
             system.file("extdata/mousecolon/filtered_feature_bc_matrix.h5", package = "semla"))
  
  mat <- LoadAndMergeMatrices(paths, verbose = FALSE)
  
  expect_true(inherits(mat, what = "dgCMatrix"))
  expect_true(all(dim(mat) == c(188, 5164)))
})

# Test if the function throws an error when provided with non-existing file paths
test_that("LoadAndMergeMatrices throws an error when provided with non-existing file paths", {
  paths <- c(system.file("extdata/mousebrain/filtered_feature_bc_matrix.h5", package = "semla"),
             "path/does/not/exist.h5")
  
  expect_error(LoadAndMergeMatrices(paths, verbose = FALSE), "Invalid path")
})

# Test if the function throws an error when provided with paths to non-10x matrix files
test_that("LoadAndMergeMatrices throws an error when provided with paths to non-10x matrix files", {
  paths <- c(system.file("extdata/mousebrain/filtered_feature_bc_matrix.h5", package = "semla"),
             system.file("extdata/mousebrain/spatial/tissue_positions_list.csv", package = "semla"))
  
  expect_error(LoadAndMergeMatrices(paths, verbose = FALSE), 
               "Invalid file format")
})
