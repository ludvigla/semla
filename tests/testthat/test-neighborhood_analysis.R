# Load required libraries
library(testthat)
library(dplyr)
library(tidyr)
library(semla)

# Define test data
se_mbrain <- readRDS(system.file("extdata", "/mousebrain/se_mbrain", package = "semla"))

se_mbrain <- se_mbrain %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(verbose = FALSE) %>%
  FindClusters(verbose = FALSE)

# Test cases
test_that("RegionNeighbors should work with default parameters", {
  res <- RegionNeighbors(se_mbrain, column_name = "seurat_clusters", column_labels = "10")
  expect_true("nb_to_10" %in% colnames(res@meta.data))
})

test_that("RegionNeighbors should work with different modes", {
  res_outer <- RegionNeighbors(se_mbrain, column_name = "seurat_clusters", column_labels = "10", mode = "outer")
  res_inner <- RegionNeighbors(se_mbrain, column_name = "seurat_clusters", column_labels = "10", mode = "inner")
  res_inner_outer <- RegionNeighbors(se_mbrain, column_name = "seurat_clusters", column_labels = "10", mode = "inner_outer")
  res_all_inner_outer <- RegionNeighbors(se_mbrain, column_name = "seurat_clusters", column_labels = "10", mode = "all_inner_outer")
  
  expect_true("nb_to_10" %in% colnames(res_outer@meta.data))
  expect_true("inner_border_10" %in% colnames(res_inner@meta.data))
  expect_true("nb_to_10" %in% colnames(res_inner_outer@meta.data))
  expect_true("nb_to_10" %in% colnames(res_all_inner_outer@meta.data))
})

test_that("RegionNeighbors should work with custom column_key", {
  res <- RegionNeighbors(se_mbrain, column_name = "seurat_clusters", column_labels = "10", column_key = "custom_key_")
  expect_true("custom_key_10" %in% colnames(res@meta.data))
})

test_that("RegionNeighbors should work with multiple column_labels", {
  res <- RegionNeighbors(se_mbrain, column_name = "seurat_clusters", column_labels = c("8", "10"))
  expect_true("nb_to_8" %in% colnames(res@meta.data))
  expect_true("nb_to_10" %in% colnames(res@meta.data))
})

test_that("RegionNeighbors should work with additional parameters for GetSpatialNetwork", {
  res <- RegionNeighbors(se_mbrain, column_name = "seurat_clusters", column_labels = "10", nNeighbors = 40, maxDist = Inf)
  expect_true("nb_to_10" %in% colnames(res@meta.data))
})

test_that("RegionNeighbors should return an error with invalid column_name", {
  expect_error(RegionNeighbors(se_mbrain, column_name = "non_existent_column", column_labels = "10"))
})

test_that("RegionNeighbors should return an error with invalid column_labels", {
  expect_error(RegionNeighbors(se_mbrain, column_name = "seurat_clusters", column_labels = "non_existent_label"))
})
 
