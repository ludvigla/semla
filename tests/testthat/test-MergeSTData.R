library(testthat)
library(semla)

# Load Seurat objects from semla package for testing
se_mbrain <- readRDS(system.file("extdata", "mousebrain/se_mbrain", package = "semla"))
se_mcolon <- readRDS(system.file("extdata", "mousecolon/se_mcolon", package = "semla"))

test_that("MergeSTData works as expected", {
  # Merge Seurat objects
  se_merged <- MergeSTData(x = se_mbrain, y = se_mcolon)
  
  # Check that the output is a Seurat object
  expect_s4_class(se_merged, "Seurat")
  
  # Check that the merged data contains cells from both input objects
  expect_equal(length(colnames(se_merged)), length(colnames(se_mbrain)) + length(colnames(se_mcolon)))
  
  # Check that the Staffli object is present and correctly merged
  merged_staffli <- GetStaffli(se_merged)
  expect_equal(length(unique(merged_staffli@meta_data$barcode)), length(colnames(se_merged)))
})
