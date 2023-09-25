library(testthat)
library(semla)

# Load Seurat objects from semla package for testing
se_mbrain <- readRDS(system.file("extdata", "mousebrain/se_mbrain", package = "semla"))

test_that("SubsetSTData works as expected", {
  # Subset Seurat object
  se_subset <- SubsetSTData(se_mbrain, spots = colnames(se_mbrain)[1:1000])
  
  # Check that the output is a Seurat object
  expect_s4_class(se_subset, "Seurat")
  
  # Check that the number of cells matches the expected number
  expect_equal(length(colnames(se_subset)), 1000)
  
  # Check that the subsetted cells are correct
  expect_equal(colnames(se_subset), colnames(se_mbrain)[1:1000])
  
  # Check that the Staffli object is present and correctly subsetted
  expect_equal(length(unique(GetStaffli(se_subset)@meta_data$barcode)), 1000)
  
  # Subset Seurat object
  se_subset <- SubsetSTData(se_mbrain, expression = nFeature_Spatial > 60, features = rownames(se_mbrain)[1:100])
  
  expect_equal(ncol(se_subset), 686)
  expect_equal(nrow(se_subset), 100)
})
