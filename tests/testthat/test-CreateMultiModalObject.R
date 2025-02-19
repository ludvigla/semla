library(testthat)
library(semla)

# Load Seurat objects from semla package for testing
se_mod1 <- readRDS(system.file("extdata", "mousebrain/se_mbrain", package = "semla"))
se_mod2 <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))

test_that("CreateMultiModalObject returns a correct object", {
  # Test function for aggregating function "mean"
  se_mmo <- CreateMultiModalObject(object_ref = se_mod1, object_map = se_mod2,
                                   agg_func = "mean",
                                   new_assay_name = "Modality2")
  
  ## Check that the output is a Seurat object
  expect_s4_class(se_mmo, "Seurat")
  
  # Test function for aggregating function "sum"
  se_mmo <- CreateMultiModalObject(object_ref = se_mod1, object_map = se_mod2,
                                   agg_func = "sum",
                                   new_assay_name = "Modality2")
  
  ## Check that the output is a Seurat object
  expect_s4_class(se_mmo, "Seurat")
})

test_that("CreateMultiModalObject returns an error with wrong input", {
  # Expected errors
  # expect_error(CreateMultiModalObject(object_ref = se_mod1, object_map = se_mod2, agg_func = "test"),
  #              '"arg" should be one of "mean", "sum"')
  
  expect_error(CreateMultiModalObject(object_ref = se_mod1, object_map = se_mod2, n_neighbors = 0), 
               "'n_neighbors' should be larger than 0")
})
