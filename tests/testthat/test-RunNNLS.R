library(testthat)
library(semla)
library(RcppML)

se_mbrain <- readRDS(system.file("extdata/mousebrain",
                                 "se_mbrain",
                                 package = "semla"))

# Random clusters
se_mbrain$clusters <- sample(c("A", "B", "C"), size = ncol(se_mbrain), replace = TRUE)
se_mbrain <- SetIdent(se_mbrain, value = "clusters")

test_that("RunNNLS returns valid results", {
  
  # Test with default arguments
  res1 <- suppressWarnings({RunNNLS(se_mbrain, singlecell_object = se_mbrain, groups = "clusters", 
                  singlecell_assay = "Spatial", spatial_assay = "Spatial")})
  expect_s4_class(res1, "Seurat")
  expect_true("celltypeprops" %in% names(res1@assays))
  
  # Test with return_as_dimred = TRUE
  res2 <- suppressWarnings({RunNNLS(se_mbrain, singlecell_object = se_mbrain, groups = "clusters", 
                  singlecell_assay = "Spatial", spatial_assay = "Spatial", return_as_dimred = TRUE)})
  expect_s4_class(res2, "Seurat")
  expect_true("nnls" %in% names(res2@reductions))
  
  # Test with idents
  res3 <- suppressWarnings({RunNNLS(se_mbrain, singlecell_object = se_mbrain,  
                                    singlecell_assay = "Spatial", spatial_assay = "Spatial")})
  expect_s4_class(res3, "Seurat")
  expect_true("celltypeprops" %in% names(res3@assays))
  
  # Test with rare cell type
  se_mbrain$clusters[1] <- "D"
  res4 <- suppressWarnings({RunNNLS(se_mbrain, singlecell_object = se_mbrain, groups = "clusters", 
          singlecell_assay = "Spatial", spatial_assay = "Spatial")})
  expect_s4_class(res4, "Seurat")
  expect_true("celltypeprops" %in% names(res4@assays))
  
  # Test with negative values
  mat1 <- LayerData(se_mbrain, layer = "data")
  mat1[1, 1] <- -1
  res5 <- suppressWarnings({RunNNLS(mat1, singlecell_matrix = mat1, groups = se_mbrain$clusters)})
  expect_type(res5, "double")
})

test_that("RunNNLS throws expected errors", {
  
  expect_error(RunNNLS("invalid"))
  expect_error(RunNNLS(se_mbrain, singlecell_object = "invalid"))
})
