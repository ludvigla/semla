library(testthat)
library(semla)

# Load Seurat objects from semla package for testing
se_mbrain <- readRDS(system.file("extdata", "mousebrain/se_mbrain", package = "semla"))
se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
se <- MergeSTData(se_mbrain, se_mcolon)

test_that("UpdateSeuratFromSemla works as expected", {
  # Update object for one section raw image
  se_upd <- UpdateSeuratFromSemla(se_mbrain, image_use = "raw", verbose = TRUE)
  
  ## Check that the output is a Seurat object
  expect_s4_class(se_upd, "Seurat")
  
  ## Check that we can plot the image
  # expect_s3_class(SpatialFeaturePlot(se_upd, features = "Nrgn"), "patchwork")
  
  # Update object for one section transformed image
  se_mbrain <- LoadImages(se_mbrain)
  rotation_angle <- get_array_rotation(se_mbrain, grid_pattern = "hexagonal")
  transforms <- generate_rigid_transform(sampleID = 1L, angle = 10)
  se_mbrain <- RigidTransformImages(se_mbrain, transforms = transforms)
  
  se_upd <- UpdateSeuratFromSemla(se_mbrain, image_use = "transformed", verbose = TRUE)
  
  ## Check that the output is a Seurat object
  expect_s4_class(se_upd, "Seurat")
  
  # Update object for multiple sections raw image
  se_upd <- UpdateSeuratFromSemla(se, image_use = "raw", verbose = TRUE)
  
  ## Check that the output is a Seurat object
  expect_s4_class(se_upd, "Seurat")
  
  # Update object for multiple sections transformed image
  se <- LoadImages(se)
  rotation_angle <- get_array_rotation(se, grid_pattern = "hexagonal")
  transforms <- rbind(generate_rigid_transform(sampleID = 1L, angle = rotation_angle$`1`),
                      generate_rigid_transform(sampleID = 2L, angle = 12))
  se <- RigidTransformImages(se, transforms = transforms)
  
  se_upd <- UpdateSeuratFromSemla(se, image_use = "transformed", verbose = TRUE)
  
  ## Check that the output is a Seurat object
  expect_s4_class(se_upd, "Seurat")
  
  # ## Check that the output is a Seurat object
  # expect_s4_class(se_upd, "Seurat")
  # 
  # ## Check that we can plot the image
  # expect_s3_class(SpatialFeaturePlot(se_upd, features = "Nrgn"), "patchwork")
  # 
  # # Update object for multiple sections 
  # se_upd <- UpdateSeuratFromSemla(se_upd, image_use = "raw", verbose = TRUE)
  # 
  # ## Check that the output is a Seurat object
  # expect_s4_class(se_upd, "Seurat")
})