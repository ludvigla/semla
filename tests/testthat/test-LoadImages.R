library(testthat)
library(semla)

# Define test
test_that("LoadImages function loads images and scales them appropriately (default)", {
  # Create example input
  imgs <- c(system.file("extdata/mousebrain", "spatial/tissue_lowres_image.jpg", package = "semla"),
            system.file("extdata/mousecolon", "spatial/tissue_lowres_image.jpg", package = "semla"))
  
  # Call function
  rsts <- LoadImages(imgs, image_height = 300)
  heights <- sapply(rsts, nrow)
  
  # Check output
  expect_true(all(heights == 300))
  expect_true(length(rsts) == 2)
  expect_true(all(sapply(rsts, class) == "raster"))
})

# Define test
test_that("LoadImages function loads images and scales them appropriately (Seurat)", {
  # Create example input
  se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
  se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
  se_merged <- MergeSTData(se_mbrain, se_mcolon)
  
  # Call function
  se_merged <- LoadImages(se_merged, image_height = 300)
  st_object <- GetStaffli(se_merged)
  imgs <- st_object@rasterlists$raw
  heights <- sapply(imgs, nrow)
  
  # Check output
  expect_true(all(heights == 300))
  expect_true(length(imgs) == 2)
  expect_true(all(sapply(imgs, class) == "raster"))
})
