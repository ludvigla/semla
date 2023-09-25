library(testthat)
library(semla)

se_mbrain <- readRDS(system.file("extdata/mousebrain",
                                 "se_mbrain",
                                 package = "semla")) |> LoadImages()

test_that("ExportDataForViewer exports data correctly", {
  
  # Set up a temporary directory for testing
  test_dir <- tempdir()
  
  # Test exporting data with default arguments
  outpath <- ExportDataForViewer(se_mbrain, outdir = test_dir, overwrite = TRUE)
  
  # Check if the output path exists
  expect_true(dir.exists(outpath), info = "Output directory exists")
  
  # Check if files were created in the output directory
  expect_true(dir.exists(file.path(outpath, "tiles1")))
  expect_true(file.exists(file.path(outpath, "image_info_1.json")))
  expect_true(file.exists(file.path(outpath, "coords_Visium_1.json")))
})

test_that("ExportDataForViewer returns errors correctly", {
  
  # Set up a temporary directory for testing
  test_dir <- tempdir()
  
  # Test exporting data with wrong sampleID
  expect_error(ExportDataForViewer(se_mbrain, sampleIDs = 2, outdir = test_dir, overwrite = TRUE))
  
  # Invalid image
  se_mbrain@tools$Staffli@imgs <- "invalid"
  expect_error(ExportDataForViewer(se_mbrain, outdir = test_dir, overwrite = TRUE))
  
  # Invalid outdir
  expect_error(ExportDataForViewer(se_mbrain, outdir = "__invalid", overwrite = TRUE))
  expect_error(ExportDataForViewer(se_mbrain, outdir = 1L, overwrite = TRUE))
  
})
