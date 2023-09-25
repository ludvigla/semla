library(testthat)
library(semla)
library(magick)
library(dplyr)

# Load example data
lowres_image_file <- system.file("extdata/mousebrain/spatial",
                                 "tissue_lowres_image.jpg",
                                 package = "semla")
im <- image_read(lowres_image_file)

coordinates_file <- system.file("extdata/mousebrain/spatial",
                                "tissue_positions_list.csv",
                                package = "semla")
xy <- LoadSpatialCoordinates(coordinatefiles = coordinates_file)

scalefactors <- system.file("extdata/mousebrain/spatial",
                            "scalefactors_json.json",
                            package = "semla") |> 
  jsonlite::read_json()

xy <- xy |> 
  mutate(across(all_of(c("pxl_row_in_fullres",
                         "pxl_col_in_fullres")), ~.x*scalefactors$tissue_lowres_scalef))

se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))

test_that("MaskImages.default returns a masked image with the correct class", {
  im_masked <- MaskImages(im, xy_coords = xy, verbose = FALSE)
  expect_true(inherits(im_masked, "magick-image"))
})

test_that("MaskImages.default throws an error when given an incorrect input class", {
  expect_error(MaskImages.default("Invalid class", xy_coords = xy, verbose = FALSE),
               "Invalid class 'character', expected a 'magick-image' object")
})

test_that("MaskImages.Seurat returns a Seurat object with masked images", {
  se_mbrain <- LoadImages(se_mbrain)
  se_mbrain_masked <- MaskImages(se_mbrain, verbose = FALSE)
  expect_true(inherits(se_mbrain_masked, "Seurat"))
  expect_true(!is.null(se_mbrain_masked@tools$Staffli@rasterlists[["raw"]]))
})

test_that("MaskImages.Seurat throws an error when given an incorrect section number", {
  se_mbrain <- LoadImages(se_mbrain)
  expect_error(MaskImages(se_mbrain, section_numbers = 99, verbose = FALSE),
               "Invalid numbers for 'section_numbers'. Sections available:")
})

