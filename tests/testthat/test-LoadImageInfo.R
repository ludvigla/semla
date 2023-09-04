library(testthat)
library(dplyr)
library(semla)

imgfiles <-
  c(system.file("extdata/mousebrain/spatial",
                "tissue_lowres_image.jpg",
                package = "semla"),
    system.file("extdata/mousecolon/spatial",
                "tissue_lowres_image.jpg",
                package = "semla"))
img_info <- LoadImageInfo(imgfiles)

test_that("LoadImageInfo returns a tibble", {
  expect_s3_class(img_info, "tbl")
})

test_that("LoadImageInfo returns correct column names", {
  expect_equal(colnames(img_info), c("format", "width", "height", "colorspace", "matte", "filesize", "density", "sampleID", "type"))
})

test_that("LoadImageInfo returns correct column classes", {
  expect_equal(sapply(img_info, class) |> as.character(), c("character","integer","integer", "character","logical","integer", "character", "character", "character"))
})

test_that("LoadImageInfo returns correct number of rows names", {
  expect_equal(nrow(img_info), 2)
})

test_that("LoadImageInfo returns correct number of sampleIDs", {
  expect_equal(table(img_info$sampleID) |> length(), 2)
})
