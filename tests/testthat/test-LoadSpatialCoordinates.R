library(testthat)
library(dplyr)
library(semla)

coordinatefiles <-
  c(system.file("extdata/mousebrain/spatial",
                "tissue_positions_list.csv",
                package = "semla"),
    system.file("extdata/mousecolon/spatial",
                "tissue_positions_list.csv",
                package = "semla"))
coordinates <- LoadSpatialCoordinates(coordinatefiles)

test_that("LoadSpatialCoordinates returns a tibble", {
  expect_s3_class(coordinates, "tbl")
})

test_that("LoadSpatialCoordinates returns correct column names", {
  expect_equal(colnames(coordinates), c("barcode", "selected", "y", "x", "pxl_row_in_fullres", "pxl_col_in_fullres", "sampleID"))
})

test_that("LoadSpatialCoordinates returns correct column classes", {
  expect_equal(sapply(coordinates, class) |> as.character(), c("character","integer","integer","integer","integer","integer","integer"))
})

test_that("LoadSpatialCoordinates returns correct number of rows names", {
  expect_equal(nrow(coordinates), 5164)
})

test_that("LoadSpatialCoordinates returns correct number of sampleIDs", {
  expect_equal(table(coordinates$sampleID) |> length(), 2)
})
