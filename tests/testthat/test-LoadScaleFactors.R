library(testthat)
library(semla)

jsonfiles <-
  c(system.file("extdata/mousebrain/spatial",
                "scalefactors_json.json",
                package = "semla"),
    system.file("extdata/mousecolon/spatial",
                "scalefactors_json.json",
                package = "semla"))
scalefactors <- LoadScaleFactors(jsonfiles)

test_that("LoadScaleFactors returns a tibble", {
  expect_s3_class(scalefactors, "tbl")
})

test_that("LoadScaleFactors returns correct column names", {
  expect_equal(colnames(scalefactors), c("spot_diameter_fullres","tissue_hires_scalef","fiducial_diameter_fullres", "tissue_lowres_scalef", "sampleID"))
})

test_that("LoadScaleFactors returns correct column classes", {
  expect_equal(sapply(scalefactors, class) |> as.character(), c("numeric", "numeric", "numeric", "numeric", "character"))
})

test_that("LoadScaleFactors returns correct number of rows names", {
  expect_equal(nrow(scalefactors), 2)
})

test_that("LoadScaleFactors returns correct number of sampleIDs", {
  expect_equal(table(scalefactors$sampleID) |> length(), 2)
})
