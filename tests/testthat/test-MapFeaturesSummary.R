library(testthat)
library(semla)

se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))

# Test Case 1: Test if MapFeaturesSummary throws an error when subplot_type is missing
test_that("MapFeaturesSummary returns a ggplot object", {
  p <- MapFeaturesSummary(se_mcolon, features = "Nrgn", subplot_type = "violin")
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapFeaturesSummary(se_mcolon, features = "Nrgn", subplot_type = "box")
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapFeaturesSummary(se_mcolon, features = "Nrgn", subplot_type = "histogram")
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapFeaturesSummary(se_mcolon, features = "Nrgn", subplot_type = "density")
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapFeaturesSummary(se_mcolon, features = "Nrgn", subplot_type = "density", crop_area = c(0.4, 0.4, 0.5, 0.5))
  expect_identical(class(p)[1], "patchwork")
})
