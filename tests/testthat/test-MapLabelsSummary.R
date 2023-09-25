library(testthat)
library(semla)

se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla")) |> 
  LoadImages()

# Test Case 1: Test if MapFeaturesSummary throws an error when subplot_type is missing
test_that("MapFeaturesSummary returns a ggplot object", {
  p <- MapLabelsSummary(se_mcolon, column_name = "orig.ident")
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapLabelsSummary(se_mcolon, column_name = "orig.ident", bar_display = "percent")
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapLabelsSummary(se_mcolon, column_name = "orig.ident", bar_display = "count")
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapLabelsSummary(se_mcolon, column_name = "orig.ident", bar_width = 2)
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapLabelsSummary(se_mcolon, column_name = "orig.ident", image_use = "raw")
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapLabelsSummary(se_mcolon, column_name = "orig.ident", hide_legend = TRUE)
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapLabelsSummary(se_mcolon, column_name = "orig.ident", label_by = "sample_id")
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapLabelsSummary(se_mcolon, column_name = "orig.ident", colors = c("red", "blue"))
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapLabelsSummary(se_mcolon, column_name = "orig.ident", crop_area = c(0.4, 0.4, 0.5, 0.5))
  expect_identical(class(p)[1], "patchwork")
  
  p <- MapLabelsSummary(se_mcolon, column_name = "orig.ident", section_number = 1)
  expect_identical(class(p)[1], "patchwork")
})
