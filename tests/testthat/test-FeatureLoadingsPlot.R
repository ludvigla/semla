library(testthat)
library(semla)

se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
se_mcolon <- se_mcolon |> ScaleData() |> RunPCA()

test_that("PlotFeatureLoadings function returns valid plots and handles errors", {
  
  # Test case 1: Check if the function returns a ggplot object
  expect_s3_class(PlotFeatureLoadings(se_mcolon), "ggplot")
  
  # Test case 2: Check if the function returns a ggplot object for different modes
  expect_s3_class(PlotFeatureLoadings(se_mcolon, mode = "barplot"), "ggplot")
  expect_s3_class(PlotFeatureLoadings(se_mcolon, mode = "dotplot"), "ggplot")
  expect_s3_class(PlotFeatureLoadings(se_mcolon, mode = "heatmap", dims = 1:5), "ggplot")
  
  # Test case 3: Check if the function handles invalid mode gracefully
  expect_error(PlotFeatureLoadings(se_mcolon, mode = "invalid_mode"), class = "error")
  
  # Test case 4: Check if the function handles invalid type gracefully
  expect_error(PlotFeatureLoadings(se_mcolon, type = "invalid_type"), class = "error")
  
  # Test case 5: Check if the function handles invalid dims input gracefully
  expect_error(PlotFeatureLoadings(se_mcolon, dims = "invalid_dims"), class = "error")
  
  # Test case 6: Check if the function handles invalid nfeatures input gracefully
  expect_error(PlotFeatureLoadings(se_mcolon, nfeatures = 0), class = "error")
  
  # Test case 7: Check if the function handles invalid fill input gracefully
  expect_error(PlotFeatureLoadings(se_mcolon, fill = 123), class = "error")
  
  # Test case 8: Check if the function handles invalid color input gracefully
  expect_error(PlotFeatureLoadings(se_mcolon, color = 123), class = "error")
  
  # Test case 9: Check if the function handles invalid bar_width input gracefully
  expect_error(PlotFeatureLoadings(se_mcolon, bar_width = "invalid_width"), class = "error")
  
  # Test case 10: Check if the function handles invalid pt_size input gracefully
  expect_error(PlotFeatureLoadings(se_mcolon, pt_size = "invalid_size"), class = "error")
  
  # Test case 11: Check if the function handles invalid pt_stroke input gracefully
  expect_error(PlotFeatureLoadings(se_mcolon, pt_stroke = "invalid_stroke"), class = "error")
  
  # Test case 12: Check if the function handles invalid color_by_loadings input gracefully
  expect_error(PlotFeatureLoadings(se_mcolon, color_by_loadings = "invalid_color"), class = "error")
  
  # Test case 13: Check if the function handles invalid gradient_colors input gracefully
  expect_error(PlotFeatureLoadings(se_mcolon, gradient_colors = 123), class = "error")
  
  # More checks
  expect_s3_class(PlotFeatureLoadings(se_mcolon, mode = "barplot", dims = 1, color_by_loadings = TRUE), "ggplot")
  
  expect_s3_class(PlotFeatureLoadings(se_mcolon, mode = "barplot", dims = 1, type = "centered"), "ggplot")
  
  expect_s3_class(PlotFeatureLoadings(se_mcolon, mode = "barplot", dims = 1, type = "positive"), "ggplot")
  
  expect_s3_class(PlotFeatureLoadings(se_mcolon, mode = "barplot", dims = 1, type = "negative"), "ggplot")
})
  