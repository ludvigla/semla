library(testthat)
library(semla)

# Test the centroid_angles_plot function
test_that("centroid_angles_plot returns a ggplot object", {
  # Call the centroid_angles_plot function
  plot <- centroid_angles_plot(nbreaks = 9, centroid_size = 8)
  
  # Check if the returned object is a ggplot object
  expect_s3_class(plot, "ggplot")
})

# Test the AnglePlot function
test_that("AnglePlot returns a patchwork object", {
  
  # Create a Seurat object for testing (replace with your own test data)
  se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
  
  # Call the AnglePlot function
  se_mcolon <- DisconnectRegions(se_mcolon, column_name = "selection", selected_groups = "GALT")
  plot <- AnglePlot(se_mcolon, column_name = "GALT_split", selected_group = "S1_region1")
  
  # Check if the returned object is a patchwork object
  expect_s3_class(plot, "patchwork")
  
  se_mcolon <- LoadImages(se_mcolon)
  plot <- AnglePlot(se_mcolon, column_name = "GALT_split", selected_group = "S1_region1", image_use = "raw")
  
  # Check if the returned object is a patchwork object
  expect_s3_class(plot, "patchwork")
  
  plot <- AnglePlot(se_mcolon, column_name = "GALT_split", selected_group = "S1_region1", crop_area = c(0.4, 0.4, 0.6, 0.6))
  
  # Check if the returned object is a patchwork object
  expect_s3_class(plot, "patchwork")
  
  plot <- AnglePlot(se_mcolon, column_name = "GALT_split", selected_group = "S1_region1", override_plot_dims = TRUE)
  
  # Check if the returned object is a patchwork object
  expect_s3_class(plot, "patchwork")
})