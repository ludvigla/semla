library(testthat)

# Load a dummy Seurat object for testing
se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
se_merged <- MergeSTData(se_mbrain, se_mcolon) |> 
  ScaleData() |> 
  LoadImages() |> 
  RigidTransformImages(transforms = generate_rigid_transform(sampleID = 2, angle = 45))

# Define a test case
test_that("MapLabels works as expected", {
  
  # Test basic functionality
  expect_s3_class(MapLabels(se_merged, column_name = "selection"), "patchwork")
  
  # Test image_use parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", image_use = "raw"), "patchwork")
  
  # Test coords_use parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", coords_use = "transformed"), "patchwork")
  
  # Test crop_area parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", crop_area = c(0.1, 0.1, 0.9, 0.9)), "patchwork")
  
  # Test pt_size parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", pt_size = 2), "patchwork")
  
  # Test pt_alpha parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", pt_alpha = 0.5), "patchwork")
  
  # Test pt_stroke parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", pt_stroke = 0.5), "patchwork")
  
  # Test scale_alpha parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", scale_alpha = TRUE), "patchwork")
  
  # Test section_number parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", section_number = 2), "patchwork")
  
  # Test label_by parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", label_by = "sample_id"), "patchwork")
  p <- MapLabels(se_merged, column_name = "selection", label_by = "sample_id")
  lbls <- c(p$labels$title, p$patches$plots[[1]]$labels$title)
  expect_true(all(lbls == c("mousecolon", "mousebrain")))
  
  # Test ncol parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", ncol = 2), "patchwork")
  
  # Test scalebar parameters
  expect_s3_class(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE, scalebar_gg = scalebar(x = 1000)), "patchwork")
  
  # Test colors parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", colors = c("red", "blue")), "patchwork")
  
  # Test that the correct samples and column_name are plotted
  p <- MapLabels(se_merged, column_name = "selection", return_plot_list = TRUE)
  
  # Test that the list contains correct classes
  expect_true(all(sapply(p, function(p) class(p)[1]) == "gg"))
  
  # Test that the correct column_name are plotted
  expect_true(all(c("1", "2") == names(p)))
})


# Test argument validation
test_that("MapLabels throws appropriate errors for invalid arguments", {
  # Test invalid feature
  expect_error(MapLabels(se_merged, column_name = "not_a_label"))
  
  # Test invalid image_use
  expect_error(MapLabels(se_merged, column_name = "selection", image_use = "test"))
  
  # Test with invalid section_number
  expect_error(MapLabels(se_merged, column_name = "selection", section_number = 3),
               "is out of range")
  
  # Test with invalid scalebar_height
  expect_error(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE,
                           scalebar_height = 2),
               "Expected a numeric of length 1 between 0 and 1")
  
  # Test with invalid scalebar_position
  expect_error(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE,
                           scalebar_position = c(1.5, 0.5)),
               "Expected values between 0 and 1")
  
  # Test with invalid scalebar_gg
  expect_error(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE,
                           scalebar_gg = NULL),
               "Expected a 'ggplot' object")
  
  # test with invalid crop_area
  expect_error(MapLabels(se_merged, column_name = "selection", crop_area = c(-1, 0.3, 0.6, 0.6)),
               "'crop_area' can only take values between 0-1")
  expect_error(MapLabels(se_merged, column_name = "selection", crop_area = c(0.3, 0.3, 0.3, 0.6)),
               "'left' value needs to be lower that 'right' value")
  expect_error(MapLabels(se_merged, column_name = "selection", crop_area = c(0.3, 0.3, 0.6, 0.3)),
               "'top' value needs to be lower that 'bottom' value")
  expect_error(MapLabels(se_merged, column_name = "selection", crop_area = c(0.3, 0.3, 0.3)),
               "Invalid length for 'crop_area', expected a 'numeric' vector of length 4")
  expect_error(MapLabels(se_merged, column_name = "selection", crop_area = TRUE),
               "Invalid class 'logical' for 'crop_area', expected 'numeric'")
})
