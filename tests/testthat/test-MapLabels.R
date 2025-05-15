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
  expect_s3_class(MapLabels(se_merged, column_name = "selection", shape = "tile"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", shape = "raster"), "patchwork")
  
  # Test image_use parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", image_use = "raw"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", image_use = "raw", shape = "tile"), "patchwork")
  
  # Test coords_use parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", coords_use = "transformed"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", coords_use = "transformed", shape = "tile"), "patchwork") # returns a plot but ignores the transformation, as it takes the xy coordiantes. if you plot with the image it works liek intended
  expect_s3_class(MapLabels(se_merged, column_name = "selection", coords_use = "transformed", shape = "raster"), "patchwork")
  
  # Test crop_area parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", crop_area = c(0.1, 0.1, 0.9, 0.9)), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", crop_area = c(0.1, 0.1, 0.9, 0.9), shape = "tile"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", crop_area = c(0.1, 0.1, 0.9, 0.9), shape = "raster"), "patchwork")
  
  # Test pt_size parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", pt_size = 2), "patchwork")
  
  # Test spot_side parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", spot_side = c(2, 4), shape = "tile"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", spot_side = c(200, 400), shape = "tile", image_use = "raw"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", spot_side = c(200, 400), shape = "tile", image_use = "transformed"), "patchwork")
  
  # Test pt_alpha parameter
  ## point
  expect_s3_class(MapLabels(se_merged, column_name = "selection", pt_alpha = 0.5), "patchwork")
  ## tile
  expect_s3_class(MapLabels(se_merged, column_name = "selection", pt_alpha = 0.5, shape = "tile"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", pt_alpha = 0.5, shape = "tile", image_use = "raw"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", pt_alpha = 0.5, shape = "tile", image_use = "transformed"), "patchwork")
  ## raster
  expect_s3_class(MapLabels(se_merged, column_name = "selection", pt_alpha = 0.5, shape = "raster"), "patchwork")
  
  # Test pt_stroke parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", pt_stroke = 0.5), "patchwork")
  
  # Test scale_alpha parameter
  ## point
  expect_s3_class(MapLabels(se_merged, column_name = "selection", scale_alpha = TRUE), "patchwork")
  ## tile
  expect_s3_class(MapLabels(se_merged, column_name = "selection", scale_alpha = TRUE, shape = "tile"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", scale_alpha = TRUE, shape = "tile", image_use = "raw"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", scale_alpha = TRUE, shape = "tile", image_use = "transformed"), "patchwork")
  ## raster
  expect_s3_class(MapLabels(se_merged, column_name = "selection", scale_alpha = TRUE, shape = "raster"), "patchwork")
  
  # Test section_number parameter
  ## point
  expect_s3_class(MapLabels(se_merged, column_name = "selection", section_number = 2), "patchwork")
  ## tile
  expect_s3_class(MapLabels(se_merged, column_name = "selection", section_number = 2, shape = "tile"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", section_number = 2, shape = "tile", image_use = "raw"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", section_number = 2, shape = "tile", image_use = "transformed"), "patchwork")
  ## raster
  expect_s3_class(MapLabels(se_merged, column_name = "selection", section_number = 2, shape = "raster"), "patchwork")
  
  # Test label_by parameter
  ## point
  p <- MapLabels(se_merged, column_name = "selection", label_by = "sample_id")
  expect_s3_class(p, "patchwork")
  lbls <- c(p$labels$title, p$patches$plots[[1]]$labels$title)
  expect_true(all(lbls == c("mousecolon", "mousebrain")))
  ## tile
  p <- MapLabels(se_merged, column_name = "selection", label_by = "sample_id", shape = "tile")
  expect_s3_class(p, "patchwork")
  lbls <- c(p$labels$title, p$patches$plots[[1]]$labels$title)
  expect_true(all(lbls == c("mousecolon", "mousebrain")))
  ## raster
  p <- MapLabels(se_merged, column_name = "selection", label_by = "sample_id", shape = "raster")
  lbls <- c(p$labels$title, p$patches$plots[[1]]$labels$title)
  expect_true(all(lbls == c("mousecolon", "mousebrain")))
  
  # Test ncol parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", ncol = 1), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", shape = "tile", ncol = 1), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", shape = "raster", ncol = 1), "patchwork")
  
  # Test scalebar parameters
  ## point
  expect_s3_class(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE, scalebar_gg = scalebar(x = 1000)), "patchwork")
  ## tile. suppress tile image rotated warning
  expect_s3_class(suppressWarnings(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE, shape = "tile", image_use = "raw")), "patchwork")
  expect_s3_class(suppressWarnings(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE, shape = "tile", image_use = "transformed")), "patchwork")
  expect_s3_class(suppressWarnings(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE, scalebar_gg = scalebar(x = 1000), shape = "tile", image_use = "raw")), "patchwork")
  expect_s3_class(suppressWarnings(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE, scalebar_gg = scalebar(x = 1000), shape = "tile", image_use = "transformed")), "patchwork")
  
  # Test colors parameter
  expect_s3_class(MapLabels(se_merged, column_name = "selection", colors = c("red", "blue")), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", colors = c("red", "blue"), shape = "tile"), "patchwork")
  expect_s3_class(MapLabels(se_merged, column_name = "selection", colors = c("red", "blue"), shape = "raster"), "patchwork")
  
  # Test that the correct samples and column_name are plotted
  p1 <- MapLabels(se_merged, column_name = "selection", return_plot_list = TRUE)
  p2 <- MapLabels(se_merged, column_name = "selection", return_plot_list = TRUE, shape = "tile")
  p3 <- MapLabels(se_merged, column_name = "selection", return_plot_list = TRUE, shape = "raster")
  
  # Test that the list contains correct classes
  expect_true(all(sapply(p1, function(p) class(p)[1]) == "gg"))
  expect_true(all(sapply(p2, function(p) class(p)[1]) == "gg"))
  expect_true(all(sapply(p3, function(p) class(p)[1]) == "gg"))
  
  # Test that the correct column_name are plotted
  expect_true(all(c("1", "2") == names(p1)))
  expect_true(all(c("1", "2") == names(p2)))
  expect_true(all(c("1", "2") == names(p3)))
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
  
  # Test HE with raster
  expect_error(MapLabels(se_merged, column_name = "selection", image_use = "raw", shape = "raster"))
  
  # Test error for scalebar with tiles/raster and no HE
  expect_error(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE, shape = "tile"))
  expect_error(MapLabels(se_merged, column_name = "selection", add_scalebar = TRUE, shape = "raster"))
})
