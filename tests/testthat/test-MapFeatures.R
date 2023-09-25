# Load required packages and data
library(testthat)
library(semla)
library(viridis)
library(ggplot2)

se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
se_merged <- MergeSTData(se_mbrain, se_mcolon) |> 
  ScaleData() |> 
  LoadImages() |> 
  RigidTransformImages(transforms = generate_rigid_transform(sampleID = 2, angle = 45))

# Define a test case
test_that("MapFeatures works as expected", {
  # Test basic functionality
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip")), "patchwork")
  
  # Test slot parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), slot = "scale.data"), "patchwork")
  
  # Test image_use parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), image_use = "raw"), "patchwork")
  
  # Test coords_use parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), coords_use = "transformed"), "patchwork")
  
  # Test crop_area parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), crop_area = c(0.1, 0.1, 0.9, 0.9)), "patchwork")
  
  # Test pt_size parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), pt_size = 2), "patchwork")
  
  # Test pt_alpha parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), pt_alpha = 0.5), "patchwork")
  
  # Test pt_stroke parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), pt_stroke = 0.5), "patchwork")
  
  # Test scale_alpha parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), scale_alpha = TRUE), "patchwork")
  
  # Test section_number parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), section_number = 2), "patchwork")
  
  # Test label_by parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), label_by = "sample_id"), "patchwork")
  p <- MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), label_by = "sample_id")
  lbls <- c(p$labels$title, p$patches$plots[[1]]$labels$title)
  expect_true(all(lbls == c("mousecolon", "mousebrain")))
  
  # Test ncol parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), ncol = 2), "patchwork")
  
  # Test blend parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), blend = TRUE), "patchwork")
  
  # Test scalebar parameters
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), add_scalebar = TRUE), "patchwork")
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), add_scalebar = TRUE, scalebar_gg = scalebar(x = 1000)), "patchwork")
  
  # Test colors parameter
  expect_s3_class(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), colors = magma(n = 11, direction = -1)), "patchwork")
  
  # Test that the correct samples and features are plotted
  p <- MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), return_plot_list = TRUE)
  
  # Test that the list contains correct classes
  expect_true(all(sapply(p, class) == "list"))
  expect_true(all(sapply(p[[1]], function(p) class(p)[1]) == "gg"))
  
  # Test that the correct features are plotted
  expect_true(all(c("1", "2") == names(p)))
  expect_true(all(c("Clu", "Slc6a3", "Vip") == names(p[[1]])))
  
  # Test that the expected number of plots are produced
  p <- MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"))
  expect_equal(p$patches$layout$ncol, 3)
})


# Test argument validation
test_that("MapFeatures throws appropriate errors for invalid arguments", {
  # Test invalid feature
  expect_error(MapFeatures(se_merged, features = "not_a_feature"))
  
  # Test invalid slot
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), slot = "not_a_slot"))
  
  # Test invalid image_use
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), image_use = "test"))
               
  # Test with invalid section_number
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), section_number = 3),
               "section_number = 3 out of range")
  
  # Test with invalid scalebar_height
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), add_scalebar = TRUE,
                           scalebar_height = 2),
               "Expected a numeric of length 1 between 0 and 1")
  
  # Test with invalid scalebar_position
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), add_scalebar = TRUE,
                           scalebar_position = c(1.5, 0.5)),
               "Expected values between 0 and 1")
  
  # Test with invalid scalebar_gg
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), add_scalebar = TRUE,
                           scalebar_gg = NULL),
               "Expected a 'ggplot' object")
  
  # Test with invalid blending
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip", "Th"), blend = TRUE))
  expect_error(MapFeatures(se_merged, features = c("Clu"), blend = TRUE))
  
  # Test with invalid blend order
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), blend = TRUE, blend_order = 1:4))
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), blend = TRUE, blend_order = 2:1))
  
  # test with invalid crop_area
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), crop_area = c(-1, 0.3, 0.6, 0.6)),
               "'crop_area' can only take values between 0-1")
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), crop_area = c(0.3, 0.3, 0.3, 0.6)),
               "'left' value needs to be lower that 'right' value")
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), crop_area = c(0.3, 0.3, 0.6, 0.3)),
               "'top' value needs to be lower that 'bottom' value")
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), crop_area = c(0.3, 0.3, 0.3)),
               "Invalid length for 'crop_area', expected a 'numeric' vector of length 4")
  expect_error(MapFeatures(se_merged, features = c("Clu", "Slc6a3", "Vip"), crop_area = TRUE),
               "Invalid class 'logical' for 'crop_area', expected 'numeric'")
})
