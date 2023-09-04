library(testthat)
library(semla)

test_that("MapMultipleFeatures works with example dataset", {
  # Load example Visium data
  se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
  se_mbrain <- LoadImages(se_mbrain)
  
  # Select features to plot
  sel_features <- c("Th", "Trh", "Calb2", "Prkcd")
  
  # Test with shared scale
  p1 <- MapMultipleFeatures(
    object = se_mbrain,
    features = sel_features,
    scale = "shared"
  )
  
  # Test with free scale
  p2 <- MapMultipleFeatures(
    object = se_mbrain,
    features = sel_features,
    scale = "free"
  )
  
  p3 <- MapMultipleFeatures(
    object = se_mbrain,
    features = sel_features, 
    return_plot_list = TRUE
  )
  
  p4 <- MapMultipleFeatures(
    object = se_mbrain,
    features = sel_features, 
    section_number = 1
  )
  
  p5 <- MapMultipleFeatures(
    object = se_mbrain,
    features = sel_features, 
    label_by = "orig.ident"
  )
  
  p6 <- MapMultipleFeatures(
    object = se_mbrain,
    features = sel_features, 
    crop_area = c(0.3, 0.3, 0.6, 0.6)
  )
  
  p7 <- MapMultipleFeatures(
    object = se_mbrain,
    features = sel_features, 
    colors = c("red")
  )
  
  p8 <- MapMultipleFeatures(
    object = se_mbrain,
    features = sel_features, 
    image_use = "raw"
  )
  
  p9 <- MapMultipleFeatures(
    object = se_mbrain,
    features = sel_features, 
    override_plot_dims = TRUE
  )
  
  expect_s3_class(p1, class = "patchwork")
  expect_s3_class(p2, class = "patchwork")
  expect_type(p3, type = "list")
  expect_s3_class(p4, class = "patchwork")
  expect_s3_class(p5, class = "patchwork")
  expect_s3_class(p6, class = "patchwork")
  expect_s3_class(p7, class = "patchwork")
  expect_s3_class(p8, class = "patchwork")
  expect_s3_class(p9, class = "patchwork")
})
