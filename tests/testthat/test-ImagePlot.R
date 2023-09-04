library(testthat)
library(semla)

# Load a test Seurat object with semla
se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))

# Test that function throws an error if the images are not loaded
test_that("ImagePlot throws error", {
  expect_error(ImagePlot(se_mbrain), "Images have not been loaded yet")
})

se_mbrain <- LoadImages(se_mbrain)

# Test that function runs without errors
test_that("ImagePlot runs without errors", {
  expect_invisible(ImagePlot(se_mbrain))
})

# Test that function returns a ggplot object when return_as_gg = TRUE
test_that("ImagePlot returns ggplot object when return_as_gg = TRUE", {
  ggplot_obj <- ImagePlot(se_mbrain, return_as_gg = TRUE)
  expect_s3_class(ggplot_obj, "ggplot")
})

# Test that function throws an error if label_by is not a character or factor
test_that("ImagePlot throws error if label_by is not character or factor", {
  expect_error(ImagePlot(se_mbrain, label_by = "nFeature_Spatial"), 
               "Invalid class 'integer' for 'label_by' column. Expected a 'character' of 'factor'.")
})

# Test that function throws an error if label_by is not present in the Seurat object
test_that("ImagePlot throws error if label_by is not present in Seurat object", {
  expect_error(ImagePlot(se_mbrain, label_by = "nonexistent_column"), "not present in the Seurat object meta data")
})

# Test that function throws an error if label_by column has multiple labels per tissue section
test_that("ImagePlot throws error if label_by column has multiple labels per tissue section", {
  se_mbrain$test <- sample(c("brain", "colon"), ncol(se_mbrain), replace = TRUE)
  expect_error(ImagePlot(se_mbrain, label_by = "test"), "Invalid 'label_by' meta data column.")
})

# Test that function throws an error if transformed images are not available in object
test_that("ImagePlot throws error if transformed images are not available in object", {
  expect_error(ImagePlot(se_mbrain, image_use = "transformed"), "Transformed images are not available in this object")
})

# Test that function throws an error if crop_area is not numeric
test_that("ImagePlot throws error if crop_area is not numeric", {
  expect_error(ImagePlot(se_mbrain, crop_area = "not_a_number"), "Invalid class 'character' for 'crop_area', expected 'numeric'")
})

# Test that function throws an error if crop_area is not of length 4
test_that("ImagePlot throws error if crop_area is not of length 4", {
  expect_error(ImagePlot(se_mbrain, crop_area = c(0.1, 0.2, 0.3)), "Invalid length for 'crop_area', expected a 'numeric' vector of length 4")
})

# Test that function throws an error if crop_area values are not between 0 and 1
test_that("ImagePlot throws error if crop_area values are not between 0 and 1", {
  expect_error(ImagePlot(se_mbrain, crop_area = c(1.1, 0.2, 0.3, 0.4)), "'crop_area' can only take values between 0-1")
})

# Test that function throws an error if sampleIDs are out of range
test_that("ImagePlot throws error if sampleIDs are out of range", {
  expect_error(ImagePlot(se_mbrain, sampleIDs = c(0, 11)), "'sampleIDs' is out of range.")
})