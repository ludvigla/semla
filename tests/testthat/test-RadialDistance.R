library(testthat)
library(semla)

se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))

test_that("RadialDistance function returns valid results", {
  
  # Test case 1: Check if the function returns a data frame
  expect_s4_class(RadialDistance(se_mcolon, column_name = "selection", selected_groups = "GALT"), "Seurat")
  
  # Test case 2: Check if the function handles angles correctly
  se_mcolon <- DisconnectRegions(se_mcolon, column_name = "selection", selected_groups = "GALT")
  expect_s4_class(RadialDistance(se_mcolon, column_name = "GALT_split", selected_groups = "S1_region1", angles = c(0, 90)), "Seurat")
  
  # Test case 3: Check if the function handles angles_nbreaks correctly
  expect_s4_class(RadialDistance(se_mcolon, column_name = "GALT_split", selected_groups = "S1_region1", angles_nbreaks = 4), "Seurat")
  
  # Test case 4: Check if the function handles remove_singletons correctly
  expect_s4_class(RadialDistance(se_mcolon, column_name = "GALT_split", selected_groups = "S1_region1", remove_singletons = FALSE), "Seurat")
  
  # Test case 5: Check if the function handles convert_to_microns correctly
  expect_s4_class(RadialDistance(se_mcolon, column_name = "GALT_split", selected_groups = "S1_region1", convert_to_microns = TRUE), "Seurat")
})

test_that("RadialDistance function returns valid errors", {
  # Test case 6: Check if the function handles invalid input  (e.g., non-existent spots)
  expect_error(RadialDistance(se_mcolon, column_name = "missing"), class = "error")
  
  # Test case 7: Check if the function handles invalid angles 
  expect_error(RadialDistance(se_mcolon, column_name = "GALT_split", selected_groups = "S1_region1", angles = c(90, 0)), class = "error")
})