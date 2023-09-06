library(testthat)
library(semla)

se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))

# Test the export_graph function
test_that("export_graph exports a JSON file", {
  
  # Define a temporary directory for testing
  test_outdir <- tempdir()
  
  # Call the export_graph function
  export_graph(se_mbrain, sampleID = 1L, outdir = test_outdir, verbose = FALSE)
  
  # Check if the JSON file was created in the specified directory
  json_file <- file.path(test_outdir, "network_Visium_1.json")
  expect_true(file.exists(json_file))
  
  # Clean up: Remove the temporary directory and file
  unlink(json_file, recursive = TRUE)
})

