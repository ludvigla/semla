library(testthat)
library(dplyr)

# Test if the function returns a Seurat object
test_that("ReadVisiumData returns a Seurat object", {
  infoTable <- tibble(samples = system.file("extdata/mousebrain/filtered_feature_bc_matrix.h5", package = "semla"),
                      imgs = system.file("extdata/mousebrain/spatial/tissue_lowres_image.jpg", package = "semla"),
                      spotfiles = system.file("extdata/mousebrain/spatial/tissue_positions_list.csv", package = "semla"),
                      json = system.file("extdata/mousebrain/spatial/scalefactors_json.json", package = "semla"))
  
  se <- ReadVisiumData(infoTable, verbose = FALSE)
  se_sf <- ReadVisiumData(infoTable |> select(-all_of("json")) |> mutate(scalefactor = 0.5), verbose = FALSE)
  
  expect_true(class(se) == "Seurat")
  expect_true(all(dim(se) == c(188, 2560)))
  
  expect_true(class(se_sf) == "Seurat")
  expect_true(all(dim(se_sf) == c(188, 2560)))
})

# Test if the function returns a Seurat object
test_that("ReadVisiumData returns a an error when the infoTable is incorrect", {
  infoTable <- tibble(samples = system.file("extdata/mousebrain/filtered_feature_bc_matrix.h5", package = "semla"),
                      imgs = system.file("extdata/mousebrain/spatial/tissue_lowres_image.jpg", package = "semla"),
                      spotfiles = system.file("extdata/mousebrain/spatial/tissue_positions_list.csv", package = "semla"),
                      json = system.file("extdata/mousebrain/spatial/scalefactors_json.json", package = "semla"))
  
  # Should throw error when samples, imgs or spotfiles is missing
  expect_error(ReadVisiumData(infoTable |> select(-samples), verbose = FALSE), 
               "One or several of 'samples', 'imgs' and 'spotfiles' are missing from infoTable")
  expect_error(ReadVisiumData(infoTable |> select(-imgs), verbose = FALSE), 
               "One or several of 'samples', 'imgs' and 'spotfiles' are missing from infoTable")
  expect_error(ReadVisiumData(infoTable |> select(spotfiles), verbose = FALSE), 
               "One or several of 'samples', 'imgs' and 'spotfiles' are missing from infoTable")
  expect_error(ReadVisiumData(infoTable |> select(-json), verbose = FALSE), 
               "One of 'json' or 'scalefactor' columns needs to be provided")
  
  # Throw error if column classes are invalid
  expect_error(ReadVisiumData(infoTable |> mutate(samples = NA), verbose = FALSE), 
               "Invalid column classes in 'infoTable'. Expected 'character' vectors")
  expect_error(ReadVisiumData(infoTable |> mutate(samples = NA), verbose = FALSE), 
               "Invalid column classes in 'infoTable'. Expected 'character' vectors")
  expect_error(ReadVisiumData(infoTable |> mutate(samples = NA), verbose = FALSE), 
               "Invalid column classes in 'infoTable'. Expected 'character' vectors")
  expect_error(ReadVisiumData(infoTable |> mutate(samples = NA), verbose = FALSE), 
               "Invalid column classes in 'infoTable'. Expected 'character' vectors")
  
  # Throw error if files are missing
  expect_error(ReadVisiumData(infoTable |> mutate(samples = NA_character_), verbose = FALSE), 
               "Missing file")
  expect_error(ReadVisiumData(infoTable |> mutate(samples = NA_character_), verbose = FALSE), 
               "Missing file")
  expect_error(ReadVisiumData(infoTable |> mutate(samples = NA_character_), verbose = FALSE), 
               "Missing file")
  expect_error(ReadVisiumData(infoTable |> mutate(samples = NA_character_), verbose = FALSE), 
               "Missing file")
})
