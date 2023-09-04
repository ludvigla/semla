library(testthat)
library(tibble)

# Load a dummy Seurat object for testing
coordfile <- system.file("extdata/mousecolon/spatial",
                         "tissue_positions_list.csv",
                         package = "semla")
coords <- read.csv(coordfile, header = FALSE) |>
  filter(V2 == 1) |>
  select(V1, V6, V5) |>
  setNames(nm = c("barcode", "x", "y")) |>
  bind_cols(sampleID = 1) |>
  as_tibble()

galt_spots_file <- system.file("extdata/mousecolon",
                               "galt_spots.csv",
                                package = "semla")
galt_spots <- read.csv(galt_spots_file) |>
  as_tibble()

spots <- galt_spots$barcode[galt_spots$selection == "GALT"]

test_that("default method works as expected", {
  
  # Check that the returned object is of class character
  expect_type(DisconnectRegions(coords, spots), "character")
  
  # it raises an error with an invalid coords object
  expect_error(DisconnectRegions(coords |> select(-sampleID), spots), 
               "Invalid number of columns '3'. Expected 4.")
  
  # it raises an error with an invalid spots object
  expect_error(DisconnectRegions(coords, spots = c("A", "B")))
  
  # it raises an error with no disconnected components
  expect_error(DisconnectRegions(coords, spots = coords$barcode), 
               "Found no disconnected components.")
  
  # it correctly labels each region and singletons
  labels <- DisconnectRegions(coords, spots) |> table()
  expect_equal(names(labels), c("S1_region1", "S1_region2", "S1_region3",
                                "S1_region4", "S1_region5", "S1_region6",
                                "S1_region7", "S1_region8", "singleton"))
  expect_equal(labels |> as.integer(), c(67, 6, 6, 5, 3, 3, 2, 2, 12))
})

# Load a dummy Seurat object for testing
se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
se_merged <- MergeSTData(se_mbrain, se_mcolon) |> 
  ScaleData() |> 
  LoadImages() |> 
  RigidTransformImages(transforms = generate_rigid_transform(sampleID = 2, angle = 45))

test_that("Seurat method works as expected", {
  
  # the returned object is of class Seurat
  expect_s4_class(DisconnectRegions(se_merged, column_name = "selection", selected_groups = "GALT"), 
              "Seurat")

  # it throws an error when the column name is invalid
  expect_error(DisconnectRegions(se_merged, column_name = "missing_column", selected_groups = "GALT"))
  
  # it throws an error when the selected group is invalid
  expect_error(DisconnectRegions(se_merged, column_name = "selection", selected_groups = "missing_group"))
  
  # it raises an error with no disconnected components
  expect_error(DisconnectRegions(SubsetSTData(se_merged, spots = setdiff(colnames(se_merged), 
                                                                         c("GGAAAGTCTTGATTGT-1",
                                                                           "TCAAAGTCACGGCGTC-1",
                                                                           "TTGTCTCGGCAAGATG-1"))), 
                                 column_name = "orig.ident"), 
               "Found no disconnected components.")
  
  # it correctly labels each region and singletons
  labels <- DisconnectRegions(se_merged, column_name = "selection", selected_groups = "GALT")[[]] |> 
    pull(GALT_split) |> 
    table()
  expect_equal(names(labels), c("S2_region1", "S2_region2", "S2_region3",
                                "S2_region4", "S2_region5", "S2_region6",
                                "S2_region7", "S2_region8", "singleton"))
  expect_equal(labels |> as.integer(), c(67, 6, 6, 5, 3, 3, 2, 2, 12))
})
