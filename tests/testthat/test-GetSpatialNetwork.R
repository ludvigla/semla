library(testthat)

# Create test data
coordfiles <- c(system.file("extdata/mousebrain/spatial",
                            "tissue_positions_list.csv",
                            package = "semla"),
                system.file("extdata/mousecolon/spatial",
                            "tissue_positions_list.csv",
                            package = "semla"))

# Load coordinate data into a tibble
xys <- do.call(rbind, lapply(seq_along(coordfiles), function(i) {
  coords <- setNames(read.csv(coordfiles[i], header = FALSE),
                     nm = c("barcode", "selection", "grid_y", "grid_x", "y", "x"))
  coords$sampleID <- i
  coords <- coords |>
    dplyr::filter(selection == 1) |>
    dplyr::select(barcode, x, y, sampleID) |>
    tibble::as_tibble()
  return(coords)
}))

# Test that function returns a list
test_that("GetSpatialNetwork returns a list", {
  expect_type(GetSpatialNetwork(xys), "list")
})

# Test that function returns a list with the same number of elements as unique sample IDs
test_that("GetSpatialNetwork returns a list with the same number of elements as unique sample IDs", {
  expect_equal(length(GetSpatialNetwork(xys)), length(unique(xys$sampleID)))
})

# Test that function throws an error when passed an object with invalid class
test_that("GetSpatialNetwork throws an error when passed an object with invalid class", {
  expect_error(GetSpatialNetwork(1:10), "Invalid class 'integer'.")
})

# Test that function throws an error when passed an object with invalid number of columns
test_that("GetSpatialNetwork throws an error when passed an object with invalid number of columns", {
  expect_error(GetSpatialNetwork(xys[, -1]), "Invalid number of columns '3'. Expected 4.")
})

# Test that function throws an error when passed an object with invalid column names
test_that("GetSpatialNetwork throws an error when passed an object with invalid column names", {
  expect_error(GetSpatialNetwork(xys |> rename(foo = barcode)), "Expected 'barcode', 'x', 'y' and 'sampleID'")
})

# Test that function throws an error when passed an object with invalid column classes
test_that("GetSpatialNetwork throws an error when passed an object with invalid column classes", {
  expect_error(GetSpatialNetwork(xys |> mutate(x = "string")), "Invalid column class.")
})

# Test that function calculates the correct number of nearest neighbors
test_that("GetSpatialNetwork calculates the correct number of nearest neighbors", {
  expect_equal(nrow(GetSpatialNetwork(xys)[[1]]), 14818)
})

# Test that function discards spots with fewer neighbors than minK
test_that("GetSpatialNetwork discards spots with fewer neighbors than minK", {
  expect_equal(nrow(GetSpatialNetwork(xys, minK = 6)[[1]]), 0)
})

se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))

# Test that function throws an error when passed a Seurat object without a Staffli object
test_that("GetSpatialNetwork throws an error when passed a Seurat object without Staffli object", {
  tmp <- se_mbrain
  tmp@tools$Staffli <- NULL
  expect_error(GetSpatialNetwork(tmp), "'Staffli' object is missing from tools slot.")
})