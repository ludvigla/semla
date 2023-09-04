library(testthat)
library(dplyr)
library(tibble)

# Helper function to compare transformed coordinates within a tolerance range
expect_coords_equal <- function(coords1, coords2, tolerance = 1e-5) {
  expect_true(all(abs(coords1 - coords2) < tolerance))
}

# Test with valid input
test_that("CoordTransform works with valid input", {
  xy <- tibble(x = 1:20, y = 1:20)
  result <- CoordTransform(xy, angle = 45)
  expect_s3_class(result, "tbl")
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 20)
})

# Test with invalid input
test_that("CoordTransform throws errors with invalid input", {
  xy <- tibble(x = 1:20, y = 1:20)
  expect_error(CoordTransform(xy, angle = "invalid"), "Invalid class 'character' for angle, expected a 'numeric'")
  expect_error(CoordTransform(xy, angle = c(10, 20)), "Invalid length '2' for angle, expected 1 value")
  expect_error(CoordTransform(xy, angle = 400), "Invalid value `400`for angle. Value should be in range \\[-360, 360\\]")
  expect_error(CoordTransform(xy, center = "invalid"), "Invalid class 'character' for center, expected a 'numeric'")
  expect_error(CoordTransform(xy, center = c(10, 20, 30)), "Invalid length '3' for angle, expected 2 values")
  expect_error(CoordTransform(xy, xy_offset = "invalid"), "Invalid class 'character' for xy_offset, expected a 'numeric'")
  expect_error(CoordTransform(xy, xy_offset = c(10, 20, 30)), "Invalid length '3' for angle, expected 2 values")
  expect_error(CoordTransform(xy, angle = 0, xy_offset = c(0, 0), scalefactor = 1), "Only default values provided which would result in no tranformation.")
})

# Test rotation
test_that("CoordTransform correctly rotates coordinates", {
  xy <- tibble(x = 1:2, y = 1:2)
  result <- CoordTransform(xy, angle = 90, center = c(0, 0))
  expected <- tibble(tr_x = c(-1, -2), tr_y = c(1, 2))
  expect_coords_equal(result, expected)
})

# Test translation
test_that("CoordTransform correctly translates coordinates", {
  xy <- tibble(x = 1:2, y = 1:2)
  result <- CoordTransform(xy, angle = 0, xy_offset = c(10, 20))
  expected <- tibble(tr_x = c(11, 12), tr_y = c(21, 22))
  expect_coords_equal(result, expected)
})

# Test scaling
test_that("CoordTransform correctly scales coordinates", {
  xy <- tibble(x = 1:2, y = 1:2)
  result <- CoordTransform(xy, angle = 0, scalefactor = 2)
  expected <- tibble(tr_x = c(0.5, 2.5), tr_y = c(0.5, 2.5))
  expect_coords_equal(result, expected)
})

# Test combination of transformations
test_that("CoordTransform correctly applies combination of transformations", {
  xy <- tibble(x = 1:2, y = 1:2)
  result <- CoordTransform(xy, angle = 90, center = c(0, 0), xy_offset = c(10, 20), scalefactor = 2)
  expected <- tibble(tr_x = c(8, 6), tr_y = c(22, 24))
  expect_coords_equal(result, expected)
})

# Test edge case with zero rotation and translation
test_that("CoordTransform handles edge case with zero rotation and translation", {
  xy <- tibble(x = 1:20, y = 1:20)
  expect_error(CoordTransform(xy, angle = 0, xy_offset = c(0, 0)), 
               "Only default values provided which would result in no tranformation.")
})

# Test edge case with negative angle
test_that("CoordTransform handles edge case with negative angle", {
  xy <- tibble(x = 1:2, y = 1:2)
  result <- CoordTransform(xy, angle = -90, center = c(0, 0))
  expected <- tibble(tr_x = c(1, 2), tr_y = c(-1, -2))
  expect_coords_equal(result, expected)
})


# Test edge case with angle over 360 degrees
test_that("CoordTransform handles edge case with angle over 360 degrees", {
  xy <- tibble(x = 1:2, y = 1:2)
  expect_error(CoordTransform(xy, angle = 450, center = c(0, 0)),
               "Value should be in range")
})  
