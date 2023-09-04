library(testthat)
library(tibble)

# Test that CoordMirror handles invalid input correctly
test_that("CoordMirror handles invalid input", {
  # Test invalid format
  expect_error(CoordMirror(1:10, mirror.x = TRUE), "Invalid class")
  
  # Test invalid number of columns
  expect_error(CoordMirror(matrix(1:6, nrow = 2), mirror.x = TRUE), "Expected 2 columns")
  
  # Test non-numeric columns
  expect_error(CoordMirror(data.frame(x = 1:5, y = letters[1:5])), "Invalid column classes.")
  
  # Test neither mirror.x nor mirror.y set
  expect_error(CoordMirror(matrix(1:10, ncol = 2)), "One of 'mirror.x' or 'mirror.y' or both need to be selected.")
})

# Test that CoordMirror mirrors along the x-axis correctly
test_that("CoordMirror mirrors along x-axis", {
  input_coords <- tibble(x = 1:5, y = 1:5)
  expected_output <- tibble(tr_x = 5:1, tr_y = 1:5)
  
  mirrored_coords <- CoordMirror(input_coords, mirror.x = TRUE)
  expect_equal(mirrored_coords, expected_output)
})

# Test that CoordMirror mirrors along the y-axis correctly
test_that("CoordMirror mirrors along y-axis", {
  input_coords <- tibble(x = 1:5, y = 1:5)
  expected_output <- tibble(tr_x = 1:5, tr_y = 5:1)
  
  mirrored_coords <- CoordMirror(input_coords, mirror.y = TRUE)
  expect_equal(mirrored_coords, expected_output)
})

# Test that CoordMirror mirrors along both x and y axes correctly
test_that("CoordMirror mirrors along x and y axes", {
  input_coords <- tibble(x = 1:5, y = 1:5)
  expected_output <- tibble(tr_x = 5:1, tr_y = 5:1)
  
  mirrored_coords <- CoordMirror(input_coords, mirror.x = TRUE, mirror.y = TRUE)
  expect_equal(mirrored_coords, expected_output)
})

# Test that CoordMirror works with specified center
test_that("CoordMirror works with custom center", {
  input_coords <- tibble(x = 1:5, y = 1:5)
  center <- c(3, 3)
  expected_output <- tibble(tr_x = c(5, 4, 3, 2, 1), tr_y = c(5, 4, 3, 2, 1))
  
  mirrored_coords <- CoordMirror(input_coords, mirror.x = TRUE, mirror.y = TRUE, center = center)
  expect_equal(mirrored_coords, expected_output)
})

# Test that CoordMirror returns a tibble
test_that("CoordMirror returns a tibble", {
  input_coords <- tibble(x = 1:5, y = 1:5)
  mirrored_coords <- CoordMirror(input_coords, mirror.x = TRUE)
  expect_s3_class(mirrored_coords, "tbl")
})

