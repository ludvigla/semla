library(testthat)
library(semla)

# Create test data
set.seed(123)
pm <- matrix(runif(100), ncol = 2)
qm <- matrix(runif(100), ncol = 2)
xy_ref <- matrix(runif(100), ncol = 2)
xy_query <- matrix(runif(100), ncol = 2)

test_that("icp returns valid results", {
  result <- icp(xy_ref, xy_query)
  
  # Check if the result is a list with "y_transf" and "rot_mat" elements
  expect_type(result, "list")
  expect_named(result, c("y_transf", "rot_mat"))
  expect_type(result$y_transf, "double")
  expect_type(result$rot_mat, "double")
})
