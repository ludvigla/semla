library(beakr)

# Test case 1: Test if the function can host a directory and return a beakr instance
test_that("file_server returns a beakr instance", {
  beakr::stopAllServers()
  fs <- file_server(hostDir = tempdir())
  expect_equal(class(fs), c("Beakr", "R6"))
  beakr::stopServer(fs)
})

# Test case 2: Test if an error is raised if the directory to host files from does not exist
test_that("an error is raised if the directory to host files from does not exist", {
  expect_error(suppressWarnings({file_server(hostDir = "invalid_directory")}))
})
