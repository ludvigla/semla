library(testthat)
library(patchwork)
library(ggplot2)
library(semla)

# Test ThemeClean function
test_that("ThemeClean returns a theme object with expected modifications", {
  theme_obj <- ThemeClean()
  expect_s3_class(theme_obj, "theme")
})

# Test ThemeLegendRight function
test_that("ThemeLegendRight returns a theme object with expected modifications", {
  theme_obj <- ThemeLegendRight()
  
  expect_s3_class(theme_obj, "theme")
})

# Test ThemeClean
test_that("ThemeClean returns a theme with blank title, subtitle, and no legend", {
  th <- ThemeClean()
  expect_true(is(th, "theme"))
  expect_equal(th$plot.title, element_blank())
  expect_equal(th$plot.subtitle, element_blank())
  expect_equal(th$legend.position, "none")
})

# Test ThemeLegendRight
test_that("ThemeLegendRight returns a theme with legend on the right and horizontal text", {
  th <- ThemeLegendRight()
  expect_true(is(th, "theme"))
  expect_equal(th$legend.position, "right")
  expect_equal(th$legend.text, element_text(angle = 0, hjust = 1))
})

# Test ModifyPatchworkTitles
test_that("ModifyPatchworkTitles modifies patchwork titles correctly", {
  # Create a patchwork object
  p <- (ggplot(mtcars, aes(mpg, disp)) + geom_point()) / (ggplot(mtcars, aes(mpg, hp)) + geom_point())
  
  # Modify titles
  p_modified <- ModifyPatchworkTitles(p, titles = c("First", "Second"))
  
  # Check that the titles have been modified correctly
  expect_equal(p_modified$patches$plots[[1]]$labels$title, "First")
  expect_equal(p_modified$labels$title, "Second")
})


# Test ModifyPatchworkTitles with invalid inputs
test_that("ModifyPatchworkTitles throws an error with invalid inputs", {
  # Not a patchwork object
  expect_error(ModifyPatchworkTitles(ggplot(mtcars, aes(mpg, disp)) + geom_point(), titles = c("First")), 
               "Expected a 'patcwork' object")
  
  # Not a character vector for titles
  p <- (ggplot(mtcars, aes(mpg, disp)) + geom_point()) / (ggplot(mtcars, aes(mpg, hp)) + geom_point())
  expect_error(ModifyPatchworkTitles(p, titles = list("First", "Second")), 
               "Expected titles to be 'character' vector")
  
  # Incorrect number of titles
  expect_error(ModifyPatchworkTitles(p, titles = c("A", "B", "C")), 
               "Invalid length of 'titles'. Expected a character vector of length 2.")
})
