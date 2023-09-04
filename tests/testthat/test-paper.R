library(testthat)

# Test paper
test_that("paper returns an htmlwidget", {
  skip_on_cran()  # Skip this test on CRAN because it requires external dependencies
  data <- list()  # Need actual data for a meaningful test
  widget <- paper(data)
  expect_s3_class(widget, "htmlwidget")
})

# Test widget_html.paper
test_that("widget_html.paper returns an HTML tag list", {
  skip_on_cran()  # Skip this test on CRAN because it requires external dependencies
  tags <- widget_html.paper("id", "style", "class")
  expect_s3_class(tags, "shiny.tag.list")
})

# Test paperOutput
test_that("paperOutput returns an HTML output element", {
  skip_on_cran()  # Skip this test on CRAN because it requires external dependencies
  output <- paperOutput("outputId")
  expect_s3_class(output, "shiny.tag.list")
})

# Test renderPaper
test_that("renderPaper returns a shiny render function", {
  skip_on_cran()  # Skip this test on CRAN because it requires external dependencies
  renderFunc <- renderPaper(expr = { })
  expect_s3_class(renderFunc, "shiny.render.function")
})
