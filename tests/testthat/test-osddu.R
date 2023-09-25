library(testthat)
library(semla)

# Test osddu
test_that("osddu returns an htmlwidget", {
  data <- list()  # Need actual data for a meaningful test
  widget <- osddu(data)
  expect_s3_class(widget, "htmlwidget")
})

# Test widget_html.osddu
test_that("widget_html.osddu returns an HTML tag list", {
  tags <- widget_html.osddu("id", "style", "class")
  expect_s3_class(tags, "shiny.tag.list")
})

# Test osdduOutput
test_that("osdduOutput returns an HTML output element", {
  output <- suppressWarnings({osdduOutput("outputId")})
  expect_s3_class(output, "shiny.tag.list")
})

# Test renderOsddu
test_that("renderOsddu returns a shiny render function", {
  renderFunc <- renderOsddu(expr = { })
  expect_s3_class(renderFunc, "shiny.render.function")
})

# Test Case 2: Test if renderOsddu generates a render function
test_that("renderOsddu generates a render function", {
  render_function <- renderOsddu({
    paper(data = list(), width = "100%", height = "400px")
  })
  expect_true(inherits(render_function, what = "shiny.render.function"))
})
