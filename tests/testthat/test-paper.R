library(testthat)
library(semla)

# Test paper
test_that("paper returns an htmlwidget", {
  data <- list()  # Need actual data for a meaningful test
  widget <- paper(data)
  expect_s3_class(widget, "htmlwidget")
})

# Test widget_html.paper
test_that("widget_html.paper returns an HTML tag list", {
  tags <- widget_html.paper("id", "style", "class")
  expect_s3_class(tags, "shiny.tag.list")
})

# Test paperOutput
test_that("paperOutput returns an HTML output element", {
  output <- suppressWarnings({paperOutput("outputId")})
  expect_s3_class(output, "shiny.tag.list")
})

# Test renderPaper
test_that("renderPaper returns a shiny render function", {
  renderFunc <- renderPaper(expr = { })
  expect_s3_class(renderFunc, "shiny.render.function")
})

# Test Case 2: Test if renderPaper generates a render function
test_that("renderPaper generates a render function", {
  render_function <- renderPaper({
    paper(data = list(), width = "100%", height = "400px")
  })
  expect_true(inherits(render_function, what = "shiny.render.function"))
})
