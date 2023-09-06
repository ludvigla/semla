library(testthat)
library(semla)

# Test ftrviewer
test_that("ftrviewer returns an htmlwidget", {
  widget <- ftrviewer(values = 1, opacities = 1, range = c(0, 1))
  expect_s3_class(widget, "htmlwidget")
})

# Test widget_html.ftrviewer
test_that("widget_html.ftrviewer returns an HTML tag list", {
  tags <- widget_html.ftrviewer("id", "style", "class")
  expect_s3_class(tags, "shiny.tag.list")
})

# Test ftrviewerOutput
test_that("ftrviewerOutput returns an HTML output element", {
  output <- suppressWarnings({ftrviewerOutput("outputId")})
  expect_s3_class(output, "shiny.tag.list")
})

# Test renderftrviewer
test_that("renderftrviewer returns a shiny render function", {
  renderFunc <- renderFtrviewer(expr = { })
  expect_s3_class(renderFunc, "shiny.render.function")
})

# Test Case 2: Test if renderftrviewer generates a render function
test_that("renderftrviewer generates a render function", {
  render_function <- renderFtrviewer({
    paper(data = list(), width = "100%", height = "400px")
  })
  expect_true(inherits(render_function, what = "shiny.render.function"))
})
