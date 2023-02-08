#' Create a react app for digital unrolling
#'
#' This function is used to start an interactive widget for
#' digital unrolling. It requires a static files server to be
#' hosted in order to find load the necessary files.
#'
#' @param sampleID A section ID
#' @param host A host address
#' @param port A valid port
#' @param width,height Width and height of container
#' @param elementId The element id of the widget
#' @param quit A logical specifying is the app should quit
#'
#' @export
osddu <- function (
    sampleID = 1,
    host = "127.0.0.1",
    port = "8080",
    width = NULL,
    height = NULL,
    elementId = NULL,
    quit = FALSE
) {

  # Require htmlwidgets
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    install.packages("htmlwidgets")
  }

  # describe a React component to send to the browser for rendering.
  content <- reactR::component(
    "Osddu",
    list(sampleID = sampleID, port = port, host = host, width = width, height = height, quit = quit)
  )
  component <- reactR::reactMarkup(content)

  # create widget
  htmlwidgets::createWidget(
    name = 'osddu',
    component,
    width = width,
    height = height,
    package = 'STUtiity2',
    elementId = elementId
  )
}

#' Called by HTMLWidgets to produce the widget's root element.
#' @noRd
widget_html.osddu <- function(id, style, class, ...) {
  htmltools::tagList(
    # Necessary for RStudio viewer version < 1.2
    reactR::html_dependency_corejs(),
    reactR::html_dependency_react(),
    reactR::html_dependency_reacttools(),
    htmltools::tags$div(id = id, class = class, style = style)
  )
}

#' Shiny bindings for osddu
#'
#' Output and render functions for using osddu within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a osddu
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name osddu-shiny
#'
#' @export
osdduOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'osddu', width, height, package = 'semla')
}

#' @rdname osddu-shiny
#' @export
renderOsddu <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, osdduOutput, env, quoted = TRUE)
}
