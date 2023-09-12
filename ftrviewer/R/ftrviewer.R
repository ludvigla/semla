#' Open a feature viewer react application
#'
#' This function will open a react application that listens can
#' be used to interactively visualize categorical or numeric features
#' as spatial maps.
#'
#' @import htmlwidgets
#'
#' @export
ftrviewer <- function (
  host = "127.0.0.1",
  port = "8080",
  sampleID = 1,
  scaleByOpacity = FALSE,
  isNumeric = TRUE,
  levels = c("A", "B", "C"),
  categories = c("cat1", "cat2", "cat3"),
  colors = RColorBrewer::brewer.pal(n = 9, name = "Spectral"),
  container_width = 800,
  container_height = 800,
  elementId = NULL
) {

  stopifnot(
    inherits(host, what = "character") & length(host) == 1,
    inherits(port, what = "character") & length(port) == 1,
    inherits(sampleID, what = c("integer", "numeric")) & length(sampleID) == 1,
    inherits(scaleByOpacity, what = "logical") & length(scaleByOpacity) == 1,
    inherits(isNumeric, what = "logical") & length(isNumeric) == 1,
    inherits(levels, what = "character"),
    inherits(categories, what = "character"),
    inherits(colors, what = "character"),
    inherits(container_width, what = c("numeric", "integer")) & length(container_width) == 1,
    inherits(container_height, what = c("numeric", "integer")) & length(container_height) == 1
  )

  # describe a React component to send to the browser for rendering.
  content <- reactR::component(
    "Ftrviewer",
    list(
      port = port,
      host = host,
      width = container_width,
      height = container_height,
      sampleID = sampleID,
      scaleByOpacity = scaleByOpacity,
      isNumeric = isNumeric,
      levels = as.list(levels),
      categories = as.list(categories),
      colors = as.list(colors)
    )
  )
  component <- reactR::reactMarkup(content)

  # create widget
  htmlwidgets::createWidget(
    name = 'ftrviewer',
    component,
    width = paste0(container_width, "px"),
    height = paste0(container_height, "px"),
    package = 'ftrviewer', # This needs to be correct
    elementId = elementId
  )
}

#' Called by HTMLWidgets to produce the widget's root element.
#' @noRd
widget_html.ftrviewer <- function(id, style, class, ...) {
  htmltools::tagList(
    # Necessary for RStudio viewer version < 1.2
    reactR::html_dependency_corejs(),
    reactR::html_dependency_react(),
    reactR::html_dependency_reacttools(),
    htmltools::tags$div(id = id, class = class, style = style)
  )
}

#' Shiny bindings for ftrviewer
#'
#' Output and render functions for using ftrviewer within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a ftrviewer
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name ftrviewer-shiny
#'
#' @export
ftrviewerOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'ftrviewer', width, height, package = 'ftrviewer') # This needs to be correct
}

#' @rdname ftrviewer-shiny
#' @export
renderFtrviewer <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, ftrviewerOutput, env, quoted = TRUE)
}
