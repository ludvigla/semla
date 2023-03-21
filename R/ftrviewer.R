#' Open a feature viewer react application
#'
#' This function will open a react application that listens can
#' be used to interactively visualize categorical or numeric features
#' as spatial maps. For the application to run properly, relevant data needs
#' to be available via a file server. This function is intended to be used
#' within a shiny app 
#'
#' @param host The host address. Defaults to localhost "127.0.0.1"
#' @param port A number for a valid port
#' @param sampleID A section number
#' @param values A vector of numeric or categorical values
#' @param opacities A numeric vector of opacity values
#' @param opacity An integer of length 1 specifying a fixed opacity value
#' @param range A numeric vector of length 2 specifying a range of values (color domain)
#' @param scaleByOpacity A logical specifying if the opacity should be set by
#' \code{opacities} or \code{opacity}
#' @param isNumeric A logical specifying if the input is numeric or not
#' @param useLasso A logical specifying if the lasso selection tool should be activated
#' @param levels Category levels for coloring of values
#' @param categories A character vector with the categories available
#' @param colors A character vector of colors
#' @param container_width,container_height Container width/height in pixels
#' @param elementId The id of the viewer element
#' 
#' @examples
#' \dontrun{
#' 
#' library(semla)
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
#' se_mbrain <- LoadImages(se_mbrain)
#' 
#' datapath <- ExportDataForViewer(se_mbrain, outdir = ".")
#' 
#' # Start file server
#' file_server(datapath)
#' 
#' # Run feature viewer widget
#' ftrviewer(values = as.numeric(se_mbrain$nFeature_Spatial), opacities = rep(1, ncol(se_mbrain)), isNumeric = TRUE,
#'           range = range(se_mbrain$nFeature_Spatial))
#' 
#' # Stop file server
#' beakr::stopAllServers()
#' 
#' }
#'
#' @export
ftrviewer <- function (
    host = "127.0.0.1",
    port = "8080",
    sampleID = 1,
    values,
    opacities,
    opacity = 1,
    range,
    scaleByOpacity = FALSE,
    isNumeric = TRUE,
    useLasso = FALSE,
    levels = character(),
    categories = character(),
    colors = RColorBrewer::brewer.pal(n = 9, name = "Spectral") |> rev(),
    container_width = 800,
    container_height = 800,
    elementId = NULL
) {

  # Require htmlwidgets
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    install.packages("htmlwidgets")
  }

  stopifnot(
    inherits(host, what = "character") & length(host) == 1,
    inherits(port, what = "character") & length(port) == 1,
    inherits(values, what = c("integer", "numeric", "character", "factor")) & length(values) > 0,
    inherits(opacity, what = c("integer", "numeric")) & length(opacity) == 1,
    inherits(opacities, what = c("integer", "numeric")) & length(opacities) > 0,
    inherits(sampleID, what = c("integer", "numeric")) & length(sampleID) == 1,
    inherits(scaleByOpacity, what = "logical") & length(scaleByOpacity) == 1,
    inherits(isNumeric, what = "logical") & length(isNumeric) == 1,
    inherits(useLasso, what = "logical") & length(useLasso) == 1,
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
      values = values,
      range = range,
      opacities = opacities,
      width = container_width,
      height = container_height,
      sampleID = sampleID,
      scaleByOpacity = scaleByOpacity,
      isNumeric = isNumeric,
      levels = as.list(levels),
      categories = as.list(categories),
      colors = as.list(colors),
      useLasso = useLasso,
      opacity = opacity
    )
  )
  component <- reactR::reactMarkup(content)

  # create widget
  htmlwidgets::createWidget(
    name = 'ftrviewer',
    component,
    width = paste0(container_width, "px"),
    height = paste0(container_height, "px"),
    package = 'semla',
    elementId = elementId
  )
}

#' Called by HTMLWidgets to produce the widget's root element.
#' @noRd
widget_html.ftrviewer <- function(id, style, class, ...) {

  # Require htmltools
  if (!requireNamespace("htmltools", quietly = TRUE)) {
    install.packages("htmltools")
  }

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
  # Require htmlwidgets
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    install.packages("htmlwidgets")
  }

  htmlwidgets::shinyWidgetOutput(outputId, 'ftrviewer', width, height, package = 'semla')
}

#' @rdname ftrviewer-shiny
#' @export
renderFtrviewer <- function(expr, env = parent.frame(), quoted = FALSE) {
  # Require htmlwidgets
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    install.packages("htmlwidgets")
  }

  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, ftrviewerOutput, env, quoted = TRUE)
}
