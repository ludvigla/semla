#' Create a react app for paper JS in R
#'
#' Provided a list of images, this function is used to open an interactive app
#' built with react and paper JS. The application is mean to be used in a shiny 
#' application and is used for the \code{\link{RunAlignment}} function provided 
#' in \code{semla}.
#'
#' @param data A list of images prepared with \code{\link{prep_image}}
#' @param width Width of component
#' @param height height of component
#' @param elementId Component element ID
#' 
#' @return A \code{htmlwidget} to be used in a \code{shiny} application
#' 
#' @examples 
#' 
#' library(semla)
#' library(magick)
#' 
#' im_mbrain <- system.file("extdata/mousebrain/spatial",
#'                          "tissue_lowres_image.jpg",
#'                          package = "semla")
#' 
#' img1 <- prep_image(im_mbrain |>
#'                      image_read(),
#'                    height = 256)
#' img2 <- prep_image(im_mbrain |>
#'                      image_read() |>
#'                      image_flip(),
#'                    height = 256)
#' 
#' \dontrun{
#' if (interactive()) {
#' # Run widget
#' paper(data = list(img1, img2), width = 600, height = 600)
#' }
#' }
#'
#' @export
paper <- function (
    data,
    width = NULL,
    height = NULL,
    elementId = NULL
) {

  # Require htmlwidgets
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    abort(glue("Package {cli::col_br_magenta('htmlwidgets')} is required. Please install it with: \n",
               "install.packages('htmlwidgets')"))
  }

  content <- reactR::component(
    "Paper",
    list(data = data, width = width, height = height)
  )

  # describe a React component to send to the browser for rendering.
  component <- reactR::reactMarkup(content)

  # create widget
  htmlwidgets::createWidget(
    name = 'paper',
    component,
    width = width,
    height = height,
    package = 'semla',
    elementId = elementId
  )
}

#' Called by HTMLWidgets to produce the widget's root element.
#' @noRd
widget_html.paper <- function (
    id,
    style,
    class,
    ...
) {

  # Require htmltools
  if (!requireNamespace("htmltools", quietly = TRUE)) {
    abort(glue("Package {cli::col_br_magenta('htmltools')} is required. Please install it with: \n",
               "install.packages('htmltools')"))
  }

  htmltools::tagList(
    # Necessary for RStudio viewer version < 1.2
    reactR::html_dependency_corejs(),
    reactR::html_dependency_react(),
    reactR::html_dependency_reacttools(),
    htmltools::tags$div(id = id, class = class, style = style)
  )
}

#' Shiny bindings for paper
#'
#' Output and render functions for using paper within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a paper
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#' is useful if you want to save an expression in a variable.
#'   
#' @return A \code{shiny.tag.list} output or a \code{shiny.render.function} 
#' function produced to be used in a \code{shiny} app. \code{paperOutput} 
#' is used in the UI and \code{renderPaper} on the server side.
#'
#' @name paper-shiny
#'
#' @export
#'
paperOutput <- function (
    outputId,
    width = '100%',
    height = '400px'
){
  # Require htmlwidgets
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    abort(glue("Package {cli::col_br_magenta('htmlwidgets')} is required. Please install it with: \n",
               "install.packages('htmlwidgets')"))
  }

  htmlwidgets::shinyWidgetOutput(outputId, 'paper', width, height, package = 'semla')
}

#' @rdname paper-shiny
#' @export
#'
renderPaper <- function (
    expr,
    env = parent.frame(),
    quoted = FALSE
) {
  # Require htmlwidgets
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    abort(glue("Package {cli::col_br_magenta('htmlwidgets')} is required. Please install it with: \n",
               "install.packages('htmlwidgets')"))
  }

  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, paperOutput, env, quoted = TRUE)
}
