#' @include generics.R
#' @include checks.R
#'
NULL


#' @param width,height Width and height of paper widget given in pixels.
#'
#' @section default method:
#' Takes a list of images prepared with the \code{\link{prep_image}} function
#' and send them to the interactive application. When the application stops
#' (after clicking "quit & save") a tibble is returned with information about
#' the rigid transformations applied to the images.
#'
#' @importFrom shiny fluidPage actionButton tableOutput uiOutput reactiveTimer
#' observe fluidRow column p h4 helpText strong code h5 observeEvent stopApp
#' runApp renderTable renderUI hr
#' @importFrom shinyjs useShinyjs extendShinyjs
#' @importFrom tibble as_tibble
#'
#' @rdname manual-transform-images
#'
#' @examples
#' \dontrun{
#'
#' library(STUtility2)
#' library(magick)
#'
#' im_mbrain <- system.file("extdata/mousebrain/spatial",
#'                          "tissue_hires_image.png",
#'                          package = "STUtility2")
#'
#' img1 <- prep_image(im_mbrain |>
#'                      image_read(),
#'                    width = 512)
#' img2 <- prep_image(im_mbrain |>
#'                      image_read() |>
#'                      image_flip(),
#'                    width = 512)
#'
#' transforms <- RunAlignment(object = list(img1, img2))
#'
#' }
#'
#' @export
#'
RunAlignment.default <- function (
    object,
    width = '800px',
    height = '650px',
    ...
) {

  # Validate object data
  stopifnot(is.list(object))
  .validate_image_data(object)

  # Open communication with react app through window
  jsCode <- "shinyjs.getWindow = function(){Shiny.setInputValue('myVal', window.MyLib)}"

  # Create UI for app
  ui <- fluidPage(

    useShinyjs(),
    extendShinyjs(text = jsCode, functions = c("getWindow")),
    actionButton("help", "Help"),
    actionButton("quit", "Quit & Save"),
    tableOutput("table"),
    paperOutput("paperWidget", height = height, width = width),
    uiOutput("HelpBox")

  )

  # Create server side
  server <- function(input, output, session) {

    # Fire event every 50ms
    autoCheck <- reactiveTimer(50)

    # Observe timer and return latest transformations
    # this is used to communicate transformations from
    # the app back to R
    observe({
      autoCheck()
      shinyjs::js$getWindow()
    })

    # Send image data to widget
    output$paperWidget <- renderPaper({
      paper(data = object)
    })

    # Fill table when values are passed from react app
    output$table <- renderTable(do.call(rbind, lapply(input$myVal, as.data.frame)))

    # Helpt text
    output$HelpBox = renderUI({
      if (input$help %% 2){
        fluidRow(column(12, helpText(h4("How to use the tool"), hr())),
                 column(12, helpText(h5(strong("NB:"), " The dashed lines highlight the borders of the canvas area. Only part of images that are inside this area will be kept."))),
                 column(12, column(3, helpText(p(strong("Select images")))),
                        column(9, p("use the blue image buttons above the viewer to select images to show"))),
                 column(12, column(3, helpText(p(strong("Move to front")))),
                        column(9, p("click on an image to move it to the foreground"))),
                 column(12, column(3, helpText(p(strong("Move")))),
                        column(9, p("click & drag to move an image"))),
                 column(12, column(3, helpText(p(strong("Scale")))),
                        column(9, p("hold & drag the top right corner to scale an image"))),
                 column(12, column(3, helpText(p(strong("Rotate")))),
                        column(9, p("hold the ", code("SHIFT"), " key + hold & drag the top right corner to rotate an image"))),
                 column(12, column(3, helpText(p(strong("Increase transparency")))),
                        column(9, p("hold the ", code("q"), " key and click on an image to increase its transparency"))),
                 column(12, column(3, helpText(p(strong("Decrease transparency")))),
                        column(9, p("hold the ", code("w"), " key and click on an image to decrease its transparency"))),
                 column(12, column(3, helpText(p(strong("Reset image")))),
                        column(9, p("hold the ", code("r"), " key and click on an image to reset it back to its original shape"))),
                 column(12, column(3, helpText(p(strong("Flip horisontally")))),
                        column(9, p("hold the ", code("f"), " key and click on an image to flip it horisontally")))
        )
      } else {
        return()
      }
    })

    # Run code on exit
    observeEvent(input$quit, {
      stopApp(returnValue = do.call(rbind, lapply(input$myVal, as.data.frame)))
    })

  }

  transformations <- runApp(list(ui = ui, server = server))

  return(transformations |> as_tibble())
}


#' Create a react app for paper JS in R
#'
#' Provided a list of images, this function is used to open an interactive app
#' built with react and paper JS.
#'
#' @param data A list of images prepared with \code{\link{prep_image}}
#' @param width Width of component
#' @param height height of component
#' @param elementId Component element ID
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
    install.packages("htmlwidgets")
  }

  # validate data
  # if (!class(data) == "list") stop(sprintf("invalid class '%s' of data object", class(data)))
  # if (!all(names(data) == c("data", "dimx", "dimy", "length"))) stop("Invalid structure of data.")
  # if (data$length == 0) stop("The length of the buffer is 0.")

  content <- reactR::component(
    "Paper",
    list(data = data)
  )

  # describe a React component to send to the browser for rendering.
  component <- reactR::reactMarkup(content)

  # create widget
  htmlwidgets::createWidget(
    name = 'paper',
    component,
    width = width,
    height = height,
    package = 'paper',
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
#'   is useful if you want to save an expression in a variable.
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
    install.packages("htmlwidgets")
  }

  htmlwidgets::shinyWidgetOutput(outputId, 'paper', width, height, package = 'paper')
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
    install.packages("htmlwidgets")
  }

  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, paperOutput, env, quoted = TRUE)
}


#' Function used to prepare images for paper JS react app
#'
#' @param input An object of class 'magick-image' or a path
#' to an image in png or jpeg format
#' @param width width of image sent to react app
#'
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @importFrom rlang abort
#' @importFrom glue glue
#' @importFrom tools file_ext
#' @importFrom magick image_read image_scale image_data
#'
#' @return a list with an array buffer, the image dimensions and
#' the length of the array buffer
#'
#' @export
#'
prep_image <- function (
    input,
    width = 256
) {

  # Set global variables to NULL
  R <- G <- B <- A <- key <- ord <- NULL

  # Check width
  if (between(x = width, left = 256, right = 512))
    if (width < 256) abort(glue("The image width ({width}) has to be between 256 and 512 pixels"))

  # Check input
  stopifnot(
    inherits(input, what = c("magick-image", "character"))
  )
  if (inherits(input, what = "magick-image")) {
    im <- input
  }
  if (inherits(input, what = "character")) {
    if (length(input) > 1) abort(glue("{length(input)} images provided, expected 1."))
    if (!file.exists(input)) stop("Invalid path to image.")
    file.ext <- file_ext(input)
    if (!file.ext %in% c("png", "jpg", "jpeg")) stop(sprintf("Invalid file type '.%s'", file.ext))
    im <- image_read(input)
  }

  # Read image, scale it to have a specified width and
  # convert it to an array
  im <- im |>
    image_scale(paste0(width)) |>
    image_data() |>
    as.integer()

  # Create an array buffer
  # ------------------------------------
  # Each pixel is defined by RGBA values and the values of d
  # are sorted by pixel, i.e. R1, G1, B1, A1, R2, G2, B2, A2, ...
  # where 1, 2, ... corresponds to the pixels
  d <- data.frame(R = as.vector(t(im[, , 1])),
                  G = as.vector(t(im[, , 2])),
                  B = as.vector(t(im[, , 3])),
                  A = 255) %>%
    mutate(ord = 1:n()) |>
    pivot_longer(c(R, G, B, A), names_to = "key", values_to = "value") |>
    mutate(key = factor(key, levels = c("R", "G", "B", "A"))) |>
    arrange(ord, key)

  # Create a list object containing the array buffer and
  # some additional properties such as the x, y dimensions and
  # the length of the buffer
  d <- list(data = d$value,
            dimx = ncol(im),
            dimy = nrow(im),
            length = nrow(im)*ncol(im)*4)

  return(d)
}


#' Validate data sent to \code{\link{RunAlignment}}
#'
#' @param data A list of images prepared with \code{\link{prep_image}}
#'
#' @noRd
.validate_image_data <- function (data) {
  for (i in seq_along(data)) {
    x <- data[[i]]
    if (!all(names(x) == c("data", "dimx", "dimy", "length")))
      abort(glue("Image {i} is invalid"))
  }
}
