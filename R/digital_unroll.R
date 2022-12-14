#' @include generics.R
#' @include checks.R
#'
NULL

#' Cut spatial networks in folded tissues
#'
#' Opens an interactive viewer where a spatial network is visualized on top
#' of an H&E image.
#'
#' @details Each spot is connected to adjacent neighbors by edges and
#' the edges can be cut by holding the SHIFT key while moving the cursor
#' across them. Cut edges can be mended by holding the CTRL key while moving
#' the cursor across them. The aim is to cut edges between spots located in separate
#' layers. The output is a `tbl_graph` object representing the spatial network
#' which can be processed further with \code{\link{AdjustTissueCoordinates}}
#' to perform "digital unrolling".
#'
#' A more detailed tutorial can be found on the `STUtility2` website.
#'
#' @param datadir A directory containing network data and image tiles
#' @param container_width,container_height Set height and width of container
#' @param overwrite_network_json Logical specifying if the JSON file
#' containing the spatial network should be overwritten after completion
#' @param verbose Print messages
#' @inheritParams file_server
#'
#' @import rlang
#' @import glue
#' @import cli
#' @importFrom shiny h2
#'
#' @rdname digital-unroll
#' @family spatial-methods
#'
#' @export
#'
CutSpatialNetwork <- function (
  datadir,
  host = "127.0.0.1",
  port = 8080L,
  container_width = '800px',
  container_height = '800px',
  overwrite_network_json = TRUE,
  verbose = TRUE
) {

  # Set global variables to NULL
  name <- x <- y <- target <- x_end <- y_end <- keep <- NULL

  # Import tidygraph
  if (!requireNamespace("tidygraph"))
    install.packages("tidygraph")
  if (!requireNamespace("beakr"))
    install.packages("beakr")

  # Check input
  stopifnot(
    inherits(datadir, what = "character"),
    dir.exists(datadir),
    inherits(host, what = "character"),
    length(host) == 1,
    inherits(port, what = c("numeric", "integer")),
    length(port) == 1
  )

  # Start file server
  if (verbose) cli_alert_info("Starting static file server")
  beakr::stopAllServers()
  fs <- try({file_server(hostDir = datadir, host = host, port = port)})
  print(fs)
  if (inherits(fs, what = "try-error"))
    abort(c(
      "Failed to start static file server",
      "*" = glue("Make sure to provide a valid host address of class 'character'"),
      "*" = glue("Make sure to provide a valid port of class 'integer'"),
      "*" = glue(
        "Make sure that a static file server is not already running. ",
        "You can stop all servers with beakr::stopAllServers()"
      )
    ))

  # Open communication with react app through window
  jsCode_getdata <- "shinyjs.getNetwork = function(){Shiny.setInputValue('network', window.savedData)}"

  # Create UI for app
  ui <- fluidPage(

    useShinyjs(),
    extendShinyjs(text = jsCode_getdata, functions = c("getNetwork")),
    h2(),
    actionButton("help", "Help", width = "100px"),
    actionButton("quit", "Quit & Save", width = "100px"),
    osdduOutput("osdduWidget", height = container_height, width = container_width),
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
      shinyjs::js$getNetwork()
    })

    # Send image data to widget
    output$osdduWidget <- renderOsddu({
      osddu(port = port, host = host, width = container_width, height = container_height, quit = (input$quit %% 2) == 1)
    })

    # Helpt text
    output$HelpBox = renderUI({
      if (input$help %% 2){
        fluidRow(column(12, helpText(h4("How to use the tool"), hr())),
                 column(12, column(3, helpText(p(strong("Cut edges")))),
                        column(9, p("Hold SHIFT and move the cursor across an edge to cut it. Cut edges are colored red."))),
                 column(12, column(3, helpText(p(strong("Glue edges")))),
                        column(9, p("Hold CTRL and move the cursor across an edge to glue it together. Intact edges are colored black."))))
      } else {
        return()
      }
    })

    observeEvent(input$network, {
      stopApp(returnValue = input$network)
    })
  }

  # Run application and return network on quit
  vals <- runApp(list(ui = ui, server = server))
  if (verbose) {
    cli_alert_info("Retrieved data from application")
    cli_alert_info("Converting data to a tidygraph object")
  }

  # Format network
  nodes <- do.call(rbind, lapply(vals$nodes, function(x) as_tibble(x)))
  nodes_filtered <- nodes |>
    dplyr::select(name, x, y)
  links <- do.call(rbind, lapply(vals$links, function(x) as_tibble(x)))
  links_filtered <- links |>
    dplyr::select(source, target, x, y, x_end, y_end, keep)
  tidygr <- tidygraph::tbl_graph(nodes = nodes_filtered, edges = links_filtered)

  # Overwrite json
  if (overwrite_network_json) {
    if (verbose) cli_alert_info("Overwriting spatial network file with new results")
    # Put nodes and edges into a list
    data <- list(nodes = nodes, links = links)
    data_json <- data |>
      write_json(auto_unbox = TRUE,
                 path = file.path(datadir, "data_Visium.json"))
  }

  # Stop static file server
  if (verbose) cli_alert_info("Stopping file server")
  beakr::stopServer(fs)

  if (verbose) cli_alert_success("Finished!")
  return(tidygr)
}


#' Create a react app for digital unrolling
#'
#' This function is used to start an interactive widget for
#' digital unrolling. It requires a static files server to be
#' hosted in order to find load the necessary files.
#'
#' @noRd
osddu <- function (
    port = "8080",
    host = "127.0.0.1",
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
    list(port = port, host = host, width = width, height = height, quit = quit)
  )
  component <- reactR::reactMarkup(content)

  # create widget
  htmlwidgets::createWidget(
    name = 'osddu',
    component,
    width = width,
    height = height,
    package = 'osddu',
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
  htmlwidgets::shinyWidgetOutput(outputId, 'osddu', width, height, package = 'osddu')
}

#' @rdname osddu-shiny
#' @export
renderOsddu <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, osdduOutput, env, quoted = TRUE)
}

