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
#' layers. The output is a \code{tbl_graph} object representing the spatial network
#' which can be processed further with \code{\link{AdjustTissueCoordinates}}
#' to perform "digital unrolling".
#'
#' @section Examples:
#' A tutorial can be found on our [package website](https://ludvigla.github.io/semla/).
#' Got to tutorials -> Digital unrolling
#'
#' @param object A \code{Seurat} object created with \code{semla}
#' @param sampleID An integer specifying a tissue section in the \code{Seurat} object
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
#' @family spatial-methods
#' @seealso export_graph
#'
#' @author Ludvig Larsson
#' 
#' @returns A \code{tbl_graph} object representing a spatial network
#'
#' @export
#'
CutSpatialNetwork <- function (
  object,
  datadir = NULL,
  sampleID = 1,
  host = "127.0.0.1",
  port = 8080L,
  container_width = '800px',
  container_height = '800px',
  overwrite_network_json = TRUE,
  verbose = TRUE
) {

  # Set global variables to NULL
  name <- x <- y <- target <- x_end <- y_end <- keep <- NULL

  # Check Seurat object
  .check_seurat_object(object)

  # Import tidygraph
  if (!requireNamespace("tidygraph"))
    abort(glue("Package {cli::col_br_magenta('tidygraph')} is required. Please install it with: \n",
               "install.packages('tidygraph')"))
  if (!requireNamespace("beakr"))
    abort(glue("Package {cli::col_br_magenta('beakr')} is required. Please install it with: \n",
               "install.packages('beakr')"))

  # Check input parameters
  stopifnot(
    inherits(sampleID, what = c("numeric", "integer")),
    length(sampleID) == 1,
    inherits(host, what = "character"),
    length(host) == 1,
    inherits(port, what = c("numeric", "integer")),
    length(port) == 1
  )

  # Validate sampleID
  available_sampleIDs <- GetStaffli(object)@image_info$sampleID |> as.integer()
  if (!sampleID %in% available_sampleIDs)
    abort(glue("Invalid sampleID: {sampleID}"))

  # Check tile paths
  if (!is.null(datadir)) {
    if (!(inherits(datadir, what = "character") & (length(datadir) == 1)))
      abort(glue("Invalid datadir {datadir}, expected a 'character' of length 1"))
    if (!dir.exists(datadir)) abort(glue("Directory {datadir} does not exist"))

    if (verbose) cli_alert_info("Got data directory {datadir}")
    if (verbose) cli_alert_info("Checking for required files ...")
    tilepath <- paste0(datadir, paste0("/tiles", sampleID))
    if (!file.exists(tilepath)) {
      abort(glue("Path {tilepath} is missing for sample {sampleID}"))
    }
    infopath <- paste0(datadir, paste0("/image_info_", sampleID, ".json"))
    if (!file.exists(infopath)) {
      abort(glue("Path {infopath} is missing for sample {sampleID}"))
    }
    networkpath <- paste0(datadir,  paste0("/network_Visium_", sampleID, ".json"))
    if (!file.exists(networkpath)) {
      abort(glue("Path {networkpath} is missing for sample {sampleID}"))
    }
    datapath <- datadir
    clean_after_close <- FALSE
  } else {
    cli_alert_info("Attempting to tile H&E image(s)")
    img <- GetStaffli(object)@imgs[sampleID]
    if (!file.exists(img)) {
      abort(glue("{img} is not a valid path. Update path to @imgs in 'Staffli' slot"))
    }
    cli_alert("  Tiling sample {sampleID} H&E image to a temporary directory")
    dirs <- TileImage(im = image_read(img), sampleID = sampleID)
    datapath <- dirs$datapath
    cli_alert("  Exporting spatial network for sample {sampleID}")
    export_graph(object = object, sampleID = sampleID, outdir = datapath, verbose = verbose)
    clean_after_close <- TRUE
  }

  # Start file server
  if (verbose) cli_alert_info("Starting static file server")
  beakr::stopAllServers()
  fs <- try({file_server(hostDir = datapath, host = host, port = port)})
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
      osddu(sampleID = sampleID,
            host = host,
            port = port,
            width = container_width,
            height = container_height,
            quit = (input$quit %% 2) == 1)
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
  vals <- runApp(list(ui = ui, server = server), launch.browser = TRUE)
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
                 path = file.path(datapath, paste0("network_Visium_", sampleID, ".json")))
  }

  # Stop static file server
  if (verbose) cli_alert_info("Stopping file server")
  beakr::stopServer(fs)

  # Remove temporary directory
  if (clean_after_close) {
    if (verbose) cli_alert_info("Removing temporary directory")
    unlink(x = datapath, recursive = TRUE)
  }

  if (verbose) cli_alert_success("Finished!")
  return(tidygr)
}
