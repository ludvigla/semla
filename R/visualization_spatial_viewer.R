#' @include generics.R
#' @include checks.R
#'
NULL

# TODO: update docs
#' Interactive spatial feature viewer
#'
#' A more detailed tutorial can be found on the `STUtility2` website.
#'
#' @param object A Seurat object
#' @param datadir A directory spatial data and image tiles
#' @param selected_features A character vector of features to select for viewer
#' @param sampleIDs A vector of section IDs to use for the viewer
#' @param container_width,container_height Set height and width of container
#' @param overwrite_network_json Logical specifying if the JSON file
#' containing the spatial network should be overwritten after completion
#' @param verbose Print messages
#' @inheritParams file_server
#'
#' @import rlang
#' @import glue
#' @import cli
#' @import shiny
#'
#' @family spatial-visualization
#'
#' @export
#'
FeatureViewer <- function (
    object,
    datadir = NULL,
    selected_features = NULL,
    sampleIDs = 1,
    host = "127.0.0.1",
    port = 8080L,
    container_width = 800,
    container_height = 800,
    verbose = TRUE
) {

  # Set global variables to NULL
  sampleID <- sampleID <- NULL

  # Check Seurat object
  .check_seurat_object(object)

  # Check tile paths
  if (!is.null(datadir)) {
    if (verbose) cli_alert_info("Got data directory {datadir}")
    if (verbose) cli_alert_info("Checking for required files ...")
    for (i in sampleIDs) {
      tilepath <- paste0(datadir, paste0("/tiles", sampleID))
      if (!file.exists(tilepath)) {
        abort(glue("Path {tilepath} is missing for sample {sampleID}"))
      }
      infopath <- paste0(datadir, paste0("/image_info", sampleID, ".json"))
      if (!file.exists(infopath)) {
        abort(glue("Path {infopath} is missing for sample {sampleID}"))
      }
      coordpath <- paste0(datadir,  paste0("/coords_Visium_", sampleID, ".json"))
      if (!file.exists(coordpath)) {
        abort(glue("Path {coordpath} is missing for sample {sampleID}"))
      }
    }
    datapath <- datadir
    clean_after_close <- FALSE
  } else {
    cli_alert_info("Attempting to tile H&E image(s)")
    imgs <- GetStaffli(object)@imgs[sampleIDs]
    for (i in seq_along(imgs)) {
      if (!file.exists(imgs[i])) {
        abort(glue("{imgs[i]} is not a valid path. Update path to @imgs in 'Staffli' slot"))
      }
      cli_alert("  Tiling sample {i} H&E image to temporary directory")
      dirs <- TileImage(im = image_read(imgs[i]), sampleID = sampleIDs[i])
      datapath <- dirs$datapath
      cli_alert("  Exporting Visium coordinates for sample {i}")
      export_coordinates(object = object, sampleID = sampleIDs[i], outdir = datapath, verbose = verbose)
    }
    clean_after_close <- TRUE
  }

  # Fetch data
  #se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "STUtility2"))
  #expr_data <- GetAssayData(object, slot = "data")
  all_features <- c(rownames(GetAssayData(object, slot = "data")), sapply(object@reductions, function(x) {
    colnames(x@cell.embeddings)
  }))

  # Fetch categorical data
  .categorical_data <- object[[]] |>
    select_if(function(col) {is.factor(col) | is.character(col)})
  categorical_features <- colnames(.categorical_data)

  # Define a list of colors for each category
  .colors_saved <- sapply(.categorical_data, function(x) {
    if (inherits(x, what =  "factor")) {
      setNames(.set_colors(length(levels(x))), nm = levels(x))
    } else {
      setNames(.set_colors(length(unique(x))), nm = unique(x))
    }
  })

  # Import beakr
  if (!requireNamespace("beakr"))
    install.packages("beakr")
  if (!requireNamespace("colourpicker"))
    install.packages("colourpicker")

  # Check input
  stopifnot(
    inherits(datapath, what = "character"),
    dir.exists(datapath),
    inherits(host, what = "character") & length(host) == 1,
    inherits(port, what = c("numeric", "integer")) & length(port) == 1
  )

  # Start file server
  if (verbose) cli_alert_info("Starting static file server")
  beakr::stopAllServers()
  fs <- try({file_server(hostDir = datapath, host = host, port = port)})
  if (verbose) print(fs)
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
  # This will make it possible to retrieve lasso selections
  jsCode_getselection <- "shinyjs.getSelection = function(){Shiny.setInputValue('selbarcodes', window.curSelection)}"

  # Adds a small trick to make an actionButton respond o ENTER
  jscode <- '
    $(function() {
      var $els = $("[data-proxy-click]");
      $.each(
        $els,
        function(idx, el) {
          var $el = $(el);
          var $proxy = $("#" + $el.data("proxyClick"));
          $el.keydown(function (e) {
            if (e.keyCode == 13) {
              $proxy.click();
            }
          });
        }
      );
    });
    '

  # Create UI for app
  ui <- fluidPage(

    # Add custom js code to retrieve lasso selection
    useShinyjs(),
    extendShinyjs(text = jsCode_getselection, functions = c("getSelection")),

    # Add custom jscode to make it possible to save text on ENTER
    shiny::tags$head(shiny::tags$script(HTML(jscode))),

    # UI for main panel including the feature viewer widget
    ftrviewerOutput("ftrviewerWidget", height = container_height, width = container_width),
    uiOutput("HelpBox"),

    # UI for floating sidebar
    absolutePanel(
      actionButton("quit", "Quit & Save"), # Quits application and save changes
      actionButton("help", "Help"), # Opens up a help menu
      fastSliderInput("opacity", "Opacity", min = 0, max = 1, value = 1, step = 0.05), # Slider for opacity values
      checkboxInput("lasso", "Lasso", value = FALSE), # Activates a lasso selection tool
      selectizeInput("feature", label = "Feature", choices = NULL), # numeric features

      # This panel only opens up if a numeric feature is selected
      conditionalPanel(
        condition = "!output.panelStatus",
        checkboxInput("scalealpha", "Scale alpha"),
        fastSliderInput("trim", "Trim", min = 0, max = 1, value = c(0, 1), step = 0.01),
        selectizeInput("colorscale", "Colorscale", choices = .color_scales(info = TRUE), options = list(create = TRUE))
      ),

      selectizeInput("category", "Select input", choices = categorical_features, options = list(create = TRUE)),

      # This panel only opens up if a categorical feature is selected
      conditionalPanel(
        condition = "output.panelStatus",
        tagAppendAttributes(
          textInput("label", "Label"),
          `data-proxy-click` = "label_save"
        ),
        colourpicker::colourInput("color", label = "Select color", value = "#FFA500")
      ),
      actionButton("label_save", "", style = "visibility: hidden;"), # Used for the actionButton ENTER press trick

      top = "100px",
      right = "20px",
      width = "150px",
      draggable = TRUE
    )

  )

  # Create server side
  server <- function(input, output, session) {

    # Select type of feature
    # Create reactive value (responds on input)
    output$panelStatus <- reactive({
      rv$lastBtn == "category"
    })
    outputOptions(output, "panelStatus", suspendWhenHidden = FALSE)
    rv <- reactiveValues(lastBtn = "feature", curFeature = "gene1")
    # # Set last button to "category" when a categorical feature is selected
    observeEvent(input$category, {
      if (input$category > 0 ) {
        rv$lastBtn = "category"
        rv$curFeature = input$category
        if (!input$category %in% colnames(.categorical_data)) {
          .categorical_data[, input$category] <<- "background"
          if (!rv$curFeature %in% names(.colors_saved)) {
            .colors_saved[[rv$curFeature]] <<- setNames(.set_colors(n = 1), nm = "background")
          }
        }
        rv$values = .categorical_data[, input$category]
        rv$opacities = rep(1, length(rv$values))
        if (inherits(rv$values, what =  "character")) {
          rv$levels <- unique(rv$values)
        } else {
          rv$levels <- levels(rv$values)
        }
        rv$colors <- .colors_saved[[rv$curFeature]]
        rv$isNumeric = FALSE
      }
    })
    # Set last button to "feature" when a numeric feature is selected
    observeEvent(input$feature, {
      if (input$feature > 0 ) {
        rv$lastBtn = "feature"
        rv$curFeature = input$feature
        rv$values = FetchData(object, vars = input$feature) |> pull(all_of(input$feature))
        rv$range = range(rv$values) |> round(digits = 2)
        rv$opacities = scales::rescale(rv$values, to = c(0, 1), from = rv$range)
        rv$colors = .color_scales(input$colorscale)
        rv$isNumeric = TRUE
      }
    })
    observeEvent(input$colorscale, {
      if (rv$lastBtn == "feature") {
        rv$colors = .color_scales(input$colorscale)
      }
    })

    # Server side updating of available all_features
    updateSelectizeInput(session, 'feature', choices = all_features, selected = all_features[1], server = TRUE)

    # Fire event every 50ms
    autoCheck <- reactiveTimer(50)

    # Observe timer and return latest transformations
    # this is used to communicate transformations from
    # the app back to R
    observe({
      autoCheck()
      shinyjs::js$getSelection()
    })


    # Observe label save
    observeEvent(input$label_save, {
      if (input$label != "") {
        updateCheckboxInput(session, "lasso", value = FALSE)
      }
    })

    observeEvent(input$selbarcodes, {
      if (length(input$selbarcodes) > 1) {
        if (inherits(.categorical_data[, input$category], what = "factor")) {
          tmp <- .categorical_data |>
            select(contains(input$category)) |>
            mutate_if(is.factor, as.character)
          tmp[input$selbarcodes, input$category] <- input$label
          tmp <- tmp |>
            mutate(across(everything(), ~factor(.x)))
          .categorical_data[, input$category] <<- tmp[, input$category]
        } else {
          .categorical_data[input$selbarcodes, input$category] <<- input$label
        }
        if (!input$label %in% rv$levels) {
          tmp_cols <- c(.colors_saved[[rv$curFeature]][rv$levels], setNames(input$color, nm = input$label))
          rv$levels <- c(rv$levels, input$label)
          rv$colors <- tmp_cols
          .colors_saved[[rv$curFeature]] <<- tmp_cols
        } else {
          tmp_cols <- .colors_saved[[rv$curFeature]]
          tmp_cols <- c(tmp_cols[-match(input$label, names(tmp_cols))], setNames(input$color, nm = input$label))
          .colors_saved[[rv$curFeature]] <<- tmp_cols
          rv$colors <- tmp_cols
        }
        rv$values <- .categorical_data[, input$category]
        updateTextInput(session, "label", value = "")
      }
    })

    # Send image data to widget
    output$ftrviewerWidget <- renderFtrviewer({
      ftrviewer(values = rv$values,
                range = .trim_range(x = rv$range, minCutoff = input$trim[1], maxCutoff = input$trim[2]),
                levels = rv$levels,
                opacities = rv$opacities,
                scaleByOpacity = input$scalealpha,
                colors = rv$colors |> unname(),
                isNumeric = rv$isNumeric,
                useLasso = input$lasso,
                opacity = input$opacity)
    })

    # Help text
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

    # Stop app when pressing quite button
    observeEvent(input$quit, {
      stopApp()
    })
  }

  # Run application and return network on quit
  runApp(list(ui = ui, server = server), launch.browser = TRUE)
  if (verbose) {
    cli_alert_info("Retrieved data from application")
  }

  # Stop static file server
  if (verbose) cli_alert_info("Stopping file server")
  beakr::stopServer(fs)

  # Save new categorical variables
  if (verbose) cli_alert_info("Saving changes to meta data")
  initial_categories <- colnames(.categorical_data)
  new_mdata <- object[[]] |>
    select(-contains(initial_categories)) |>
    bind_cols(.categorical_data)
  object@meta.data <- new_mdata

  # Remove temporary directory
  if (clean_after_close) {
    if (verbose) cli_alert_info("Removing temporary directory")
    unlink(x = datapath, recursive = TRUE)
  }

  if (verbose) cli_alert_success("Finished!")
  return(object)
}


#' Trim values to quantiles
#'
#' @noRd
.trim_range <- function (
  x,
  minCutoff = 0,
  maxCutoff = 1
) {
  if (!is.null(x)) {
    x_d <- x - min(x)
    low <- x_d[2]*minCutoff + min(x)
    high <- x_d[2]*maxCutoff + min(x)
    return(c(low, high))
  }
}


#' Select or show color palettes
#'
#' @importFrom grDevices colorRampPalette
#'
#'
#' @noRd
.color_scales <- function (
  colorscale = "viridis",
  info = FALSE
) {
  colscales <- list(
    viridis = viridis::viridis(n = 50),
    magma = viridis::magma(n = 50, direction = -1),
    heat = colorRampPalette(c("darkblue", "cyan", "yellow", "red", "darkred"))(50),
    spectral = colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "Spectral") |> rev())(50),
    reds = colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Reds"))(50),
    blues = colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Blues"))(50)
  )
  if (info) {
    return(names(colscales))
  } else {
    return(colscales[[colorscale]])
  }
}

#' Set categorical colors
#'
#' @noRd
.set_colors <- function (
  n
) {
  if (n == 1) {
    return("#C8C8C8")
  } else {
    .gg_color_hue(n = n)
  }
}


#' Fast sliderInput
#'
#' Modified version of sliderInput
#'
#' @importFrom shiny sliderInput
#'
#' @noRd
fastSliderInput <- function(...) {
  sld <- sliderInput(...)
  sld$children[[2]]$attribs$`data-immediate` <- "true"
  return(sld)
}


# TODO: fix categories in src js
#' Open a feature viewer react application
#'
#' This function will open a react application that listens can
#' be used to interactively visualize categorical or numeric features
#' as spatial maps.
#'
#' @param host The host address. Defaults to localhost "127.0.0.1"
#' @param port A number for a valid port
#' @param sampleID A section number
#' @param values A vector of numeric or categorical values
#' @param opacities A numeric vector of opacity values
#' @param opacity An integer of length 1 specifying a fixed opacity value
#' @param range A numeric vector of length 2 specifying a range of values (color domain)
#' @param scaleByOpacity A logical specifying if the opacity should be set by
#' `opacities` or `opacity`
#' @param isNumeric A logical specifying if the input is numeric or not
#' @param useLasso A logical specifying if the lasso selection tool should be activated
#' @param levels category levels for coloring of values
#' @param categories not yet implemented
#' @param colors A vector of colors
#' @param container_width,container_height Container width/height in pixels
#' @param elementId The id of the viewer element
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
    levels,
    categories = c("cat1", "cat2", "cat3"),
    colors = RColorBrewer::brewer.pal(n = 9, name = "Spectral"),
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
    package = 'ftrviewer',
    elementId = elementId
  )
}

#' Called by HTMLWidgets to produce the widget's root element.
#' @noRd
widget_html.ftrviewer <- function(id, style, class, ...) {

  # Require htmlwidgets
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    install.packages("htmlwidgets")
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

  htmlwidgets::shinyWidgetOutput(outputId, 'ftrviewer', width, height, package = 'ftrviewer')
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
