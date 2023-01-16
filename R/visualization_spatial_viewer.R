#' @include generics.R
#' @include checks.R
#'
NULL

# TODO: fix bug when only 1 spot is selected
# TODO: fix feature range when 1 sample is selected
#' Interactive spatial feature viewer
#'
#' `FeatureViewer` opens up an interactive shiny application where
#' one can zoom and pan the H&E image while overlaying and color selected features.
#' It is also possible to add new categorical features or modify existing
#' categorical features using a lasso tool.
#'
#' The viewer requires a tiled H&E image along with some additional image data to work.
#'
#' A more detailed tutorial can be found on the `STUtility2` website. You can also
#' get more detailed instructions by pressing the help icon in the app.
#'
#' @param object A `Seurat` object
#' @param slot A slot to use for Assay data
#' @param datadir A directory spatial data and image tiles
#' @param selected_features A character vector of features to select for viewer
#' @param sampleIDs A vector of section IDs to use for the viewer
#' @param container_width,container_height Set height and width of container
#' @param verbose Print messages
#'
#' @return A `Seurat` object
#'
#' @author Ludvig Larsson
#'
#' @inheritParams file_server
#'
#' @seealso TileImage
#'
#' @import rlang
#' @import glue
#' @import cli
#' @import shiny
#'
#' @family spatial-visualization-methods
#'
#' @export
#'
FeatureViewer <- function (
    object,
    slot = "data",
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
  sampleID <- barcode <- NULL

  if (!requireNamespace("shinyBS"))
    install.packages("shinyBS")

  # Check Seurat object
  .check_seurat_object(object)

  # Fetch all features
  all_features <- c(rownames(GetAssayData(object, slot = slot)),
                    sapply(object@reductions, function(x) {
                        colnames(x@cell.embeddings)
                      })) |> unlist() |> unique()
  mdata_num_features <- object[[]] |>
    select_if(is.numeric) |>
    colnames()
  all_features <- c(all_features, mdata_num_features)
  if (sum(duplicated(all_features)) > 0) {
    cli_alert_danger("Some features were found in multiple places")
    all_features <- unique(all_features)
  }

  # Validate input
  all_features <- .check_features(selected_features = selected_features, all_features = all_features)
  available_sampleIDs <- GetStaffli(object)@image_info$sampleID |> as.integer()
  if (any(!sampleIDs %in% available_sampleIDs))
    abort(glue("Invalid sampleIDs: {paste(sampleIDs[!sampleIDs %in% available_sampleIDs], collapse = ', ')}"))


  # Check tile paths
  params <- .check_tile_paths(object = object, datadir = datadir, sampleIDs = sampleIDs, verbose = verbose)
  datapath <- params$datapath
  clean_after_close <- params$clean_after_close

  # Get sample IDs barcodes
  barcodes <- GetStaffli(object)@meta_data |>
    select(barcode, sampleID)
  barcodes <- split(barcodes$barcode, barcodes$sampleID)

  # Fetch categorical data
  .categorical_data <- object[[]] |>
    select_if(function(col) {is.factor(col) | is.character(col)})
  .categorical_data[is.na(.categorical_data)] <- "NA"
  categorical_features <- colnames(.categorical_data)

  # Define a list of colors for each category
  .colors_saved <- lapply(.categorical_data, function(x) {
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
    conditionalPanel(
      condition = "output.panelStatus",
      absolutePanel(
        plotOutput("legend", height = "200px"),
        width = "60px",
        top = "50px",
        left = paste0(container_width - 50, "px"),
        draggable = TRUE
      )
    ),
    conditionalPanel(
      condition = "!output.panelStatus",
      uiOutput("catlegend")
    ),
    shinyBS::bsModal("HelpBox", trigger = "help", title = "Usage instructions", size = "large",
            column(12, column(12, p())),
            column(12, column(12, p())),
            column(12, column(4, p(strong("Zoom and pan"))),
                   column(8, p("Scroll to zoom in and out or drag the viewer by holding and dragging the cursor. ",
                               "This feature is inactivated when the lasso tool is in use."))),
            column(12, column(4, p(strong("Numeric features"))),
                   column(8, p("Select a numeric feature from the ", code("Feature"), " drop-down list to color spots with a continuous color scale."))),
            column(12, column(4, p(strong("Categorical features"))),
                   column(8, p("Select a categorical feature from the ", code("Category"), " drop-down list to color spots by category labels. ",
                               "Update the viewer to show the selected category by clicking on the ", code("Category"),
                               " drop-down list and press ", code("ENTER"), "."))),
            column(12, column(4, p(strong("Add new categories"))),
                   column(8, p("Click on the ", code("Category"), " drop-down list, write a new name and press ", code("ENTER"), "."))),
            column(12, column(4, p(strong("Lasso selection"))),
                   column(8, p("Click on ", icon("fa-solid fa-draw-polygon", verify_fa = FALSE),
                               " and select spots with the cursor. ",
                               "Hold ", code("SHIFT"), "to deselect spots.",
                               "To save the selection, make sure that a ", code("Category"),
                               " is selected (or add a new one), write a new ", code("Label"),
                               " name and press ", code("ENTER"), ". ",
                               "Use the color picker tool (", code("Select color"), ") to select a color for the ",
                               code("Label"), " before pressing ", code("ENTER"), "."))),
            column(12, column(4, p(strong("Scale alpha"))),
                   column(8, p("Click on the ", code("Scale alpha"), " tick box to add opacity to spots ",
                               "proportional to the numeric feature values."))),
            column(12, column(4, p(strong("Center zero"))),
                   column(8, p("Click on the ", code("Center zero"), " to center the color scale at 0. ",
                               "If this option is used, use a divergent color palette from ",
                               code("Colorscale"), " such as 'RdBu'."))),
            column(12, column(4, p(strong("Trim values"))),
                   column(8, p("Adjust the  ", code("Trim"), " sliders to trim the lower/upper bounds of numeric feature values "))),
            column(12, column(4, p(strong("Scale alpha"))),
                   column(8, p("Click on ", icon("fa-regular fa-floppy-disk", verify_fa = FALSE),
                               " to save changes and quit the application. If any new selections have been added, these will be ",
                               " stored in the 'meta.data' slot of the returned ", code("Seurat"), " object. ",
                               " The results will NOT be returned if the app is closed from the R session.")))),

    # UI for floating sidebar
    absolutePanel(
      shinyBS::bsButton("lasso", "", icon = icon("fa-solid fa-draw-polygon", verify_fa = FALSE)),
      shinyBS::bsTooltip(id = "lasso", title = "lasso selection", placement = "right", trigger = "hover"),
      shinyBS::bsButton("help", "", icon = icon("fa-solid fa-question", verify_fa = FALSE)), # Opens up a help menu
      shinyBS::bsButton("quit", "", icon = icon("fa-regular fa-floppy-disk", verify_fa = FALSE)), # Quits application and save changes
      shinyBS::bsTooltip(id = "quit", title = "save & quit", placement = "right", trigger = "hover"),
      fluidRow(),
      selectizeInput("sample", "Sample", selected = sampleIDs[1], choices = sampleIDs),
      fastSliderInput("opacity", "Opacity", min = 0, max = 1, value = 1, step = 0.05), # Slider for opacity values
      selectizeInput("feature", label = "Feature", selected = selected_features[1], choices = NULL), # numeric features

      # This panel only opens up if a numeric feature is selected
      conditionalPanel(
        condition = "output.panelStatus",
        checkboxInput("scalealpha", "Scale alpha"),
        checkboxInput("centerzero", "Center zero"),
        sliderInput("trim", "Trim", min = 0, max = 1, value = c(0, 1), step = 0.01),
        selectizeInput("colorscale", "Colorscale", choices = .color_scales(info = TRUE), options = list(create = TRUE))
      ),

      # selectizeInput("category", "Select input", choices = categorical_features, options = list(create = TRUE)),
      tagAppendAttributes(
        selectizeInput("category", "Category", choices = categorical_features, options = list(create = TRUE)),
        `data-proxy-click` = "category_save",
      ),
      actionButton("category_save", "", style = "visibility: hidden;"),

      # This panel only opens up if a categorical feature is selected
      conditionalPanel(
        condition = "!output.panelStatus",
        tagAppendAttributes(
          textInput("label", "Label"),
          `data-proxy-click` = "label_save"
        ),
        colourpicker::colourInput("color", label = "Select color", value = "#FFA500"),
        actionButton("label_save", "", style = "visibility: hidden;") # Used for the actionButton ENTER press trick
      ),

      top = "60px",
      left = paste0(container_width + 20, "px"),
      width = "150px",
      draggable = TRUE
    )

  )

  # Create server side
  server <- function(input, output, session) {

    # Select type of feature
    # Create reactive value (responds on input)
    rv <- reactiveValues(lastBtn = "feature", curFeature = selected_features[1],
                         category = 1, sample = sampleIDs[1], lasso = FALSE)
    output$panelStatus <- reactive({
      rv$lastBtn == "feature"
    })
    outputOptions(output, "panelStatus", suspendWhenHidden = FALSE)

    # Switch sample ID
    observeEvent(input$sample, {
      rv$curbarcodes <- barcodes[[input$sample]]
    })

    # # Set last button to "category" when a categorical feature is selected
    observeEvent(toListenCategory(), {
      if (input$category > 0) {
        rv$lastBtn = "category"
        rv$curFeature = input$category
        if (!input$category %in% colnames(.categorical_data)) {
          .categorical_data[, input$category] <<- "NA"
          if (!rv$curFeature %in% names(.colors_saved)) {
            .colors_saved[[rv$curFeature]] <<- setNames(.set_colors(n = 1), nm = "NA")
          }
        }
        rv$values = .categorical_data[rv$curbarcodes, input$category]
        rv$opacities = rep(1, length(rv$values))
        if (inherits(rv$values, what =  "character")) {
          rv$levels <- unique(rv$values)
        } else {
          rv$levels <- levels(rv$values)
        }
        rv$colors <- .colors_saved[[rv$curFeature]][rv$levels]
        rv$isNumeric = FALSE
      }
    })

    # Set last button to "feature" when a numeric feature is selected
    toListenFeatureInput <- reactive({
      list(input$centerzero, input$feature, input$trim)
    })
    observeEvent(toListenFeatureInput(), {
      if (!is.null(input$feature)) {
        if (input$feature %in% all_features) {
          rv$lastBtn = "feature"
          rv$curFeature = input$feature
          rv$values = FetchData(object, cells = rv$curbarcodes, vars = input$feature) |> pull(all_of(input$feature))
          if (input$centerzero) {
            maxAbsVal <- max(abs(rv$values))
            rv$range = c(-maxAbsVal, maxAbsVal) |> round(digits = 2)
          } else {
            rv$range = range(rv$values) |> round(digits = 2)
          }
          if (all(c(0, 1) == input$trim)) {
            trimmed_range <- rv$range
            tmp_vals <- rv$values
          } else {
            trimmed_range <- .trim_range(rv$range, minCutoff = input$trim[1], input$trim[2])
            tmp_vals <- rv$values
            tmp_vals[tmp_vals < trimmed_range[1]] <- trimmed_range[1]
            tmp_vals[tmp_vals > trimmed_range[2]] <- trimmed_range[2]
          }
          rv$opacities = scales::rescale(tmp_vals, to = c(0, 1), from = trimmed_range)
          rv$colors = .color_scales(input$colorscale)
          rv$isNumeric = TRUE
          rv$allClear <- TRUE
        }
      }
    })

    # listen to changes in colorscale and update
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
        #updateCheckboxInput(session, "lasso", value = FALSE)
        rv$lasso <- !rv$lasso
        rv$label <- input$label
        shinyBS::updateButton(session, "lasso", label = "", block = F, style = "default")
      }
    })

    # Trigger category change on ENTER
    observeEvent(input$category_save, {
      rv$category <- rv$category + 1
    })
    toListenCategory <- reactive({
      list(input$category, rv$category)
    })

    # Listen to change in lasso
    observeEvent(input$lasso, {
      rv$lasso <- !rv$lasso
      if (rv$lasso) {
        shinyBS::updateButton(session, "lasso", label = "", block = F, style = "success")
      } else {
        shinyBS::updateButton(session, "lasso", label = "", block = F, style = "default")
      }
    })

    # Listen for changes in selbarcodes and add selection if
    # a label is submitted
    observeEvent(input$selbarcodes, {
      if ((length(input$selbarcodes) > 0) & (input$label != "")) {
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
        if (!rv$label %in% rv$levels) {
          tmp_cols <- c(.colors_saved[[input$category]][rv$levels], setNames(input$color, nm = rv$label))
          rv$levels <- c(rv$levels, rv$label)
          rv$colors <- tmp_cols
          .colors_saved[[input$category]] <<- tmp_cols
        } else {
          tmp_cols <- .colors_saved[[input$category]]
          tmp_cols <- c(tmp_cols[setdiff(names(tmp_cols), rv$label)], setNames(input$color, nm = rv$label))
          .colors_saved[[input$category]] <<- tmp_cols
          rv$colors <- tmp_cols[rv$levels]
        }
        rv$values <- .categorical_data[rv$curbarcodes, input$category]
        updateTextInput(session, "label", value = "")
      }
    })

    # Send image data to widget
    output$ftrviewerWidget <- renderFtrviewer({
      if (length(rv$isNumeric) > 0) {
        ftrviewer(values = rv$values,
                  sampleID = input$sample |> as.integer(),
                  range = .trim_range(x = rv$range, minCutoff = input$trim[1], maxCutoff = input$trim[2]),
                  levels = rv$levels,
                  opacities = rv$opacities,
                  scaleByOpacity = input$scalealpha,
                  colors = rv$colors |> unname(),
                  isNumeric = rv$isNumeric,
                  useLasso = rv$lasso,
                  opacity = input$opacity)
      }
    })

    # Render categorical feature legend on client side
    output$catlegend <- renderUI({
      panel_width = max(sapply(rv$levels, nchar) |> max(), 15)
      absolutePanel(
        plotOutput("legend_cat", height = paste0(length(rv$levels)*20 + 30, "px")),
        width = paste0(panel_width*9 + 10, "px"),
        top = "50px",
        left = paste0(container_width - (panel_width*9 + 10), "px"), #paste0(container_width - (min(panel_width, 15)*10 + 10), "px"),
        draggable = TRUE
      )
    })

    # Render color legend for numerical features
    output$legend <- renderPlot({
      if (!is.null(rv$allClear) & rv$isNumeric) {
        cur_range <- .trim_range(x = rv$range, minCutoff = input$trim[1], maxCutoff = input$trim[2])
        p <- .create_legend(minVal = cur_range[1],
                            maxVal = cur_range[2],
                            levels = rv$levels,
                            colors = rv$colors,
                            isNumeric = rv$isNumeric,
                            title = input$feature)
        return(p)
      }
    }, bg = "transparent")

    # Render color legend for categorical features
    output$legend_cat <- renderPlot({
      if (!is.null(rv$allClear) & !rv$isNumeric) {
        p <- .create_legend(minVal = rv$range[1],
                            maxVal = rv$range[2],
                            levels = rv$levels,
                            colors = rv$colors,
                            isNumeric = rv$isNumeric,
                            title = input$category)
        return(p)
      }
    }, bg = "transparent")

    # Stop app when pressing quit button
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
  new_mdata[new_mdata == "NA"] <- NA
  object@meta.data <- new_mdata

  # Remove temporary directory
  if (clean_after_close) {
    if (verbose) cli_alert_info("Removing temporary directory")
    unlink(x = datapath, recursive = TRUE)
  }

  return(object)
}


#' Check input features
#'
#' @noRd
.check_features <- function (
    selected_features,
    all_features
) {
  if (!is.null(selected_features)) {
    if (!all(selected_features %in% all_features)) {
      # Get intersect
      features_left <- intersect(selected_features, all_features)
      if (length(features_left) == 0) abort("selected_features is invalid, foun 0 features in Seurat object")
      cli_alert_info(col_br_red("{length(selected_features) - length(features_left)} features found in Seruat object"))
      all_features <- features_left
    } else {
      all_features <- selected_features
    }
  }
  return(all_features)
}


#' Check tile paths
#'
#' @param object A `Seurat` object
#' @param datadir A path to a directory with requried input files
#' @param sampleIDs A vector with sampleIDs
#' @param verbose Print messages
#'
#' @noRd
.check_tile_paths <- function (
  object,
  datadir,
  sampleIDs,
  verbose
) {

  # Set global variables to NULL
  sampleID <- NULL

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
    if (verbose) cli_alert_info("Attempting to tile H&E image(s)")
    imgs <- GetStaffli(object)@imgs[sampleIDs]
    for (i in seq_along(imgs)) {
      if (!file.exists(imgs[i])) {
        abort(glue("{imgs[i]} is not a valid path. Update path to @imgs in 'Staffli' slot"))
      }
      if (verbose) cli_alert("  Tiling sample {i} H&E image to temporary directory")
      dirs <- TileImage(im = image_read(imgs[i]), sampleID = sampleIDs[i], verbose = verbose)
      datapath <- dirs$datapath
      if (verbose) cli_alert("  Exporting Visium coordinates for sample {sampleIDs[i]}")
      export_coordinates(object = object, sampleNumber = sampleIDs[i], outdir = datapath, verbose = verbose)
    }
    clean_after_close <- TRUE
  }
  return(list(datapath = datapath, clean_after_close = clean_after_close))
}


#' Create a custom legend
#'
#' @param minVal,maxVal Range for numeric values
#' @param levels Levels for categorical values
#' @param ntiles Number of tiles for color bar
#' @param colors A chaacter vector of colors for color bar
#' @param isNumeric Logical speicifying if numeric data is passed
#' @param title Legend title
#'
#' @import ggplot2
#'
#' @noRd
.create_legend <- function (
  minVal = 0,
  maxVal = 1,
  levels = NULL,
  ntiles = 100,
  colors = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
  isNumeric = TRUE,
  title = "category"
) {

  # Set global variables to NULL
  y <- NULL

  if (!isNumeric) {
    colors <- colors[levels]
    levels <- sapply(levels, function(x) substr(x, start = 1, stop = 20))
    colors <- setNames(colors, nm = levels)

    df <- data.frame(yc = seq(0, length(levels) - 1) + 0.5, row.names = levels)
    new_labels <- sapply(rownames(df), function(x) {
      ifelse(nchar(x) < 15, paste0(x, paste(rep(" ", 15 - nchar(x)), collapse = "")), x)
    })
    p <- ggplot() +
      theme_void() +
      scale_y_continuous(breaks = df$yc, labels = new_labels, position = "right") +
      theme(axis.ticks.y = element_line(),
            plot.background = element_rect(fill = "#FFFFFFAA", colour = NA),
            axis.text.y = element_text(size = 14, hjust = 0.1, colour = "black"),
            plot.title = element_text(size = 16, colour = "black", margin = margin(b = 5)),
            plot.margin = margin(t = 3, r = 3, b = 3, l = 3)) +
      coord_fixed() +
      ggtitle(title)
    for (lvl in rownames(df)) {
      p <- p +
        gg_circle(r = 0.4, xc = 0.5, yc = df[lvl, "yc", drop = TRUE], fill = colors[lvl])
    }
  } else {
  df <- data.frame(y = seq(minVal, maxVal, length.out = ntiles))
  p <- ggplot(df, aes("", y, fill = y)) +
    geom_raster() +
    theme(legend.position = "none") +
    scale_y_continuous(position = "right") +
    theme(axis.text.y = element_text(size = 14, colour = "black"), axis.title = element_blank(),
          plot.background = element_rect(fill = "#FFFFFFAA", colour = NA),
          panel.background = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(size = 16, colour = "black"),
          plot.margin = margin(3, 3, 3, 3),
          panel.grid = element_blank()) +
    scale_fill_gradientn(colours = colors) +
    ggtitle(title)
  }
  return(p)
}

#' Create a circle in ggplot2
#'
#' @param r Circle radius
#' @param xc,yc Circle center
#' @param color Circle stroke color
#' @param fille Circle fill color
#' @param ... parameters passed to `annotate`
#'
#' @import ggplot2
#'
#' @noRd
gg_circle <- function (
  r,
  xc,
  yc,
  color = "black",
  fill = NA,
  ...
) {
  x <- xc + r*cos(seq(0, pi, length.out = 100))
  ymax <- yc + r*sin(seq(0, pi, length.out = 100))
  ymin <- yc + r*sin(seq(0, -pi, length.out = 100))
  annotate("ribbon", x = x, ymin = ymin, ymax = ymax, color = color, fill = fill, ...)
}


#' Trim values to quantiles
#'
#' @param x Value range
#' @param minCutoff,maxCutoff Minimum and maximum cutoff given as
#' values beetween 0-1
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
#' @param colorscale Name of color scale
#' @param info Should names of color palettes be returned?
#'
#' @importFrom grDevices colorRampPalette
#'
#' @noRd
.color_scales <- function (
  colorscale = "viridis",
  info = FALSE
) {
  if (!requireNamespace("viridis"))
    install.packages("viridis")
  colscales <- list(
    viridis = viridis::viridis(n = 50),
    magma = viridis::magma(n = 50, direction = -1),
    heat = colorRampPalette(c("darkblue", "cyan", "yellow", "red", "darkred"))(50),
    spectral = colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "Spectral") |> rev())(50),
    reds = colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Reds"))(50),
    blues = colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Blues"))(50),
    RdBu = colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev())(50)
  )
  if (info) {
    return(names(colscales))
  } else {
    return(colscales[[colorscale]])
  }
}

#' Set categorical colors
#'
#' @param n Number of colors
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
    categories = character(),
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
