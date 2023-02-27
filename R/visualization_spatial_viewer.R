#' @include generics.R
#' @include checks.R
#'
NULL


#' Interactive spatial feature viewer
#'
#' \code{FeatureViewer} opens up an interactive shiny application where
#' one can zoom and pan the H&E image while overlaying and color selected features.
#' It is also possible to add new categorical features or modify existing
#' categorical features with a lasso tool.
#'
#' The viewer requires a tiled H&E image along with some additional image data to work. By default,
#' the function will try to create these files and export them to a temporary directory, but it
#' is also possible to export the files before running the app and provide the path to the
#' data directory with \code{datadir} (see \code{\link{ExportDataForViewer}}).
#'
#' A more detailed tutorial can be found on the \code{semla} website. You can also
#' get more detailed instructions by pressing the help icon in the app.
#'
#' @param object A \code{Seurat} object
#' @param slot A slot to use for Assay data
#' @param datadir A directory spatial data and image tiles
#' @param selected_features A character vector of features to select for viewer
#' @param sampleIDs An integer vector of section IDs to use for the viewer. All sections
#' will be used by default.
#' @param custom_color_palettes A names list of color vectors to use as custom color palettes
#' @param categorical_colors A named list of character vectors with color names. The name of
#' each character vector should correspond to a categorical variable in the meta.data slot of the
#' \code{Seurat} object. Each character vector should be named where each name corresponds to a
#' label of the category.
#' @param container_width,container_height Set height and width of container
#' @param verbose Print messages
#'
#' @family feature-viewer-methods
#'
#' @return A \code{Seurat} object
#'
#' @author Ludvig Larsson
#'
#' @inheritParams file_server
#'
#' @seealso \code{\link{TileImage}},\code{\link{export_coordinates}}
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
    sampleIDs = NULL,
    host = "127.0.0.1",
    port = 8080L,
    custom_color_palettes = NULL,
    categorical_colors = NULL,
    container_width = 800,
    container_height = 800,
    verbose = TRUE
) {

  # Set global variables to NULL
  sampleID <- barcode <- pxl_col_in_fullres <- pxl_row_in_fullres <- NULL

  if (!requireNamespace("shinyBS"))
    install.packages("shinyBS")

  # Create link to javascript and css files for package
  # This has to be loaded before running the app to make sure
  # that the BS UI elements work properly
  shiny::addResourcePath("sbs", system.file("www", package = "shinyBS"))

  # Check Seurat object
  .check_seurat_object(object)

  # Check custom_color_palettes
  .color_palettes <- sapply(.color_scales(info = TRUE), function(nm) {
    .color_scales(colorscale = nm)
  }) |> setNames(nm = .color_scales(info = TRUE))
  if (!is.null(custom_color_palettes)) {
    .check_custom_palettes(custom_color_palettes)
    .color_palettes <- c(custom_color_palettes, .color_palettes)
  }

  # Fetch all features
  # This will include all feature names of the DefaultAssay + slot,
  # the dimensionality reduction vector names
  all_features <- c(rownames(GetAssayData(object, slot = slot)),
                    sapply(object@reductions, function(x) {
                        colnames(x@cell.embeddings)
                      })) |> unlist() |> unique()
  # Select all numeric features from the meta.data slot
  mdata_num_features <- object[[]] |>
    select_if(is.numeric) |>
    colnames()
  all_features <- c(all_features, mdata_num_features)
  # If there are multiple features with the same names,
  # throw a warning and make the feature names unique
  if (sum(duplicated(all_features)) > 0) {
    cli_alert_danger("Some features were found in multiple places")
    all_features <- unique(all_features)
  }

  # Validate input
  all_features <- .check_features(selected_features = selected_features, all_features = all_features)
  available_sampleIDs <- GetStaffli(object)@image_info$sampleID |> as.integer()
  sampleIDs <- sampleIDs %||% available_sampleIDs # Select all sampleIDs is sampleIDs=NULL
  if (any(!sampleIDs %in% available_sampleIDs))
    abort(glue("Invalid sampleIDs: {paste(sampleIDs[!sampleIDs %in% available_sampleIDs], collapse = ', ')}"))

  # Fetch categorical data
  # Only select character or factor columns from meta.data slot
  .categorical_data <- object[[]] |>
    select_if(function(col) {is.factor(col) | is.character(col)})
  .categorical_data[is.na(.categorical_data)] <- "NA"
  categorical_features <- colnames(.categorical_data)

  # Check custom_color_palettes
  if (!is.null(categorical_colors)) {
    .feature_viewer_colors <- .check_categorical_colors(categorical_colors, categorical_features, .categorical_data)
  } else {
    # Define a list of colors for each category
    # .feature_viewer_colors will work as a dictionary for the app to select colors from
    .feature_viewer_colors <- lapply(.categorical_data, function(x) {
      if (inherits(x, what =  "factor")) {
        setNames(.set_colors(length(levels(x))), nm = levels(x))
      } else {
        setNames(.set_colors(length(unique(x))), nm = unique(x))
      }
    })
  }

  # Check tile paths
  params <- .check_tile_paths(object = object, datadir = datadir, sampleIDs = sampleIDs, verbose = verbose)
  datapath <- params$datapath
  
  # Load subset barcodes
  image_info <- GetStaffli(object)@image_info |>
    filter(sampleID %in% sampleIDs)
  
  # Check if image has been padded
  if ("pad" %in% colnames(image_info)) {
    pad <- strsplit(image_info[image_info$sampleID %in% sampleIDs, ]$pad, "x") |> unlist() |> as.integer() |> split(f = rep(sampleIDs, 4))
    spatial_coords <- GetStaffli(object)@meta_data |> 
      group_by(sampleID) |> 
      group_split()
    spatial_coords <- lapply(seq_along(spatial_coords), function(i) {
      max_imwidth <- image_info[image_info$sampleID == i, ]$full_width
      max_imheight <- image_info[image_info$sampleID == i, ]$full_height
      spatial_coords[[i]] |> mutate(pxl_col_in_fullres = pxl_col_in_fullres - pad[[i]][1]) |> 
        mutate(pxl_row_in_fullres = pxl_row_in_fullres - pad[[i]][3]) |> 
        filter(between(x = pxl_col_in_fullres, left = 0, right = max_imwidth - (pad[[i]][1] + pad[[i]][2]))) |> 
        filter(between(x = pxl_row_in_fullres, left = 0, right = max_imheight - (pad[[i]][3] + pad[[i]][4])))
    })
    spatial_coords <- do.call(bind_rows, spatial_coords)
    .categorical_data <- .categorical_data[spatial_coords$barcode, ]
  } else {
    spatial_coords <- GetStaffli(object)@meta_data
  }
  
  # This will make sure the temporary directory is removed when the app is closed
  clean_after_close <- params$clean_after_close

  # Get sample IDs barcodes
  barcodes <- spatial_coords |>
    select(barcode, sampleID)
  # Split barcodes by sampleIDs
  barcodes <- split(barcodes$barcode, barcodes$sampleID)

  # Import beakr
  if (!requireNamespace("beakr"))
    install.packages("beakr")
  # Import colourpicker
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
  if (verbose) {
    cli_alert_info("Hosting file server at http://{host}:{port}")
  }
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
        width = "90px",
        top = "50px",
        left = paste0(container_width - 90, "px"),
        draggable = TRUE
      )
    ),
    conditionalPanel(
      condition = "!output.panelStatus",
      uiOutput("catlegend")
    ),
    # Create the help menu BS modal
    .helpMenu(),

    # UI for floating sidebar
    absolutePanel(
      shinyBS::bsButton("lasso", "", icon = icon("draw-polygon", verify_fa = FALSE)),
      shinyBS::bsTooltip(id = "lasso", title = "lasso selection", placement = "right", trigger = "hover"),
      shinyBS::bsButton("help", "", icon = icon("question", verify_fa = FALSE)), # Opens up a help menu
      shinyBS::bsButton("quit", "", icon = icon("floppy-disk", verify_fa = FALSE)), # Quits application and save changes
      shinyBS::bsTooltip(id = "quit", title = "save & quit", placement = "right", trigger = "hover"),
      fluidRow(),
      selectizeInput("sample", "Sample", selected = sampleIDs[1], choices = sampleIDs),
      fastSliderInput("opacity", "Opacity", min = 0, max = 1, value = 1, step = 0.05), # Slider for opacity values
      selectizeInput("feature", label = "Feature", selected = selected_features[1], choices = NULL), # numeric features

      # This panel only opens up if a numeric feature is selected
      conditionalPanel(
        condition = "output.panelStatus",
        shinyBS::bsButton("scalealpha", "", icon = icon("signal", verify_fa = FALSE), type = "toggle"),
        shinyBS::bsTooltip(id = "scalealpha", title = "Add transparency", placement = "right", trigger = "hover"),
        shinyBS::bsButton("centerzero", "", icon = icon("arrows-up-down", verify_fa = FALSE), type = "toggle"),
        shinyBS::bsTooltip(id = "centerzero", title = "Center data", placement = "right", trigger = "hover"),
        shinyBS::bsButton("revpal", "", icon = icon("backward", verify_fa = FALSE), type = "toggle"),
        shinyBS::bsTooltip(id = "revpal", title = "Invert colors", placement = "right", trigger = "hover"),
        sliderInput("trim", "Trim", min = 0, max = 1, value = c(0, 1), step = 0.01),
        selectizeInput("colorscale", "Colorscale", choices = names(.color_palettes), options = list(create = TRUE))
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
                         category = 1, sample = sampleIDs[1], lasso = FALSE, allClear = FALSE)

    # panelStatus is used to tell whether a numeric or categorical feature has been selected
    # If a categorical value is selected, the Label input field and the color picker will appear
    output$panelStatus <- reactive({rv$lastBtn == "feature"})
    outputOptions(output, "panelStatus", suspendWhenHidden = FALSE)

    # Switch sample ID when selecting a section ID from the Sample drop-down list
    observeEvent(input$sample, {
      rv$curbarcodes <- barcodes[[input$sample]]
      # Update reactive value rv$allClear to make sure that the
      # color legends aren't rendered until the feature values
      # have been updated
      rv$allClear <- FALSE
    })

    # listen to changes in input$category and reactive value rv$category
    toListenCategory <- reactive({
      list(input$category, rv$category, rv$curbarcodes, rv$levels)
    })

    # # Set last button to "category" when a categorical feature is selected
    observeEvent(toListenCategory(), {
      rv$allClear <- FALSE
      if (input$category > 0) { # Trigger only when input$category exists
        rv$lastBtn = "category"
        rv$curFeature = input$category
        # If the selected category doesn't exists, create a new one
        if (!input$category %in% colnames(.categorical_data)) {
          .categorical_data[, input$category] <<- "NA"
          # Update color dictionary for new category
          if (!rv$curFeature %in% names(.feature_viewer_colors)) {
            .feature_viewer_colors[[rv$curFeature]] <<- setNames(.set_colors(n = 1), nm = "NA")
          }
        }
        # Select categorical data and save to reactive variable rv$values
        rv$values = .categorical_data[rv$curbarcodes, input$category]
        # Generate opacity values for selected category data
        rv$opacities = rep(1, length(rv$values))
        # Create levels for character vectors or select levels for factors
        if (inherits(rv$values, what =  "character")) {
          rv$levels <- unique(.categorical_data |> pull(all_of(input$category)))
        } else {
          rv$levels <- levels(rv$values)
        }
        # Define colors based on color dictionary
        rv$colors <- .feature_viewer_colors[[rv$curFeature]][rv$levels]
        # Switch isNumeric to FALSE to tell app that the data should be treated as categorical
        rv$isNumeric = FALSE
        rv$allClear <- TRUE
      }
    })

    # Make the following inputs reactive: centerzero, feature and trim
    # Any changes to these inputs will rerender the app
    toListenFeatureInput <- reactive({
      list(input$centerzero, input$feature, input$trim, rv$curbarcodes, input$revpal)
    })

    # Listen to changes in input values: centerzero, feature and trim
    observeEvent(toListenFeatureInput(), {
      rv$allClear <- FALSE
      if (!is.null(input$feature)) { # Trigger only when input$feature exists
        if (input$feature %in% all_features) {
          rv$lastBtn = "feature"
          rv$curFeature = input$feature
          # fetch numeric data with FetchData and pull out the vector
          rv$values = FetchData(object, cells = rv$curbarcodes, vars = input$feature) |> pull(all_of(input$feature))
          # Redefine the range of the color bar when centerzero is active
          if (input$centerzero) {
            maxAbsVal <- max(abs(rv$values))
            rv$range = c(-maxAbsVal, maxAbsVal)# |> round(digits = 2)
          } else {
            # TODO: handle data with low values
            rv$range = range(rv$values)# |> round(digits = 2)
          }
          # Trim data if the trim sliders have been changed or keep as is
          # if the trim sliders are default
          if (all(c(0, 1) == input$trim)) {
            trimmed_range <- rv$range
            tmp_vals <- rv$values
          } else {
            trimmed_range <- .trim_range(rv$range, minCutoff = input$trim[1], input$trim[2])
            tmp_vals <- rv$values
            tmp_vals[tmp_vals < trimmed_range[1]] <- trimmed_range[1]
            tmp_vals[tmp_vals > trimmed_range[2]] <- trimmed_range[2]
          }
          # Define opacity values based on trimmed values
          rv$opacities = scales::rescale(tmp_vals, to = c(0, 1), from = trimmed_range)
          # Select input color scale
          rv$colors = .color_palettes[[input$colorscale]]
          if (input$revpal) {
            rv$colors <- rev(rv$colors)
          }
          # Switch isNumeric to TRUE to tell app that the data should be treated as numeric
          # This will trigger showing the scale alpha, center zero, trim and colorscale UI elements
          rv$isNumeric <- TRUE
          rv$allClear <- TRUE
        }
      }
    })

    # listen to changes in colorscale and update
    observeEvent(input$colorscale, {
      if (rv$lastBtn == "feature") {
        rv$colors = .color_palettes[[input$colorscale]]
      }
    })

    # Server side updating of available all_features
    # This is a faster alternative to selectizeInput, where the drop-down list is rendered server side
    updateSelectizeInput(session, 'feature', choices = all_features, selected = all_features[1], server = TRUE)

    # Fire event every 50ms to listen for changes in lasso selection
    autoCheck <- reactiveTimer(50)

    # Observe timer and return latest transformations
    # this is used to communicate transformations from
    # the app back to R
    observe({
      autoCheck()
      shinyjs::js$getSelection()
    })

    # Observe ENTER press when a label has been added
    # If a new label is added and the user presses ENTER, the label will be saved
    # and the status of the lasso selection tool will switch
    observeEvent(input$label_save, {
      if (input$label != "") {
        #updateCheckboxInput(session, "lasso", value = FALSE)
        rv$lasso <- !rv$lasso
        rv$label <- input$label
        # Update style of lasso icon to deselected
        shinyBS::updateButton(session, "lasso", label = "", block = FALSE, style = "default")
      }
    })

    # Update UI buttons for scale alpha, center zero and invert colors
    observeEvent(input$scalealpha, {
      if (input$scalealpha) {
        shinyBS::updateButton(session, "scalealpha", label = "", block = FALSE, style = "success")
      } else {
        shinyBS::updateButton(session, "scalealpha", label = "", block = FALSE, style = "default")
      }
    })
    observeEvent(input$centerzero, {
      if (input$centerzero) {
        shinyBS::updateButton(session, "centerzero", label = "", block = FALSE, style = "success")
      } else {
        shinyBS::updateButton(session, "centerzero", label = "", block = FALSE, style = "default")
      }
    })
    observeEvent(input$revpal, {
      if (input$revpal) {
        shinyBS::updateButton(session, "revpal", label = "", block = FALSE, style = "success")
      } else {
        shinyBS::updateButton(session, "revpal", label = "", block = FALSE, style = "default")
      }
    })

    # Trigger category change on ENTER
    observeEvent(input$category_save, {
      rv$category <- rv$category + 1
    })

    # Listen to change in lasso
    observeEvent(input$lasso, {
      rv$lasso <- !rv$lasso
      if (rv$lasso) {
        # Update style of lasso icon when it is selected
        shinyBS::updateButton(session, "lasso", label = "", block = FALSE, style = "success")
      } else {
        # Update style of lasso icon when it is deselected
        shinyBS::updateButton(session, "lasso", label = "", block = FALSE, style = "default")
      }
    })

    # Listen for changes in selbarcodes and add selection if
    # a label is submitted
    observeEvent(input$selbarcodes, {
      # Only trigger if at least 1 barcode has been selected and if the label name is not ""
      if ((length(input$selbarcodes) > 0) & (input$label != "")) {
        # Handle factor data
        if (inherits(.categorical_data[, input$category], what = "factor")) {
          tmp <- .categorical_data |>
            select(contains(input$category)) |>
            mutate_if(is.factor, as.character)
          tmp[input$selbarcodes, input$category] <- input$label
          tmp <- tmp |>
            mutate(across(everything(), ~factor(.x)))
          .categorical_data[, input$category] <<- tmp[, input$category, drop = TRUE]
        } else {
          .categorical_data[input$selbarcodes, input$category] <<- input$label
        }
        # Add new color if the label doesn't exist yet
        if (!rv$label %in% rv$levels) {
          tmp_cols <- c(.feature_viewer_colors[[input$category]][rv$levels], setNames(input$color, nm = rv$label))
          rv$levels <- c(rv$levels, rv$label)
          rv$colors <- tmp_cols
          .feature_viewer_colors[[input$category]] <<- tmp_cols
        } else {
          tmp_cols <- .feature_viewer_colors[[input$category]]
          tmp_cols <- c(tmp_cols[setdiff(names(tmp_cols), rv$label)], setNames(input$color, nm = rv$label))
          .feature_viewer_colors[[input$category]] <<- tmp_cols
          rv$colors <- tmp_cols[rv$levels]
        }
        rv$values <- .categorical_data[rv$curbarcodes, input$category]
        updateTextInput(session, "label", value = "")
      }
    })


    # Send data to ftrviewer widget
    output$ftrviewerWidget <- renderFtrviewer({
      if (length(rv$isNumeric) > 0) {
        ftrviewer(host = host,
                  port = paste0(port),
                  values = rv$values,
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

    # Render color legend for numerical features as a separate, floating UI element
    output$legend <- renderPlot({
      if (rv$allClear & rv$isNumeric) {
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

    # Render color legend for categorical features as a separate, floating UI element
    output$legend_cat <- renderPlot({
      if (rv$allClear & !rv$isNumeric) {
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
  
  # Fetch new categories if available
  initial_categories <- colnames(object[[]])
  new_categories <- setdiff(colnames(.categorical_data), initial_categories)

  # Create a copy of the original meta data object
  mdata_copy <- object[[]] |> 
    mutate_if(is.factor, as.character)
  if (length(new_categories) > 0) {
    # Create empty columns for new data
    mdata_copy[new_categories] <- NA_character_
  }

  # Fetch factor levels
  factor_levels <- lapply(.categorical_data, function(x) {
    if (inherits(x, what = "factor")) {
      levels(x)
    }
  })
  factor_levels <- factor_levels[!sapply(factor_levels, is.null)]
  
  # Convert all factors to character vectors to enable conversion 
  .categorical_data <- .categorical_data |> 
    mutate_if(is.factor, as.character)
  
  # Convert all "NA" to NA
  .categorical_data[.categorical_data == "NA"] <- NA_character_
  
  # Place modified data in meta data copy
  for (ctgry in colnames(.categorical_data)) {
    mdata_copy[rownames(.categorical_data), ctgry] <- .categorical_data[, ctgry, drop = TRUE] |> as.character()
  }
  
  # Restore factor levels
  if (length(factor_levels) > 0) {
    for (ctgry in names(factor_levels)) {
      mdata_copy[, ctgry] <- factor(mdata_copy[, ctgry, drop = TRUE], levels = factor_levels[[ctgry]])
    }
  }
  
  # Place new meta data in Seurat object
  object@meta.data <- mdata_copy

  # Remove temporary directory
  if (clean_after_close) {
    if (verbose) cli_alert_info("Removing temporary directory")
    unlink(x = datapath, recursive = TRUE)
  }

  return(object)
}


#' Add a help menu to FeatureViewer app
#'
#' @noRd
.helpMenu <- function() {
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
                          column(8, p("Click on ", icon("draw-polygon", verify_fa = FALSE),
                                      " and select spots with the cursor. ",
                                      "Hold ", code("SHIFT"), "to deselect spots.",
                                      "To save the selection, make sure that a ", code("Category"),
                                      " is selected (or add a new one), write a new ", code("Label"),
                                      " name and press ", code("ENTER"), ". ",
                                      "Use the color picker tool (", code("Select color"), ") to select a color for the ",
                                      code("Label"), " before pressing ", code("ENTER"), "."))),
                   column(12, column(4, p(strong("Add transparency"))),
                          column(8, p("Click on the ", icon("signal", verify_fa = FALSE), " icon to add opacity to spots ",
                                      "proportional to the numeric feature values."))),
                   column(12, column(4, p(strong("Center numeric data"))),
                          column(8, p("Click on the ", icon("arrows-up-down", verify_fa = FALSE), " icon to center the color scale at 0. ",
                                      "If this option is used, use a divergent color palette from ",
                                      code("Colorscale"), " such as 'RdBu'."))),
                   column(12, column(4, p(strong("Invert colors"))),
                          column(8, p("Click on the ", icon("backward", verify_fa = FALSE), " icon to flip the color palette ",
                                      " selected from ", code("Colorscale")))),
                   column(12, column(4, p(strong("Trim values"))),
                          column(8, p("Adjust the  ", code("Trim"), " sliders to trim the lower/upper bounds of numeric feature values "))),
                   column(12, column(4, p(strong("Save changes and quit"))),
                          column(8, p("Click on ", icon("floppy-disk", verify_fa = FALSE),
                                      " to save changes and quit the application. If any new selections have been added, these will be ",
                                      " stored in the 'meta.data' slot of the returned ", code("Seurat"), " object. ",
                                      " The results will NOT be returned if the app is closed from the R session."))))
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

#' Check if colors are valid
#'
#' @importFrom grDevices col2rgb
#'
#' @noRd
.areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}


#' Check custom color palettes
#'
#' @noRd
.check_custom_palettes <- function (
  custom_color_palettes
) {
  stopifnot(inherits(custom_color_palettes, what = "list"),
            length(custom_color_palettes) > 0,
            inherits(names(custom_color_palettes), what = "character"))
  if (!all(sapply(custom_color_palettes, class) == "character")) {
    abort("Invalid custom_color_palettes. Expected a named list of character vectors.")
  }
  for (cols in custom_color_palettes) {
    check_cols <- .areColors(cols)
    if (!all(check_cols)) {
      abort(glue("Invalid color(s) {paste(cols[!check_cols], collapse = ', ')}"))
    }
  }
  if (any(names(custom_color_palettes) %in% .color_scales(info = TRUE))) {
    abort(glue("Color palette {intersect(names(custom_color_palettes), .color_scales(info = TRUE))} already exists"))
  }
}


#' Check categoircal colors
#'
#' @import rlang
#' @import glue
#'
#' @noRd
.check_categorical_colors <- function (
  categorical_colors,
  categorical_features,
  .categorical_data
) {
  stopifnot(inherits(categorical_colors, what = "list"),
            length(categorical_colors) > 0,
            inherits(names(categorical_colors), what = "character"))
  if (!all(sapply(categorical_colors, class) == "character")) {
    abort("Invalid categorical_colors. Expected a named list of named character vectors.")
  }
  for (cols in categorical_colors) {
    if (!inherits(names(cols), what = "character")) {
      abort("Invalid categorical_colors. Expected named character vectors in every element of the list.")
    }
    check_cols <- .areColors(cols)
    if (!all(check_cols)) {
      abort(glue("Invalid color(s) {paste(cols[!check_cols], collapse = ', ')}"))
    }
  }
  if (!all(categorical_features %in% names(categorical_colors))) {
    abort("Invalid categorical_colors. Names do not match with available categorical features.")
  }
  for (ftr_name in categorical_features) {
    x <- .categorical_data[, ftr_name, drop = TRUE]
    if (inherits(x, what = "factor")) {
      x_unique <- levels(x)
      if (!all(x_unique %in% names(categorical_colors[[ftr_name]]))) {
        abort(glue("Invalid categorical_colors. Missing colors for categorical feature {ftr_name}."))
      }
    } else {
      x_unique <- unique(x)
      if (!all(x_unique %in% names(categorical_colors[[ftr_name]]))) {
        abort(glue("Invalid categorical_colors. Missing colors for categorical feature {ftr_name}."))
      }
    }
    categorical_colors[[ftr_name]] <- categorical_colors[[ftr_name]][x_unique]
  }
  return(categorical_colors)
}


#' Check tile paths
#'
#' @param object A \code{Seurat} object
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
      tilepath <- paste0(datadir, paste0("/tiles", i))
      if (!file.exists(tilepath)) {
        abort(glue("Path {tilepath} is missing for sample {i}"))
      }
      infopath <- paste0(datadir, paste0("/image_info_", i, ".json"))
      if (!file.exists(infopath)) {
        abort(glue("Path {infopath} is missing for sample {i}"))
      }
      coordpath <- paste0(datadir,  paste0("/coords_Visium_", i, ".json"))
      if (!file.exists(coordpath)) {
        abort(glue("Path {coordpath} is missing for sample {i}"))
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
      dirs <- TileImage(im = image_read(imgs[i]), sampleID = sampleIDs[i], overwrite = TRUE, verbose = verbose)
      datapath <- dirs$datapath
      if (verbose) cli_alert("  Exporting Visium coordinates for sample {sampleIDs[i]}")
      export_coordinates(object = object, sampleNumber = sampleIDs[i], outdir = datapath, overwrite = TRUE, verbose = verbose)
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
        .gg_circle(r = 0.4, xc = 0.5, yc = df[lvl, "yc", drop = TRUE], fill = colors[lvl])
    }
  } else {
  # Try to shorten title
  if (nchar(title) > 10) {
    title <- gsub(pattern = "\\.|\ |_", replacement = "\n", x = title)
    if (nchar(title) > 20) {
      title <- paste0(substr(title, start = 1, stop = 17), "...")
    }
    title <- strwrap(title, width = 12) |> paste(collapse = "\n")
  }
  df <- data.frame(y = seq(minVal, maxVal, length.out = ntiles))
  p <- ggplot(df, aes("", y, fill = y)) +
    geom_raster() +
    theme(legend.position = "none") +
    scale_y_continuous(position = "right", labels = .format_numeric_label_text) +
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


#' Format label numbers to have an exact width
#'
#' @noRd
.format_numeric_label_text <- function(x) {
  x <- format(x = x, digits = 3)
  pad <- 8 - nchar(x)
  sapply(seq_along(x), function(i) {
    paste0(x[i], paste(rep(" ", pad[i]), collapse = ""))
  })
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
.gg_circle <- function (
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
  if (requireNamespace("scico", quietly = TRUE)) {
    scico_paletets <- lapply(scico::scico_palette_names(), function(nm) {
      scico::scico(n = 11, palette = nm)
    }) |> setNames(scico::scico_palette_names())
    colscales <- c(colscales, scico_paletets)
  }
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
