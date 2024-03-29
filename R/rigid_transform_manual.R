#' @include generics.R
#' @include checks.R
#'
NULL


#' @param container_width,container_height Width and height of paper widget given in pixels.
#' @param launch_browser Should the app be opened in the default browser? Default is TRUE.
#' Set \code{launch_browser = getOption("shiny.launch.browser", interactive())} to use the RStudio 
#' built-in browser.
#'
#' @section default method:
#' Takes a list of images prepared with the \code{\link{prep_image}} function
#' and send them to the interactive application. When the application stops
#' (after clicking "quit & save") a tibble is returned with information about
#' the rigid transformations applied to the images.
#' 
#' @section Seurat method:
#' Takes a \code{Seurat} object with at least 2 tissue sections and opens an interactive
#' shiny application where rigid transformations can be applied to the images.
#' When the application stops (after clicking "quit & save"), the transformations 
#' are applied to the images using \code{\link{RigidTransformImages}} and a 
#' \code{Seurat} object is returned. If no transformations are supplied, the function
#' will return the input \code{Seurat} object unmodified.
#'
#' @importFrom shiny fluidPage actionButton tableOutput uiOutput reactiveTimer
#' observe fluidRow column p h4 helpText strong code h5 observeEvent stopApp
#' runApp renderTable renderUI hr
#' @importFrom shinyjs useShinyjs extendShinyjs
#' @importFrom tibble as_tibble
#'
#' @rdname manual-transform-images
#' @family transforms
#'
#' @author Ludvig Larsson
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
#' # Align images and find transform
#' \dontrun{
#' if (interactive()) {
#'   transforms <- RunAlignment(object = list(img1, img2))
#' }
#' }
#'
#' @export
RunAlignment.default <- function (
    object,
    container_width = '800px',
    container_height = '650px',
    launch_browser = TRUE,
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
    paperOutput("paperWidget", height = container_height, width = container_width),
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
      paper(data = object, width = container_width, height = container_height)
    })

    # Fill table when values are passed from react app
    output$table <- renderTable(do.call(rbind, lapply(input$myVal, as.data.frame)))

    # Helpt text
    output$HelpBox = renderUI({
      if (input$help %% 2){
        fluidRow(column(12, helpText(h4("How to use the tool"), hr())),
                 column(12, helpText(h5(strong("NB:"), " The dashed lines highlight the borders of the canvas area defined by the first image in the stack. ",
                                        "Only parts of the first image that are inside this area will be kept. Every image will be transformed independently, ",
                                        "meaning that each image will keep their original dimensions."))),
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

  transformations <- runApp(list(ui = ui, server = server), launch.browser = launch_browser)

  # Add image dimensions to results
  transformations <- transformations |>
    bind_cols(do.call(bind_rows, lapply(object, function(x) {tibble(dimx = x$dimx, dimy = x$dimy)})))

  return(transformations |> as_tibble())
}

#' @param image_height An integer used for rescaling images
#' @param verbose Print messages
#' 
#' @importFrom magick image_read
#' @import rlang
#' @import cli
#' @import glue
#' 
#' @rdname manual-transform-images
#' 
#' @examples 
#' 
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
#' se_merged <- MergeSTData(se_mbrain, se_mcolon) |> 
#'   LoadImages()
#' 
#' # Run alignment application on a Seurat object
#' \dontrun{
#' if (interactive()) {
#'   se_merged <- RunAlignment(se_merged)
#' }
#' }
#' 
#' @export
RunAlignment.Seurat <- function (
  object,
  image_height = 400,
  container_width = '800px',
  container_height = '650px',
  launch_browser = TRUE,
  verbose = TRUE,
  ...
) {

  # Check Seurat object
  .check_seurat_object(object)

  # Check that images are present
  .check_seurat_images(object)

  # Check if images have sufficient resolution
  rstrs <- GetStaffli(object)@rasterlists[["raw"]]
  rstrs <- setNames(rstrs, nm = paste0(seq_along(rstrs)))
  for (nm in names(rstrs)) {
    rst <- rstrs[[nm]]
    if (nrow(rst) < image_height) abort(c(glue("Pre-processing failed for sample {nm}: ",
                                             "image_height ({image_height}) has to be smaller than ",
                                             "the height of the loaded image ({nrow(rst)})."),
                                          glue("Either decrease the value for image_height or rerun ",
                                             "LoadImages() with a higher image_height.")))
  }

  images <- lapply(rstrs, function(rst) {
    im <- prep_image(rst |>
                       image_read(),
                       height = image_height)
  }) |>
    setNames(nm = names(rstrs))

  # Run app
  transforms <- RunAlignment(images, 
                             container_width = container_width, 
                             container_height = container_height, 
                             launch_browser = launch_browser,
                             ...)
  transforms$shift_x <- -transforms$shift_x
  transforms$shift_y <- -transforms$shift_y

  # Apply transformations to images
  transforms_tibble <- do.call(bind_rows, lapply(1:nrow(transforms), function(i) {
    dimx <- transforms[i, "dimx", drop = TRUE]
    dimy <- transforms[i, "dimy", drop = TRUE]
    res <- try({
      generate_rigid_transform(sampleID = i,
                               mirror_y = transforms[i, "flip", drop = TRUE],
                               angle = transforms[i, "angle", drop = TRUE],
                               tr_x = transforms[i, "shift_x", drop = TRUE]/dimx,
                               tr_y = transforms[i, "shift_y", drop = TRUE]/dimy,
                               scalefactor = transforms[i, "scalefactor", drop = TRUE])
    }, silent = TRUE)

    if (inherits(res, what = "try-error")) {
      if (verbose) cli_alert_info(col_br_magenta("No transformation applied to sample {i}"))
      return(NULL)
    } else {
      return(res)
    }
  }))
  
  if (length(transforms_tibble) == 0) {
    cli_alert_warning("Found no transformations. Returning unmodified {col_br_magenta('Seurat')} object")
    return(object)
  }

  object <- RigidTransformImages(object, transforms = transforms_tibble, verbose = verbose)

  # Return object
  return(object)
}


#' Prepare images for paper JS react app
#'
#' @param input An object of class \code{magick-image} or a path
#' to an image in png or jpeg format
#' @param height Height of image sent to react app
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
#' @examples
#' library(semla)
#' library(magick)
#' 
#' im <- system.file("extdata/mousebrain", "spatial/tissue_lowres_image.jpg", package = "semla")
#' 
#' # Prep image
#' im_prepped <- prep_image(im)
#'
#' @export
#'
prep_image <- function (
    input,
    height = 256
) {

  # Set global variables to NULL
  R <- G <- B <- A <- key <- ord <- NULL

  # Check height
  if (between(x = height, left = 256, right = 512))
    if (height < 256) abort(glue("The image height ({height}) has to be between 256 and 512 pixels"))

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

  # Read image, scale it to have a specified height and
  # convert it to an array
  im <- im |>
    image_scale(paste0("x", height)) |>
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
