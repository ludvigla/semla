#' @include checks.R
#'
NULL

#' Plot H&E images
#'
#' If images are loaded into the \code{Seurat} object with \code{LoadImages},
#' this function can be used to quickly plot the images in a grid. If
#' you have applied transformations to the images, e.g. with
#' \code{RigidTransformImages}, you can specify \code{image_use = "transformed"}
#' to plot the transformed images instead.
#'
#' @param object A Seurat object
#' @param label_by A string specifying a meta data column to label plots by.
#' This needs to be a \code{character} or \code{factor} and multiple labels for each section
#' are not allowed.
#' @param image_use String specifying image type to use, either 'raw' or
#' 'transformed'
#' @param crop_area A numeric vector of length 4 specifying a rectangular area to crop
#' the plots by. These numbers should be within 0-1. The x-axis is goes from left=0 to
#' right=1 and the y axis is goes from top=0 to bottom=1. The order of the values are
#' specified as follows: \code{crop_area = c(left, top, right, bottom)}. The crop area
#' will be used on all tissue sections and cannot be set for each section individually.
#' using crop areas of different sizes on different sections can lead to unwanted side
#' effects as the point sizes will remain constant. In this case it is better to generate
#' separate plots for different tissue sections.
#' @param sampleIDs An integer vector with section numbers to plot
#' @param ncol An integer value specifying the number of columns in the plot
#' grid
#' @param mar Margins around each plot. See \code{\link{par}} for details
#' @param return_as_gg Should the plot be returned as a \code{ggplot} object?
#'
#' @importFrom rlang abort %||%
#' @importFrom graphics layout
#' @importFrom dplyr select bind_cols
#' @importFrom graphics par title
#'
#' @family spatial-visualization-methods
#'
#' @author Ludvig Larsson
#'
#' @return draws a plot of the H&E images
#'
#' @examples
#'
#' library(semla)
#'
#' # Load example Visium data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
#' se_merged <- MergeSTData(se_mbrain, se_mcolon)
#'
#' # ImagePlot will throw an error if no images are loaded
#' \dontrun{
#' ImagePlot(se_merged)
#' }
#'
#' # Load images
#' se_merged <- LoadImages(se_merged)
#' ImagePlot(se_merged)
#'
#' # Plot only selected tissue sections
#' ImagePlot(se_merged, sampleIDs = 1)
#'
#' # Change order of plot
#' ImagePlot(se_merged, sampleIDs = 2:1)
#'
#' # Add a sample_id column and label plots
#' se_merged$sample_id <- ifelse(GetStaffli(se_merged)@meta_data$sampleID == 1, "brain", "colon")
#' ImagePlot(se_merged, label_by = "sample_id")
#'
#' # Reload images in higher resolution, crop image and remove margins
#' se_merged <- LoadImages(se_merged, image_height = 1e3)
#' se_merged <- LoadImages(se_merged, image_height = 1.5e3)
#' ImagePlot(se_merged, crop_area = c(0.4, 0.4, 0.7, 0.7), sampleIDs = 1, mar = c(0, 0, 0, 0))
#' ImagePlot(se_merged, crop_area = c(0.45, 0.55, 0.65, 0.7), sampleIDs = 2, mar = c(0, 0, 0, 0))
#'
#' @export
#'
ImagePlot <- function (
    object,
    label_by = NULL,
    image_use = c("raw", "transformed"),
    crop_area = NULL,
    sampleIDs = NULL,
    ncol = NULL,
    mar = c(1, 1, 1, 1),
    return_as_gg = FALSE
) {

  # Run checks
  .check_seurat_object(object)
  .check_seurat_images(object)
  image_use <- match.arg(image_use, choices = c("raw", "transformed"))

  # obtain Staffli object
  st_object <- GetStaffli(object)

  # Check label_by
  if (!is.null(label_by)) {
    if (!is.character(label_by)) abort("'label_by' should be a string specifying a meta data column in the Seurat object.")
    if (!label_by %in% colnames(object[[]])) abort(glue("'{label_by}' is not present in the Seurat object meta data."))
    label_by_vec <- object[[]] |>
      pull(all_of(label_by))
    if (!class(label_by_vec) %in% c("character", "factor")) abort(glue("Invalid class '{class(label_by_vec)}' for 'label_by' column. ",
                                                                       "Expects a 'character' of 'factor'."))
    label_by_vec  <- sapply(split(label_by_vec, st_object@meta_data$sampleID), function(x) {
      y <- unique(x)
      if (length(y) > 1) abort("Invalid 'label_by' meta data column. ",
                               "Cannot have multiple labels per tissue section.")
      return(y)
    })
  }

  # fetch images
  if (!image_use %in% names(st_object@rasterlists))
    abort(glue("Transformed images are not available in this object."))
  images <- st_object@rasterlists[[image_use]]

  # Use all images if indices are not specified
  sampleIDs <- sampleIDs %||% {
    seq_along(images)
  }

  # Check if indices are OK
  if (any(!sampleIDs %in% seq_along(images))) abort(glue("'sampleIDs' is out of range. Section numbers ",
                                                               "{paste(seq_along(images), collapse = ', ')} are available."))

  # Select sections
  images <- images[sampleIDs]
  if (!is.null(label_by)) {
    label_by_vec <- label_by_vec[sampleIDs]
  }

  # Validate crop_area
  if (!is.null(crop_area)) {
    if (!is.numeric(crop_area)) abort(glue("Invalid class '{class(crop_area)}' for 'crop_area', expected 'numeric'"))
    if (length(crop_area) != 4) abort(glue("Invalid length for 'crop_area', expected a 'numeric' vector of length 4"))
    if (!all(between(x = crop_area, left = 0, right = 1))) abort("'crop_area' can only take values between 0-1")
    images <- lapply(images, function(im) {
      x_start <- round(crop_area[1]*ncol(im))
      x_end <- round(crop_area[3]*ncol(im))
      y_start <- round(crop_area[2]*nrow(im))
      y_end <- round(crop_area[4]*nrow(im))
      im <- im[y_start:y_end, x_start:x_end]
      return(im)
    })
  }

  # Define dimensions of plot grid
  ncols <- ncol %||% ceiling(sqrt(length(x = images)))
  nrows <- ceiling(length(x = images)/ncols)

  if (!return_as_gg) {
    # Create a plot layout
    layout.matrix <- t(matrix(c(1:length(images), rep(0, nrows*ncols - length(images))), nrow = ncols, ncol = nrows))
    layout(mat = layout.matrix)
    
    # plot images
    for (i in seq_along(images)) {
      rst <- images[[i]]
      par(mar = mar)
      plot(rst)
      if (!is.null(label_by)) {
        title(main = label_by_vec[i])
      }
    }
  } else {
    plots <- lapply(seq_along(images), function(i) {
      p <- ggplot()
      if (!is.null(label_by)) {
        p <- p + ggtitle(label = label_by_vec[i])
      }
      p <- p +
        theme(plot.margin = margin(t = mar[1], r = mar[2], b = mar[3], l = mar[4])) +
        inset_element(p = images[[i]], left = 0, bottom = 0, right = 1, top = 1)
    })
    p <- wrap_plots(plots, ncol = ncols)
    return(p)
  }
}
