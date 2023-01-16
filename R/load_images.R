#' @include checks.R
#' @include generics.R
#'
NULL

#' @param image_height an integer value specifying the height of the down-scaled images
#' @param verbose print messages
#'
#' @importFrom rlang abort
#' @import cli
#' @importFrom glue glue
#' @importFrom magick image_read image_scale geometry_area image_info
#' @importFrom tools file_ext
#' @importFrom grDevices as.raster
#'
#' @section Seurat:
#' If a Seurat object is provided, the images will be loaded as `raster` objects
#' and stored inside the \code{Staffli} object that is located in the \code{tools}
#' slot.
#'
#' @section default method:
#' If a character vector of image paths are provided, the images will be loaded,
#' then down-scaled based on \code{image_height} and returned as a list of `raster`
#' objects. Only JPEG and PNG images are supported.
#'
#' @rdname load-images
#'
#' @export
#'
LoadImages.default <- function (
  object,
  image_height = 400,
  verbose = TRUE,
  ...
) {

  if (verbose) cli_h2("Load H&E images")

  # Check input
  if (!is.character(object)) abort("Invalid class '{class(object)}', expected a 'character' vector.")

  # Include images from example data if object is either "mousebrain" or "mousecolon"
  for (i in seq_along(object)) {
    if (object[i] %in% c("mousebrain", "mousecolon")) {
      object[i] <- .load_ext_images(object[i])
    }
  }

  # Check files
  for (f in object) {
    if (!file_ext(f) %in% c("png", "jpg", "jpeg"))
      abort(glue("Only PNG and JPEG images are supported, got file extension .{file_ext(f)}"))
    if (!file.exists(f)) abort(glue("File {f} doesn't exist."))
  }

  # Load images
  raw_rasters <- lapply(object, function(f) {
    if (verbose) cli_alert_info("Loading image from {f}")
    im <- f |>
      image_read()
    info <- im |>
      image_info()
    if ((!inherits(image_height, what = "numeric")) | (length(image_height) != 1))
      abort("image_height should be a numeric of length 1")
    if (info$height < image_height) {
      abort(glue("image_height has to be smaller than or equal to {info$height}px"))
    }
    if (info$height != image_height) {
      im <- im |>
        image_scale(geometry = geometry_area(height = image_height))
      if (verbose) cli_alert_info("Scaled image from {info$height}x{info$width} to {image_height}x{ncol(rst)} pixels")
    } else {
      if (verbose) cli_alert_info("Loaded H&E image in full resolution")
    }
    rst <- im |> as.raster()
    return(rst)
  })

  if (verbose) cli_alert_info("Saving loaded H&E images as 'rasters' in Seurat object")
  return(raw_rasters)
}

#' @importFrom rlang %||%
#'
#' @rdname load-images
#'
#' @examples
#'
#' #' library(STUtility2)
#'
#' # Load example Visium data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "STUtility2"))
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "STUtility2"))
#' se_merged <- MergeSTData(se_mbrain, se_mcolon)
#'
#' # Load images
#' se_merged <- LoadImages(se_merged)
#'
#' @export
#'
LoadImages.Seurat <- function (
    object,
    image_height = 400,
    verbose = TRUE,
    ...
) {

  # validate Seurat object
  .check_seurat_object(object)

  # Get the Staffli object
  st_object <- GetStaffli(object)

  # Load images
  raw_rasters <- LoadImages(st_object@imgs, image_height = image_height, verbose = verbose)

  # Save loaded images in object
  st_object@rasterlists$raw <- raw_rasters
  st_object@image_height <- image_height

  # Set new images if running examples
  for (i in seq_along(st_object@imgs)) {
    if (st_object@imgs[i] %in% c("mousebrain", "mousecolon")) {
      st_object@imgs[i] <- .load_ext_images(st_object@imgs[i])
    }
  }

  # Place Staffli object in Seurat object
  object@tools$Staffli <- st_object

  return(object)
}


#' Utility function for setting images in package Seurat objects
#'
#' @param dataset One of "mousebrain" or "mousecolon"
#'
#' @return Path to "tissue_hires_image.jpg"
#'
#' @noRd
.load_ext_images <- function (
  dataset = "mousebrain"
) {
  hiresresimagefile <- system.file(paste0("extdata/", dataset, "/spatial"),
                                   "tissue_hires_image.jpg",
                                   package = "STUtility2")
  return(hiresresimagefile)
}
