#' @include checks.R
#'
NULL

#' @param image_height an integer value specifying the height of the down-scaled images
#' @param verbose print messages
#'
#' @importFrom rlang abort inform
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

  if (verbose) cli::cli_h2("Load H&E images")

  # Check files
  if (!is.character(object)) abort("Invalid class '{class(object)}', expected a 'character' vector.")
  for (f in object) {
    if (!file_ext(f) %in% c("png", "jpg", "jpeg")) abort(glue("Only PNG and JPEG images are supported, got file extension .{file_ext(f)}"))
    if (!file.exists(f)) abort(glue("File {f} doesn't exist."))
  }

  # Load images
  raw_rasters <- lapply(object, function(f) {
    if (verbose) inform(c("i" = glue("Loading image from {f}")))
    im <- f |>
      image_read()
    info <- im |>
      image_info()
    rst <- im |>
      image_scale(geometry = geometry_area(height = image_height)) |>
      as.raster()

    if (verbose) inform(c("v" = glue("Scaled image from {info$height}x{info$width} to {image_height}x{ncol(rst)} pixels")))
    return(rst)
  })

  if (verbose) inform(c("i" = glue("Saving loaded H&E images as 'rasters' in Seurat object")))
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

  # Place Staffli object in Seurat object
  object@tools$Staffli <- st_object

  return(object)
}
