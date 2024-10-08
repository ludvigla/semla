#' @include checks.R
#' @include generics.R
#' @include staffli_object.R
NULL

#' @param image_height An integer specifying the height of the down-scaled images
#' @param pad_info A tibble with information about how the image should be padded
#' @param verbose print messages
#' 
#' @rdname load-images
#' 
#' @importFrom rlang abort
#' @import cli
#' @importFrom glue glue
#' @importFrom magick image_read image_scale geometry_area image_info
#' @importFrom tools file_ext
#' @importFrom grDevices as.raster
#' 
#' @section default method:
#' If a character vector of image paths are provided, the images will be loaded,
#' then down-scaled based on \code{image_height} and returned as a list of \code{raster}
#' objects. Only JPEG and PNG images are supported.
#' 
#' @section Seurat:
#' If a Seurat object is provided, the images will be loaded as \code{raster} objects
#' and stored inside the \code{Staffli} object that is located in the \code{tools}
#' slot.
#' 
#' @examples 
#' 
#' library(semla)
#' 
#' # Get paths for example images
#' mousebrain_jpg <- system.file("extdata/mousebrain", 
#'                               "spatial/tissue_lowres_image.jpg", 
#'                               package = "semla")
#' mousecolon_jpg <- system.file("extdata/mousecolon", 
#'                               "spatial/tissue_lowres_image.jpg", 
#'                               package = "semla")
#' 
#' rasters <- LoadImages(c(mousebrain_jpg, mousecolon_jpg))
#' 
#' # Save graphical parameters
#' oldpar <- par(no.readonly = TRUE)
#' 
#' # plot images
#' par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))
#' for (rst in rasters) {
#'   plot(rst)
#' }
#' 
#' # Reset graphical parameters
#' par(oldpar)
#' 
#' @export
#' 
LoadImages.default <- function (
  object,
  image_height = 400,
  pad_info = NULL,
  verbose = TRUE,
  ...
) {

  if (verbose) cli_h2("Loading H&E images")

  # Check input
  if (!is.character(object)) abort("Invalid class '{class(object)}', expected a 'character' vector.")

  # Check padding
  if (!is.null(pad_info)) {
    if (nrow(pad_info) != length(object)) abort("pad_info needs to have the same length as object")
    stopifnot(inherits(pad_info, what = "tbl"))
  }

  # Include images from example data if object is either "mousebrain" or "mousecolon"
  for (i in seq_along(object)) {
    if (object[i] %in% c("mousebrain", "mousecolon")) {
      object[i] <- .load_ext_images(object[i])
    }
  }

  # Load images
  raw_rasters <- lapply(seq_along(object), function(i) {
    f <- object[i]
    
    # Check file extension
    if (file.exists(f)) {
      if (!file_ext(f) %in% c("png", "jpg", "jpeg"))
        abort(glue("Only PNG and JPEG images are supported, got file extension .{file_ext(f)}"))
    }
    
    if (verbose) cli_alert_info("Loading image from {f}")
    im <- try({f |> image_read()}, silent = TRUE)
    if (inherits(im, "try-error"))
      abort(glue("Invalid path/url {f}. Make sure to have ",
                 "valid image paths before running {col_br_magenta('LoadImages')}"))
    info <- im |>
      image_info()

    # Check padding
    if (!is.null(pad_info)) {
      im_blank <- image_blank(width = info$width + pad_info[i, "before_x", drop = TRUE] + pad_info[i, "after_x", drop = TRUE],
                  height = info$height + pad_info[i, "before_y", drop = TRUE] + pad_info[i, "after_y", drop = TRUE],
                  color = "white")
      if (any(pad_info[i, c("before_x", "before_y"), drop = TRUE] |> as.numeric() > 0)) {
        im <- image_composite(im_blank, composite_image = im,
                            offset = geometry_area(x_off = pad_info[i, "before_x", drop = TRUE],
                                                   y_off = pad_info[i, "before_y", drop = TRUE]))
      } else {
        im <- image_composite(im_blank, composite_image = im)
      }
    }

    if ((!inherits(image_height, what = "numeric")) | (length(image_height) != 1))
      abort("image_height should be a numeric of length 1")
    if (info$height < image_height) {
      abort(glue("image_height has to be smaller than or equal to {info$height}px"))
    }
    if (info$height != image_height) {
      im <- im |>
        image_scale(geometry = geometry_area(height = image_height))
      im_info_scaled <- image_info(im)
      if (verbose) cli_alert_info("Scaled image from {info$height}x{info$width} to {image_height}x{im_info_scaled$width} pixels")
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
#' @family pre-process
#' @rdname load-images
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

  # Set global variables to NULL
  pad <- sampleID <- width <- height <- full_width <- full_height <- NULL
  before_x <- after_x <- before_y <- after_y <- NULL

  # validate Seurat object
  .check_seurat_object(object)

  # Get the Staffli object
  st_object <- GetStaffli(object)

  # Check if images needs to be padded
  if ("pad" %in% colnames(st_object@image_info)) {
    if (verbose) {
      cli_alert_warning("White space will be added to H&E images to fit all spots on the images")
    }
    pad_info <- st_object@image_info |>
      select(pad, sampleID, width, height, full_width, full_height) |>
      separate(col = "pad", into = c("before_x", "after_x", "before_y", "after_y"), sep = "x") |>
      mutate_if(is.character, as.integer) |>
      mutate(before_x = (before_x/full_width)*width,
             after_x = (after_x/full_width)*width,
             before_y = (before_y/full_height)*height,
             after_y = (after_y/full_height)*height) |>
      mutate(across(before_x:after_y, ~ceiling(.x)))
  } else {
    pad_info <- NULL
  }

  # Load images
  raw_rasters <- LoadImages(st_object@imgs, image_height = image_height, pad_info = pad_info, verbose = verbose)

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


#' Utility function for setting images
#'
#' @param dataset One of "mousebrain" or "mousecolon"
#'
#' @return Path to "tissue_lowres_image.jpg"
#'
#' @noRd
.load_ext_images <- function (
  dataset = "mousebrain"
) {
  hiresresimagefile <- system.file(paste0("extdata/", dataset, "/spatial"),
                                   "tissue_lowres_image.jpg",
                                   package = "semla")
  return(hiresresimagefile)
}
