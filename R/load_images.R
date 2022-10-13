#' @include checks.R
#'
NULL

#' @importFrom rlang abort
#' @importFrom glue glue
#' @importFrom magick image_read image_scale geometry_area
#' @importFrom tools file_ext
#'
#' @rdname load-images
#'
#' @export
#'
LoadImages.default <- function (
  image_paths,
  image_height = 400
) {

  # Check files
  if (!is.character(image_paths)) abort("Invalid class '{class(image_paths)}', expected a 'character' vector.")
  for (f in image_paths) {
    if (!file_ext(f) %in% c("png", "jpg", "jpeg")) abort(glue("Only PNG and JPEG images are supported, got file extension .{file_ext(f)}"))
    if (!file.exists(f)) abort(glue("File {f} doesn't exist."))
  }

  # Load images
  raw_rasters <- lapply(image_paths, function(f) {
    f |>
      image_read() |>
      image_scale(geometry = geometry_area(height = image_height)) |>
      as.raster()
  })

  return(raw_rasters)
}


#' @param object An object
#'
#' @rdname load-images
#'
#' @export
#'
LoadImages.Seurat <- function (
    object,
    image_height = 400
) {

  # validate Seurat object
  .check_seurat_object(object)

  # Get the Staffli object
  st_object <- GetStaffli(object)

  # Load images
  raw_rasters <- LoadImages(st_object@imgs, image_height = image_height)

  # Save loaded images in object
  st_object@rasterlists$raw <- raw_rasters
  st_object@image_height <- image_height

  # Place Staffli object in Seurat object
  object@tools$Staffli <- st_object

  return(object)
}
