#' @include generics.R
#'
NULL


#' Utility function to generate a tibble with image transformations
#'
#' @param sampleID an integer specifying a sample ID
#' @inheritParams CoordTransform
#' @inheritParams CoordAndImageTransform
#'
#' @importFrom tibble tibble
#' @importFrom glue glue
#' @importFrom rlang abort
#'
#' @return a tibble with image transformations
#'
#' @export
#'
generate_rigid_transform <- function (
  sampleID = 1,
  mirror_x = FALSE,
  mirror_y = FALSE,
  angle = 0,
  tr_x = 0,
  tr_y = 0
) {

  # Check input parameters
  stopifnot(is.logical(mirror_x) & (length(mirror_x) == 1),
            is.logical(mirror_y) & (length(mirror_y) == 1),
            is.numeric(angle) & (length(angle) == 1),
            is.numeric(tr_x) & (length(tr_x) == 1),
            is.numeric(tr_y) & (length(tr_y) == 1))

  transforms <- tibble(sampleID, mirror_x, mirror_y, angle, tr_x, tr_y)

  # Check transforms
  if (.check_transforms(transforms)) abort(glue("'transforms' cannot take default values. ",
                                          "This would result in no transformation."))

  return(transforms)
}

#' @param image an image of class `magick-image`, `raster`, `StoredSpatialImage`
#' or a path to an image in PNG or JPEG format
#' @param xy_coords spot coordinates that can be mapped to \code{image}
#' @param verbose print messages
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#' @importFrom dplyr select between
#'
#' @rdname transform-images
#'
#' @section default method:
#' Object is a tibble containing information about the image dimensions that are
#' matched with \code{xy_coords}. The coordinates \code{xy_coords} do not have to
#' match the input image as long as the image dimensions in \code{object} are matched.
#' This is to ensure that the coordinates can be transformed correctly regardless of the
#' dimensions of the \code{image}. For example, the coordinates provided by spaceranger
#' in "tissue_positions.csv" can be mapped to the original H&E image and therefore the
#' original H&E image dimensions are required to correctly specify the dimensions
#' when using these coordinates for plotting. Doing so, we can map the coordinates
#' correctly on the H&E image regardless of its size.
#'
#' The required columns are:
#'
#' \itemize{
#'   \item{\code{full_width, full_height}: the dimensions of the image that \code{xy_coords} map to}
#'   \item{\code{sampleID}: an integer specifying a sample ID}
#'   \item{\code{
#'      mirror_x, mirror_y}: TRUE/FALSE specifying if the image should be mirrored along
#'      the x- and/or y-axis
#'      }
#'   \item{\code{angle}: a numeric specifying an angle to rotate the image by in degrees}
#'   \item{
#'      \code{tr_x, tr_y}: numeric values specifying translations along the x- and/or y-axis.
#'      \code{tr_x, tr_y} have to be values between -1 and 1 where 0 means no translation and
#'      1 is equal to the image width or height. Negative values will shift the image along the axis in the
#'      opposite direction. For example, setting \code{tr_x = 0.5} will move the image 50% of the image width
#'      to the right and setting \code{tr_x = -0.5} will move the image 50% of the image width to the left
#'      }
#' }
#'
#' @return A list with two elements:
#' \itemize{
#'   \item{"im_transf": An object of class `magick-image` representing the transformed image}
#'   \item{"xy_transf": An object of class `tble` representing the transformed coordinates}
#' }
#'
#' @examples
#'
#' library(STUtility2)
#' library(tibble)
#'
#' transforms <- generate_rigid_transform(mirror_x = TRUE, angle = 30, tr_x = 0.2, tr_y = -0.2)
#'
#' # Combine image dimensions with transforms.
#' # These are the dimensions of the H&E image used as input for
#' # spaceranger count.
#' transforms <- tibble(full_width = 18107, full_height = 19242) |>
#'   bind_cols(transforms)
#' transforms
#' # get example coordinate file
#' coordinatefile <- system.file("extdata/mousebrain/spatial",
#'                               "tissue_positions_list.csv",
#'                               package = "STUtility2")
#'
#' # Load coordinates
#' # These coordinates are defined on the H&E image used as input for
#' # spaceranger count.
#' xy <- LoadSpatialCoordinates(coordinatefiles = coordinatefile, verbose = T)
#' xy
#'
#' # Load image
#' lowresimagefile <- system.file("extdata/mousebrain/spatial",
#'                                "tissue_lowres_image.png",
#'                                package = "STUtility2")
#' im <- image_read(lowresimagefile)
#'
#' # Transform image and coordinates
#' transf_res <- RigidTransformImages(transforms, image = im, xy_coords = xy)
#'
#' ggplot(transf_res$xy_transf, aes(tr_x, tr_y)) +
#'   geom_point(color = "red", alpha = 0.5) +
#'   geom_segment(aes(x = 0, y = 0.5*19242, xend = 0.2*18107, yend = 0.5*19242),
#'                arrow = arrow(length = unit(0.5, "cm"))) +
#'   geom_vline(xintercept = 0.2*18107, linetype = "dashed") +
#'   geom_segment(aes(x = 0.5*18107, y = 19242, xend = 0.5*18107, yend = 0.8*19242),
#'                arrow = arrow(length = unit(0.5, "cm"))) +
#'   geom_hline(yintercept = 0.8*19242, linetype = "dashed") +
#'   scale_x_continuous(limits = c(0, 18107), expand = c(0, 0)) +
#'   scale_y_reverse(limits = c(19242, 0), expand = c(0, 0)) +
#'   labs(x = expression("x"["transformed"]),
#'        y = expression("y"["transformed"]),
#'        title = "20% right, 20% up, 30 deg rotation") +
#'   theme_void() +
#'   theme(axis.text = element_text(),
#'         axis.title.x = element_text(),
#'         axis.title.y = element_text(angle = 90)) +
#'   coord_fixed() +
#'   # Insert H&E image
#'   inset_element(p = as.raster(transf_res$im_transf),
#'                 left = 0, bottom = 0,
#'                 right = 1, top = 1,
#'                 on_top = FALSE)
#'
#' @export
#'
RigidTransformImages.default <- function (
  object,
  image,
  xy_coords,
  verbose = TRUE,
  ...
) {

  # Run checks for object
  if (!inherits(object, what = "tbl")) abort(glue("Invalid class {class(object)} of 'object', extected a 'tbl'"))
  if (nrow(object) != 1) abort(glue("Invalid dimensions for 'object', expected 1 row but got {nrow(object)}"))
  if (!all(c("full_width", "full_height", "sampleID",
         "mirror_x", "mirror_y", "angle", "tr_x", "tr_y") %in% colnames(object))) {
    abort("Missing columns in input tibble 'object'")
  }
  if (verbose) inform(c("i" = glue("Transforming image and coordinates for sample {object$sampleID}")))

  # Check image
  if (!inherits(image, what = c("StoredSpatialImage", "raster", "character", "magick-image")))
    abort(glue("Invalid class {class(image)} of 'image', expected one of ",
               "'StoredSpatialImage', 'raster', 'character', 'magick-image'"))

  # Check spatial coordinates
  if (!inherits(xy_coords, what = "tbl"))
    abort(glue("Invalid class {class(xy_coords)} of 'xy_coords', extected a 'tbl'"))
  if (!all(c("pxl_col_in_fullres", "pxl_row_in_fullres") %in% colnames(xy_coords)))
    abort("Couldn't find coordinates in 'xy_coords', expected 'pxl_col_in_fullres' and 'pxl_row_in_fullres'")
  xy_coords <- xy_coords |>
    select(all_of(c("pxl_col_in_fullres", "pxl_row_in_fullres")))
  if (verbose) inform(c("i" = glue("Fetched spot coordinates")))

  # Check transformations
  if (.check_transforms(object)) abort(glue("Transformations cannot take default values. ",
                                           "This would result in no transformation."))
  if (verbose) inform(c("i" = glue("Supplied transformations are valid")))

  # Make sure that translations are between -1-1
  checks <- sapply(object |> select(tr_x, tr_y), function(x) {
    between(x, left = -1, right = 1)
  })
  if (!all(checks)) abort("'tr_x' and 'tr_y' have to be between 0 and 1.")

  # Convert image to a "magick-image"
  if (inherits(image, what = "StoredSpatialImage")) image <- image_read(image@path)
  if (inherits(image, what = "raster")) image <- image_read(image)
  if (inherits(image, what = "character")) {
    if (!file.exists(image)) abort("File doesn't exist.")
    image <- image_read(image)
  }

  # make sure that image is of class "magick-image"
  if (!inherits(image, what = "magick-image")) abort("Invalid image format.")

  # get image info
  iminfo <- image_info(image)

  if (verbose)
    inform(c(
      ">" = glue("Mirror along x-axis: {object$mirror_x}"),
      ">" = glue("Mirror along x-axis: {object$mirror_y}"),
      ">" = glue("Rotation angle: {object$angle}"),
      ">" = glue("Translation along x axis: {object$tr_x*100}%"),
      ">" = glue("Translation along y axis: {object$tr_y*100}%")
    )
    )

  # Run transformation
  transf_res <- CoordAndImageTransform(
    im = image,
    xy_coords = xy_coords,
    mirror_x = object$mirror_x,
    mirror_y = object$mirror_y,
    angle = object$angle,
    xy_offset_image = c(round(object$tr_x*iminfo$width),
                        round(object$tr_y*iminfo$height)),
    xy_offset_spots = c(round(object$tr_x*object$full_width),
                        round(object$tr_y*object$full_height)),
    imcenter = c(object$full_width/2, object$full_height/2))

  if (verbose) inform(c("v" = "Returning transformed images"))

  # return results
  return(transf_res)

}

#' Check image transformation tibble
#'
#' @param object a `tibble` with transformation parameters
#'
#' @return a logical telling if transformations other
#' than default have been set
#'
.check_transforms <- function (
    object
) {
  checks <- c(!object$mirror_x, !object$mirror_y, object$angle == 0, object$tr_x == 0, object$tr_y == 0)
  return(all(checks))
}
