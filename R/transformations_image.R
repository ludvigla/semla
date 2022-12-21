#' @title Function used to apply translations to an image object of class `magick-image`
#'
#' @description This function takes an `magick-image` object as input and applies translations
#'  defined by the \code{xy_offset} argument. The output image dimensions will remain the same
#'  as the input image dimensions, meaning that the transformation might result in cropping
#'  the image.
#'
#' @param im An image of class `magick-image`
#' @param xy_offset A numeric vector of length 2 providing the offsets along
#' the x- and y-axis given as pixels. These values should not exceed the image dimensions.
#'
#' @family transforms
#'
#' @author Ludvig Larsson
#'
#' @return An object of class `magick-image`
#'
#' @importFrom magick image_info image_blank image_crop image_append
#' @importFrom dplyr mutate case_when
#' @importFrom rlang %||%
#'
#' @examples
#'
#' library(magick)
#' library(STUtility2)
#' lowresimagefile <- system.file("extdata/mousebrain/spatial",
#'                                "tissue_lowres_image.jpg",
#'                                package = "STUtility2")
#' im <- image_read(lowresimagefile)
#'
#' # move image 100 pixels to the right and 100 pixels down
#' im_transformed <- ImageTranslate(im, xy_offset = c(100, 100))
#' im_transformed
#'
#' # move image 100 pixels to the left and 20 pixels up
#' im_transformed <- ImageTranslate(im, xy_offset = c(-100, 20))
#' im_transformed
#'
#' @export
ImageTranslate <- function (
    im,
    xy_offset = NULL
) {

  # get image info (width and height)
  info <- image_info(im)

  # Set xy_offset to (0, 0) if NULL
  xy_offset <- xy_offset %||% c(0, 0)

  # Define crop window depending on x, y offsets
  # x_crop & y_crop defines a "crop window" to trim image along the x or y axis
  # example: "100x120+20+30" means "cut a rectangle with the width of 100 pixels,
  # height of 120 pixels, offset by 20 pixels to the right and 30 pixels down
  # from the top left corner".
  crop_info <- data.frame(x_off = xy_offset[1], y_off = xy_offset[2]) |>
    mutate(x_crop = case_when(x_off < 0 ~ paste0(info$width - abs(xy_offset[1]), "x", info$height, "+", abs(xy_offset[1]), "+0"),
                              x_off > 0 ~ paste0(info$width - xy_offset[1], "x", info$height, "+0+0")),
           y_crop = case_when(y_off < 0 ~ paste0(info$width, "x", info$height - abs(xy_offset[2]), "+0+", abs(xy_offset[2])),
                              y_off > 0 ~ paste0(info$width, "x", info$height - xy_offset[2], "+0+0")))

  # If an x offset is provided, crop image and prepend/append white space (left or right)
  if (crop_info$x_off != 0) {
    im_margin <- image_blank(width = abs(xy_offset[1]), height = info$height) # Create empty space
    img_crop <- image_crop(im, geometry = crop_info$x_crop) # Crop based on geometry
    # Create an image stack
    if (crop_info$x_off < 0) {
      img_stack <- c(img_crop, im_margin) # prepend empty space
    } else {
      img_stack <- c(im_margin, img_crop) # append empty space
    }
    # Stack FALSE for merging along x axis
    im <- image_append(img_stack, stack = FALSE)
  }

  # If an y offset is provided, crop image and prepend/append white space (top or bottom)
  if (crop_info$y_off != 0) {
    im_margin <- image_blank(height = abs(xy_offset[2]), width = info$width) # Create empty space
    img_crop <- image_crop(im, geometry = crop_info$y_crop) # Crop based on geometry
    # Create an image stack
    if (crop_info$y_off < 0) {
      img_stack <- c(img_crop, im_margin) # prepend empty space
    } else {
      img_stack <- c(im_margin, img_crop) # append empty space
    }
    # Stack TRUE for merging along y axis
    im <- image_append(img_stack, stack = TRUE)
  }

  return(im)
}

#' @title Apply rotations and translations to an image
#'
#' @description This function takes an `magick-image` object as input and applies translations
#' defined by the \code{angle} and \code{xy_offset} arguments. The output image dimensions
#' will remain the same as the input image dimensions, meaning that the transformation might
#' result in cropping the image. If you don't want this behavior, you should use `image_rotate`
#' from the `magick` R package instead.
#'
#' @param im An image of class `magick-image`
#' @param angle An integer value specifying the rotation angle [-360, 360]
#' @param xy_offset A numeric vector of length 2 providing the offsets along
#' the x- and y-axis given as pixels. These values should not exceed the image dimensions.
#' @param scalefactor A numeric value specifying a scaling factor between [0, 3]
#'
#' @family tranforms
#'
#' @author Ludvig Larsson
#'
#' @importFrom magick image_info image_rotate image_blank image_extent image_composite geometry_area
#' @import dplyr
#' @import rlang
#' @import glue
#'
#' @return An object of class `magick-image`
#'
#' @examples
#'
#' library(magick)
#' library(STUtility2)
#' lowresimagefile <- system.file("extdata/mousebrain/spatial",
#'                                "tissue_lowres_image.jpg",
#'                                package = "STUtility2")
#' im <- image_read(lowresimagefile)
#'
#' # rotate image 45 degrees clockwise, move image 100 pixels to the right and 100 pixels down
#' im_transformed <- ImageTransform(im, angle = 45, xy_offset = c(100, 100))
#' im_transformed
#'
#' # rotate image 45 degrees counter-clockwise, move image 100 pixels to the left and 20 pixels up
#' im_transformed <- ImageTransform(im, angle = -45, xy_offset = c(-100, 20))
#' im_transformed
#'
#' # shrink image to helf width/height
#' im_transformed <- ImageTransform(im, scalefactor = 0.5)
#' im_transformed
#'
#' @export
ImageTransform <- function (
    im,
    angle = 0,
    xy_offset = c(0, 0),
    scalefactor = 1
) {
  # Check input parameters
  if (!inherits(im, what = "magick-image")) abort(glue("Invalid class '{class(im)}' for im, expected a 'magick-image'"))
  if (!inherits(angle, what = c("numeric", "integer"))) abort(glue("Invalid class '{class(angle)}' for angle, expected a 'numeric'"))
  if (length(angle) != 1) abort(glue("Invalid length '{length(angle)}' for angle, expected 1 value"))
  if (!between(x = angle, left = -360, right = 360)) abort(glue("Invalid value `{angle}`for angle. Value should be in range [-360, 360]"))
  if (!inherits(xy_offset, what = c("numeric", "integer"))) abort(glue("Invalid class '{class(xy_offset)}' for xy_offset, expected a 'numeric'"))
  if (length(xy_offset) != 2) abort(glue("Invalid length '{length(xy_offset)}' for angle, expected 2 values"))
  if (!inherits(scalefactor, what = c("numeric", "integer"))) abort(glue("Invalid class '{class(scalefactor)}' for scalefactor, expected a 'numeric'"))
  if (length(scalefactor) != 1) abort(glue("Invalid length '{length(scalefactor)}' for scalefactor, expected 1 value"))
  if (!between(x = scalefactor, left = 0, right = 3)) abort(glue("Invalid value `{scalefactor}`for scalefactor. Value should be in range [0, 3]"))

  if (angle == 0 & all(xy_offset == c(0, 0)) & scalefactor == 1) abort("Only default values provided which would result in no tranformation.")

  # Get width/height of input image
  info <- image_info(im)
  sf <- 1

  # Rotate image and get width/height of rotated image
  if (angle != 0) {
    im_rot <- image_rotate(im, degrees = angle)
    rot_info <- image_info(im_rot)
    # Create a blank image with the same dimensions as the rotated image
    im_blank <- image_blank(width = rot_info$width, rot_info$height, color = "#FFFFFFFF")
    # Place rotated image on top of blank image and apply translations
    # The first step is required to be able to crop the image later
    im <- image_composite(im_blank, composite_image = im_rot)
  }

  # Rescale image
  if (scalefactor != 1) {
    im_cur_info <- image_info(im)
    im_scaled <- im |> image_scale(geometry = geometry_area(width = im_cur_info$width*scalefactor,
                                                            height = im_cur_info$height*scalefactor))
    im_scaled_info <-  image_info(im_scaled)
    if (scalefactor < 1) {
      im_blank <- image_blank(width = im_cur_info$width, im_cur_info$height, color = "#FFFFFFFF")
      offset_x <- (im_cur_info$width - im_scaled_info$width)/2
      offset_y <- (im_cur_info$height - im_scaled_info$height)/2
      im_scaled <- image_composite(im_blank, composite_image = im_scaled,
                                   offset = geometry_area(x_off = offset_x, y_off = offset_y))
    } else {
      offset_x <- (im_scaled_info$width - im_cur_info$width)/2
      offset_y <- (im_scaled_info$height - im_cur_info$height)/2
      im_scaled <- im_scaled |> image_crop(geometry = geometry_area(width = im_cur_info$width,
                                                                    height = im_cur_info$height,
                                                                    x_off = offset_x,
                                                                    y_off = offset_y))
    }
    im <- im_scaled
    sf <- im_cur_info$height/info$height
  }

  # Apply translations
  if (any(xy_offset != c(0, 0))) {
    im <- im |>
      ImageTranslate(xy_offset)
  }

  if (angle != 0) {
    # Crop rotated image at center so that it get the same dimensions as the input image
    im <- image_extent(im, geometry = geometry_area(width = info$width, height = info$height))
  }

  return(im)
}
