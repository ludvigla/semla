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
#' frink <- image_read("https://jeroen.github.io/images/frink.png")
#' frink
#'
#' # move image 100 pixels to the right and 100 pixels down
#' frink_transformed <- ImageTranslate(frink, xy_offset = c(100, 100))
#' frink_transformed
#'
#' # move image 100 pixels to the left and 20 pixels up
#' frink_transformed <- ImageTranslate(frink, xy_offset = c(-100, 20))
#' frink_transformed
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

#' @title Function used to apply rotations and translations to an image object of class `magick-image`
#'
#' @description This function takes an `magick-image` object as input and applies translations
#'  defined by the \code{angle} and \code{xy_offset} arguments. The output image dimensions
#'  will remain the same as the input image dimensions, meaning that the transformation might
#'  result in cropping the image. If you don't want this behavior, you should use `image_rotate`
#'  from the `magick` R package instead.
#'
#' @param im An image of class `magick-image`
#' @param angle An integer value specifying the rotation angle (0-360)
#' @param xy_offset A numeric vector of length 2 providing the offsets along
#' the x- and y-axis given as pixels. These values should not exceed the image dimensions.
#'
#' @family tranforms
#'
#' @importFrom magick image_info image_rotate image_blank image_extent image_composite
#'
#' @return An object of class `magick-image`
#'
#' @examples
#'
#' library(magick)
#' library(STUtility2)
#' frink <- image_read("https://jeroen.github.io/images/frink.png")
#' frink
#'
#' # rotate image 45 degrees clockwise, move image 100 pixels to the right and 100 pixels down
#' frink_transformed <- ImageTransform(frink, angle = 45, xy_offset = c(100, 100))
#' frink_transformed
#'
#' # rotate image 45 degrees counter-clockwise, move image 100 pixels to the left and 20 pixels up
#' frink_transformed <- ImageTransform(frink, angle = -45, xy_offset = c(-100, 20))
#' frink_transformed
#'
#' @export
ImageTransform <- function (
    im,
    angle,
    xy_offset
) {
  # Get width/height of input image
  info <- image_info(im)
  # Rotate image and get width/height of roaated image
  im_rot <- image_rotate(im, degrees = angle)
  rot_info <- image_info(im_rot)
  # Create a blank image with the same dimensions as the rotated image
  im_blank <- image_blank(width = rot_info$width, rot_info$height, color = "#FFFFFFFF")
  # Place rotated image on top of blank image and apply translations
  # The first step is required to be able to crop the image later
  im_combined <- image_composite(im_blank, composite_image = im_rot) |>
    ImageTranslate(xy_offset)
  # Crop rotated image at center so that it get the same dimensions as the input image
  gmtry <- paste0(info$width, "x", info$height)
  im_crop <- image_extent(im_combined, geometry = gmtry)
  return(im_crop)
}
