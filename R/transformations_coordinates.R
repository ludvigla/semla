#' @include extdata.R
#'
NULL
#' Mirror coordinates
#'
#' @details The coordinate system for `xy_coords` should match the dimensions of the
#' image. In other words, the coordinates should map spots to the tissue section on H&E image.
#' A 3x3 transformation matrix is constructed by combining the following matrices:
#'
#' \itemize{
#'    \item{\eqn{T(-x, -y)}: Translate coordinates to origin, i.e. (0, 0) becomes the new center \cr
#'
#'    | 1 | 0 | \eqn{-center_{x}} |
#'    | - | - | :-: |
#'    | 0 | 1 | \eqn{-center_{y}} |
#'    | 0 | 0 | 1 |
#'
#'    }
#'    \item{\eqn{M_{x}}: Reflect coordinates along x-axis
#'
#'    | -1 | 0 | 0 |
#'    | - | - | - |
#'    | 0 | 1 | 0 |
#'    | 0 | 0 | 1 |
#'
#'    }
#'    \item{\eqn{M_{y}}: Reflect coordinates along y-axis \cr
#'
#'    | 1 | 0 | 0 |
#'    | - | - | - |
#'    | 0 | -1 | 0 |
#'    | 0 | 0 | 1 |
#'
#'    }
#'    \item{\eqn{T(x, y)}: Translate coordinates back to center \cr
#'
#'    | 1 | 0 | \eqn{center_{x}} |
#'    | - | - | :-: |
#'    | 0 | 1 | \eqn{center_{y}} |
#'    | 0 | 0 | 1 |
#'
#'    }
#' }
#'
#' Then, these matrices are combined to form the final transformation matrix:
#'
#' \eqn{T_{final} = T(x, y)*M_{y}*M_{x}*T(-x, -y)}
#'
#' Which can be used to transform our input coordinates:
#'
#' \eqn{xy_{out} = T_{final}*xy_{in}}
#'
#' @param xy_coords A `matrix`, `data.frame` or `tibble` object with numeric x, y coordinates.
#' @param mirror.x,mirror.y Logical specifying whether the coordinates should be reflected along
#' the x-axis, y-axis or both.
#' @param center Optional point (x, y) specifying the center of reflection.
#'
#' @import rlang
#'
#' @family transforms
#'
#' @author Ludvig Larsson
#'
#' @return A `tbl` object with transformed coordinates
#'
#' @examples
#' library(ggplot2)
#'
#' # Create a data.frame with x, y coordinates
#' xy <- data.frame(x = 1:20, y = 1:20)
#'
#' # Reflect coordinates along x axis
#' xy_mx <- CoordMirror(xy, mirror.x = TRUE) |> setNames(nm = c("x", "y"))
#'
#' # Reflect along both x, and y axes
#' xy_mxy <- CoordMirror(xy, mirror.x = TRUE, mirror.y = TRUE) |> setNames(nm = c("x", "y"))
#'
#' # Combine all coordinates
#' xy_all <- do.call(rbind, list(cbind(xy, type = "original", ord = 1:20),
#'                     cbind(xy_mx, type = "mirror_x", ord = 1:20),
#'                     cbind(xy_mxy, type = "mirror_x_and_y", ord = 1:20)))
#' xy_all$type <- factor(xy_all$type, levels = c("original", "mirror_x", "mirror_x_and_y"))
#'
#' # Now we can see the effects of mirroring
#' # Mirror x flips the coordinates along the x axis while mirror_x and mirror_y
#' # effectively inverts the coordinates, thus changing the order of the points
#' ggplot(xy_all, aes(x, y)) +
#'   geom_point(color = "steelblue", size = 7, alpha = 0.5) +
#'   geom_text(aes(label = ord)) +
#'   facet_grid(~type) +
#'   geom_vline(xintercept = 10.5, linetype = "dashed", color = "red") +
#'   geom_hline(yintercept = 10.5, linetype = "dashed", color = "green")
#'
#' @export
#'
#' @md
CoordMirror <- function (
    xy_coords,
    mirror.x = FALSE,
    mirror.y = FALSE,
    center = NULL
) {

  # Check xy_coords
  if (!any(c("tbl", "data.frame", "matrix") %in% class(xy_coords))) abort(glue::glue("Invalid format '{class(xy_coords)}' of xy_coords"))
  if (!ncol(xy_coords) == 2) abort(glue::glue("Expected 2 columns, got '{ncol(xy_coords)}'"))
  if (!all(sapply(xy_coords, is.numeric))) abort(glue::glue("Invalid column classes."))

  # Define center if NULL
  center <- center %||% colMeans(xy_coords)
  stopifnot(is.numeric(center), length(center) == 2)

  # Check mirror.x, mirror.y
  if (all(!c(mirror.x, mirror.y))) abort("One of 'mirror.x' or 'mirror.y' or both need to be selected.")

  # Create translation matrix to origin
  backTr <- matrix(c(1, 0, -center[1],
                     0, 1, -center[2],
                     0, 0, 1), byrow = TRUE, ncol = 3)

  # Create translation matrix back to center + optional offset
  forwardTr <- matrix(c(1, 0, center[1],
                        0, 1, center[2],
                        0, 0, 1), byrow = TRUE, ncol = 3)

  # Define mirror matrices
  mirrX <- matrix(c(ifelse(mirror.x, -1, 1), 0, 0,
                    0, 1, 0,
                    0, 0, 1),
                  byrow = TRUE, ncol = 3)
  mirrY <- matrix(c(1, 0, 0,
                    0, ifelse(mirror.y, -1, 1), 0,
                    0, 0, 1),
                  byrow = TRUE, ncol = 3)

  # Combine transformations
  combTr <- forwardTr%*%mirrY%*%mirrX%*%backTr

  # Apply transformations to input coordinates
  xy_coords_transformed <- (combTr%*%rbind(xy_coords |> t(), 1) |> t())[, 1:2] |>
    tibble::as_tibble(.name_repair = "minimal") |>
    setNames(nm  = c("tr_x", "tr_y"))

  return(xy_coords_transformed)
}


#' Apply transformations to a set of x, y coordinates
#'
#' @details The coordinate system for `xy_coords` should match the dimensions of the
#' image. In other words, the coordinates should map spots to the tissue section on H&E image.
#' Translations and rotations are done by multiplying `xy_coords` with the following
#' transformation matrix \eqn{T_{final}} as described below:
#'
#' \itemize{
#'    \item{\eqn{T(-x, -y)}: Translate coordinates to origin, i.e. (0, 0) becomes the new center \cr
#'
#'    | 1 | 0 | \eqn{-center_{x}} |
#'    | - | - | :-: |
#'    | 0 | 1 | \eqn{-center_{y}} |
#'    | 0 | 0 | 1 |
#'
#'    }
#'    \item{\eqn{R}: Rotate coordinates around origin \cr
#'
#'    | \eqn{cos(\alpha)} | \eqn{-sin(\alpha)} | 0 |
#'    | :-: | :-: | :-: |
#'    | \eqn{sin(\alpha)} | \eqn{cos(\alpha)} | 0 |
#'    | 0 | 0 | 1 |
#'
#'    }
#'    \item{\eqn{T(x, y)}: Translate coordinates back to center, or optionally to center + an `xy_offset` \cr
#'
#'    | 1 | 0 | \eqn{center_{x} + offset_{x}} |
#'    | - | - | :-: |
#'    | 0 | 1 | \eqn{center_{y} + offset_{y}} |
#'    | 0 | 0 | 1 |
#'
#'    }
#' }
#'
#' Then, these matrices are combined to form the final transformation matrix:
#'
#' \eqn{T_{final} = T(x, y)*R*T(-x, -y)}
#'
#' Which can be used to transform our input coordinates:
#'
#' \eqn{xy_{out} = T_{final}*xy_{in}}
#'
#' The scaling is handled separated after the translations and rotations.
#'
#' @param xy_coords A `matrix`, `data.frame` or `tibble` object with numeric x, y coordinates.
#' @param angle Numeric value specifying the degree of rotation. Use negative angles
#' for counter-clockwise rotation. The value needs to be in the range (-360, 360)
#' @param center Optional point (x, y) specifying the center of rotation.
#' @param xy_offset Optional point (x, y) specifying the translation.
#' @param scalefactor A numeric value specifying a scaling factor between (0, 3)
#'
#' @import rlang
#' @import glue
#'
#' @family transforms
#'
#' @author Ludvig Larsson
#'
#' @return A `tbl` object with transformed coordinates
#'
#' @examples
#' # Create a data.frame with x, y coordinates
#' xy <- data.frame(x = 1:20, y = 1:20)
#'
#' # Rotate coordinates 45 degrees clockwise around the center
#' xy_rotated <- CoordTransform(xy, angle = 45)
#' plot(xy)
#' points(xy_rotated, col = "red")
#'
#' @export
#'
#' @md
CoordTransform <- function (
  xy_coords,
  angle = 0,
  center = NULL,
  xy_offset = c(0, 0),
  scalefactor = 1
) {

  # Check input parameters
  if (!inherits(angle, what = c("numeric", "integer"))) abort(glue("Invalid class '{class(angle)}' for angle, expected a 'numeric'"))
  if (length(angle) != 1) abort(glue("Invalid length '{length(angle)}' for angle, expected 1 value"))
  if (!between(x = angle, left = -360, right = 360)) abort(glue("Invalid value `{angle}`for angle. Value should be in range [-360, 360]"))
  if (!is.null(center)) {
    if (!inherits(center, what = c("numeric", "integer"))) abort(glue("Invalid class '{class(center)}' for center, expected a 'numeric'"))
    if (length(center) != 2) abort(glue("Invalid length '{length(center)}' for angle, expected 2 values"))
  }
  if (!inherits(xy_offset, what = c("numeric", "integer"))) abort(glue("Invalid class '{class(xy_offset)}' for xy_offset, expected a 'numeric'"))
  if (length(xy_offset) != 2) abort(glue("Invalid length '{length(xy_offset)}' for angle, expected 2 values"))

  # Check xy_coords
  if (!any(c("tbl", "data.frame", "matrix") %in% class(xy_coords))) abort(glue::glue("Invalid format '{class(xy_coords)}' of xy_coords"))
  if (!ncol(xy_coords) == 2) abort(glue::glue("Expected 2 columns, got '{ncol(xy_coords)}'"))
  if (!all(sapply(xy_coords, is.numeric))) abort(glue::glue("Invalid column classes."))

  if (angle == 0 & all(xy_offset == c(0, 0)) & scalefactor == 1) abort("Only default values provided which would result in no tranformation.")

  # Define center if NULL
  center <- center %||% colMeans(xy_coords)

  # Define conversion function from degrees to radians
  deg2rad <- function(deg) {(deg * pi) / (180)}

  # Convert angle to radians
  rad <- deg2rad(angle)

  # Create rotation matrix
  rotTr <- matrix(c(cos(rad), -sin(rad), 0,
                     sin(rad), cos(rad), 0,
                     0, 0, 1), byrow = TRUE, ncol = 3)

  # Create translation matrix to origin
  backTr <- matrix(c(1, 0, -center[1],
                     0, 1, -center[2],
                     0, 0, 1), byrow = TRUE, ncol = 3)

  # Create translation matrix back to center + optional offset
  forwardTr <- matrix(c(1, 0, center[1] + xy_offset[1],
                     0, 1, center[2] + xy_offset[2],
                     0, 0, 1), byrow = TRUE, ncol = 3)

  # Combine transformations
  combTr <- forwardTr%*%rotTr%*%backTr

  # Apply transformations to input coordinates
  xy_coords_transformed <- (combTr%*%rbind(xy_coords |> t(), 1) |> t())[, 1:2] |>
    tibble::as_tibble(.name_repair = "minimal") |>
    setNames(nm  = c("tr_x", "tr_y"))

  # Scale coordinates
  if (scalefactor != 1) {
    xy_centered <- t(xy_coords_transformed) - (center + xy_offset)
    xy_scaled <- xy_centered*scalefactor
    xy_coords_transformed <- t(xy_scaled + (center + xy_offset)) |>
      tibble::as_tibble(.name_repair = "minimal")
  }

  return(xy_coords_transformed)
}

# TODO fix bug with CytAssist data when spots are located outside of H&E
#' Apply transformation to paired image and coordinates
#'
#' SRT data generated with Visium consists of an H&E image and a gene expression
#' matrix. The columns of the gene expression matrix correspond to spots that can
#' be mapped to the H&E image using a set of pixel coordinates.
#' This function takes an image and its corresponding set of spot coordinates as
#' input and applies the same transformation to the image and spot coordinates
#' simultaneously.
#'
#' @details
#' Mirroring is prioritized and will be applied to the image before applying rotations
#' and translations.
#'
#' @param mirror_x,mirror_y Logical specifying if the image and spots should be mirrored
#' along the x- and/or y-axis
#' @param im Image of class `magick-image`, `StoredSpatialImage`, `raster` or a
#' path to an external image file.
#' @param imcenter A numeric vector of length 2 specifying the image center. Not required
#' if the spot coordinates are already aligned to the H&E image.
#' @param xy_offset_image A numeric vector of length 2 specifying the translations
#' along the x-a dn y-axes for the image
#' @param xy_offset_spots A numeric vector of length 2 specifying the translations
#' along the x- and y-axes for the spots. If \code{xy_offset_image} is \code{NULL},
#' \code{xy_offset_spots} will be set to \code{xy_offset_image} as it is assumed that
#' the spots are matched with the image.
#'
#' @inheritParams CoordTransform
#'
#' @importFrom magick image_flop image_flip
#' @importFrom rlang %||%
#'
#' @family transforms
#'
#' @author Ludvig Larsson
#'
#' @return A list with two elements:
#' \itemize{
#'   \item{"im_transf": An object of class `magick-image` representing the transformed image}
#'   \item{"xy_transf": An object of class `tbl` representing the transformed coordinates}
#' }
#'
#' @examples
#'
#' library(STUtility2)
#' library(ggplot2)
#' library(patchwork)
#' library(magick)
#' library(dplyr)
#'
#' # get example coordinate file
#' coordinatefile <- system.file("extdata",
#'                               "mousebrain/spatial/tissue_positions_list.csv",
#'                               package = "STUtility2")
#'
#' # Load coordinates
#' xy <- LoadSpatialCoordinates(coordinatefiles = coordinatefile)
#' xy
#'
#' # Load image
#' lowresimagefile <- system.file("extdata/mousebrain/spatial",
#'                                "tissue_lowres_image.jpg",
#'                                package = "STUtility2")
#' im <- image_read(lowresimagefile)
#'
#' # read scalefactors
#' scalefactorfile <- system.file("extdata/mousebrain/spatial",
#'                                "scalefactors_json.json",
#'                                package = "STUtility2")
#' scalefactors <- jsonlite::read_json(scalefactorfile)
#' scalefactors
#'
#' # Convert coordinates using appropriate scalefactor, in this
#' # case the scalefactor for tissue_image_lowres
#' xy <- xy |>
#'   mutate(across(pxl_col_in_fullres:pxl_row_in_fullres,
#'                 ~ .x*scalefactors$tissue_lowres_scalef))
#'
#' # Note that the y axis needs to be reversed and you need to
#' # specify the axis limits uisng the dimensions of the image
#' ggplot(xy, aes(pxl_col_in_fullres, pxl_row_in_fullres)) +
#'   geom_point(color = "red", alpha = 0.5) +
#'   scale_x_continuous(limits = c(0, image_info(im)$width), expand = c(0, 0)) +
#'   scale_y_reverse(limits = c(image_info(im)$height, 0), expand = c(0, 0)) +
#'   labs(x = expression("x"["original"]),
#'        y = expression("y"["original"]),
#'        title = "Original image and coordinates") +
#'   theme_void() +
#'   theme(axis.text = element_text(),
#'         axis.title.x = element_text(),
#'         axis.title.y = element_text(angle = 90)) +
#'   coord_fixed() +
#'   # Insert H&E image
#'   inset_element(p = as.raster(im),
#'                 left = 0, bottom = 0,
#'                 right = 1, top = 1,
#'                 on_top = FALSE)
#'
#' # Select coordinates for transformation
#' xy_coords <- xy |>
#'   select(pxl_col_in_fullres, pxl_row_in_fullres)
#'
#' # Apply transformations
#' transf_res <- CoordAndImageTransform(im, xy_coords, angle = 45, xy_offset_image = c(100, 100))
#'
#' # Add selected to transf_res$xy_transf
#' transf_res$xy_transf$selected <- xy$selected
#'
#' # Note that the y axis needs to be reversed and you need to
#' # specify the axis limits uisng the dimensions of the image
#' ggplot(transf_res$xy_transf, aes(tr_x, tr_y)) +
#'   geom_point(color = "red", alpha = 0.5) +
#'   scale_x_continuous(limits = c(0, image_info(im)$width), expand = c(0, 0)) +
#'   scale_y_reverse(limits = c(image_info(im)$height, 0), expand = c(0, 0)) +
#'   labs(x = expression("x"["transformed"]),
#'        y = expression("y"["transformed"]),
#'        title = "Transformed image and coordinates") +
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
CoordAndImageTransform <- function (
  im,
  xy_coords,
  mirror_x = FALSE,
  mirror_y = FALSE,
  angle = 0,
  scalefactor = 1,
  imcenter = NULL,
  xy_offset_image = c(0, 0),
  xy_offset_spots = NULL
) {

  # Set global variables to NULL
  tr_x <- tr_y <- NULL

  # Set spots offset if NULL
  xy_offset_spots <- xy_offset_spots %||% xy_offset_image

  # Read image
  if (inherits(im, what = "StoredSpatialImage")) im <- image_read(im@path)
  if (inherits(im, what = "raster")) im <- image_read(im)
  if (inherits(im, what = "character")) {
    if (!file.exists(im)) abort("File doesn't exist.")
    im <- image_read(im)
  }

  # make sure that image is of class "magick-image"
  if (!inherits(im, what = "magick-image")) abort("Invalid image format.")

  # Get image width and height and calculate the center of the image
  info <- image_info(im)
  imcenter <- imcenter %||% c(info$width/2, info$height/2)

  # Mirror image
  if (mirror_x | mirror_y) {
    if (mirror_x) {
      im <- image_flop(im)
    }
    if (mirror_y) {
      im <- image_flip(im)
    }
    xy_coords_mirror <- xy_coords |>
      CoordMirror(mirror.x = mirror_x,
                  mirror.y = mirror_y,
                  center = imcenter)

    # Set new im and xy_coords
    xy_coords <- xy_coords_mirror
  }


  # Apply transformations
  if (!(angle == 0 & all(xy_offset_image == c(0, 0)) & scalefactor == 1)) {
    im_transformed <- ImageTransform(im = im, angle = angle, xy_offset = xy_offset_image, scalefactor = scalefactor)
    xy_coords_transformed <- CoordTransform(xy_coords, angle, center = imcenter, xy_offset = xy_offset_spots, scalefactor = scalefactor)
    xy_coords_transformed <- xy_coords_transformed |>
      select(tr_x, tr_y)
  } else {
    im_transformed <- im
    xy_coords_transformed <- xy_coords
  }

  # Return list of results
  return(list(im_transf = im_transformed, xy_transf = xy_coords_transformed))
}
