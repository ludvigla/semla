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
#' @importFrom rlang %||%
#'
#' @family transforms
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
#' # Now we can sse the effects of mirroring
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


#' Rotate a set of x, y coordinates around a center point
#'
#' @details The coordinate system for `xy_coords` should match the dimensions of the
#' image. In other words, the coordinates should map spots to the tissue section on H&E image.
#' The transformation is conducted by creating a transformation matrix:
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
#' Then, these matrcies are combined to form the final transformation matrix:
#'
#' \eqn{T_{final} = T(x, y)*R*T(-x, -y)}
#'
#' Which can be used to transform our input coordinates:
#'
#' \eqn{xy_{out} = T_{final}*xy_{in}}
#'
#' @param xy_coords A `matrix`, `data.frame` or `tibble` object with numeric x, y coordinates.
#' @param angle Numeric value specifying the degree of rotation. Use negative angles
#' for  counter-clockwise rotation.
#' @param center Optional point (x, y) specifying the center of rotation.
#' @param xy_offset Optional point (x, y) specifying the translation.
#'
#' @importFrom rlang %||%
#'
#' @family transforms
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
  xy_offset = NULL
) {

  # Check xy_coords
  if (!any(c("tbl", "data.frame", "matrix") %in% class(xy_coords))) abort(glue::glue("Invalid format '{class(xy_coords)}' of xy_coords"))
  if (!ncol(xy_coords) == 2) abort(glue::glue("Expected 2 columns, got '{ncol(xy_coords)}'"))
  if (!all(sapply(xy_coords, is.numeric))) abort(glue::glue("Invalid column classes."))

  # Define center if NULL
  center <- center %||% colMeans(xy_coords)
  stopifnot(is.numeric(center), length(center) == 2)

  # Define translation if NULL
  xy_offset <- xy_offset %||% c(0, 0)
  stopifnot(is.numeric(xy_offset), length(xy_offset) == 2)

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
    setNames(nm  = c("tr_x", "tr_y")) #|>
    #select(tr_x, tr_y)

  return(xy_coords_transformed)
}


#' Apply transformation to paired image and coordinates
#'
#' SRT data generated with Visium consists of an H&E image and a gene expression
#' matrix. The columns of the gene expression matrix correspond to spots that can
#' be mapped to the H&E image using a set of pixel coordinates.
#' This function takes an image and its corresponding set of spot coordinates as
#' input and applies the same transformation to the image and spot coordinates
#' simultaneously.
#'
#' @param im Image of class `magick-image`, `StoredSpatialImage`, `raster` or a
#' path to an external image file.
#'
#' @inheritParams CoordTransform
#'
#' @family transforms
#'
#' @return A list with two elements:
#' \itemize{
#'   \item{"im_transf": An object of class `magick-image` representing the transformed image}
#'   \item{"xy_transf": An object of class `tble` representing the transformed coordinates}
#' }
#'
#' @export
CoordAndImageTransform <- function (
  im,
  xy_coords,
  angle = 0,
  xy_offset = NULL
) {

  # Set global variables to NULL
  tr_x <- tr_y <- NULL

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
  imcenter <- c(info$height - info$height/2, info$width/2)

  # Apply transformations
  im_transformed <- ImageTransform(im, angle, xy_offset)
  xy_coords_transformed <- CoordTransform(xy_coords, -angle, center = imcenter, xy_offset = xy_offset)
  xy_coords_transformed <- xy_coords_transformed |>
    dplyr::rename(tr_x = tr_y, tr_y = tr_x) |>
    select(tr_x, tr_y)

  # Return list of results
  return(list(im_transf = im_transformed, xy_transf = xy_coords_transformed))
}
