#' @include checks.R
#'
NULL


#' Find array transformation for Visium coordinates
#' 
#' The Visium coordinates are either defined by their array position in a 
#' grid or as pixel coordinates which map to a registered H&E image. The
#' Space Ranger pre-processing pipeline registers the coordinates to the
#' H&E image which can introduce rotations. \code{get_array_rotation} helps
#' calculating these rotation. 
#' 
#' @section Algorithm:
#' The first step is to take the array coordinates and the pixel coordinates 
#' and run center and scaling on both by subtracting the mean and dividing by 
#' the maximum radius. This will center the coordinates at (0, 0) and make sure 
#' that the coordinates share the same scale. As the two sets of points are paired, 
#' we can then use the Kabsch algorithm to find the rotation. 
#' 
#' The coordinates might be flipped along the x- or y-axis in which case you can 
#' use the \code{flip_x}, \code{flip_y} arguments to flip it back.
#' 
#' @param object An object of class \code{Seurat} created with \code{semla}
#' @param flip_x,flip_y Logical specifying if the x- or y-axis is flipped
#' @param grid_pattern One of "hexagonal" or "square". Use "hexagonal" for VisiumV1 data
#' @param return_aligned_coords If TRUE, a list is returned with the rotation angle and 
#' the aligned coordinates
#' @param verbose Print messages
#' 
#' @import rlang
#' @import dplyr
#' @import cli
#' 
#' @examples 
#' 
#' library(semla)
#' 
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
#' 
#' # Compare coordinates
#' # Note that the array coordinates needs an additional transformation 
#' # for VisiumV1 data since the grid is arranged in a hexagonal pattern
#' coords <- GetCoordinates(se_mcolon) |> 
#'   mutate(y = y*sqrt(3)) # Scale y values for hexaognal grids
#' p1 <- ggplot(coords, aes(x, y)) + geom_point() + 
#'    scale_y_reverse() + coord_fixed()
#' p2 <- ggplot(coords, aes(pxl_col_in_fullres, pxl_row_in_fullres)) + 
#'    geom_point() + scale_y_reverse() + coord_fixed()
#' p1 + p2
#' 
#' # The rotation between the two sets of coordinates is sutble,
#' # but we can find it with get_array_rotation
#' rotation_angle <- get_array_rotation(se_mcolon, grid_pattern = "hexagonal")
#' 
#' # generate transforms
#' transforms <- generate_rigid_transform(sampleID = 1L, angle = rotation_angle$`1`)
#' 
#' # Load H&E images
#' se_mcolon <- LoadImages(se_mcolon)
#' 
#' # Apply rotations
#' se_mcolon <- RigidTransformImages(se_mcolon, transforms = transforms)
#' 
#' # Plot coordinates
#' coords <- GetCoordinates(se_mcolon) |> 
#'   mutate(alpha = case_when(between(x = pxl_row_in_fullres_transformed, 
#'                                    left = 4280, right = 4360) ~ 1, TRUE ~ 0.2))
#' 
#' # The difference is quite subtle, but if you look closely you can see
#' # that the transformed coordinates (black) arranged in a square grid, 
#' # whereas the original coordinates (red) are slighlty rotated
#' ggplot(coords) +
#'   geom_point(aes(pxl_col_in_fullres, pxl_row_in_fullres, color = "original"), 
#'              size = 3, alpha = coords$alpha) +
#'   geom_point(aes(pxl_col_in_fullres_transformed, pxl_row_in_fullres_transformed, 
#'                  color = "transformed"), 
#'              alpha = coords$alpha) +
#'   coord_fixed() +
#'   scale_color_manual(values = c("original" = "red", "transformed" = "black")) +
#'   scale_y_reverse() +
#'   geom_segment(data = tibble(x = 2e3, xend = 8e3, y = 4320, yend = 4320),
#'                aes(x = x, xend = xend, y = y, yend = yend))
#' 
#' @export
#' 
get_array_rotation <- function (
  object,
  flip_x = FALSE,
  flip_y = FALSE,
  grid_pattern = c("hexagonal", "square"),
  return_aligned_coords = FALSE,
  verbose = TRUE
) {
  
  # Check Seurat object
  .check_seurat_object(object)
  
  # Check grid_pattern
  grid_pattern <- match.arg(grid_pattern, choices = c("hexagonal", "square"))
  
  # Check flip args
  if (flip_x & flip_y) abort("Only one of flip_x or flip_y can be set to TRUE")
  
  # Get pixel coordinates
  coords <- GetCoordinates(object)
  
  # get sampleIDs
  sampleIDs <- unique(coords$sampleID)
  
  if (!all(c("x", "y") %in% colnames(coords))) abort("Missing array coordinates 'x', 'y' from Staffli object `meta_data`")
  coords <- coords |> 
    group_by(.data[["sampleID"]]) |> 
    group_split()
  
  results <- lapply(seq_along(coords), function(i) {
    
    if (verbose) cli_alert_info("Finding rotation angle for sample {i}")
    xy <- coords[[i]]
    
    # Get basic array coordinates and scale
    if (verbose) cli_alert("   Scaling array coordinates")
    xy_arr <- xy |> select(all_of(c("x", "y"))) |> as.matrix() |> scale(center = TRUE, scale = FALSE)
    if (grid_pattern == "hexagonal") {
      xy_arr[, 2] <- xy_arr[, 2]*sqrt(3)
    }
    radius <- sqrt(rowSums(xy_arr^2))
    max_radius <- max(radius)
    xy_arr <- xy_arr/max_radius
    
    # Flip x axis
    if (flip_x) {
      xy_arr[, 1] <- -xy_arr[, 1]
    }
    if (flip_y) {
      xy_arr[, 2] <- -xy_arr[, 2]
    }
    
    # Scale pixel coordinates
    if (verbose) cli_alert("   Scaling pixel coordinates")
    xy_pixel <- xy |> select(all_of(c("pxl_col_in_fullres", "pxl_row_in_fullres"))) |> as.matrix() |> scale(center = TRUE, scale = FALSE)
    radius_pixel <- sqrt(rowSums(xy_pixel^2))
    max_radius_pixel <- max(radius_pixel)
    xy_pixel <- xy_pixel/max_radius_pixel
    
    # Run Kabsch algorithm to align points
    if (verbose) cli_alert("   Running Kabsch algorithm")
    kb <- kabsch(xy_pixel, xy_arr)
    
    # Get rotation anlge from rotation matrix
    angle <- ifelse(flip_x | flip_y, -(180/pi)*acos(kb$um[1, 1]), (180/pi)*acos(kb$um[1, 1]))
    if (verbose) cli_alert("   The rotation angle between the array and pixel coordinates is {round(angle, digits = 3)} degrees")
    
    if (return_aligned_coords) {
      aligned_coords <- kb$qm
      aligned_coords <- t(t(aligned_coords)*max_radius_pixel + attr(xy_pixel, "scaled:center"))
      return(list(angle = angle, aligned_coords = aligned_coords |> 
                    as_tibble() |> 
                    setNames(nm = c("pxl_col_in_fullres_transformed", "pxl_row_in_fullres_transformed")) |> 
                    mutate(barcode = xy$barcode) |> 
                    select(all_of(c("barcode", "pxl_col_in_fullres_transformed", "pxl_row_in_fullres_transformed")))))
    } else {
      return(setNames(angle, nm = "angle"))
    }
  }) |> setNames(nm = paste0(sampleIDs))
  
  return(results)
}
