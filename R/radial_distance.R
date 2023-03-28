#' @include generics.R
#' @include checks.R
#' @include spatial_utils.R
#' @include disconnect_regions.R
#'
NULL

#' @description Calculates the radial distances to all spots from the borders of
#' a selected defined region.
#'
#' @section Scenario:
#' The region of interest could for example be an isolated tumor in the tissue
#' section surrounded by stroma. If we are interested in expressional changes
#' from the tumor core and outwards, we can use radial distances to model such
#' changes. Below are a few examples for how radial distances can be used to in
#' address certain question in this scenario:
#'
#' \itemize{
#'    \item{Identify genes that vary inside the tumor, e.g. tumor edge vs tumor core.}
#'    \item{Identify cell types located at the tumor edge}
#'    \item{Characterize the surrounding tumor microenvironment just outside the tumor edge}
#'    \item{Identify cell types whose abundances change with distance to tumor}
#' }
#'
#' @section Algorithm:
#' First, the border spots of the selected region is identified based on the
#' spatial network of nearest neighbors identified with \code{\link{GetSpatialNetwork}}.
#' For each spot outside of this border, the distance is calculated to its nearest border spot.
#' Spots located inside the selected region will have negative distances and
#' spots located outside of the selected region will have positive distances.
#'
#' @section Search interval:
#' The microenvironment of the region of interest might be extremely heterogeneous
#' depending on the direction from its center. For this reason, it can be useful to narrow
#' down the search area by defining a smaller angle interval with \code{angles}. Alternatively,
#' you can split the radial distances into an even number of slices with \code{angles_nbreaks}.
#' When using a predefined search interval, the region of interest (e.g. manual annotation)
#' should not contain multiple spatially disconnected regions. Angles are calculated from
#' center of the region of interest so it only makes sense to investigate one region at the time.
#' You can use \code{\link{DisconnectRegions}} to split a categorical variable that contains
#' multiple spatially disconnected regions.
#'
#' @section default method:
#' \code{object} should be a tibble with four columns:
#'
#' \itemize{
#'  \item{'barcode' : character vector with spot IDs}
#'  \item{'x', 'y' : numeric vectors with pixel coordinates}
#'  \item{'sampleID' : numeric vector with sample IDs}
#' }
#'
#' If \code{angles} and/or \code{angles_nbreaks} are set, the function will return a
#' tibble with spot IDS, sampleIDs, the angle between the region center and
#' spots and the radial distances. If  \code{angles_nbreaks} is provided, a
#' fifth column will be provided that groups spots into even intervals based on angles.
#' Otherwise, the default is to return numeric vector with radial distances
#' for all spots in \code{object}.
#' 
#' @section Seurat method:
#' If \code{object} is is a \code{Seurat} object created with \code{semla},
#' the results are returned to the \code{meta.data} slot.
#'
#' @param spots A character vector with spot IDs present \code{object}. These spots typically
#' represent one particular tissue structure identified either by data-driven clustering
#' or by the tissue histology.
#' @param angles A numeric vector of length 2 specifying a "search interval" of angles
#' to compute the radial distances for. Values between 0 and 360 are accepted where
#' \code{angles[1] < angles[2]}. The angles are defined in a clockwise manner, where right=0,
#' down=90, left=180 and up=270. The angles are calculated relative to the region
#' center and can therefore only be used when a single connected region is present. If
#' there are multiple, spatially disconnected regions present, use \code{\link{DisconnectRegions}}
#' to split the spatially disconnected regions first.
#' @param angles_nbreaks An integer specifying a number of intervals to cut the "search interval"
#' into. This can be useful if you want to group radial distances into different directions from
#' the region center.
#' @param remove_singletons Logical specifying if 'singletons' should be excluded. Spatially disconnected
#' regions are not allowed when a "search interval" is defined, but single spots without neighbors are
#' not detected as disconnected components. Single spots will most likely not interfere when calculating
#' the centroid of the region of interest and can therefore be kept.
#' @param convert_to_microns Logical specifying if pixel distances should be converted to microns. 
#' This requires the \code{dbscan} R package to be installed. When this option sis set to TRUE, the method
#' will first attempt to estimate the center to center distance between adjacent spots which 
#' corresponds to 100 microns in Visium. The center to center distance will then be used to 
#' convert the radial distances. Note that if no spots are adjacent in the dataset or if any other data 
#' type than Visium is used, the distances will not correspond to micrometers.
#'
#' @import dplyr
#' @importFrom rlang abort warn
#' @import cli
#' @importFrom glue glue
#' @importFrom dbscan kNN
#' @importFrom zeallot %<-%
#'
#' @rdname radial-distance
#'
#' @examples
#'
#' library(semla)
#' library(ggplot2)
#' library(patchwork)
#' library(RColorBrewer)
#'
#' # Get coordinates
#' galt_spots_file <- "~/semla/repo/semla/inst/extdata/mousecolon/galt_spots.csv"
#' galt_spots <- read.csv(galt_spots_file) |>
#'   as_tibble()
#'
#' # read coordinates
#' coordfile <- system.file("extdata/mousecolon/spatial",
#'                          "tissue_positions_list.csv",
#'                          package = "semla")
#' coords <- read.csv(coordfile, header = FALSE) |>
#'   filter(V2 == 1) |>
#'   select(V1, V6, V5) |>
#'   setNames(nm = c("barcode", "x", "y")) |>
#'   bind_cols(sampleID = 1) |>
#'   as_tibble()
#'
#' # Select spots
#' spots <- galt_spots$barcode[galt_spots$selection == "GALT"]
#' head(spots)
#'
#' # Calculate radial distances
#' radial_distances <- RadialDistance(coords, spots)
#' gg <- bind_cols(coords, r_dist =  radial_distances) |>
#'   left_join(y = galt_spots, by = "barcode")
#'
#' # Convert to sqrt scale
#' gg$r_dist_sqrt <- sign(gg$r_dist)*sqrt(abs(gg$r_dist))
#'
#' # Make plot
#' p1 <- ggplot(gg, aes(x, y, color = r_dist_sqrt)) +
#'   geom_point() +
#'   scale_y_reverse() +
#'   scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
#' p2 <- ggplot(gg, aes(x, y, color = selection)) +
#'   geom_point() +
#'   scale_y_reverse()
#'
#' # Wrap plots
#' wrap_plots(p2, p1, ncol = 2) &
#'   coord_fixed() &
#'   theme_void()
#'
#' \donttest{
#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' # Calculate radial distances for fixed angle interval
#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' # NB: This should only be run on a single region! In
#' # this example, the disconnected regions have to be split first
#' disconnected_regions <- DisconnectRegions(coords, spots)
#' spots_keep <- names(disconnected_regions[disconnected_regions == "S1_region1"])
#'
#' # Calculate radial distances between 200-300 degrees
#' radial_distances <- RadialDistance(coords, spots_keep, angles = c(200, 300))
#' gg <- coords |>
#'   select(-sampleID) |>
#'   left_join(y = radial_distances, by = "barcode")
#' gg$r_dist_sqrt <- sign(gg$r_dist)*sqrt(abs(gg$r_dist))
#'
#' # Plot radial distances
#' ggplot(gg, aes(x, y, color = r_dist_sqrt)) +
#'   geom_point() +
#'   scale_y_reverse() +
#'   coord_fixed() +
#'   theme_void() +
#'   scale_color_gradientn(colours = viridis::magma(n = 9))
#'
#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' # Calculate radial distances for multiple angle intervals
#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' # NB: This should only be run on a single region! In
#' # this example, the disconnected regions have to be split first
#' disconnected_regions <- DisconnectRegions(coords, spots)
#' spots_keep <- names(disconnected_regions[disconnected_regions == "S1_region1"])
#'
#' # Calculate radial distances between 0-360 degrees and split
#' # these into 8 slices
#' radial_distances <- RadialDistance(coords, spots_keep,
#'                                    angles = c(0, 360), angles_nbreaks = 8)
#' gg <- coords |>
#'   select(-sampleID) |>
#'   left_join(y = radial_distances, by = "barcode")
#' gg$r_dist_sqrt <- sign(gg$r_dist)*sqrt(abs(gg$r_dist))
#'
#' # Color slices
#' p1 <- ggplot(gg, aes(x, y, color = intervals)) +
#'   geom_point() +
#'   scale_y_reverse() +
#'   coord_fixed() +
#'   theme_void() +
#'   scale_color_manual(values = RColorBrewer::brewer.pal(n = 8, name = "Spectral"))
#'
#' # Plot distances
#' p2 <- ggplot(gg, aes(intervals, r_dist)) +
#'   geom_jitter()
#'
#' # Now we can group our radial distances by slice
#' p1 + p2
#' }
#'
#' @export
#'
RadialDistance.default <- function (
  object,
  spots,
  angles = NULL,
  angles_nbreaks = NULL,
  remove_singletons = TRUE,
  convert_to_microns = FALSE,
  verbose = TRUE,
  ...
) {

  # Set global variables to NULL
  barcode <- x <- y <- sampleID <- from  <- to <- centroids <- diff_x <- diff_y <- angle <- intervals <- r_dist <- nn <- NULL

  # Check object class
  if (!any(class(object) %in% c("data.frame", "matrix", "tbl")))
    abort(glue("Invalid class '{class(object)}'."))
  if (ncol(object) != 4)
    abort(glue("Invalid number of columns '{ncol(object)}'. Expected 4."))
  if (!all(colnames(object) == c("barcode", "x", "y", "sampleID")))
    abort("Required columns are: 'barcode', 'x', 'y', and 'sampleID'")
  if (!all(
    object |> summarize(
      check_barcode = is.character(barcode),
      check_x = is.numeric(x),
      check_y = is.numeric(y),
      check_sample = is.numeric(sampleID)
    ) |>
    unlist()
  )) {
    abort(glue("Invalid column class(es)."))
  }

  # Check that spots are in object
  stopifnot(
    inherits(spots, what = "character"),
    length(spots) > 0,
    all(spots %in% object$barcode)
  )

  # Get spatial network
  spatnet <- GetSpatialNetwork(object)
  if (length(spatnet) > 1) {
    abort(glue("Default method can only handle 1 tissue section at the time, got {length(spatnet)}"))
  }
  spatnet_region <- spatnet[[1]] |>
    filter(from %in% spots, to %in% spots) |>
    group_by(from) |>
    mutate(nn = n())

  # Check for singletons
  if (remove_singletons) {
    spatnet_region2 <- spatnet_region |> filter(nn > 1)
    spots_filtered <- spatnet_region2$from |> unlist() |> unique()
    if (verbose) cli_alert_info("Removing {length(spots) - length(spots_filtered)} spots with 0 neighbors.")
    spots <- spots_filtered
  }

  if (verbose) cli_alert_info("Extracting border spots from a region with {length(spots)} spots")

  # Find border
  border_spots <- RegionNeighbors(spatnet, spots = spots, outer_border = FALSE)
  inside_spots <- setdiff(spots, border_spots)

  # Get indices for spot groups
  border_spots_indices <- match(border_spots, object$barcode)
  outside_spots_indices <- -match(spots, object$barcode)
  inside_spots_indices <- match(inside_spots, object$barcode)
  if (verbose) {
    cli_alert("  Detected {length(border_spots_indices)} spots on borders")
    cli_alert("  Detected {length(inside_spots_indices)} spots inside borders")
    cli_alert("  Detected {nrow(object) - length(spots)} spots outside borders")
  }

  # Find nearest neighbor at selected spots border for
  # spot outside and inside the border
  knn_spatial_outside <- kNN(x = object[border_spots_indices, c("x", "y")] |> as.matrix(),
                             query = object[outside_spots_indices, c("x", "y")] |> as.matrix(),
                             k = 1)
  knn_spatial_inside <- kNN(x = object[border_spots_indices, c("x", "y")] |> as.matrix(),
                            query = object[inside_spots_indices, c("x", "y")] |> as.matrix(),
                            k = 1)
  if (verbose) cli_alert_success("Returning radial distances")

  # Get radial distances
  radial_dists <- setNames(numeric(length = nrow(object)), nm = object$barcode)
  radial_dists_outside <- setNames(knn_spatial_outside$dist[, 1, drop = TRUE], nm = object$barcode[outside_spots_indices])
  radial_dists_inside <- setNames(knn_spatial_inside$dist[, 1, drop = TRUE], nm = object$barcode[inside_spots_indices])
  radial_dists[names(radial_dists_outside)] <- radial_dists_outside
  radial_dists[names(radial_dists_inside)] <- -radial_dists_inside
  
  # Convert to microns
  if (convert_to_microns) {
    if (!requireNamespace("dbscan")) {
      abort(glue("Package {cli::col_br_magenta('dbscan')} is required. Please install it with: \n",
                 "install.packages('dbscan')"))
    }
    center_to_center_pixel_distances <- sapply(object |> group_by(sampleID) |> group_split(), function(xy) {
      kNN(x = xy |> select(x, y), k = 1)$dist |> min()
    })
    radial_dists <- radial_dists/(center_to_center_pixel_distances/100)
  }

  # Split radial distances by angle if angles are specified
  if (!is.null(angles_nbreaks)) {
    stopifnot(is.numeric(angles_nbreaks),
              length(angles_nbreaks) == 1)
    angles <- angles %||% c(0, 360)
  }
  if (!is.null(angles)) {
    # Check angles
    c(angles, centroids) %<-% .check_angles(list(spatnet_region) |> setNames(nm = "1"), object, spots, angles, angles_nbreaks, verbose)
    # Calculate angles from center point
    dist_angle <- object |>
      group_by(sampleID) |>
      group_split() |>
      setNames(nm = names(spatnet))
    dist_angle <- do.call(bind_rows, lapply(names(dist_angle), function(nm) {
      dist_angle[[nm]] |>
        mutate(diff_x = x - centroids[[nm]][1],
             diff_y = y -  centroids[[nm]][2]) |>
        mutate(angle = atan2(diff_y, diff_x)*(180/pi)) |>
        mutate(angle = case_when(angle <= 0 ~ angle + 360, TRUE ~ angle)) |>
        select(-diff_x, -diff_y) |>
        filter(between(x = angle, left = angles[1], right = angles[2]))
    }))
    # Split data if angled_length_out is set
    if (!is.null(angles_nbreaks)) {
      dist_angle <- dist_angle |>
        mutate(intervals = cut(angle, breaks = seq(angles[1], angles[2], length.out = angles_nbreaks + 1))) |>
        arrange(intervals) |>
        mutate(intervals = factor(intervals, levels = unique(intervals)))
      dist_angle$r_dist <- radial_dists[dist_angle$barcode]
    } else {
      dist_angle$r_dist <- radial_dists[dist_angle$barcode]
    }
    radial_dists <- dist_angle |> select(barcode, sampleID, angle, r_dist, contains("intervals"))
  }

  # Return radial distances
  return(radial_dists)

}


#' @param column_name A character specifying the name of a column in your meta data that contains
#'  categorical data, e.g. clusters or manual selections
#' @param selected_groups A character vector to select specific groups in \code{column_name} with.
#' All groups are selected by default, but the common use case is to select a region of interest.
#' @param verbose Print messages
#'
#' @import dplyr
#' @import cli
#' @importFrom rlang abort
#' @importFrom glue glue
#'
#' @rdname radial-distance
#' @family spatial-methods
#'
#' @author Ludvig Larsson
#'
#' @examples
#'
#' library(semla)
#' library(ggplot2)
#' library(patchwork)
#' library(tidyr)
#' library(RColorBrewer)
#'
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
#' se_mcolon <- RadialDistance(se_mcolon, column_name = "selection", selected_groups = "GALT")
#'
#' # Plot results
#' p1 <- MapLabels(se_mcolon, column_name = "selection")
#' p2 <- MapFeatures(se_mcolon, features = "r_dist_GALT", colors = c("lightgray", "black"))
#' p1 | p2
#'
#' #  Plot expression as function of distance
#' sel_genes <- c("Clu", "Tagln")
#' gg <- FetchData(se_mcolon, vars = c(sel_genes, "r_dist_GALT")) |>
#'   pivot_longer(all_of(sel_genes), names_to = "variable", values_to = "value")
#'
#' # Plot features
#' p1 <- MapFeatures(se_mcolon, features = sel_genes)
#'
#' # Plot expression as a function of distance
#' p2 <- ggplot(gg, aes(r_dist_GALT, value, color = variable)) +
#'   geom_smooth(method = "loess", span = 0.2, formula = y ~ x) +
#'   geom_vline(xintercept = 0, linetype = "dashed") +
#'   theme_minimal()
#'
#' # Combine plots
#' p1/p2
#'
#' # It can also be useful to apply transformations to the distances
#' se_mcolon$r_dist_GALT_sqrt <- sign(se_mcolon$r_dist_GALT)*sqrt(abs(se_mcolon$r_dist_GALT))
#' MapFeatures(se_mcolon, features = "r_dist_GALT_sqrt", pt_size = 2,
#'             colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev())
#'
#' @export
#'
RadialDistance.Seurat <- function (
    object,
    column_name,
    selected_groups = NULL,
    angles = NULL,
    angles_nbreaks = NULL,
    remove_singletons = TRUE,
    convert_to_microns = FALSE,
    verbose = TRUE,
    ...
) {

  # Set global variables to NULL
  barcode <- pxl_col_in_fullres <- pxl_row_in_fullres <- sampleID <- angle <- r_dist <- NULL

  # validate input
  .check_seurat_object(object)
  .validate_column_name(object, column_name)
  selected_groups <- .validate_selected_labels(object, selected_groups, column_name)

  # Select spots
  spots_list <- .get_spots_list(object, selected_groups, column_name)

  # Get coordinates
  coords_list <- .get_coords_list(object)

  # Calculate radial distances
  radial_distances <- do.call(bind_rows, lapply(names(coords_list), function(nm) {
    if (verbose) cli_alert_info("Running calculations for sample {nm}")
    sample_radial_distances <- lapply(names(spots_list), function(lbl) {
      if (verbose) cli_alert_info("Calculating radial distances for group '{lbl}'")
      if (is.null(spots_list[[lbl]][[nm]])) {
        if (verbose) cli_alert_warning("Found no spots for groups '{lbl}' in section {nm}. Returning NA values for section {nm}")
        res <- rep(NA_real_, nrow(coords_list[[nm]])) |> 
          setNames(nm = coords_list[[nm]]$barcode)
      } else {
        res <- RadialDistance(object = coords_list[[nm]],
                             spots = spots_list[[lbl]][[nm]],
                             verbose = verbose,
                             angles = angles,
                             angles_nbreaks = angles_nbreaks,
                             remove_singletons = remove_singletons,
                             convert_to_microns = convert_to_microns,
          ...)
      }
      if (inherits(res, what = "numeric")) {
        res <- tibble(barcode = names(res), res) |>
          setNames(nm = c("barcode", paste0("r_dist_", lbl)))
      }
      if (ncol(res) > 2) {
        res <- res |>
          select(barcode, angle, r_dist, contains("intervals"))
        colnames(res) <- c("barcode", paste0(colnames(res)[2:ncol(res)], "_", lbl))
      }
      return(res)
    })
    sample_radial_distances <- Reduce(\(x, y) left_join(x, y, by = "barcode"), sample_radial_distances)
    return(sample_radial_distances)
  }))

  # Return data to Seurat object
  object@meta.data <- object@meta.data |>
    select(-contains(colnames(radial_distances))) |>
    bind_cols(radial_distances[match(colnames(object), radial_distances$barcode), ] |> select(-barcode))

  return(object)
}


#' Check that selected angles are valid
#'
#' @param spatnet A list of tibbles containing spatial networks generated with
#' \code{\link{GetSpatialNetwork}}
#' @param coords A tibble with spot coordinates
#' @param spots A character vector with spot IDs
#' @param angles A numeric vector of length 2 defining the search area
#' @param angles_nbreaks Number of intervals to split search area into
#' @param verbose Print messages
#'
#' @import dplyr
#' @import cli
#' @importFrom rlang warn abort %||%
#' @importFrom stats median
#'
#' @noRd
.check_angles <- function (
  spatnet,
  coords,
  spots,
  angles,
  angles_nbreaks,
  verbose
) {

  # Set global variables to NULL
  x <- y <- sampleID <- dist <- NULL

  # Check if spatnet is connected
  if (.is_disconnected(spatnet)) {
    abort("Detected diconnected components. Cannot compute radial distances for fixed angles when there are multiple diconnected components present in data.")
  }

  # Check angles
  if (!is.null(angles)) {
    stopifnot(is.numeric(angles),
              length(angles) == 2,
              between(x = angles, left = 0, right = 360),
              angles[2] > angles[1])
    # Check that graph is connected
    if (.is_disconnected(spatnet))
      abort("Dividing radial distances into angular intervals only works for regions with one center.")
    # Find centroid for selected spots
    centroids <- coords[match(spots, coords$barcode), ] |>
      group_by(sampleID) |>
      group_split() |>
      lapply(function(x) {
        x |>  summarize(x = median(x),
                        y = median(y)) |>
          as.numeric()
      }) |>
      setNames(names(spatnet))
    # Check that centroids are within selected region
    check_centroids <- sapply(names(centroids), function(nm) {
      sample_coords <- coords |>
        filter(sampleID == nm) |>
        mutate(dist = sqrt(rowSums((centroids[[nm]] - cbind(x, y))^2))) |>
        slice_min(order_by = dist)
      return(sample_coords$barcode %in% spots)
    })
    if (!all(check_centroids))
      warn(glue("Center outside selected region. Make sure that the region center is valid."))
    if (verbose) cli_alert_info("Setting search area between {angles[1]} and {angles[2]} degrees from region center")
    if (!is.null(angles_nbreaks)) {
      if (verbose) cli_alert_info("Splitting search area into {angles_nbreaks} interval(s)")
    }
  }

  return(list(angles, centroids))
}
