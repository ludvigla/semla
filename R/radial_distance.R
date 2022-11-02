#' @include generics.R
#' @include checks.R
#' @include spatial_utils.R
#'
NULL

#' @description Calculates the radial distances to surrounding spots from a selected
#' group of spots covering a defined regions. The region could for example
#' represent an isolated tumor in the tissue section surrounded by stroma. First,
#' the border of the selected region is defined and for each spot outside
#' of this border, the distance is calculated to its nearest border spot.
#' Spots located inside the selected region will have negative distances and
#' spots located outside of the selected region will have positive distances.
#' Having access to the radial distances can be useful when inspecting changes
#' in gene expression as a function of distance to a region of interest.
#'
#' @param spots A character vector with spot IDs present `object`
#'
#' @import dplyr
#' @importFrom rlang abort inform
#' @importFrom cli cli_h3
#' @importFrom glue glue
#' @importFrom dbscan kNN
#'
#' @rdname radial-distance
#'
#' @return A numeric vector with radial distances. If the input object is of class
#' `Seurat`, the radial distances will be returned in the meta data slot.
#'
#'
#' @examples
#'
#' \dontrun{
#' library(STUtility2)
#' library(ggplot2)
#' library(patchwork)
#'
#' # Get coordinates
#' galt_spots_file <- "~/STUtility2/repo/STUtility2/inst/extdata/mousecolon/galt_spots.csv"
#' galt_spots <- read.csv(galt_spots_file) |>
#'   as_tibble()
#'
#' # read coordinates
#' coordfile <- system.file("extdata/mousecolon/spatial",
#'                          "tissue_positions_list.csv",
#'                          package = "STUtility2")
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
#' gg <- gg |>
#'   mutate(r_dist_sqrt = case_when(r_dist < 0 ~ -sqrt(abs(r_dist)),
#'                                  r_dist >= 0 ~ sqrt(r_dist)))
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
#' }
#'
#' @export
#'
RadialDistance.default <- function (
  object,
  spots,
  verbose = TRUE,
  ...
) {

  # Set global variables to NULL
  barcode <- x <- y <- sampleID <- NULL

  # Check object object class
  if (!any(class(object) %in% c("data.frame", "matrix", "tbl")))
    abort(glue("Invalid class '{class(object)}'."))
  if (ncol(object) != 4)
    abort(glue("Invalid number of columns '{ncol(object)}'. Expected 4."))
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

  if (verbose) inform(c("i" = glue("Extracting border spots from a region with {length(spots)} spots")))

  # Find border
  spatnet <- GetSpatialNetwork(object)
  border_spots <- RegionNeighbors(spatnet, spots = spots, outer_border = FALSE)
  inside_spots <- setdiff(spots, border_spots)

  # Get indices for spot groups
  border_spots_indices <- match(border_spots, object$barcode)
  outside_spots_indices <- -match(spots, object$barcode)
  inside_spots_indices <- match(inside_spots, object$barcode)
  if (verbose) inform(c(">" = glue("Detected {length(border_spots_indices)} spots on borders"),
                        ">" = glue("Detected {length(inside_spots_indices)} spots inside borders"),
                        ">" = glue("Detected {nrow(object) - length(spots)} spots outside borders")))

  # Find nearest neighbor at selected spots border for
  # spot outside and inside the border
  knn_spatial_outside <- kNN(x = object[border_spots_indices, c("x", "y")] |> as.matrix(),
                             query = object[outside_spots_indices, c("x", "y")] |> as.matrix(),
                             k = 1)
  knn_spatial_inside <- kNN(x = object[border_spots_indices, c("x", "y")] |> as.matrix(),
                            query = object[inside_spots_indices, c("x", "y")] |> as.matrix(),
                            k = 1)
  if (verbose) inform(c("v" = glue("Returning radial distances")))

  # Get radial distances
  radial_dists <- setNames(numeric(length = nrow(object)), nm = object$barcode)
  radial_dists_outside <- setNames(knn_spatial_outside$dist[, 1, drop = TRUE], nm = object$barcode[outside_spots_indices])
  radial_dists_inside <- setNames(knn_spatial_inside$dist[, 1, drop = TRUE], nm = object$barcode[inside_spots_indices])
  radial_dists[names(radial_dists_outside)] <- radial_dists_outside
  radial_dists[names(radial_dists_inside)] <- -radial_dists_inside

  # Return radial distances
  return(radial_dists)

}


#' @param column_name A character specifying the name of a column in your meta data that contains
#'  categorical data, e.g. clusters or manual selections
#' @param sel_groups A character vector to select specific groups in \code{column_name} with.
#' All groups are selected by default, but the common use case is to select a region of interest.
#' @param verbose Print messages
#'
#' @import dplyr
#' @importFrom rlang inform abort
#' @importFrom glue glue
#'
#' @rdname radial-distance
#'
#' @examples
#'
#' library(STUtility2)
#' library(ggplot2)
#' library(patchwork)
#' library(tidyr)
#'
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "STUtility2"))
#' se_mcolon <- RadialDistance(se_mcolon, column_name = "selection", sel_groups = "GALT")
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
#' MapFeatures(se_mcolon, features = "r_dist_GALT_sqrt", pt_size = 3,
#'             colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev())
#'
#' @export
#'
RadialDistance.Seurat <- function (
    object,
    column_name,
    sel_groups = NULL,
    verbose = TRUE,
    ...
) {

  # Set global variables to NULL
  barcode <- pxl_col_in_fullres <- pxl_row_in_fullres <- sampleID <- NULL

  # validate input
  .check_seurat_object(object)
  .validate_column_name(object, column_name)
  sel_groups <- .validate_selected_labels(object, sel_groups, column_name)

  # Select spots
  spots_list <- .get_spots_list(object, sel_groups, column_name)

  # Get coordinates
  coords_list <- .get_coords_list(object)

  # Calculate radial distances
  distances <- do.call(bind_rows, lapply(names(coords_list), function(nm) {
    if (verbose) inform(glue("Running calculations for sample {nm}"))
    radial_distances <- do.call(bind_cols, lapply(names(spots_list), function(lbl) {
      if (verbose) inform(glue("Calculating radial distances for group '{lbl}'"))
      tibble(RadialDistance(coords_list[[nm]], spots = spots_list[[lbl]][[nm]], verbose = verbose, ...)) |>
        setNames(nm = paste0("r_dist_", lbl))
    }))
  }))

  # Return data to Seurat object
  object@meta.data <- object@meta.data |>
    select(-contains(colnames(distances))) |>
    bind_cols(distances)

  return(object)
}
