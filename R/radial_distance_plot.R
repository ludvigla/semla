#' @include generics.R
#' @include checks.R
#'
NULL

#' Draw a centroid angle plot
#'
#' @param nbreaks Number of intervals to cut the angles into
#' @param centroid_size Size of the centroid spot
#'
#' @import ggplot2
#' 
#' @returns A \code{ggplot} object
#'
#' @examples
#' # Draw a plot
#' centroid_angles_plot(9)
#'
#' @export
#'
centroid_angles_plot <- function (
  nbreaks = 9,
  centroid_size = 8
) {

  # Set global variables to NULL
  x <- y <- x_lab <- y_lab <- angle <- NULL

  # Plot selection, center and directions
  arrows <- tibble(x = cos(seq(0, 2*pi, length.out = nbreaks)),
                   y = sin(seq(0, 2*pi, length.out = nbreaks)),
                   x_lab = 1.2*cos(seq(0, 2*pi, length.out = nbreaks)),
                   y_lab = 1.2*sin(seq(0, 2*pi, length.out = nbreaks)),
                   angle = round(seq(0, 360, length.out = nbreaks)))
  p <- ggplot() +
    scale_y_reverse() +
    coord_fixed() +
    annotate("path", lty = "dashed",
             x = cos(seq(0, 2*pi, length.out = 100)),
             y = sin(seq(0, 2*pi, length.out = 100))) +
    geom_segment(data = arrows, aes(x = 0, xend = x, y = 0, yend = y),
                 arrow = arrow(type = "closed", angle = 20, length = unit(0.05, "npc")), size = 0.5) +
    geom_text(data = arrows[-nbreaks, ], aes(x_lab, y_lab, label = paste0(angle, "\u00b0")),
              size = 4, color = "black", fontface = "bold") +
    geom_point(aes(0, 0, fill = "center"),
               size = centroid_size, shape = 21, fill = "lightgray") +
    theme_void()

  return(p)

}

#' Angle plot
#'
#' Draws on angle plot on top of a selected region. The plot is meant to help
#' defining angle intervals for computing radial distances with \code{\link{RadialDistance}}.
#'
#' @details
#' The region of interest is selected with \code{selected_group} and has to be a group
#' of a categorical variable selected with \code{column_name} which is stored in the
#' meta data of the input Seurat object. The selected region has to be spatially connected.
#'
#' @param object An object of class \code{Seurat}
#' @param selected_group A label defining a group of spots found in a column of
#' the meta data slot specified by \code{column_name}
#' @param radius A numeric value between 0.1 and 1 specifying the size of the
#' "angle plot" overlaid on the spatial plot.
#'
#' @inheritParams MapLabels
#' @inheritParams centroid_angles_plot
#'
#' @family spatial-visualization-methods
#'
#' @author Ludvig Larsson
#'
#' @import dplyr
#' @importFrom Seurat AddMetaData
#'
#' @examples
#' library(semla)
#'
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
#' se_mcolon <- DisconnectRegions(se_mcolon, column_name = "selection", selected_groups = "GALT")
#' se_mcolon <- LoadImages(se_mcolon)
#' AnglePlot(se_mcolon, column_name = "GALT_split", selected_group = "S1_region1", pt_size = 2,
#'         image_use = "raw", crop_area = c(0.4, 0.5, 0.7, 0.8), radius = 0.4, nbreaks = 12)
#'
#' @return An object of class \code{patchwork}
#'
#' @export
#'
AnglePlot <- function (
  object,
  column_name,
  selected_group,
  radius = 0.3,
  nbreaks = 9,
  centroid_size = 8,
  image_use = NULL,
  coords_use = "raw",
  crop_area = NULL,
  pt_size = 1,
  pt_alpha = 1,
  pt_stroke = 0,
  section_number = NULL,
  label_by = NULL,
  ncol = NULL,
  colors = NULL,
  override_plot_dims = FALSE,
  drop_na = FALSE,
  ...
) {

  # Set global variables to NULL
  from <- to <- barcode <- sampleID <- cx <- cy <- NULL

  # validate input
  .check_seurat_object(object)
  .validate_column_name(object, column_name)
  if (is.null(selected_group)) abort("'selected_group' cannot be NULL")
  select_groups <- .validate_selected_labels(object, selected_group, column_name)
  stopifnot(inherits(radius, what = "numeric"), length(radius) == 1)
  if (!between(x = radius, left = 0.1, right = 1))
    abort(glue("'radius' must take a value between 0.1 and 1, got {radius}"))

  # Get spots
  spots <- FetchData(object, vars = column_name) |>
    filter(!! sym(column_name) %in% select_groups) |>
    rownames()

  # Check if spatnet is connected
  spatnet <- GetSpatialNetwork(object)
  spatnet_region <- lapply(spatnet, function(x) {
    x |>
      filter(from %in% spots, to %in% spots)
  })
  if (.is_disconnected(spatnet_region)) {
    abort("Detected diconnected components. Cannot compute radial distances for fixed angles when there are multiple diconnected components present in data.")
  }

  # Create a new variable with selected group
  tmp_data <- object[[]] |>
    pull(all_of(column_name))
  tmp_data <- ifelse(tmp_data %in% selected_group, tmp_data, NA)
  object <- AddMetaData(object, metadata = tmp_data, col.name = paste0(column_name, "_"))

  # Create plots
  plots <- MapLabels(object = object,
                     column_name = paste0(column_name, "_"),
                     image_use = image_use,
                     coords_use = coords_use,
                     crop_area = crop_area,
                     pt_size = pt_size,
                     pt_alpha = pt_alpha,
                     pt_stroke = pt_stroke,
                     section_number = section_number,
                     label_by = label_by,
                     split_labels = FALSE,
                     ncol = NULL,
                     colors = colors,
                     override_plot_dims = override_plot_dims,
                     return_plot_list = TRUE,
                     drop_na = drop_na,
                     ... = ...)

  # Prep data
  coords_columns <- .get_coords_column(image_use, coords_use)
  coords <- GetStaffli(object)@meta_data |>
    bind_cols(object[[]] |> select(all_of(paste0(column_name, "_")))) |>
    na.omit() |>
    select(barcode, all_of(coords_columns), sampleID)
  image_dims <- GetStaffli(object)@image_info
  centroids <- coords |>
    group_by(sampleID) |>
    summarize(cx = median(!! sym(coords_columns[1])),
              cy = image_dims$full_height - median(!! sym(coords_columns[2]))) |>
    mutate(cxs = cx/image_dims$full_width, cys = cy/image_dims$full_height)

  # Get radius plots
  rad_plot <- centroid_angles_plot(nbreaks = nbreaks,
                                   centroid_size = centroid_size)

  # Create crop area if override_plot_dims = TRUE
  if (override_plot_dims) {
    # Split data by sampleID
    new_dims <- .get_limits(GetStaffli(object)@meta_data, coords_columns)
    crop_area <- c(min(new_dims$x_start/image_dims$full_width),
                   min(new_dims$y_start/image_dims$full_height),
                   max(new_dims$full_width/image_dims$full_width),
                   max(new_dims$full_height/image_dims$full_height))
  }

  # Add radius plots to label plots
  plots_with_radius <- lapply(names(plots), function(nm) {
    label_plot <- plots[[nm]]
    cxy <- centroids[centroids$sampleID == nm, 4:5] |> as.numeric()
    if (!is.null(crop_area)) {
      cxy[1] <- rescale(x = cxy[1], from = crop_area[c(1, 3)], to = c(0, 1))
      cxy[2] <- rescale(x = 1 - cxy[2], from = crop_area[c(4, 2)], to = c(0, 1))
    }
    label_and_rad_plot <- label_plot + inset_element(rad_plot,
                                                     left = cxy[1] - radius,
                                                     bottom = cxy[2] - radius,
                                                     right = cxy[1] + radius,
                                                     top = cxy[2] + radius)
    return(label_and_rad_plot)
  })

  # Arrange plots
  ncol <- ncol %||% ceiling(sqrt(length(plots)))
  p_wrapped <- wrap_plots(plots_with_radius, ncol = ncol)

  return(p_wrapped)

}
