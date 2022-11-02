#' @include checks.R
#'
NULL

# TODO: make compatible with scale?
#' Map features spatially and add a summary plot next to it
#'
#' This function is a wrapped for \code{\link{MapFeatures}} which allows you to
#' add a boxplot, histogram, violin plot or a density histogram showing the
#' distribution of the selected feature next to the spatial feature plot.
#'
#' Note that currently, only 1 feature can be selected
#'
#' @param fill_color Fill color for summary plot
#' @param subplot_type Select a summary plot to place next to the spatial plot:
#' \itemize{
#'    \item{"box" : boxplot}
#'    \item{"violin" : violin plot}
#'    \item{"histogram" : histogram}
#'    \item{"density" : density histogram}
#' }
#' @inheritParams MapFeatures
#'
#' @importFrom patchwork plot_layout wrap_plots area plot_spacer
#' @import rlang
#' @import dplyr
#' @import glue
#'
#' @family visualization
#'
#' @author Lovisa Franzen
#'
#' @return a `patchwork` object
#'
#' @examples
#'
#' library(tibble)
#'
#' samples <- Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"),
#'                                 "/*/filtered_feature_bc_matrix.h5"))
#' imgs <- Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"),
#'                                 "/*/spatial/tissue_hires_image.png"))
#' spotfiles <- Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"),
#'                                      "/*/spatial/tissue_positions_list.csv"))
#' json <- Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"),
#'                                 "/*/spatial/scalefactors_json.json"))
#'
#' # Create a tibble/data.frame with file paths
#' infoTable <- tibble(samples, imgs, spotfiles, json, sample_id = c("mousebrain", "mousecolon"))
#'
#' # Create Seurat object
#' se <- ReadVisiumData(infoTable = infoTable) |>
#'   NormalizeData()
#'
#' # Add boxplot
#' MapFeaturesSummary(se, features = "Nrgn", subplot_type = "box")
#'
#' # Add violin plot
#' MapFeaturesSummary(se, features = "Nrgn", subplot_type = "violin")
#'
#' # Add histogram
#' MapFeaturesSummary(se, features = "Nrgn", subplot_type = "histogram")
#'
#' # Add density histogram
#' MapFeaturesSummary(se, features = "Nrgn", subplot_type = "density")
#'
#' @export
#'
MapFeaturesSummary <- function (
  object,
  features,
  subplot_type = c("box", "violin", "histogram", "density"),
  image_use = NULL,
  coords_use = "raw",
  crop_area = NULL,
  pt_size = 1,
  pt_alpha = 1,
  pt_stroke = 0,
  scale_alpha = FALSE,
  section_number = NULL,
  label_by = NULL,
  ncol = NULL,
  colors = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
  fill_color = NULL,
  scale = c("shared", "free"),
  override_plot_dims = FALSE,
  max_cutoff = NULL,
  min_cutoff = NULL,
  ...
) {
  # TODO: How to deal with multiple features?

  # Set global variables to NULL
  sampleID <- NULL

  # Check Seurat object
  .check_seurat_object(object)

  # validate subplot_type and fill_color
  stopifnot(
    is.character(subplot_type),
    length(subplot_type) == 1
  )
  subplot_type <- match.arg(subplot_type, choices = c("box", "violin", "histogram", "density"))

  # Check features
  stopifnot(
    is.character(features),
    length(features) == 1
  )

  # Check color - pick a fill color based on mid scale color
  colors_scale <- colors %||% RColorBrewer::brewer.pal(8, "Reds")
  fill_color <- fill_color %||% colors_scale[round(length(colors_scale)/2)]
  if (!is.character(fill_color)) abort(glue("Invalid class '{class(fill_color)}' for 'fill_color', expected a 'character'"))
  if (!length(fill_color) == 1) abort(glue("Only 1 fill color can be provided"))

  # Get data (code copied from MapFeatures())
  data_use <- GetStaffli(object)@meta_data |>
    bind_cols(FetchData(object, vars = features) |> as_tibble())

  # Split data by sampleID (code copied from MapFeatures())
  data <- data_use |>
    group_by(sampleID) |>
    group_split() |>
    setNames(nm = unique(data_use$sampleID))

  # Plot MapFeatures
  p_list <- MapFeatures(
    object,
    features = features,
    image_use = image_use,
    coords_use = coords_use,
    crop_area = crop_area,
    pt_size = pt_size,
    pt_alpha = pt_alpha,
    pt_stroke = pt_stroke,
    scale_alpha = scale_alpha,
    section_number = section_number,
    label_by = label_by,
    ncol = ncol,
    scale = scale,
    colors = colors_scale,
    override_plot_dims = override_plot_dims,
    max_cutoff = max_cutoff,
    min_cutoff = min_cutoff,
    return_plot_list = TRUE
  )

  # Add names to list
  list_names <- (section_number %||% 1:length(p_list)) |> paste0()
  p_list <- setNames(p_list, nm = list_names)

  # Set up base plot
  p_stat_base <- .plot_base()

  # fkns
  plot_stat_fkn <- switch(subplot_type,
                "box" = get(".plot_box"),
                "violin" = get(".plot_violin"),
                "histogram" = get(".plot_histogram"),
                "density" = get(".plot_density"))

  # Crop data if crop_area is set
  if (!is.null(crop_area)) {
    if (!is.numeric(crop_area)) abort(glue("Invalid class '{class(crop_area)}' for 'crop_area', expected 'numeric'"))
    if (length(crop_area) != 4) abort(glue("Invalid length for 'crop_area', expected a 'numeric' vector of length 4"))
    if (!all(between(x = crop_area, left = 0, right = 1))) abort("'crop_area' can only take values between 0-1")
    # Get coords columns
    coords_columns <- .get_coords_column(image_use, coords_use)
    # Get image dimensions
    dims <- GetStaffli(object)@image_info
    # Crop data
    data <- setNames(lapply(names(data), function(nm) {
      crop_limits <- c(round(crop_area[1]*dims[dims$sampleID == nm, "full_width", drop = TRUE]),
                       round(crop_area[2]*dims[dims$sampleID == nm, "full_height", drop = TRUE]),
                       round(crop_area[3]*dims[dims$sampleID == nm, "full_width", drop = TRUE]),
                       round(crop_area[4]*dims[dims$sampleID == nm, "full_height", drop = TRUE]))
      x <- data[[nm]] |>
        filter(between(x = !! sym(coords_columns[1]), left = crop_limits[1], right = crop_limits[3])) |>
        filter(between(x = !! sym(coords_columns[2]), left = crop_limits[2], right = crop_limits[4]))
      return(x)
    }), nm = names(data))
  }


  # Create stats plots
  map_feat_stats <- lapply(data, function(x){
    plot_stat_fkn(p_stat_base, x, features, fill_color)
  })

  # Patchwork plots together
  p_patch <- lapply(list_names, function(i){
    p <- p_list[[i]][[features]] #& theme(panel.background = element_rect(fill="grey90"))
    p_stat <- map_feat_stats[[i]]
    design <- c(area(t = 1, l = 1, b = 6, r = 5),
                area(t = 1, l = 6, b = 1, r = 6),
                area(t = 2, l = 6, b = 5, r = 6),
                area(t = 6, l = 6, b = 6, r = 6))
    p <- p +
      plot_spacer() +
      p_stat +
      plot_spacer() + plot_layout(design = design)
  })

  # Return wrapped plot grid
  ncol <- ncol %||% ceiling(sqrt(length(p_patch)))
  p_out <- wrap_plots(p_patch, ncol = ncol)

  return(p_out)
}


#' Create an empty plot to add geoms to
#'
#' @importFrom ggplot2 ggplot theme_linedraw theme element_text element_blank element_line unit
#'
#' @return a `ggplot` object
#'
#' @noRd
.plot_base <- function (
) {
  p_stat_base <- ggplot() +
    theme_linedraw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.y.right = element_line(size = 0.25),
          plot.margin = unit(c(0, 40, 0, 0), "pt"))
  return(p_stat_base)
}


#' Add a boxplot
#'
#' @param p_stat_base An empty `ggplot`
#' @param x Data to use for plot
#' @param features Selected features
#' @param fill_color Fill color for boxplot
#'
#' @importFrom ggplot2 geom_boxplot aes_string xlim scale_y_continuous
#'
#' @return a `ggplot` object
#'
#' @noRd
.plot_box <- function (
  p_stat_base,
  x,
  features,
  fill_color
) {
  p_stat <- p_stat_base +
    geom_boxplot(data = x, mapping = aes_string(x = "1", y = features),
                 outlier.size = 0.5,
                 fill = fill_color) +
    xlim(0.5, 1.75) +
    scale_y_continuous(position = "right")
  return(p_stat)
}


#' Add a violin plot
#'
#' @param p_stat_base An empty `ggplot`
#' @param x Data to use for plot
#' @param features Selected features
#' @param fill_color Fill color for violin plot
#'
#' @importFrom ggplot2 geom_violin aes_string xlim scale_y_continuous
#'
#' @return a `ggplot` object
#'
#' @noRd
.plot_violin <- function (
    p_stat_base,
    x,
    features,
    fill_color
) {
  p_stat <- p_stat_base +
    geom_violin(data = x, mapping = aes_string(x = "1", y = features),
                fill = fill_color,
                draw_quantiles = 0.5) +
    xlim(0.25,1.75) +
    scale_y_continuous(position = "right")
  return(p_stat)
}


#' Add a histogram
#'
#' @param p_stat_base An empty `ggplot`
#' @param x Data to use for plot
#' @param features Selected features
#' @param fill_color Fill color for histogram
#'
#' @importFrom ggplot2 geom_histogram aes_string coord_flip scale_x_continuous scale_y_reverse
#'
#' @return a `ggplot` object
#'
#' @noRd
.plot_histogram <- function (
    p_stat_base,
    x,
    features,
    fill_color
) {
  p_stat <- p_stat_base +
    geom_histogram(data = x, mapping = aes_string(x = features),
                 bins = 50,
                 fill = fill_color) +
    coord_flip() +
    scale_y_reverse() +
    scale_x_continuous(position = "top")
  return(p_stat)
}

#' Add a density histogram
#'
#' @param p_stat_base An empty `ggplot`
#' @param x Data to use for plot
#' @param features Selected features
#' @param fill_color Fill color for density histogram
#'
#' @importFrom ggplot2 geom_density aes_string coord_flip scale_x_continuous scale_y_reverse
#'
#' @return a `ggplot` object
#'
#' @noRd
.plot_density <- function (
    p_stat_base,
    x,
    features,
    fill_color
) {
  p_stat <- p_stat_base +
    geom_density(data = x, mapping = aes_string(x = features),
                 color = fill_color,
                 fill = fill_color,
                 alpha = 0.6) +
    coord_flip() +
    scale_y_reverse() +
    scale_x_continuous(position = "top")
  return(p_stat)
}
