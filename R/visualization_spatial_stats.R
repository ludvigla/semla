#' @include checks.R
#'
NULL


#' Map features spatially + summary statistic
#'
#' @param fill_color fill color for stats geom
#'
#' @importFrom patchwork plot_layout wrap_plots
#' @import dplyr
#'
#' @return a `patchwork` object
#'
#' @examples
#'
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
  arrange_features = c("col", "row"),
  #override_plot_dims = FALSE,
  max_cutoff = NULL,
  min_cutoff = NULL,
  ...
) {
  # TODO: How to deal with multiple features?

  # Check Seurat object
  .check_seurat_object(object)

  # validate subplot_type
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
    # Currently only one feature is possible
    # override_plot_dims = F,
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
    arrange_features = arrange_features,
    override_plot_dims = override_plot_dims,
    max_cutoff = max_cutoff,
    min_cutoff = min_cutoff,
    return_plot_list = TRUE # !Important
  )

  # Set up base plot
  p_stat_base <- .plot_base()

  # fkns
  plot_stat_fkn <- switch(subplot_type,
                "box" = get(".plot_box"),
                "violin" = get(".plot_violin"),
                "histogram" = get(".plot_histogram"),
                "density" = get(".plot_density"))

  map_feat_stats <- lapply(data, function(x){
    plot_stat_fkn(p_stat_base, x, features)
  })

  # Patchwork plots together
  p_patch <- lapply(1:length(p_list), function(i){
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

# TODO: document internal functions for RMD check

.plot_box <- function (
  p_stat_base,
  x,
  features
) {
  p_stat <- p_stat_base +
    geom_boxplot(data = x, mapping = aes_string(x = "1", y = features),
                 outlier.size = 0.5,
                 fill = fill_color) +
    xlim(0.5, 1.75) +
    scale_y_continuous(position = "right")
  return(p_stat)
}

.plot_violin <- function (
    p_stat_base,
    x,
    features
) {
  p_stat <- p_stat_base +
    geom_violin(data = x, mapping = aes_string(x = "1", y = features),
                fill = fill_color,
                draw_quantiles = 0.5) +
    xlim(0.25,1.75) +
    scale_y_continuous(position = "right")
  return(p_stat)
}

.plot_histogram <- function (
    p_stat_base,
    x,
    features
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

.plot_density <- function (
    p_stat_base,
    x,
    features
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
