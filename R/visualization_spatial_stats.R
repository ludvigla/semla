#' @include checks.R
#'
NULL

# TODO: make compatible with scale?
# TODO: Add functionality to plot multiple features (not prio)
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
#' @family spatial-visualization-methods
#'
#' @author Lovisa Franzen
#'
#' @return a `patchwork` object
#'
#' @examples
#' # Prepare Seurat object
#' se <- readRDS(system.file("extdata/mousebrain",
#'                           "se_mbrain",
#'                            package = "semla"))
#'
#' se <- se |>
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
  slot = "data",
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

  # Set global variables to NULL
  sampleID <- NULL

  # Check Seurat object
  .check_seurat_object(object)

  # Check subplot_type
  subplot_type_options <- c("box", "violin", "histogram", "density")
  subplot_type_options_coll <- glue::glue_collapse(subplot_type_options, sep = ', ', last = ' or ')
  if (missing(subplot_type)) {abort(glue("No subplot type specified. Please provide either {subplot_type_options_coll}."))}
  if (!missing(subplot_type)) {
    if (!is.character(subplot_type)) abort(glue("Invalid class '{class(subplot_type)}' of subplot_type."))
    if (length(subplot_type) > 1) abort(glue("Only 1 subplot type can be provided at the time."))
    if(!subplot_type %in% subplot_type_options) abort(glue("{subplot_type} is not a valid choice.  Please provide either {subplot_type_options_coll}"))
  }

  # Check features
  if (missing(features)) {abort(glue("No feature specified. Please provide one feature to plot"))}
  if (!missing(features)) {
    if (!is.character(features)) abort(glue("Invalid class '{class(features)}' of the provided feature"))
    if (length(features) != 1) abort(glue("Only 1 feature can be provided at the time."))
  }

  # Check color - pick a fill color based on mid scale color
  colors_scale <- colors %||% RColorBrewer::brewer.pal(8, "Reds")
  fill_color <- fill_color %||% colors_scale[round(length(colors_scale)/2)]
  if (!is.character(fill_color)) abort(glue("Invalid class '{class(fill_color)}' for 'fill_color', expected a 'character'"))
  if (!length(fill_color) == 1) abort(glue("Only 1 fill color can be provided"))

  # Get data (code copied from MapFeatures())
  data_use <- GetStaffli(object)@meta_data |>
    bind_cols(FetchData(object, vars = features, slot = slot) |> as_tibble())

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


# TODO: Deal with split_labels=T
#' Map features spatially and add a summary plot next to it
#'
#' This function is a wrapped for \code{\link{MapLabels}} which adds a
#' stacked bar plot showing the sample's proportions of each category in
#' the selected column next to the spatial feature plot.
#'
#' Note that currently, only 1 label can be selected
#'
#' @param hide_legend logical specifying whether to hide the label legend for the spatial plot. Set to (\code{TRUE}) by default.
#' @param bar_display a character vector of length 1 specifying one of "percent" or "count" for the bar plot to display. Default set to "percent".
#' @param bar_width a numeric value specifying width of the bar plot. Default set to 1.2.
#' @param bar_label_size a numeric value specifying text size of the bar plot labels. Default set to 3.
#' @inheritParams MapLabels
#'
#' @importFrom patchwork plot_layout wrap_plots area plot_spacer
#' @import rlang
#' @import dplyr
#' @import glue
#'
#' @family spatial-visualization-methods
#'
#' @author Lovisa Franzen
#'
#' @return a `patchwork` object
#'
#' @examples
#' # Prepare Seurat object
#' se <- readRDS(system.file("extdata/mousebrain",
#'                           "se_mbrain",
#'                            package = "semla"))
#'
#' se <- se |>
#'   NormalizeData()  |>
#'   ScaleData() |>
#'   FindVariableFeatures() |>
#'   RunPCA() |>
#'   FindNeighbors(reduction = "pca", dims = 1:10) |>
#'   FindClusters(resolution = 0.2)
#'
#' # Plot clusters
#' MapLabelsSummary(se, column_name = "seurat_clusters", override_plot_dims = TRUE)
#'
#'
#' @export
#'
MapLabelsSummary <- function (
    object,
    column_name,
    bar_display = "percent", # c("percent", "count"),
    bar_width = 1.2,
    bar_label_size = 3,
    image_use = NULL,
    coords_use = "raw",
    crop_area = NULL,
    pt_size = 1,
    pt_alpha = 1,
    pt_stroke = 0,
    hide_legend = TRUE,
    section_number = NULL,
    label_by = NULL,
    # split_labels = FALSE,
    ncol = NULL,
    colors = NULL,
    override_plot_dims = FALSE,
    return_plot_list = FALSE,
    drop_na = FALSE,
    ...
) {

  # Set global variables to NULL
  sampleID <- pos <- pct <- pct_round <- NULL

  # Check Seurat object
  .check_seurat_object(object)

  # Check bar_display arg
  bar_display <- match.arg(bar_display, choices = c("percent", "count"), several.ok = FALSE)

  # fetch data from Seurat object
  data_use <- GetStaffli(object)@meta_data |>
    bind_cols(FetchData(object, vars = column_name) |> as_tibble())

  # Convert label column to factor
  data_use <- data_use |>
    mutate(across(all_of(column_name), ~ factor(.x)))

  # Add label_by column if present
  if (!is.null(label_by)) {
    data_use <- data_use |>
      bind_cols(FetchData(object, vars = label_by))
  }
  label_levels <- levels(object[[column_name]][,1])

  # Split data by sampleID (code copied from MapLabels())
  data <- data_use |>
    group_by(sampleID) |>
    group_split() |>
    setNames(nm = unique(data_use$sampleID))

  # Check color, or create a new color palette for levels
  color_labels <- colors %||% .gg_color_hue(n = length(label_levels))
  if (length(color_labels) < length(label_levels)){
    rlang::warn("Too few colors provided. Picking new default colors.")
    color_labels <- .gg_color_hue(n = length(label_levels))
  } else if (length(color_labels) > length(label_levels)) {
    color_labels <- color_labels[1:length(label_levels)]
  }
  if (is.null(colors) | is.null(names(colors))) {  # Add color names if missing
    names(color_labels) <- label_levels
  }

  # Plot MapLabels and return list
  p_list <- MapLabels(object = object,
                      column_name = column_name,
                      colors = color_labels,
                      image_use = image_use,
                      coords_use = coords_use,
                      crop_area = crop_area,
                      pt_size = pt_size,
                      pt_alpha = pt_alpha,
                      pt_stroke = pt_stroke,
                      section_number = section_number,
                      label_by = label_by,
                      ncol = ncol,
                      split_labels = FALSE,  # ! Important
                      override_plot_dims = override_plot_dims,
                      return_plot_list = TRUE  # ! Important
  )

  # Add names to list
  list_names <- (section_number %||% 1:length(p_list)) |> paste0()
  p_list <- setNames(p_list, nm = list_names)

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

  # Subset data if section_number is set
  if (!is.null(section_number)) {
    data <- data[section_number]
  }

  # Group count props per cluster in each sample and make bar plots
  map_labs_bar <- lapply(data, function(x){
    names(x)[names(x) == column_name] <- "labels"
    x <- x |>
      dplyr::group_by(labels) |>
      dplyr::summarise(n = n()) |>
      dplyr::mutate(pct = round(100 * n/sum(n), 5),
                    pct_round = round(100 * n/sum(n), 1),
                    pos = 100-(cumsum(pct) - (0.5 * pct)))
    x <- x[x$pct>0,] # Remove 0% groups from bar - Make optional?

    p_bar <- ggplot(x, aes(x=1, y=pct, fill=labels)) +
      geom_bar(
        stat = "identity",
        position = "stack",
        width = bar_width
      ) +
      scale_fill_manual(values = color_labels) +
      coord_cartesian(
        ylim = c(0,100),
        xlim = c(0, 5),
        clip = 'off') +
      theme_void() +
      theme(legend.position = "none",
            plot.margin = unit(c(0, 50, 0, -20), "pt"))
  })
  map_labs_bar <- setNames(map_labs_bar, nm = list_names)

  p_patch <- lapply(list_names, function(i){
    p <- p_list[[i]]
    if(hide_legend){
      p <- p & theme(legend.position = "none")
    }

    legend_spacer <- "     "
    if(bar_width>=1.5) {
      legend_spacer <- paste(legend_spacer, rep(" ", round(bar_width/2)), collapse = "")
    }

    if (bar_display == "count") {
      p_bar <- map_labs_bar[[i]] &
        geom_text(
          aes(label = paste0(legend_spacer, labels, ": ", n),
              y = pos,
              vjust = "inward",
              hjust = 0
          ),
          size = bar_label_size)
    } else {
      p_bar <- map_labs_bar[[i]] &
        geom_text(
          aes(label = paste0(legend_spacer, labels, ": ", pct_round, "%"),
              y = pos,
              vjust = "inward",
              hjust = 0
          ),
          size = bar_label_size)
    }

    design <- c(area(t = 1, l = 1, b = 6, r = 5),
                area(t = 1, l = 6, b = 1, r = 6),
                area(t = 2, l = 6, b = 5, r = 6),
                area(t = 6, l = 6, b = 6, r = 6))
    p_out <- p +
      plot_spacer() +
      p_bar +
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
          axis.line.y.right = element_line(linewidth = 0.25),
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
#' @importFrom ggplot2 geom_boxplot xlim scale_y_continuous
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
    geom_boxplot(data = x, mapping = aes(x = 1, y = .data[[features]]),
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
#' @importFrom ggplot2 geom_violin xlim scale_y_continuous
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
    geom_violin(data = x, mapping = aes(x = 1, y = .data[[features]]),
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
#' @importFrom ggplot2 geom_histogram coord_flip scale_x_continuous scale_y_reverse
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
    geom_histogram(data = x, mapping = aes(x = .data[[features]]),
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
#' @importFrom ggplot2 geom_density coord_flip scale_x_continuous scale_y_reverse
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
    geom_density(data = x, mapping = aes(x = .data[[features]]),
                 color = fill_color,
                 fill = fill_color,
                 alpha = 0.6) +
    coord_flip() +
    scale_y_reverse() +
    scale_x_continuous(position = "top")
  return(p_stat)
}
