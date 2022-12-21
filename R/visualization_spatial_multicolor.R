#' @include checks.R
#'
NULL

#' @param scale A string specifying a scaling mode. See details below
#' @param add_colorscale_text A Logical specifying if labels should be added to
#' the color scales. The labels will be values from 0 to 1 ("low" to "high") and
#' the meaning of these values depend on the choice of \code{scale}
#' @inheritParams MapFeatures
#'
#' @details The values of the selected features will be rescaled for the visualization
#' depending on the choice of \code{scale}. If \code{scale="free"}, each numeric feature
#' will be scaled independently across the entire dataset to have values between 0 and 1.
#' This is most useful if the selected features have different scales and you want to
#' avoid some features dominating.
#'
#' If you want to keep the relative differences between different feature, you can set
#' \code{scale="shared"}. In this case, the feature values will be scaled to have values
#' between 0 and 1 where 0 and 1 corresponds to the minimum and maximum value across
#' all feature values. This could for example be useful if you are visualizing cell type
#' proportions in which case the color intensity should reflect the relative differences.
#'
#' This type of visualization should be used with some caution. Because only 1 color
#' is selected for each spot, it is biased to focus on dominant features.
#' This makes it less meaningful to focus on the exact values and makes the interpretation
#' more qualitative. It can however be useful as a way to summarize the spatial distribution
#' of features that are expressed in different regions. It is best to map such features with
#' \code{\link{MapFeatures}} first, before attempting to summarize them in one plot.
#'
#' @importFrom patchwork wrap_plots
#' @importFrom rlang %||%
#' @import dplyr
#'
#' @rdname visualize-multiple-features
#'
#' @author Ludvig Larsson
#'
#' @return A patchwork object
#'
#' @export
#'
MapMultipleFeatures.default <- function (
    object,
    scale = c("shared", "free"),
    crop_area = NULL,
    pt_size = 1,
    pt_stroke = 0,
    section_number = NULL,
    label_by = NULL,
    ncol = NULL,
    colors = NULL,
    dims = NULL,
    coords_columns = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    return_plot_list = FALSE,
    drop_na = FALSE,
    add_colorscale_text = FALSE,
    ...
) {

  # Set global variables to NULL
  barcode <- sampleID <- NULL

  # Check input parameters
  scale <- match.arg(scale, choices = c("shared", "free"))

  # Check data
  .prep_data_for_plotting(object = object, label_by = label_by, scale = "shared", arrange_features = "col",
                          coords_columns = coords_columns, multi_color = TRUE)

  # get features
  features <- object |>
    select(-barcode, -all_of(coords_columns), -sampleID, -all_of(label_by)) |>
    colnames()

  # Check features
  if (!length(features) > 1) abort(glue("Expected at least 2 features, got {length(features)}"))

  # Expand colors if missing
  colors <- colors %||% {
    .gg_color_hue(n = length(features))
  }
  if (length(features) != length(colors)) {
    cli_alert_danger(paste0("The number of colors provided does not match the number of selected features.",
                     " Selecting {length(features)} unique colors."))
    colors <- .gg_color_hue(n = length(features))
  }

  # Check section number and subset data
  if (!is.null(section_number)) {
    if (!is.numeric(section_number)) abort(glue("Invalid class '{class(section_number)}' for",
                                                " {cli::col_br_green('section_number')}, expected 'numeric'"))
    if (!section_number %in% object$sampleID) abort(glue("'section_number' = {section_number} is out of range. ",
                                                         "Select one of {paste(unique(object$sampleID), collapse = ', ')}"))
    object <- object |>
      filter(sampleID == section_number)
  }

  # Convert feature values
  if (scale == "shared") {
    # Set range of feature values
    rangeFeatureVal <- object |>
      select(all_of(features)) |>
      range()
    d_features <- object |>
      select(all_of(features)) |>
      as.matrix() |>
      as.vector() |>
      scales::rescale(to = c(0, 1), from = rangeFeatureVal) |>
      matrix(ncol = length(features), byrow = FALSE) |>
      apply(1, function(x) {
        x[-which.max(x)] <- NA
        return(x)
      }) |> t()
  }
  if (scale == "free") {
    d_features <- object |>
      select(all_of(features)) |>
      as.matrix() |>
      apply(2, scales::rescale) |>
      apply(1, function(x) {
        x[-which.max(x)] <- NA
        return(x)
      }) |> t()
  }
  colnames(d_features) <- features
  data <- object |> select(-all_of(features)) |>
    bind_cols(d_features)

  # Split data by sampleID
  data <- data |>
    group_by(sampleID) |>
    group_split() |>
    setNames(nm = unique(object$sampleID))

  # get image dimensions
  dims <- .get_dims(dims)

  # Edit dims of a crop area is provided
  if (!is.null(crop_area)) {
    c(dims, data) %<-% .crop_dims(dims, crop_area, data, coords_columns)
  }

  # Plot features on spatial coordinates for each sample
  sample_plots <- setNames(lapply(names(data), function(nm) {

    # Get data for plotting
    gg <- data[[nm]]

    # Create an appropriate plot title
    if (!is.null(label_by)) {
      cur_label <- unique(gg |> pull(all_of(label_by)))
    } else {
      cur_label <- paste0("section ", nm)
    }

    p <- .spatial_feature_plot_multiple(
      gg = gg,
      nm = nm,
      features = features,
      colors = colors,
      dims = dims,
      pt_size = pt_size,
      pt_stroke = pt_stroke,
      cur_label = cur_label,
      coords_columns = coords_columns,
      drop_na = drop_na,
      use_text = add_colorscale_text
    )
    return(p)
  }), nm = names(data))

  # Create final patchwork
  if (!return_plot_list) {
    ncol <- ncol %||% ceiling(sqrt(length(data)))
    wrapped_plots <- wrap_plots(sample_plots, ncol = ncol)
  } else {
    # return list of ggplot objects if return_plot_list = TRUE
    wrapped_plots <- sample_plots
  }

  return(wrapped_plots)
}


#'
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom Seurat FetchData
#' @importFrom rlang warn
#'
#' @rdname visualize-multiple-features
#' @family spatial-visualization-methods
#'
#' @author Ludvig Larsson
#'
#' @examples
#'
#' library(STUtility2)
#'
#' # Load example Visium data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain",
#'                                  "se_mbrain",
#'                                  package = "STUtility2"))
#'
#' # Select features to plot
#' sel_features = c("Th", "Trh", "Calb2", "Prkcd")
#'
#' # Map selected features
#' MapFeatures(se_mbrain, features = sel_features)
#'
#' # Summarize in one plot
#' MapMultipleFeatures(se_mbrain,
#'                     features = sel_features,
#'                     pt_size = 2)
#'
#' # Plot with H&E image
#' se_mbrain <- se_mbrain |> LoadImages()
#' MapMultipleFeatures(se_mbrain, image_use = "raw",
#'                     features = sel_features,
#'                     pt_size = 2)
#'
#'
#' @export
#'
MapMultipleFeatures.Seurat <- function (
    object,
    features,
    scale = c("free", "shared"),
    slot = "data",
    image_use = NULL,
    coords_use = "raw",
    crop_area = NULL,
    pt_size = 1,
    pt_stroke = 0,
    section_number = NULL,
    label_by = NULL,
    ncol = NULL,
    colors = NULL,
    override_plot_dims = FALSE,
    max_cutoff = NULL,
    min_cutoff = NULL,
    return_plot_list = FALSE,
    add_colorscale_text = FALSE,
    ...
) {

  # Set global variables to NULL
  sampleID <- NULL

  # Match args
  scale <- match.arg(scale, choices = c("free", "shared"))
  coords_use <- match.arg(coords_use, choices = c("raw", "transformed"))

  # Check Seurat object
  .check_seurat_object(object)

  # fetch data from Seurat object
  data_use <- GetStaffli(object)@meta_data |>
    bind_cols(FetchData(object, vars = features, slot = slot) |> as_tibble())

  # Add label_by column if present
  if (!is.null(label_by)) {
    data_use <- data_use |>
      bind_cols(FetchData(object, vars = label_by))
  }

  # Subset by section number
  if (!is.null(section_number)) {
    if (!is.numeric(section_number)) abort(glue("Invalid class '{class(section_number)}' for 'section_number, expected an integer"))
    if (length(section_number) != 1) abort(glue("Invalid length {length(section_number)} for 'section_number, ",
                                                "expected an integer vector of length 1"))
    if (!section_number %in% unique(data_use$sampleID)) abort(glue("{section_number} out of range. ",
                                                                   "Select a number between {paste(range(data_use$sampleID), collapse = '-')}"))
    data_use <- data_use |>
      filter(sampleID == section_number)
  }

  # Get images if image_use is provided
  if (!is.null(image_use)) {
    image_use <- match.arg(image_use, choices = c("raw", "transformed"))
    images <- .get_images(object, GetStaffli(object), image_use, section_number)
    if (image_use == "transformed") {
      coords_use <- "transformed"
    }
  }

  # Set coords_columns
  coords_columns <- .get_coords_column(image_use, coords_use)

  # Filter data to remove all redundant meta data columns
  data_use <- data_use |>
    select(all_of("barcode"),
           all_of(coords_columns),
           all_of("sampleID"),
           all_of(features),
           contains(label_by %||% character(0)))

  # Create crop area if override_plot_dims = TRUE
  if (override_plot_dims) {
    # Split data by sampleID
    image_dims <- GetStaffli(object)@image_info
    if (!is.null(section_number)) {
      image_dims <- image_dims[image_dims$sampleID == section_number, ]
    }
    new_dims <- .get_limits(data_use, coords_columns)
    crop_area <- c(min(new_dims$x_start/image_dims$full_width),
                   min(new_dims$y_start/image_dims$full_height),
                   max(new_dims$full_width/image_dims$full_width),
                   max(new_dims$full_height/image_dims$full_height))
  }

  # crop images if a crop_area is set
  if (!is.null(image_use)) {
    images <- .crop_images(crop_area, images)
  }

  # rescale data if max/min_cuttoffs have been set
  if (!is.null(min_cutoff) | !is.null(max_cutoff)) {
    data_use <- .trim_data(data_use, features, min_cutoff, max_cutoff)
  }

  # Set dims
  dims <- GetStaffli(object)@image_info

  # generate plots
  wrapped_plots <- MapMultipleFeatures(
    object = data_use,
    scale = scale,
    crop_area = crop_area,
    pt_size = pt_size,
    pt_stroke = pt_stroke,
    section_number = NULL,
    label_by = label_by,
    ncol = ncol,
    colors = colors,
    dims = dims,
    coords_columns = coords_columns,
    return_plot_list = (!is.null(image_use)) | return_plot_list,
    add_colorscale_text = add_colorscale_text
  )

  # Inject images if image_use is provided
  if (!is.null(image_use)) {
    wrapped_plots <- .inject_images(image_use, features, arrange_features = "col", wrapped_plots, images, NULL, return_plot_list, FALSE)
    if (!return_plot_list) {
      # Create final patchwork
      ncol <- ncol %||% ceiling(sqrt(nrow(dims)))
      wrapped_plots <- wrap_plots(wrapped_plots, ncol = ncol)
    }
  }

  return(wrapped_plots)

}


#' Plot multiple numeric features in 2D
#'
#' @param gg tibble with spatial coordinates and feature values
#' @param nm sample ID
#' @param features A caharcter vector with feature IDs
#' @param colors a character vector of colors to use for scale bar
#' @param dims tibble containing information about the dimensions
#' of the plotting area
#' @param pt_size point size passed to geom_point
#' @param pt_stroke point stroke width
#' @param cur_label string to use as title
#' @param coords_columns a character vector of length 2 specifying names of
#' columns in which spatial coordinates are located
#' @param drop_na should NA values be dropped from the data?
#' @param use_text Logical specifying if values should be added to the color scales
#'
#' @import ggplot2
#' @import dplyr
#'
#' @return a `ggplot` object with a spatial plot
#'
#' @noRd
.spatial_feature_plot_multiple <- function (
    gg,
    nm,
    features = NULL,
    colors,
    dims,
    pt_size = 1,
    pt_stroke = 0,
    cur_label = NULL,
    coords_columns,
    drop_na = FALSE,
    use_text = FALSE
) {

  # Set global variables to NULL
  value <- alpha <- x <- y <- color <- NULL

  # require ggnewscale
  if (!requireNamespace("ggnewscale")) {
    install.packages("ggnewscale")
  }

  # Should NA values be dropped?
  # This will only filter out rows where all features are NA
  if (drop_na) {
    gg <- gg |> filter(if_any(all_of(features), ~ !is.na(.x)))
  }

  # Create input data.frame for ggplot
  gg <- gg |>
    select(all_of(c(coords_columns, features)))

  # Get opacity values
  alpha_values <- gg |>
    select(all_of(features)) |>
    apply(1, na.omit)
  gg$alpha <- alpha_values

  # Split data by feature
  gg_split <- lapply(features, function(ftr) {
    gg |>
      filter(!is.na(!! sym(ftr)))
  }) |>
    setNames(nm = features)

  # Draw plot
  p <-
    ggplot() +
        geom_point(data = gg_split[[features[1]]], aes(
          x = !!sym(coords_columns[1]),
          y = !!sym(coords_columns[2]),
          fill = !!sym(features[1]), # If blended colors are provided, add color outside aesthetic
          alpha = alpha
        ),
        size = pt_size,
        shape = 21,
        stroke = pt_stroke
    ) +
    # Set plot dimensions (reverse y axis)
    scale_x_continuous(limits = c(dims[dims$sampleID == nm, "x_start", drop = TRUE],
                                  dims[dims$sampleID == nm, "full_width", drop = TRUE]),
                       expand = c(0, 0),
                       breaks = seq(0, dims[dims$sampleID == nm, "full_width", drop = TRUE], length.out = 11),
                       labels = seq(0, 1, length.out = 11) |> paste0()) +
    scale_y_reverse(limits = c(dims[dims$sampleID == nm, "full_height", drop = TRUE],
                               dims[dims$sampleID == nm, "y_start", drop = TRUE]),
                    expand = c(0, 0),
                    breaks = seq(0, dims[dims$sampleID == nm, "full_height", drop = TRUE], length.out = 11),
                    labels = seq(0, 1, length.out = 11) |> paste0()) +
    # Add themes
    theme_void() +
    guides(alpha = "none") +
    # Add color gradient for first feature
    scale_fill_gradientn(colours = c("white", colors[1]), limits = c(0, 1),
                         guide = guide_colourbar(title.position = "right",
                                                 frame.colour = "black",
                                                 frame.linewidth = 1,
                                                 draw.ulim = FALSE,
                                                 draw.llim = FALSE,
                                                 label = use_text)) +
    theme(legend.position = "right",
          legend.direction = "horizontal",
          legend.title.align = 0,
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(0, 10, 20, 10),
          legend.title = element_text(vjust = 0.8)) +
    # Create a title
    labs(title = ifelse(!is.null(cur_label), cur_label, NA)) +
    # Fix coordinates so that plot cannot be stretched
    coord_fixed()

  # Add new color scales
  for (i in 2:length(features)) {
    p <- p +
      ggnewscale::new_scale_fill() +
      geom_point(data = gg_split[[features[i]]], aes(
          x = !!sym(coords_columns[1]),
          y = !!sym(coords_columns[2]),
          fill = !!sym(features[i]), # If blended colors are provided, add color outside aesthetic
          alpha = alpha
        ),
        size = pt_size,
        shape = 21,
        stroke = pt_stroke
      ) +
      guides(alpha = "none") +
      scale_fill_gradientn(colours = c("white", colors[i]), limits = c(0, 1),
                           guide = guide_colourbar(title.position = "right",
                                                   frame.colour = "black",
                                                   frame.linewidth = 1,
                                                   draw.ulim = FALSE,
                                                   draw.llim = FALSE,
                                                   label = use_text))
  }

  return(p)
}
