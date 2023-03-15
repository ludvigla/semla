#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param crop_area A numeric vector of length 4 specifying a rectangular area to crop
#' the plots by. These numbers should be within 0-1. The x-axis is goes from left=0 to
#' right=1 and the y axis is goes from top=0 to bottom=1. The order of the values are
#' specified as follows: \code{crop_area = c(left, top, right, bottom)}. The crop area
#' will be used on all tissue sections and cannot be set for each section individually.
#' using crop areas of different sizes on different sections can lead to unwanted side
#' effects as the point sizes will remain constant. In this case it is better to generate
#' separate plots for different tissue sections.
#' @param pt_size A numeric value specifying the point size passed to \code{geom_point}
#' @param pt_alpha A numeric value between 0 and 1 specifying the point opacity passed
#' to \code{geom_point}. A value of 0 will make the points completely transparent
#' and a value of 1 will make the points completely opaque.
#' @param pt_stroke A numeric specifying the point stroke width
#' @param scale_alpha Logical specifying if the spot colors should be scaled together with
#' the feature values. This can be useful when you want to highlight regions with higher
#' feature values while making the background tissue visible.
#' @param label_by Character of length 1 providing a column name in \code{object} with
#' labels that can be used to provide a title for each subplot. This column should have
#' 1 label per tissue section. This can be useful when you need to provide more detailed
#' information about your tissue sections.
#' @param ncol Integer value specifying the number of columns in the output patchwork.
#' This parameter will only have an effect when the number of features provided is 1.
#' Otherwise, the patchwork will be arranged based on the \code{arrange_features} parameter.
#' @param colors A character vector of colors to use for the color scale. The colors should
#' preferably consist of a set of colors from a scientific color palette designed for sequential
#' data. Some useful palettes are available in the \code{RColorBrewer}, \code{viridis} and
#' \code{scico} R packages.
#' @param scale A character vector of length 1 specifying one of "shared" or "free" which will
#' determine how the color bars are structured. If scale is set to "shared", the color bars for
#' feature values will be shared across samples. If scale is set to "free", the color bars
#' will be independent.
#' @param arrange_features One of "row" or "col". "col" will put the features in columns
#' and samples in rows and "row" will transpose the arrangement
#' @param dims A tibble with information about the tissue sections. This information is used to
#' determine the limits of the plot area. If \code{dims} is not provided, the limits will be
#' computed directly from the spatial coordinates provided in \code{object}.
#' @param coords_columns a character vector of length 2 specifying the names of the columns in
#' \code{object} holding the spatial coordinates
#' @param return_plot_list A logical specifying if a `patchwork` or a list of `ggplot` objects should be returned.
#' By default, a `patchwork` is returned, but it can sometimes be useful to obtain the list of `ggplot` objects
#' if you want to manipulate each sub plot independently.
#' @param blend a logical specifying whether blending should be used. See the section about color 
#' blending below for more information.
#' @param blend_order An integer vector of length 2-3 specifying the order to blend features by. Only
#' active when \code{blend = TRUE}. See color blending for more information.
#' @param drop_na A logical specifying if NA values should be dropped
#' @param center_zero A logical specifying whether the color scale should be centered at 0
#' @param add_scalebar A logical specifying if a scale bar should be added to the plots
#' @param scalebar_gg A 'ggplot' object generated with \code{\link{scalebar}}. The appearance of the scale
#' bar is styled by passing parameters to \code{\link{scalebar}}.
#' @param scalebar_height A numeric value specifying the height of the scale bar relative to the
#' height of the full plot area. Has to be a value between 0 and 1. The title of the scale bar is
#' scaled with the plot and might disappear if the down-scaled text size is too small. Increasing
#' the height of the scale bar can sometimes be useful to increase the text size to make it more
#' visible in small plots.
#' @param scalebar_position A numeric vector of length 2 specifying the position of the scale bar
#' relative to the plot area. Default is to place it in the top right corner.
#'
#' @importFrom patchwork wrap_plots
#' @import dplyr
#'
#' @section color blending:
#' Color blending can \strong{only} be used with 2 or three features. If blending is activated,
#' the feature values will be rescaled and encoded as RGB colors. RGB allows for three channels
#' to be included, hence the reason why you can only use 2-3 features. When blending feature
#' values you will only get 1 plot per tissue section.
#'
#' Colors can be mixed to produce new colors; for example, if two features have similar values
#' in one spot and are encoded as blue and red, the mixed color will be purple.
#'
#' \code{blend_order} allows you to flip the order of the features so that each feature
#' is provided with a color of choice. By default, the order is 1, 2, 3 which means that the
#' first feature is "red", the second is "green" and the last feature is "blue".
#'
#' @rdname visualize-features
#'
#' @return A patchwork object
#'
#' @export
#'
MapFeatures.default <- function (
  object,
  crop_area = NULL,
  pt_size = 1,
  pt_alpha = 1,
  pt_stroke = 0,
  scale_alpha = FALSE,
  section_number = NULL,
  label_by = NULL,
  ncol = NULL,
  colors = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
  center_zero = FALSE,
  scale = c("shared", "free"),
  arrange_features = c("col", "row"),
  dims = NULL,
  coords_columns = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  return_plot_list = FALSE,
  drop_na = FALSE,
  blend = FALSE,
  blend_order = 1:3,
  add_scalebar = FALSE,
  scalebar_gg = NULL,
  scalebar_height = 0.05,
  scalebar_position = c(0.8, 0.8),
  ...
) {

  # Set global variables to NULL
  barcode <- sampleID <- pxl_col_in_fullres <- pxl_row_in_fullres <- NULL

  # Check data
  .prep_data_for_plotting(object, colors, label_by, scale, arrange_features, coords_columns)

  # Expand colors if length is 1
  if (length(colors) == 1) {
    colors <- c("lightgray", colors)
  }

  # get features
  features <- object |>
    select(-barcode, -all_of(coords_columns), -sampleID, -all_of(label_by)) |>
    colnames()

  # Split data by sampleID
  data <- object |>
    group_by(sampleID) |>
    group_split() |>
    setNames(nm = unique(object$sampleID))

  # Check section number and subset data
  if (!is.null(section_number)) {
    if (!is.numeric(section_number)) abort(glue("Invalid class '{class(section_number)}' for",
                                                " {cli::col_br_green('section_number')}, expected 'numeric'"))
    if (!section_number %in% seq_along(data)) abort(glue("'section_number' = {section_number} is out of range. ",
                                                         "Select one of {paste(seq_along(data), collapse = ', ')}"))
    data <- data[section_number]
  }

  # get image dimensions
  dims <- .get_dims(dims)

  # get feature limits
  feature_limits <- .get_feature_limits(data, coords_columns, scale = ifelse(blend, "shared", scale))

  # add blend colors if blend=TRUE
  if (blend) {
    if (!requireNamespace("farver", quietly = TRUE)) {
      install.packages("farver")
    }
    data <- .color_blender(data, features, blend_order, feature_limits, scale_alpha)
    extreme_colors <- farver::encode_colour(diag(ncol = 3, nrow = 3)*255, from = "rgb")
    extreme_colors <- extreme_colors[blend_order[1:length(features)]]
    #features <- features[blend_order[1:length(features)]]
  }

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

    # Default plotting for each feature when blend = FALSE
    if (!blend) {
      feature_plots <- lapply(features, function(ftr) {
        .spatial_feature_plot(
          gg = gg,
          nm = nm,
          ftr = ftr,
          feature_limits = feature_limits,
          colors = colors,
          dims = dims,
          pt_size = pt_size,
          pt_alpha = pt_alpha,
          pt_stroke = pt_stroke,
          scale_alpha = scale_alpha,
          coords_columns = coords_columns,
          cur_label = cur_label,
          drop_na = drop_na,
          center_zero = center_zero
        )
      })
      # Add names to feature_plots
      p <- setNames(feature_plots, nm = features)
    } else {
      p <- .spatial_feature_plot(
        gg = gg,
        nm = nm,
        colors = colors,
        dims = dims,
        all_features = features,
        extreme_colors = extreme_colors,
        pt_size = pt_size,
        pt_alpha = pt_alpha,
        pt_stroke = pt_stroke,
        scale_alpha = scale_alpha,
        cur_label = cur_label,
        coords_columns = coords_columns,
        drop_na = drop_na,
        center_zero = center_zero
      )
    }
    return(p)
  }), nm = names(data))

  # Add scalebar
  if (add_scalebar) {
    if (!requireNamespace("dbscan")) {
      install.packages("dbscan")
    }
    scalebar_width <- scalebar_gg$labels$scalebar_width
    sample_plots <- lapply(names(sample_plots), function(nm) {
      gg <- data[[nm]]
      nn_dist <- dbscan::kNN(gg |> select(all_of(coords_columns)), k = 1)$dist[, 1] |> min()
      plots <- sample_plots[[nm]]
      d <- dims |> filter(sampleID == nm)
      sf <- scalebar_width/100
      prop_width <- (nn_dist*sf)/(d$full_width - d$x_start)
      scalebar_pos <- scalebar_position %||% c(0.8, 0.8)
      scalebar_pos[1] <- ifelse((1 - prop_width) > scalebar_pos[1], scalebar_pos[1], (1 - prop_width))
      scalebar_pos[2] <- ifelse((1 - scalebar_height) > scalebar_pos[2], scalebar_pos[2], (1 - scalebar_height))
      if (blend) {
        plots <- plots +
          inset_element(p = scalebar_gg, left = scalebar_pos[1], bottom = scalebar_pos[2], align_to = "full",
                        right = scalebar_pos[1] + prop_width, top = scalebar_pos[2] + scalebar_height, on_top = TRUE)
      } else {
        plots <- lapply(names(plots), function(ftr_nm) {
          plots[[ftr_nm]] +
            inset_element(p = scalebar_gg, left = scalebar_pos[1], bottom = scalebar_pos[2], align_to = "full",
                          right = scalebar_pos[1] + prop_width, top = scalebar_pos[2] + scalebar_height, on_top = TRUE)
        }) |> setNames(nm = names(plots))
      }
    }) |> setNames(nm = names(sample_plots))
  }

  # Create final patchwork
  if (!return_plot_list) {
    wrapped_plots <- .arrange_plots(sample_plots, features, blend, arrange_features, ncol)
  } else {
    # return list of ggplot objects if return_plot_list = TRUE
    wrapped_plots <- sample_plots
  }

  return(wrapped_plots)
}

# TODO: label.by broken, e.g. orig.ident
#' @param features A character vector of features to plot. These features need to be
#' fetchable with \code{link{FetchData}}
#' @param slot Slot to pull features values from
#' @param image_use A character specifying image type to use
#' @param coords_use A character specifying coordinate type to use
#' @param section_number An integer select a tissue section number to subset data by
#' @param override_plot_dims A logical specifying whether the image dimensions should
#' be used to define the plot area. Setting \code{override_plot_dims} can be useful
#' in situations where the tissue section only covers a small fraction of the capture
#' area, which will create a lot of white space in the plots.
#' @param min_cutoff,max_cutoff A numeric value between 0-1 specifying either a lower
#' (\code{min_cutoff}) or upper (\code{max_cutoff}) limit for the data using \code{\link{quantile}}.
#' These arguments can be useful to make sure that the color map doesn't get dominated by outliers.
#'
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom Seurat FetchData
#' @importFrom rlang warn
#'
#' @rdname visualize-features
#' @family spatial-visualization-methods
#'
#' @author Ludvig Larsson
#'
#' @examples
#'
#' library(semla)
#' if (!requireNamespace("viridis"))
#'   install.packages("viridis")
#' library(viridis)
#' library(ggplot2)
#'
#' # Load example Visium data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
#' se_merged <- MergeSTData(se_mbrain, se_mcolon)
#'
#' \dontrun{
#'
#' # Select features
#' selected_features <- c("Clu", "Slc6a3", "Vip")
#'
#' # Plot selected features with custom colors
#' MapFeatures(se_merged, features = selected_features, colors = magma(n = 11, direction = -1))
#'
#' # Plot selected features with color bars scaled individually for each feature and sample
#' MapFeatures(se_merged, features = selected_features, scale = "free")
#'
#' # Plot selected features and add custom labels for subplots
#' se_merged$sample_id <- ifelse(GetStaffli(se_merged)@meta_data$sampleID == "1",
#'                               "mouse brain",
#'                               "mouse colon")
#' MapFeatures(se_merged, features = selected_features, label_by = "sample_id")
#'
#' # Plot selected features arranged by rows instead of columns
#' MapFeatures(se_merged, features = selected_features, arrange_features = "row")
#'
#' # Blend features
#' MapFeatures(se_merged, features = selected_features, blend = TRUE)
#'
#' }
#'
#' # The output is a patchwork object which is easy to manipulate
#' selected_feature <- "Th"
#'
#' # Move legend to right side of plot
#' MapFeatures(se_mbrain, features = selected_feature, pt_size = 2) &
#'   theme(legend.position = "right",
#'         legend.text = element_text(angle = 0, hjust = 0))
#'
#' # Remove titles
#' MapFeatures(se_mbrain, features = selected_feature, pt_size = 2) &
#'   theme(plot.title = element_blank(),
#'         plot.subtitle = element_blank())
#'
#' # Create a dark theme
#' MapFeatures(se_mbrain, features = selected_feature, pt_size = 2,
#'             colors = viridis(n = 11)) &
#'   theme(plot.background = element_rect(fill = "black"),
#'         panel.background = element_rect(fill = "black"),
#'         plot.title = element_text(colour = "white"),
#'         plot.subtitle = element_text(colour = "white"),
#'         legend.text = element_text(colour = "white"),
#'         legend.title = element_text(colour = "white"))
#'
#' # Add a background and axes to the plot
#' MapFeatures(se_mbrain, features = selected_feature, pt_size = 2) &
#'   theme(panel.background = element_rect(fill = "lightgray"),
#'         axis.text = element_text())
#'
#' # Move legend
#' MapFeatures(se_mbrain, features = selected_feature, pt_size = 2) &
#'   theme(panel.background = element_rect(fill = "lightgray"),
#'         legend.justification = 1)
#'
#' # Change title of fill aesthetic
#' MapFeatures(se_mbrain, features = selected_feature, pt_size = 2) &
#'   labs(title = "Gene expression")
#'
#' @export
#'
MapFeatures.Seurat <- function (
    object,
    features,
    slot = "data",
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
    center_zero = FALSE,
    scale = c("shared", "free"),
    arrange_features = c("col", "row"),
    drop_na = FALSE,
    blend = FALSE,
    blend_order = 1:3,
    override_plot_dims = FALSE,
    max_cutoff = NULL,
    min_cutoff = NULL,
    return_plot_list = FALSE,
    add_scalebar = FALSE,
    scalebar_height = 0.05,
    scalebar_gg = scalebar(x = 500, text_height = 5),
    scalebar_position = c(0.8, 0.7),
    ...
) {

  # Set global variables to NULL
  sampleID <- NULL

  # Match args
  scale <- match.arg(scale, choices = c("shared", "free"))
  arrange_features <- match.arg(arrange_features, choices = c("col", "row"))
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

  if (add_scalebar) {
    # Check input param
    if (!inherits(scalebar_gg, what = "gg")) {
      abort(glue("Invalid {col_br_magenta('scalebar_gg')}. Expected a 'ggplot' object"))
    }
    if (!between(x = scalebar_height, left = 0, right = 1)) {
      abort(glue("Invalid {col_br_magenta('scalebar_height')}. Expected a numeric of length 1 between 0 and 1"))
    }
    if (length(scalebar_position) != 2) {
      abort(glue("Invalid {col_br_magenta('scalebar_position')}. Expected a numeric of length 2"))
    }
   if (!all(between(x = scalebar_position, left = 0, right = 1))) {
     abort(glue("Invalid {col_br_magenta('scalebar_position')}. Expected values between 0 and 1"))
   }
  } else {
    scalebar_height <- NULL
    scalebar_gg <- NULL
    scalebar_position <- NULL
  }

  # generate plots
  wrapped_plots <- MapFeatures(
    object = data_use,
    crop_area = crop_area,
    pt_size = pt_size,
    pt_alpha = pt_alpha,
    pt_stroke = pt_stroke,
    scale_alpha = scale_alpha,
    section_number = NULL,
    label_by = label_by,
    ncol = ncol,
    colors = colors,
    center_zero = center_zero,
    scale = scale,
    arrange_features = arrange_features,
    dims = dims,
    coords_columns = coords_columns,
    drop_na = drop_na,
    blend = blend,
    blend_order = blend_order,
    return_plot_list = (!is.null(image_use)) | return_plot_list,
    add_scalebar = add_scalebar,
    scalebar_height = scalebar_height,
    scalebar_gg = scalebar_gg,
    scalebar_position = scalebar_position
  )

  # Inject images if image_use is provided
  if (!is.null(image_use)) {
    wrapped_plots <- .inject_images(image_use, features, arrange_features, wrapped_plots, images, NULL, return_plot_list, blend)
    if (!return_plot_list) {
      # Create final patchwork
      wrapped_plots <- .arrange_plots(wrapped_plots, features, blend, arrange_features, ncol)
    }
  }

  return(wrapped_plots)

}


#' @param crop_area A numeric vector of length 4 specifying a rectangular area to crop
#' the plots by. These numbers should be within 0-1. The x-axis is goes from left=0 to
#' right=1 and the y axis is goes from top=0 to bottom=1. The order of the values are
#' specified as follows: \code{crop_area = c(left, top, right, bottom)}. The crop area
#' will be used on all tissue sections and cannot be set for each section individually.
#' using crop areas of different sizes on different sections can lead to unwanted side
#' effects as the point sizes will remain constant. In this case it is better to generate
#' separate plots for different tissue sections.
#' @param pt_size A numeric value specifying the point size passed to \code{geom_point}
#' @param pt_alpha A numeric value between 0 and 1 specifying the point opacity passed
#' to \code{geom_point}. A value of 0 will make the points completely transparent
#' and a value of 1 will make the points completely opaque.
#' @param pt_stroke A numeric specifying the point stroke width
#' @param section_number An integer select a tissue section number to subset data by
#' @param label_by A character specifying a column name in \code{object} with
#' labels that can be used to provide a title for each subplot. This column should have
#' 1 label per tissue section. This can be useful when you need to provide more detailed
#' information about your tissue sections.
#' @param ncol An integer value specifying the number of columns in the output patchwork.
#' @param colors A character vector of colors to use for the color scale. The number of
#' colors should match the number of labels present.
#' @param dims A tibble with information about the tissue sections. This information is used to
#' determine the limits of the plot area. If \code{dims} is not provided, the limits will be
#' computed directly from the spatial coordinates provided in \code{object}.
#' @param coords_columns A character vector of length 2 specifying the names of the columns in
#' \code{object} holding the spatial coordinates
#' @param return_plot_list logical specifying whether the plots should be return as a list of
#' `ggplot` objects. If \code{return_plot_list = FALSE} (default), the plots will be arranged
#' into a `patchwork`
#' @param drop_na A logical specifying if NA values should be dropped
#' @param add_scalebar A logical specifying if a scale bar should be added to the plots
#' @param scalebar_gg A 'ggplot' object generated with \code{\link{scalebar}}. The appearance of the scale
#' bar is styled by passing parameters to \code{\link{scalebar}}.
#' @param scalebar_height A numeric value specifying the height of the scale bar relative to the
#' height of the full plot area. Has to be a value between 0 and 1. The title of the scale bar is
#' scaled with the plot and might disappear if the down-scaled text size is too small. Increasing
#' the height of the scale bar can sometimes be useful to increase the text size to make it more
#' visible in small plots.
#' @param scalebar_position A numeric vector of length 2 specifying the position of the scale bar
#' relative to the plot area. Default is to place it in the top right corner.
#'
#' @importFrom patchwork wrap_plots
#' @import dplyr
#' @importFrom zeallot %<-%
#' @importFrom rlang %||% warn
#'
#' @rdname visualize-labels
#' @family spatial-visualization-methods
#'
#' @return A `patchwork` object or a list of `ggplot` objects
#'
#' @export
#'
MapLabels.default <- function (
  object,
  crop_area = NULL,
  pt_size = 1,
  pt_alpha = 1,
  pt_stroke = 0,
  section_number = NULL,
  label_by = NULL,
  split_labels = FALSE,
  ncol = NULL,
  colors = NULL,
  dims = NULL,
  coords_columns = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  return_plot_list = FALSE,
  drop_na = FALSE,
  add_scalebar = FALSE,
  scalebar_gg = NULL,
  scalebar_height = 0.05,
  scalebar_position = c(0.8, 0.8),
  ...
) {

  # Set global variables to NULL
  barcode <- sampleID <- NULL

  # Check data
  .prep_data_for_plotting(object = object, label_by = label_by, coords_columns = coords_columns)

  # get label column
  label <- object |>
    select(-barcode, -all_of(coords_columns), -sampleID, -all_of(label_by)) |>
    colnames()

  # Check if label column only contains NA values
  if (object |> pull(all_of(label)) |> is.na() |> all()) abort(glue("Selected feature only contains NA values."))

  # Convert label column to factor
  object <- object |>
    mutate(across(all_of(label), ~ factor(.x)))

  # Split data by sampleID
  data <- object |>
    group_by(sampleID) |>
    group_split() |>
    setNames(nm = unique(object$sampleID))

  # Check section number and subset data
  if (!is.null(section_number)) {
    if (!is.numeric(section_number)) abort(glue("Invalid class '{class(section_number)}' for",
                                                " {cli::col_br_green('section_number')}, expected 'numeric'"))
    if (!section_number %in% seq_along(data)) abort(glue("'section_number' = {section_number} is out of range. ",
                                                         "Select one of {paste(seq_along(data), collapse = ', ')}"))
  }
  if (split_labels) {
    section_number <- section_number %||% {
      warn("No section_number selected. Selecting section 1.")
      1L
    }
    data <- data[section_number]
    c(data, dims, colors) %<-% .split_data_by_label(data, dims[dims$sampleID == section_number, ], label, colors, drop_na)
  }
  if (!is.null(section_number) & !split_labels) {
    data <- data[section_number]
  }

  # Obtain colors
  colors <- colors %||% .gg_color_hue(length(levels(object |> select(all_of(label)) |> pull(all_of(label)))))

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

    # Set label if available
    if (!is.null(label_by)) {
      cur_label <- unique(gg |> pull(all_of(label_by)) |> as.character())
    } else {
      cur_label <- paste0("section ", nm)
    }

    # Overwrite label if split_labels = TRUE
    if (split_labels) {
      cur_label <- paste0("label: ", nm)
    }

    p <- .spatial_label_plot(
      gg = gg,
      nm = nm,
      lbl = label,
      colors = colors,
      dims = dims,
      pt_size = pt_size,
      pt_alpha = pt_alpha,
      pt_stroke = pt_stroke,
      coords_columns = coords_columns,
      cur_label = cur_label,
      drop_na = drop_na
    )

    return(p)
  }), nm = names(data))
  
  # Add scalebar
  if (add_scalebar) {
    if (!requireNamespace("dbscan")) {
      install.packages("dbscan")
    }
    scalebar_width <- scalebar_gg$labels$scalebar_width
    sample_plots <- lapply(names(sample_plots), function(nm) {
      gg <- data[[nm]]
      nn_dist <- dbscan::kNN(gg |> select(all_of(coords_columns)), k = 1)$dist[, 1] |> min()
      plots <- sample_plots[[nm]]
      d <- dims |> filter(sampleID == nm)
      sf <- scalebar_width/100
      prop_width <- (nn_dist*sf)/(d$full_width - d$x_start)
      scalebar_pos <- scalebar_position %||% c(0.8, 0.8)
      scalebar_pos[1] <- ifelse((1 - prop_width) > scalebar_pos[1], scalebar_pos[1], (1 - prop_width))
      scalebar_pos[2] <- ifelse((1 - scalebar_height) > scalebar_pos[2], scalebar_pos[2], (1 - scalebar_height))
      plots <- plots +
        inset_element(p = scalebar_gg, left = scalebar_pos[1], bottom = scalebar_pos[2], align_to = "full",
                      right = scalebar_pos[1] + prop_width, top = scalebar_pos[2] + scalebar_height, on_top = TRUE)
    }) |> setNames(nm = names(sample_plots))
  }

  if (!return_plot_list) {
    # Create final patchwork
    ncol <- ncol %||% ceiling(sqrt(length(data)))
    wrapped_plots <- wrap_plots(sample_plots, ncol = ncol)

    return(wrapped_plots)
  } else {
    return(sample_plots)
  }
}

#' @param column_name A character specifying a meta data column holding the categorical
#' feature vector.
#' @param image_use A character specifying image type to use.
#' @param coords_use A character specifying coordinate type to use.
#' @param split_labels A logical specifying if labels should be split.
#' @param override_plot_dims A logical specifying whether the image dimensions should
#' be used to define the plot area. Setting \code{override_plot_dims} can be useful
#' in situations where the tissue section only covers a small fraction of the capture
#' area, which will create a lot of white space in the plots. The same effect can be
#' achieved with the \code{crop_area} but the crop area is instead determined directly
#' from the data.
#'
#' @importFrom rlang warn
#' @import dplyr
#' @importFrom Seurat FetchData
#' @importFrom tibble as_tibble
#' @importFrom glue glue
#' @importFrom patchwork wrap_plots
#'
#' @rdname visualize-labels
#'
#' @author Ludvig Larsson
#'
#' @examples
#'
#' library(semla)
#' library(ggplot2)
#'
#' # Load Seurat object
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
#'
#' # Run PCA and data-driven clustering
#' se_mbrain <- se_mbrain |>
#'   ScaleData() |>
#'   RunPCA() |>
#'   FindNeighbors(reduction = "pca", dims = 1:10) |>
#'   FindClusters(resolution = 0.2) |>
#'   FindClusters(resolution = 0.3)
#'
#' # Plot clusters
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2", pt_size = 2)
#'
#' # Plot clusters in split view
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2", pt_size = 0.5,
#'           section_number = 1, split_labels = TRUE, ncol = 4)
#'
#' # Combine plots with different labels
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2") |
#'   MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.3")
#'
#' # Move legend to the right side of the plot
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2", pt_size = 2) &
#'   theme(legend.position = "right")
#'
#' # Override fill aesthetic to increase point sixe in legend
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2", pt_size = 2) &
#'   guides(fill = guide_legend(override.aes = list(size = 4)))
#'
#' # Use custom colors
#' cols <- c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2", pt_size = 2, colors = cols)
#'
#' # Factor are to used to determine the color order. If you change the
#' # levels of your label of interest, the labels and colors will change order
#' se_mbrain$Spatial_snn_res.0.2 <- factor(se_mbrain$Spatial_snn_res.0.2,
#'                                         levels = sample(paste0(0:7), size = 8))
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2", pt_size = 2, colors = cols)
#'
#' # Control what group label colors by naming the color vector
#' # this way you can make sure that each group gets a desired color
#' # regardless of the factor levels
#' cols <- setNames(cols, nm = paste0(0:5))
#' cols
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2", pt_size = 2, colors = cols)
#'
#' @export
#'
MapLabels.Seurat <- function (
  object,
  column_name,
  image_use = NULL,
  coords_use = "raw",
  crop_area = NULL,
  pt_size = 1,
  pt_alpha = 1,
  pt_stroke = 0,
  section_number = NULL,
  label_by = NULL,
  split_labels = FALSE,
  ncol = NULL,
  colors = NULL,
  override_plot_dims = FALSE,
  return_plot_list = FALSE,
  drop_na = FALSE,
  add_scalebar = FALSE,
  scalebar_height = 0.05,
  scalebar_gg = scalebar(x = 500, text_height = 5),
  scalebar_position = c(0.8, 0.7),
  ...
) {

  # Match args
  coords_use <- match.arg(coords_use, choices = c("raw", "transformed"))

  # Check Seurat object
  .check_seurat_object(object)

  # fetch data from Seurat object
  data_use <- GetStaffli(object)@meta_data |>
    bind_cols(FetchData(object, vars = column_name) |> as_tibble())

  # Convert label column to factor
  data_use <- data_use |>
    mutate(across(all_of(column_name), ~ factor(.x)))

  # Add label_by column if present
  if (!is.null(label_by)) {
    if (label_by == column_name) {
      abort(glue("{col_br_magenta('column_name')} should not be the same as {col_br_magenta('label_by')}"))
    }
    data_use <- data_use |>
      bind_cols(FetchData(object, vars = label_by))
  }

  # Get images if image_use is provided
  if (!is.null(image_use)) {
    image_use <- match.arg(image_use, choices = c("raw", "transformed"))
    images <- .get_images(data_use, GetStaffli(object), image_use, section_number, column_name, split_labels)
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
           all_of(column_name),
           contains(label_by %||% character(0)))

  # Create crop area if override_plot_dims = TRUE
  if (override_plot_dims) {
    # Split data by sampleID
    image_dims <- GetStaffli(object)@image_info
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

  # Set dims
  dims <- GetStaffli(object)@image_info
  
  if (add_scalebar) {
    # Check input param
    if (!inherits(scalebar_gg, what = "gg")) {
      abort(glue("Invalid {col_br_magenta('scalebar_gg')}. Expected a 'ggplot' object"))
    }
    if (!between(x = scalebar_height, left = 0, right = 1)) {
      abort(glue("Invalid {col_br_magenta('scalebar_height')}. Expected a numeric of length 1 between 0 and 1"))
    }
    if (length(scalebar_position) != 2) {
      abort(glue("Invalid {col_br_magenta('scalebar_position')}. Expected a numeric of length 2"))
    }
    if (!all(between(x = scalebar_position, left = 0, right = 1))) {
      abort(glue("Invalid {col_br_magenta('scalebar_position')}. Expected values between 0 and 1"))
    }
  } else {
    scalebar_height <- NULL
    scalebar_gg <- NULL
    scalebar_position <- NULL
  }

  # generate plots
  wrapped_plots <- MapLabels(
    object = data_use,
    crop_area = crop_area,
    pt_size = pt_size,
    pt_alpha = pt_alpha,
    pt_stroke = pt_stroke,
    section_number = section_number,
    label_by = label_by,
    split_labels = split_labels,
    ncol = ncol,
    colors = colors,
    dims = dims,
    coords_columns = coords_columns,
    return_plot_list = (!is.null(image_use)) | return_plot_list,
    drop_na = drop_na,
    add_scalebar = add_scalebar,
    scalebar_height = scalebar_height,
    scalebar_gg = scalebar_gg,
    scalebar_position = scalebar_position
  )

  # Inject images if image_use is provided
  if (!is.null(image_use)) {
    wrapped_plots <- .inject_images(
      image_use = image_use,
      features = NULL,
      arrange_features = NULL,
      wrapped_plots = wrapped_plots,
      images = images,
      ncol = ncol,
      return_plot_list = return_plot_list,
      blend = FALSE
    )
    # Create final patchwork
    if (!return_plot_list) {
      ncol <- ncol %||% ceiling(sqrt(length(wrapped_plots)))
      wrapped_plots <- wrap_plots(wrapped_plots, ncol = ncol)
    }
  }

  return(wrapped_plots)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plotting utils
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Plot numeric features in 2D
#'
#' @param gg A tibble with spatial coordinates and feature values
#' @param nm A sample ID
#' @param ftr A feature name
#' @param feature_limits A list of tibbles containing information about
#' the feature value range
#' @param colors A character vector of colors to use for scale bar
#' @param dims A tibble containing information about the dimensions
#' of the plotting area
#' @param all_features A character vector with all features, only used
#' for blending features
#' @param extreme_colors A character vector with the hex colors, only
#' used for blending features
#' @param pt_size Point size passed to geom_point
#' @param pt_alpha Point opacity ranging from 0 to 1 passed to geom_point.
#' 0 = fully transparent, 1 = fully opaque
#' @param pt_stroke Point stroke width
#' @param scale_alpha Should the spot opacity be scaled along with the feature values?
#' @param cur_label A string to use as title
#' @param coords_columns A character vector of length 2 specifying names of
#' columns in which spatial coordinates are located
#' @param drop_na Should NA values be dropped from the data?
#' @param center_zero A logical specifying whether color scale should be centered at 0
#'
#' @import ggplot2
#' @import dplyr
#'
#' @return a `ggplot` object with a spatial plot
#'
#' @noRd
.spatial_feature_plot <- function (
    gg,
    nm,
    ftr = NULL,
    feature_limits = NULL,
    colors,
    dims,
    all_features = NULL,
    extreme_colors = NULL,
    pt_size = 1,
    pt_alpha = 1,
    pt_stroke = 0,
    scale_alpha = FALSE,
    cur_label = NULL,
    coords_columns,
    drop_na = FALSE,
    center_zero = FALSE
) {

  # Set global variables to NULL
  encoded_cols <- value <- alpha <- x <- y <- color <- NULL

  # Should NA values be dropped?
  if (drop_na) {
    gg <- gg |> filter(if_all(all_of(ftr), ~ !is.na(.x)))
  }

  # Check if encoded colors are present after setting blend = TRUE
  encoded_cols_present <- "encoded_cols" %in% colnames(gg)
  color_vec <- switch(encoded_cols_present + 1, NULL, gg |> pull(encoded_cols))

  # Get opacity values if scale_alpha = TRUE and encoded colors are present
  if (scale_alpha & encoded_cols_present) {
    alpha_values <- gg$alpha
  }

  # Create input data.frame for ggplot
  if (!encoded_cols_present) {
    # Select coordinates and features
    gg <- gg |>
      select(all_of(c(coords_columns, ftr))) |>
      setNames(nm = c(coords_columns, "value"))
  } else {
    # Remove features if blend = TRUE
    gg <- gg |>
      select(all_of(coords_columns))
  }

  # Get opacity values if scale_alpha = TRUE and encoded colors are not present
  if (scale_alpha & !encoded_cols_present) {
    alpha_values <- gg |>
      select(value) |>
      mutate(alpha = scales::rescale(
        x = value,
        to = c(0, 1),
        from = c(feature_limits[[nm]][1, ftr, drop = TRUE], feature_limits[[nm]][2, ftr, drop = TRUE])
      )) |>
      pull(alpha)
  }

  # Draw plot
  p <-
    ggplot() +
    {
      if (encoded_cols_present) {
        geom_point(data = gg, aes(
          x = !!sym(coords_columns[1]),
          y = !!sym(coords_columns[2])
          ),
          fill = color_vec, # If blended colors are provided, add custom colors outside aesthetic
          size = pt_size,
          alpha = switch(scale_alpha + 1, pt_alpha, alpha_values),
          shape = 21,
          stroke = pt_stroke
        )
      } else {
        geom_point(data = gg, aes(
          x = !!sym(coords_columns[1]),
          y = !!sym(coords_columns[2]),
          fill = value # If blended colors are provided, add color outside aesthetic
          ),
          size = pt_size,
          alpha = switch(scale_alpha + 1, pt_alpha, alpha_values),
          shape = 21,
          stroke = pt_stroke
        )
      }
    } +
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
    theme(legend.position = "top",
          legend.text = element_text(angle = 60, hjust = 1),
          legend.justification = "left",
          legend.margin = margin(t = 5, 0, 10, 0),
          plot.margin = margin(0, 10, 20, 10),
          legend.title = element_text(vjust = 0.8)) +
    # Add color gradient if blend = FALSE
    {
      if (!encoded_cols_present) {
        scale_fill_gradientn(colours = colors,
                              limits = c(ifelse(!center_zero,
                                                feature_limits[[nm]][1, ftr, drop = TRUE],
                                                -max(abs(feature_limits[[nm]][1:2, ftr, drop = TRUE]))),
                                         ifelse(!center_zero,
                                                feature_limits[[nm]][2, ftr, drop = TRUE],
                                                max(abs(feature_limits[[nm]][1:2, ftr, drop = TRUE])))))
      }
    } +
    # Create a title
    {
      if (!encoded_cols_present) {
        labs(title = ifelse(!is.null(cur_label), cur_label, NA),
             subtitle = paste("feature: ", ftr))
      }
    } +
    # Fix coordinates so that plot cannot be stretched
    coord_fixed()

  # Add a color legend when blend = TRUE
  if (encoded_cols_present) {
    p <- p +
      geom_point(data = data.frame(x = rep(Inf, length(all_features)),
                                   y = rep(Inf, length(all_features)),
                                   color = all_features),
                 aes(x, y, color = color),
                 size = 0) +
      guides(color = guide_legend(override.aes = list(size = 4))) +
      theme(legend.text = element_text(angle = 0)) +
      scale_color_manual(values = extreme_colors, labels = all_features) +
      labs(title = ifelse(!is.null(cur_label), cur_label, NA))
  }

  return(p)
}

#' Plot labels in 2D
#'
#' @param gg tibble with spatial coordinates and a label column
#' @param nm sample ID
#' @param lbl label column name
#' @param colors a character vector of colors IDs
#' @param dims tibble containing information about the dimensions
#' of the plotting area
#' @param pt_size point size passed to geom_point
#' @param pt_alpha point opacity ranging from 0 to 1 passed to geom_point.
#' 0 = fully transparent, 1 = fully opaque
#' @param pt_stroke point stroke width
#' @param coords_columns character vector of length 2 specifying names of
#' columns in which the spatial coordinates are stored
#' @param cur_label string with a title
#' @param drop_na logical specifying whether NA values should be dropped
#'
#' @import ggplot2
#' @importFrom stats na.omit
#' @import dplyr
#' @importFrom rlang !!
#'
#' @return a `ggplot` object with a spatial plot
#'
#' @noRd
.spatial_label_plot <- function (
    gg,
    nm,
    lbl,
    colors,
    dims,
    pt_size = 1,
    pt_alpha = 1,
    pt_stroke = 0,
    coords_columns,
    cur_label,
    drop_na = FALSE
) {

  # Set global variables to NULL
  variable <- NULL

  # Create input data.frame for ggplot
  gg <- gg |>
    select(all_of(c(coords_columns, lbl))) |>
    setNames(nm = c(coords_columns, "variable"))

  # Drop NA values
  if (drop_na) {
    gg <- gg |>
      na.omit()
  }

  # Rearrange colors by factor level
  if (!is.null(names(colors))) {
    colors <- colors[levels(gg$variable)]
  } else {
    colors <- setNames(colors, levels(gg$variable))
  }

  # Check colors
  if (length(colors) < length(levels(gg$variable))) {
    abort(glue("The number of colors ({length(colors)})",
               " does not match the number of labels ({length(levels(gg$variable))})."))
  } else {
    colors <- colors[levels(gg$variable)]
  }
  
  # Draw plot
  p <-
    ggplot(
      data = gg, aes(
      x = !! sym(coords_columns[1]),
      y = !! sym(coords_columns[2]),
      fill = variable)
    ) +
    geom_point(
      size = pt_size,
      alpha = pt_alpha,
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
    theme(legend.position = "top",
          legend.justification = "left",
          legend.margin = margin(t = 5, 0, 10, 0),
          plot.margin = margin(0, 10, 20, 10),
          legend.title = element_text(vjust = 0.8)) +
    # Add colors
    scale_fill_manual(values = colors) +
    # Create a title
    labs(title = ifelse(!is.null(cur_label), cur_label, NA), fill = lbl) +
    # Fix coordinates so that plot cannot be stretched
    coord_fixed()

  return(p)
}


#' Check input for compatibility
#'
#' This function is used internally by \code{MapLabels} and \code{MapFeatures} to validate
#' input parameters.
#'
#' @param object a tibble with spatial coordinates and feature values
#' @param colors a character vector of colors
#' @param label_by a column name from which to fetch labels from
#' @param scale one of "free" or "shared". "free" will fix the value limits for each
#' feature separately and "shared" will fix the value limits for all plots to be identical
#' @param arrange_features one of "row" or "col". "col" will put the features in columns
#' and samples in rows and "row" will transpose the arrangement
#' @param coords_columns a character vector with column names for spot coordinates
#' @param multi_color Logical specifying if multi-color plot is used. Necessary for
#' compatibility with \code{\link{MapMultipleFeatures}}.
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#' @import dplyr
#' @importFrom tidyr unite
#'
#' @return nothing
#'
#' @noRd
.prep_data_for_plotting <- function (
    object,
    colors,
    label_by,
    scale,
    arrange_features,
    coords_columns,
    multi_color = FALSE
) {

  # Set global variables to NULL
  merged_cols <- barcode <- sampleID <- NULL

  # check colors
  if (!missing(colors)) {
    if (!is.character(colors)) abort(glue("Invalid class '{class(colors)}' of colors"))
    if (length(colors) < 1) abort(glue("At least 1 color1 needs to be provided, got 0"))
  }

  # Check label by
  if (!is.null(label_by)) {
    if (!is.character(label_by)) abort(glue("Invalid class for 'label_by', got class '{class(label_by)}' but expected 'character'"))
    if (length(label_by) != 1) abort(glue("Invalid length of 'label_by', got {length(label_by)} elements but expected 1"))
    if (!label_by %in% colnames(object)) abort("'label_by' could not be found in the input data.")
    check_labels <- object |>
      unite(col = "merged_cols", all_of(c("sampleID", label_by))) |>
      dplyr::count(merged_cols)
    if (nrow(check_labels) > length(unique(object$sampleID))) abort("The 'label_by' column does not match the samples.")
  }

  # Check scale
  if (!missing(scale)) scale <- match.arg(scale, choices = c("shared", "free"))

  # Check arrange_features  argument
  if (!missing(arrange_features)) arrange_features <- match.arg(arrange_features , choices =  c("col", "row"))

  # Check objects
  if (!any(class(object) %in% c("tbl", "data.frame"))) abort(glue("Invalid format of feature matrix: '{class(object)}'"))
  if (any(dim(object) == 0)) abort(glue("Invalid dimensions of input: '{paste(dim(object), collapse = 'x')}'"))


  # Check that spatial coordinates are present
  if (!all(coords_columns %in% colnames(object))) abort(glue("Coordinates are missing.'"))

  # Check that barcodes and sampleID are present
  if (!"barcode" %in% colnames(object)) abort(glue("Barcodes are missing.'"))
  if (!"sampleID" %in% colnames(object)) abort(glue("Sample IDs are missing.'"))

  # Check that features are valid
  if (!missing(colors)) {
    checks <- object |>
      select(-barcode, -all_of(coords_columns), -sampleID, -all_of(label_by)) |>
      sapply(is.numeric)
    if (any(!checks)) abort(glue("Features have to be numeric/integer. \n",
                                 "The following features are not valid: \n {paste(names(checks[!checks]))}"))
  } else {
    if (!multi_color) {
      checks <- object |>
        select(-barcode, -all_of(coords_columns), -sampleID, -all_of(label_by)) |>
        sapply(function(x) {
          is.character(x) | is.factor(x)
        })
      if (length(checks) > 1) abort("Only 1 meta.data column for a categorical variable is allowed.")
      if (!checks) abort(glue("Categorical variables have to be a character/factor."))
    }
  }
}

#' Trim numeric features based on predefined cutoffs
#'
#' @param data a `tibble` containing spatial coordinates and feature values
#' @param features a character vector with feature names present in \code{data}
#' @param min_cutoff,max_cutoff numeric specifying a cutoff threshold (quantile)
#'
#' @importFrom rlang %||%
#' @importFrom dplyr between
#'
#' @return a `tibble` with trimmed feature values
#'
#' @noRd
.trim_data <- function (
    data,
    features,
    min_cutoff,
    max_cutoff
) {
  min_cutoff <- min_cutoff %||% 0
  max_cutoff <- max_cutoff %||% 1

  # Check cutoffs
  if (!all(between(x = c(min_cutoff, max_cutoff), left = 0, right = 1))) abort("min/max cutoffs cannot be outside the range 0-1")

  # Cut data
  data <- data |>
    mutate(across(
      all_of(features),
      ~ case_when(.x < quantile(.x, probs = min_cutoff, na.rm = TRUE) ~ quantile(.x, probs = min_cutoff, na.rm = TRUE),
                  .x > quantile(.x, probs = max_cutoff, na.rm = TRUE) ~ quantile(.x, probs = max_cutoff, na.rm = TRUE),
                  TRUE ~ .x)
    ))
  return(data)
}


#' Split data by label
#'
#' Takes a list of tibbles holding spot coordinates and a label column,
#' then selects one of these based on \code{section_number} to use for
#' restructuring the data. Each group in the label column will get its
#' own element in the output list of tibbles, where the label column is
#' set to the group label or "background". This prepares the data for
#' highlighting 1 group label at the time to use the \code{split_labels}
#' option in \code{\link{MapLabels}}.
#'
#' @param data a list of tibbles containing spot coordinates and label columns
#' @param dims a tibble with information about plot dimensions to use
#' @param label a string specifying the name of the label column
#' @param colors a character vector of colors that should match the number of
#' groups in the label column
#' @param drop_na should NA values be dropped?
#'
#' @import dplyr
#' @importFrom forcats fct_drop
#' @importFrom rlang %||%
#'
#' @return a list of tibbles holding spot coordinates, a tibble with information
#' about plot dimensions and a character vector of colors stored together in a list
#'
#' @noRd
.split_data_by_label <- function (
  data,
  dims,
  label,
  colors,
  drop_na = FALSE
) {

  data <- lapply(data, function(x) {
    x |>
      mutate(across(all_of(label), ~ fct_drop(.x)))
  })

  # restructure data so that each level in the selected factor
  # gets its own tibble and this way we can arrange 1 plot per level
  lvls <- data[[1]] |> pull(all_of(label)) |> levels()
  if (length(lvls) == 0) {
    abort(glue("Only NA values found for selected column variable"))
  }
  bg_name <- ifelse("background" %in% (data[[1]] |> pull(all_of(label))), "background__", "background")
  if (length(lvls) == 1) abort(glue("Only 1 group '{lvls}' present in '{label}'. Cannot split data."))
  data <- lapply(lvls, function(lvl) {
    # mutate data: if select factor values match lbl (current group)
    # set new column values to background, then set all other values
    # to NA or background
    data[[1]] |> mutate(across(all_of(label), ~ case_when(
      .x == lvl ~ lvl,
      TRUE ~ ifelse(drop_na, NA_character_, bg_name)
    ))) |>
      mutate(across(all_of(label), ~ factor(.x, levels = c(
        lvls, ifelse(drop_na, NA_character_, bg_name)
      ))))
  }) |>
    setNames(nm = lvls)

  # restructure dims
  dims <- do.call(bind_rows, lapply(lvls, function(lvl) {
    dims |> mutate(sampleID = lvl)
  }))

  # If colors is NULL, we create as many colors as the
  # number of levels + lightgrey for the background
  colors <- colors %||% setNames(c(.gg_color_hue(length(lvls)), "lightgray"),
                                 nm = c(lvls, ifelse(drop_na, NA_character_, bg_name)))
  # If colors are provided, add another color for background
  if (length(colors) == length(lvls)) {
    if (!is.null(names(colors))) {
      colors <- setNames(c(colors, "lightgray"),
                         nm = c(names(colors), ifelse(drop_na, NA_character_, bg_name)))
    } else {
      colors <- setNames(c(colors, "lightgray"),
                         nm = c(lvls, ifelse(drop_na, NA_character_, bg_name)))
    }
  }

  return(list(data, dims, colors))
}


#' Get limits for pixel coordinates
#'
#' @param x a tibble with spot coordinates
#' @param coords_columns a character vector specifying the columns
#' names for the spot coordinate vectors
#'
#' @import dplyr
#'
#' @return a tibble with limits
#'
#' @noRd
.get_limits <- function (
    x,
    coords_columns
) {

  # Get min/max values for the coordinates
  x |>
    group_by(!! sym("sampleID")) |>
    summarize(
      x_start = min(!! sym(coords_columns[1])),
      y_start = min(!! sym(coords_columns[2])),
      full_width = max(!! sym(coords_columns[1])),
      full_height = max(!! sym(coords_columns[2]))
    )
}

#' Get plot dimensions
#'
#' @param dims either NULL or a tibble with the plots dimensions
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#' @import dplyr
#' @importFrom tibble tibble
#'
#' @return a tibble with plot dimension information
#'
#' @noRd
.get_dims <- function (
    dims
) {

  # Set global variables to NULL
  full_width <- full_height <- sampleID <- x_start <- y_start <- NULL

  if (!any(class(dims) %in% "tbl")) abort(glue("Invalid class of dims object '{class(dims)[1]}'"))
  if (!any(c("full_width", "full_height") %in% colnames(dims))) abort(glue("Couldn't find 'full_width' or 'full_height' in dims object."))
  dims <- dims |>
    select(full_width, full_height, sampleID) |>
    mutate(x_start = 0, y_start = 0) |>
    select(x_start, y_start, full_width, full_height, sampleID)

  return(dims)

}

#' Get feature limits
#'
#' @param data a list of tibbles containing coordinates and feature values
#' @param coords_columns character vector specifying the column names of the coordinates
#' @param scale one of "free" or "shared". "free" will fix the value limits for each
#' feature separately and "shared" will fix the value limits for all plots to be identical
#'
#' @import dplyr
#'
#' @return a list of tibbles with feature value limits
#'
#' @noRd
.get_feature_limits <- function (
  data,
  coords_columns,
  scale = c("shared", "free")
) {

  # Set global variables to NULL
  barcode <- sampleID <- NULL

  if (scale == "free") {
    feature_limits <- lapply(data, function(x) {
      y <- x |>
        select(-barcode, -all_of(coords_columns), -contains("encoded_cols"), -sampleID)
      y <- sapply(y, range) |> as_tibble()
      return(y)
    })
  } else if (scale == "shared") {
    y <- do.call(bind_rows, data) |>
      select(-barcode, -all_of(coords_columns), -contains("encoded_cols"), -sampleID)
    y <- sapply(y, range) |> as_tibble()
    feature_limits <- setNames(lapply(names(data), function(nm) {
      y
    }), nm = names(data))
  }

  return(feature_limits)
}

#' Blend values
#'
#' @param data A list of tibbles containing coordinates and feature values
#' @param features A character vector of feature names
#' @param blend_order An integer vector specifying the order to blend values by.
#' This vector can only take values 1, 2 or 3 and at least two values needs to
#' be provided.
#' @param feature_limits A list of tibbles containing information about the
#' plot dimensions
#' @param scale_alpha A logical specifying if transparency should be added to spot colors
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#' @importFrom scales rescale
#' @import dplyr
#'
#' @return a list of tibbles similar to input data but in which the feature
#' columns have been replaced by a color vector with blended colors called
#' \code{encoded_cols}
#'
#' @return a list of tibbles containing coordinates and a color vector
#'
#' @noRd
.color_blender <- function (
    data,
    features,
    blend_order,
    feature_limits,
    scale_alpha = FALSE
) {

  if (!length(features) %in% 2:3) abort(glue("Feature blending only works with 2 or 3 features,",
                                             " but {length(features)} features were provided."))
  if (length(blend_order) > length(features)) blend_order <- blend_order[1:2]
  if (!length(features) == length(blend_order)) abort("'features' and 'blend_order' need to have the same length")
  if (!all(blend_order %in% 1:3)) abort(glue("'blend_order' should be a vector with values 1, 2 or 3."))
  if (sum(duplicated(blend_order)) > 0) abort(glue("'blend_order' cannot have repeated values."))

  data <- setNames(lapply(names(data), function(nm) {
    x <- data[[nm]]
    feature_values <- x |>
      select(all_of(features)) |>
      mutate(across(everything(), ~ rescale(.x, from = c(feature_limits[[nm]][1, cur_column(), drop = TRUE],
                                                                 feature_limits[[nm]][2, cur_column(), drop = TRUE]),
                                                    to = c(0, 1))))
    mat <- matrix(0, ncol = 3, nrow = nrow(feature_values))
    mat[, blend_order[1:length(features)]] <- (feature_values*255) |> as.matrix()
    if (!requireNamespace("farver", quietly = TRUE)) {
      install.packages("farver")
    }
    encoded_cols <- mat |>
      farver::encode_colour(from = "rgb")
    x <- bind_cols(x, encoded_cols = encoded_cols)
    x <- x |> select(-all_of(features))
    # Set opacity
    if (!scale_alpha) {
      x$alpha_values <- 1
    } else {
      alpha_values <- feature_values |>
        as.matrix() |>
        apply(1, max)
      x$alpha <- alpha_values
    }
    return(x)
  }), nm = names(data))

  return(data)
}

#' Validate and get images for plotting
#'
#' @param object a \code{Seurat} object
#' @param st_object a \code{Staffli} object
#' @param image_use string specifying images to use
#' @param section_number an integer value specifying a section number to subset
#' images by
#' @param column_name a string specifying a column in the meta data slot of the
#' Seurat object containing a factor with labels to color spots by
#' @param split_labels logical specifying if labels should be split
#'
#' @noRd
.get_images <- function (
    object,
    st_object,
    image_use = NULL,
    section_number = NULL,
    column_name = NULL,
    split_labels = FALSE
) {
  # Get images if image_use is provided
  if (!is.null(image_use)) {
    # validate image_use
    if (!is.character(image_use) & (length(image_use) == 1))
      abort(glue("'image_use' is invalid. Expected a character of length 1."))
    image_use <- match.arg(image_use, choices = c("raw", "transformed"))
    images <- st_object@rasterlists[[image_use]]
    if (image_use == "raw") {
      if (is.null(images)) abort("Images have not yet been loaded. Did you run LoadImages()?")
    }
    if (image_use == "transformed") {
      if (is.null(images)) abort("Images have not yet been transformed.")
    }
    images <- setNames(images, paste0(seq_along(images)))
    if (!is.null(section_number)) {
      images <- images[section_number]
    }
    if (split_labels) {
      lvls <- levels(object |> pull(all_of(column_name)))
      images <- setNames(lapply(lvls, function(lvl) {
        images[[1]]
      }), nm = lvls)
    }
  }
  return(images)
}

#' Crop images
#'
#' @param crop_area a numeric vector of length 4 specifying a rectangle to crop data by
#' @param images a list of `raster` objects (images)
#'
#' @return a  list of cropped `raster` objects
#'
#' @noRd
.crop_images <- function (
    crop_area,
    images
) {
  if (!is.null(crop_area)) {
    if (!is.numeric(crop_area)) abort(glue("Invalid class '{class(crop_area)}' for 'crop_area', expected 'numeric'"))
    if (length(crop_area) != 4) abort(glue("Invalid length for 'crop_area', expected a 'numeric' vector of length 4"))
    if (!all(between(x = crop_area, left = 0, right = 1))) abort("'crop_area' can only take values between 0-1")
    images <- setNames(lapply(images, function(im) {
      image_dims <- dim(im)
      x_start <- round(crop_area[1]*image_dims[2])
      x_end <- round(crop_area[3]*image_dims[2])
      y_start <- round(crop_area[2]*image_dims[1])
      y_end <- round(crop_area[4]*image_dims[1])
      im <- im[y_start:y_end, x_start:x_end]
      return(im)
    }), nm = names(images))
  }
  return(images)
}


#' Crop dims
#'
#' @param dims tibble specifying original image dimensions
#' @param crop_area a numeric vector of length 4 specifying a rectangle to crop data by
#' @param data list of `tibble` objects with spot coordinates
#' @param coords_columns character vector specifying column names for spot coordinates
#'
#' @import dplyr
#'
#' @return a tibble with modified image dimensions \code{dims} and list of `tibble` objects
#' with cropped spot coordinates
#'
#' @noRd
.crop_dims <- function (
    dims,
    crop_area,
    data,
    coords_columns
) {

  # Set global variables to NULL
  full_width <- full_height <- NULL

  dims <- dims |>
    mutate(x_start = full_width*crop_area[1],
           y_start = full_height*crop_area[2],
           full_width = full_width*crop_area[3],
           full_height = full_height*crop_area[4])
  # Filter data
  data <- setNames(lapply(names(data), function(nm) {
    x <- data[[nm]]
    x |>
      filter(if_any(
        all_of(coords_columns[1]),
        ~ between(x = .x,
                  left = dims[dims$sampleID == nm, "x_start", drop = TRUE],
                  right = dims[dims$sampleID == nm, "full_width", drop = TRUE])
      )) |>
      filter(if_any(
        all_of(coords_columns[2]),
        ~ between(x = .x,
                  left = dims[dims$sampleID == nm, "y_start", drop = TRUE],
                  right = dims[dims$sampleID == nm, "full_height", drop = TRUE])
      ))
  }), nm = names(data))
  return(list(dims, data))
}

#' Inject images to list of `ggplot` objects
#'
#' @param image_use string specifying image type to use
#' @param features a character vector of features names
#' @param arrange_features a string specifying how the `patchwork` should be arranged
#' @param wrapped_plots list of `ggplot` objects
#' @param images list of `raster` images
#' @param ncol number of columns to arrange plots by to produce a `patchwork`
#' @param return_plot_list logical specifying if the result should be returned
#' as a `patchwork` or a list of `ggplot` objects
#' @param blend logical specifying if blend is activated
#'
#' @importFrom patchwork inset_element wrap_plots
#' @importFrom rlang %||% warn
#' @importFrom ggplot2 is.ggplot
#'
#' @return `patchwork` or a list of `ggplot` objects
#'
#' @noRd
.inject_images <- function (
    image_use,
    features,
    arrange_features,
    wrapped_plots,
    images,
    ncol,
    return_plot_list,
    blend = FALSE
) {
  if (!is.null(image_use)) {
    # Inject images
    wrapped_plots <- setNames(lapply(names(wrapped_plots), function(nm) {
      plot <- wrapped_plots[[nm]]
      if (!is.ggplot(plot)) {
        plot <- setNames(lapply(plot, function(p) {
          p <- p + inset_element(
            images[[nm]],
            left = 0,
            bottom = 0,
            right = 1,
            top = 1,
            on_top = FALSE
          )
          return(p)
        }), nm = names(plot))
      } else {
        plot <- plot + inset_element(
          images[[nm]],
          left = 0,
          bottom = 0,
          right = 1,
          top = 1,
          on_top = FALSE
        )
      }
      return(plot)
    }), nm = names(wrapped_plots))
  }
  return(wrapped_plots)
}


#' Arrange plots
#'
#' @param wrapped_plots A list of `ggplot` objects
#' @param features A character vector of feature names
#' @param blend A logical specifying if colors are blended or not
#' @param arrange_features A string specifying how features should
#' be arranged
#' @param ncol An integer specifying the number of columns to arrange the
#' `patchwork` by
#'
#' @importFrom patchwork wrap_plots
#'
#' @noRd
.arrange_plots <- function (
    wrapped_plots,
    features = NULL,
    blend = FALSE,
    arrange_features = "col",
    ncol = NULL
) {

  # Set default value for nSamples
  nSamples <- 1

  # When blend = FALSE, wrapped_plots will be a nested list where the
  # first layer contains samples and the second layer contains features.
  # In order to create the final patchwork, the nested list is first unlisted
  # and reordered
  if (!blend) {
    nSamples <- length(wrapped_plots)
    wrapped_plots <- Reduce(c, wrapped_plots)
    # rearrange plots
    new_order <- match(names(wrapped_plots), features) |> order()
    wrapped_plots <- wrapped_plots[new_order]
  }
  # If there's only 1 feature 1 sample or if blend is active,
  # a patchwork will be created that ignores arrange_features
  if ((length(features) == 1) | (nSamples == 1) | blend) {
    ncol <- ncol %||% {
      #if (!blend) {
      #  ceiling(sqrt(length(wrapped_plots)))
      #} else {
      #  ceiling(sqrt(length(images)))
      #}
      ceiling(sqrt(length(wrapped_plots)))
    }
    wrapped_plots <- wrap_plots(wrapped_plots, ncol = ncol)
  } else {
    # arrange_features is only active when there are multiple samples
    # and features so that a grid can be created
    if (arrange_features == "col") {
      wrapped_plots <- wrap_plots(wrapped_plots, ncol = length(features), byrow = FALSE)
    } else {
      wrapped_plots <- wrap_plots(wrapped_plots, nrow = length(features), byrow = TRUE)
    }
  }
}

#' Get coords_columns
#'
#' util function to get appropriate column names for spot coordinates
#'
#' @param image_use string specifying image type to use
#' @param coords_use string specifying coordinate type to use
#'
#' @return a character vector with column names for spot coordinates
#'
#' @noRd
.get_coords_column <- function (
    image_use,
    coords_use
) {
  if (!is.null(image_use)) {
    coords_use <- image_use
  }
  if (coords_use == "raw") {
    coords_columns <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
  } else if (coords_use == "transformed") {
    coords_columns <- c("pxl_col_in_fullres_transformed", "pxl_row_in_fullres_transformed")
  }
  return(coords_columns)
}

#' Generate colors
#'
#' @param n integer specifying the number of colors to generate
#'
#' @importFrom grDevices hcl
#'
#' @return a character vector with color ids
#'
#' @noRd
.gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
