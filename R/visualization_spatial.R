#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param pt_size numeric value specifying the point size passed to \code{geom_point}
#' @param pt_alpha numeric value between 0 and 1 specifying the point opacity passed
#' to \code{geom_point}. A value of 0 will make the points completely transparent
#' and a value of 1 will make the points completely opaque.
#' @param pt_stroke numeric specifying the point stroke width
#' @param label_by character of length 1 providing a column name in \code{object} with
#' labels that can be used to provide a title for each subplot. This column should have
#' 1 label per tissue section. This can be useful when you need to provide more detailed
#' information about your tissue sections.
#' @param ncol integer value specifying the number of columns in the output patchwork.
#' This parameter will only have an effect when the number of features provided is 1.
#' Otherwise, the patchwork will be arranged based on the \code{arrange_features} parameter.
#' @param colors a character vector of colors to use for the color scale. The colors should
#' preferably consist of a set of colors from a scientific color palette designed for sequential
#' data. Some useful palettes are available in the \code{RColorBrewer}, \code{viridis} and
#' \code{scico} R packages.
#' @param scale a character vector of length 1 specifying one of "shared" or "free" which will
#' determine how the color bars are structured. If scale is set to "shared", the color bars for
#' feature values will be shared across samples. If scale is set to "free", the color bars
#' will be independent.
#' @param dims a tibble with information about the tissue sections. This information is used to
#' determine the limits of the plot area. If \code{dims} is not provided, the limits will be
#' computed directly from the spatial coordinates provided in \code{object}.
#' @param coords_columns a character vector of length 2 specifying the names of the columns in
#' \code{object} holding the spatial coordinates
#' @param blend a logical specifying whether blending should be used. See color blending for more information.
#' @param blend_order an integer vector of length 2-3 specifying the order to blend features by. Only
#' active when \code{blend = TRUE}. See color blending for more information.
#'
#' @importFrom patchwork wrap_plots
#' @importFrom dplyr select group_by group_split all_of pull
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
  pt_size = 1,
  pt_alpha = 1,
  pt_stroke = 0,
  label_by = NULL,
  ncol = NULL,
  colors = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
  scale = c("shared", "free"),
  arrange_features = c("col", "row"),
  dims = NULL,
  coords_columns = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  blend = FALSE,
  blend_order = 1:3
) {

  # Check data
  .prep_data_for_plotting(object, colors, label_by, scale, arrange_features, coords_columns)

  # Expand colors if length is 1
  if (length(colors) == 1) {
    colors <- c("lightgray", colors)
  }

  # get features
  features <- object |>
    select(-barcode, -pxl_col_in_fullres, -pxl_row_in_fullres, -sampleID, -all_of(label_by)) |>
    colnames()

  # Split data by sampleID
  data <- object |>
    group_by(sampleID) |>
    group_split() |>
    setNames(nm = unique(object$sampleID))

  # get image dimensions
  dims <- .get_dims(data, dims)

  # get feature limits
  feature_limits <- .get_feature_limits(data, coords_columns, scale = ifelse(blend, "shared", scale))

  # add blend colors if blend=TRUE
  if (blend) {
    data <- .color_blender(data, features, blend_order, feature_limits)
    extreme_colors <- encode_colour(diag(ncol = 3, nrow = 3)*255, from = "rgb")
    features <- features[blend_order[1:length(features)]]
  }

  # Plot features on spatial coordinates for each sample
  sample_plots <- setNames(lapply(names(data), function(nm) {

    # Get data for plotting
    gg <- data[[nm]]

    # Set label if available
    if (!is.null(label_by)) {
      cur_label <- unique(gg |> pull(all_of(label_by)))
    } else {
      cur_label <- paste0("section ", nm)
    }

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
          cur_label = cur_label,
          coords_columns = coords_columns
        )
      })
      # Arrange features by col or row
      p <- .arrange_plots(feature_plots, by.col = arrange_features == "col")
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
        cur_label = cur_label,
        coords_columns = coords_columns
      )
    }
    return(p)
  }), nm = names(data))

  # Create final patchwork
  if (length(features) == 1 | length(data) == 1 | blend) {
    ncol <- ncol %||% ceiling(sqrt(length(features)))
    wrapped_plots <- wrap_plots(sample_plots, ncol = ncol)
  } else {
    if (!is.null(ncol)) warn("'ncol' will not be used when more than 1 feature is provided")
    wrapped_plots <- .arrange_plots(sample_plots, by.col = arrange_features != "col")
  }

  return(wrapped_plots)
}

#' Map numeric features in 2D using a Seurat object
#'
#' @param features a character vector of features to plot. These features need to be
#' fetchable with \code{link{FetchData}}
#' @param override_plot_dims a logical specifying whether the image dimensions should
#' be used to define the plot area. Setting \code{override_plot_dims} can be useful
#' in situations where the tissue section only covers a small fraction of the capture
#' area, which will create a lot of white space in the plots.
#' @param min_cutoff,max_cutoff a numeric value between 0-1 specifying either a lower
#' (\code{min_cutoff}) or upper (\code{max_cutoff}) limit for the data using \code{\link{quantile}}.
#' These arguments can be useful to make sure that the color map doesn't get dominated by outliers.
#'
#' @importFrom dplyr bind_cols
#' @importFrom tibble as_tibble
#' @importFrom Seurat FetchData
#'
#' @rdname visualize-features
#'
#' @examples
#'
#' library(STUtility2)
#'
#' # Load example Visium data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "STUtility2"))
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "STUtility2"))
#' se_merged <- MergeSTData(se_mbrain, se_mcolon)
#'
#' \dontrun{
#'
#' # Select features
#' selected_features <- c("Clu", "Slc6a3", "Vip")
#'
#' # Plot selected features with custom colors
#' MapFeatures(se_merged, features = selected_features, colors = viridis::magma(n = 11, direction = -1))
#'
#' # Plot selected features with color bars scaled individually for each feature and sample
#' MapFeatures(se_merged, features = selected_features, scale = "free")
#'
#' # Plot selected features and add custom labels for subplots
#' se_merged$sample_id <- ifelse(GetStaffli(se_merged)@meta_data$sampleID == "1", "mouse brain", "mouse colon")
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
#' MapFeatures(se_mbrain, features = selected_feature, pt_size = 2, colors = viridis::viridis(n = 11)) &
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
    pt_size = 1,
    pt_alpha = 1,
    pt_stroke = 0,
    label_by = NULL,
    ncol = NULL,
    colors = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
    scale = c("shared", "free"),
    arrange_features = c("col", "row"),
    blend = FALSE,
    blend_order = 1:3,
    override_plot_dims = FALSE,
    max_cutoff = NULL,
    min_cutoff = NULL
) {

  # Match args
  scale <- match.arg(scale, choices = c("shared", "free"))
  arrange_features <- match.arg(arrange_features, choices = c("col", "row"))

  # Check Seurat object
  .check_seurat_object(object)

  # fetch data from Seurat object
  data_use <- GetStaffli(object)@meta_data |>
    bind_cols(FetchData(object, vars = features) |> as_tibble())

  # Add label_by column if present
  if (!is.null(label_by)) {
    data_use <- data_use |>
      bind_cols(FetchData(object, vars = label_by))
  }

  # rescale data if max/min_cuttoffs have been set
  if (!is.null(min_cutoff) | !is.null(max_cutoff)) {
    min_cutoff <- min_cutoff %||% 0
    max_cutoff <- max_cutoff %||% 1

    # Check cutoffs
    if (!all(between(x = c(min_cutoff, max_cutoff), left = 0, right = 1))) abort("min/max cutoffs cannot be outside the range 0-1")

    # Cut data
    data_use <- data_use |>
      mutate(across(
        all_of(features),
        ~ case_when(.x < quantile(.x, probs = min_cutoff) ~ quantile(.x, probs = min_cutoff),
                    .x > quantile(.x, probs = max_cutoff) ~ quantile(.x, probs = max_cutoff),
                    TRUE ~ .x)
      ))
  }

  # generate plots
  wrapped_plots <- MapFeatures(
    object = data_use,
    pt_size = pt_size,
    pt_alpha = pt_alpha,
    pt_stroke = pt_stroke,
    label_by = label_by,
    ncol = ncol,
    colors = colors,
    scale = scale,
    arrange_features = arrange_features,
    dims = switch(override_plot_dims + 1, GetStaffli(object)@image_info, NULL),
    coords_columns = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    blend = blend,
    blend_order = blend_order
  )

  return(wrapped_plots)

}


#' @param pt_size numeric value specifying the point size passed to \code{geom_point}
#' @param pt_alpha numeric value between 0 and 1 specifying the point opacity passed
#' to \code{geom_point}. A value of 0 will make the points completely transparent
#' and a value of 1 will make the points completely opaque.
#' @param pt_stroke numeric specifying the point stroke width
#' @param label_by character of length 1 providing a column name in \code{object} with
#' labels that can be used to provide a title for each subplot. This column should have
#' 1 label per tissue section. This can be useful when you need to provide more detailed
#' information about your tissue sections.
#' @param ncol integer value specifying the number of columns in the output patchwork.
#' @param colors a character vector of colors to use for the color scale. The number of
#' colors should match the number of labels present.
#' @param dims a tibble with information about the tissue sections. This information is used to
#' determine the limits of the plot area. If \code{dims} is not provided, the limits will be
#' computed directly from the spatial coordinates provided in \code{object}.
#' @param coords_columns a character vector of length 2 specifying the names of the columns in
#' \code{object} holding the spatial coordinates
#'
#' @importFrom patchwork wrap_plots
#' @importFrom dplyr select group_by group_split all_of pull
#' @importFrom zeallot %<-%
#'
#' @rdname visualize-labels
#'
#' @return A patchwork object
#'
#' @export
#'
MapLabels.default <- function (
    object,
    pt_size = 1,
    pt_alpha = 1,
    pt_stroke = 0,
    section_number = NULL,
    label_by = NULL,
    split_labels = FALSE,
    ncol = NULL,
    colors = NULL,
    dims = NULL,
    coords_columns = c("pxl_col_in_fullres", "pxl_row_in_fullres")
) {

  # Check data
  .prep_data_for_plotting(object = object, label_by = label_by, coords_columns = coords_columns)

  # get features
  label <- object |>
    select(-barcode, -pxl_col_in_fullres, -pxl_row_in_fullres, -sampleID, -all_of(label_by)) |>
    colnames()

  # Convert label column to factor
  object <- object |>
    mutate(across(all_of(label), ~ factor(.x)))

  # Split data by sampleID
  data <- object |>
    group_by(sampleID) |>
    group_split() |>
    setNames(nm = unique(object$sampleID))

  # Split data by label if split_labels = TRUE
  if (split_labels) {
    c(data, dims, colors) %<-% .split_data_by_label(data, dims, section_number, label, colors)
  } else if (!is.null(section_number)) {
    data <- data[section_number]
  }

  # Obtain colors
  colors <- colors %||% .gg_color_hue(length(unique(object |> select(all_of(label)) |> pull(all_of(label)))))

  # get image dimensions
  dims <- .get_dims(data, dims)

  # Plot features on spatial coordinates for each sample
  sample_plots <- setNames(lapply(names(data), function(nm) {

    # Get data for plotting
    gg <- data[[nm]]

    # Set label if available
    if (!is.null(label_by)) {
      cur_label <- unique(gg |> pull(all_of(label_by)))
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
      cur_label = cur_label
    )

    return(p)
  }), nm = names(data))

  # Create final patchwork
  ncol <- ncol %||% ceiling(sqrt(length(data)))
  wrapped_plots <- wrap_plots(sample_plots, ncol = ncol)

  return(wrapped_plots)
}


#' @param column_name a string specifying a meta data column holding the categorical
#' feature vector
#' @param override_plot_dims a logical specifying whether the image dimensions should
#' be used to define the plot area. Setting \code{override_plot_dims} can be useful
#' in situations where the tissue section only covers a small fraction of the capture
#' area, which will create a lot of white space in the plots.
#'
#' @rdname visualize-labels
#'
#' @examples
#'
#' library(STUtility2)
#' library(ggplot2)
#'
#' # Load Seurat object
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "STUtility2"))
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
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2", pt_size = 0.5, section_number = 1, split_labels = TRUE, ncol = 4)
#'
#' # Combine plots with different labels
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2") | MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.3")
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
#' cols <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2", pt_size = 2, colors = cols)
#'
#' # Factor are to used to determine the color order. If you change the
#' # levels of your label of interest, the labels and colors will change order
#' se_mbrain$Spatial_snn_res.0.2 <- factor(se_mbrain$Spatial_snn_res.0.2, levels = sample(paste0(0:7), size = 8))
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2", pt_size = 2, colors = cols)
#'
#' # Control what group label colors by naming the color vector
#' # this way you can make sure that each group gets a desired color
#' # regardless of the factor levels
#' cols <- setNames(cols, nm = paste0(0:7))
#' cols
#' MapLabels(se_mbrain, column_name = "Spatial_snn_res.0.2", pt_size = 2, colors = cols)
#'
#' @export
#'
MapLabels.Seurat <- function (
    object,
    column_name,
    pt_size = 1,
    pt_alpha = 1,
    pt_stroke = 0,
    section_number = NULL,
    label_by = NULL,
    split_labels = FALSE,
    ncol = NULL,
    colors = NULL,
    override_plot_dims = FALSE
) {

  # Check Seurat object
  .check_seurat_object(object)

  # fetch data from Seurat object
  data_use <- GetStaffli(object)@meta_data |>
    bind_cols(FetchData(object, vars = column_name) |> as_tibble())

  # Add label_by column if present
  if (!is.null(label_by)) {
    data_use <- data_use |>
      bind_cols(FetchData(object, vars = label_by))
  }

  # generate plots
  wrapped_plots <- MapLabels(
    object = data_use,
    pt_size = pt_size,
    pt_alpha = pt_alpha,
    pt_stroke = pt_stroke,
    section_number = section_number,
    label_by = label_by,
    split_labels = split_labels,
    ncol = ncol,
    colors = colors,
    dims = switch(override_plot_dims + 1, GetStaffli(object)@image_info, NULL),
    coords_columns = c("pxl_col_in_fullres", "pxl_row_in_fullres")
  )

  return(wrapped_plots)

}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plotting utils
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Plot numeric features in 2D
#'
#' @param gg tibble with spatial coordinates and feature values
#' @param nm sample ID
#' @param ftr feature name
#' @param feature_limits list of tibbles containing information about
#' the feature value range
#' @param dims tibble containing information about the dimensions
#' of the plotting area
#' @param all_features a character vector with all features, only used
#' for blending features
#' @param extreme_colors a character vector with the hex colors, only
#' used for blending features
#' @param pt_size point size passed to geom_point
#' @param pt_alpha point opacity ranging from 0 to 1 passed to geom_point.
#' 0 = fully transparent, 1 = fully opaque
#' @param pt-stroke point stroke width
#'
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous scale_y_reverse theme_void theme scale_color_gradientn labs coord_fixed aes
#'
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
    cur_label = NULL,
    coords_columns
) {

  # Check if encoded colors are present after setting blend = TRUE
  encoded_cols_present <- "encoded_cols" %in% colnames(gg)
  color_vec <- switch(encoded_cols_present + 1, NULL, gg |> pull(encoded_cols))

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
          alpha = pt_alpha,
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
          alpha = pt_alpha,
          shape = 21,
          stroke = pt_stroke
        )
      }
    } +
    # Set plot dimensions (reverse y axis)
    scale_x_continuous(limits = c(dims[dims$sampleID == nm, "x_start", drop = TRUE],
                                  dims[dims$sampleID == nm, "full_width", drop = TRUE]),
                       expand = c(0, 0)) +
    scale_y_reverse(limits = c(dims[dims$sampleID == nm, "full_height", drop = TRUE],
                               dims[dims$sampleID == nm, "y_start", drop = TRUE]),
                    expand = c(0, 0)) +
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
                              limits = c(feature_limits[[nm]][1, ftr, drop = TRUE],
                                         feature_limits[[nm]][2, ftr, drop = TRUE]))
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
#' @param dims tibble containing information about the dimensions
#' of the plotting area
#' @param all_features a character vector with all features, only used
#' @param pt_size point size passed to geom_point
#' @param pt_alpha point opacity ranging from 0 to 1 passed to geom_point.
#' 0 = fully transparent, 1 = fully opaque
#' @param pt_stroke point stroke width
#'
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous scale_y_reverse theme_void theme scale_color_gradientn labs coord_fixed aes
#'
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
    cur_label
) {

  # Create input data.frame for ggplot
  gg <- gg |>
    select(all_of(c(coords_columns, lbl))) |>
    setNames(nm = c(coords_columns, "variable"))

  # Rearrange colors by factor level
  if (!is.null(names(colors))) {
    colors <- colors[levels(gg$variable)]
  } else {
    colors <- setNames(colors, levels(gg$variable))
  }

  # Check colors
  if (length(colors) != length(levels(gg$variable))) abort(glue("The number of colors ({length(colors)})",
                                                           " does not match the number of labels ({length(levels(gg$variable))})."))


  # Draw plot
  p <-
    ggplot() +
    geom_point(data = gg, aes(
        x = !! sym(coords_columns[1]),
        y = !! sym(coords_columns[2]),
        fill = variable),
      size = pt_size,
      alpha = pt_alpha,
      shape = 21,
      stroke = pt_stroke
    ) +
    # Set plot dimensions (reverse y axis)
    scale_x_continuous(limits = c(dims[dims$sampleID == nm, "x_start", drop = TRUE],
                                  dims[dims$sampleID == nm, "full_width", drop = TRUE]),
                       expand = c(0, 0)) +
    scale_y_reverse(limits = c(dims[dims$sampleID == nm, "full_height", drop = TRUE],
                               dims[dims$sampleID == nm, "y_start", drop = TRUE]),
                    expand = c(0, 0)) +
    # Add themes
    theme_void() +
    theme(legend.position = "top",
          legend.justification = "left",
          legend.margin = margin(t = 5, 0, 10, 0),
          plot.margin = margin(0, 10, 20, 10),
          legend.title = element_text(vjust = 0.8)) +
    # Add colors
    scale_fill_manual(values = colors) +
    #scale_fill_manual(values = c("1" = "red", "0" = "blue", "2" = "grey", "3" = "orange",
    #                             "4" = "green", "5" = "yellow", "6" = "lightgray", "7" = "magenta")) +
    # Create a title
    labs(title = ifelse(!is.null(cur_label), cur_label, NA)) +
    # Fix coordinates so that plot cannot be stretched
    coord_fixed()

  return(p)
}


#' Check input for compatibility
#'
#' @param object a tibble with spatial coordinates and feature values
#' @param colors a character vector of colors
#' @param scale one of "free" or "shared". "free" will fix the value limits for each
#' feature separately and "shared" will fix the value limits for all plots to be identical
#' @param arrange_features one of "row" or "col". "col" will put the features in columns
#' and samples in rows and "row" will transpose the arrangement
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#'
.prep_data_for_plotting <- function (
    object,
    colors,
    label_by,
    scale,
    arrange_features,
    coords_columns
) {

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
      tidyr::unite(col = "merged_cols", all_of(c("sampleID", label_by))) |>
      count(merged_cols)
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
    checks <- object |>
      select(-barcode, -all_of(coords_columns), -sampleID, -all_of(label_by)) |>
      sapply(function(x) {
        is.character(x) | is.factor(x)
      })
    if (length(checks) > 1) abort("Only 1 label column is allowed.")
    if (!checks) abort(glue("Label columns has to be a character/factor."))
  }
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
#' @param section_number an integer specifying the section to keep for splitting
#' @param label a string specifying the name of the label column
#' @param colors a character vector of colors that should match the number of
#' groups in the label column
#'
#' @importFrom dplyr mutate across all_of case_when filter
#'
#' @return a list of tibbles holding spot coordinates, a tibble with information
#' about plot dimensions and a character vector of colors stored together in a list
#'
.split_data_by_label <- function (
  data,
  dims,
  section_number,
  label,
  colors
) {
  section_number <- section_number %||% {
    warn("No section_number selected. Selecting the section 1.")
    1
  }
  if (!is.numeric(section_number)) abort(glue("Invalid class '{class(section_number)}' for",
                                              " {cli::col_br_green('section_number')}, expected 'numeric'"))
  if (!section_number %in% seq_along(data)) abort(glue("'section_number' = {section_number} is out of range. ",
                                                       "Select one of {paste(seq_along(data), collapse = ', ')}"))
  data <- data[section_number]
  data <- lapply(data, function(x) {
    x |>
      mutate(across(all_of(label), ~ fct_drop(.x)))
  })

  # restructure data
  lvls <- data[[1]] |> pull(all_of(label)) |> levels()
  data <- lapply(lvls, function(lvl) {
    data[[1]] |> mutate(across(all_of(label), ~ case_when(.x != lvl ~ "background",
                                                          TRUE ~ lvl))) |>
      mutate(across(all_of(label), ~ factor(.x, levels = c(lvls, "background"))))
  }) |>
    setNames(nm = lvls)

  # restructure dims
  dims <- do.call(bind_rows, lapply(lvls, function(lvl) {
    dims |> filter(sampleID == section_number) |> mutate(sampleID = lvl)
  }))

  # Overwrite colors
  colors <- colors %||% setNames(c(.gg_color_hue(length(lvls)), "lightgray"), nm = c(lvls, "background"))
  if (length(colors) == length(lvls)) {
    if (!is.null(names(colors))) {
      colors <- setNames(c(colors, "lightgray"), nm = c(names(colors), "background"))
    } else {
      colors <- setNames(c(colors, "lightgray"), nm = c(lvls, "background"))
    }
  }

  return(list(data, dims, colors))
}


#' Get plot dimensions
#'
#' @param data a list of tibbles containing coordinates and feature values
#' @param dims either NULL or a tibble from \code{\link{GetSfaffli(se)@image_info}}
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#' @importFrom dplyr mutate select
#' @importFrom tibble tibble
#'
#' @return a tibble with plot dimension information
#'
.get_dims <- function (
    data,
    dims
) {

  # Get image dimensions
  if (is.null(dims)) {
    dims <- do.call(bind_rows, lapply(names(data), function(nm) {
      x <- data[[nm]]
      tibble(x_start = min(x$pxl_col_in_fullres),
             y_start = min(x$pxl_row_in_fullres),
             full_width = max(x$pxl_col_in_fullres),
             full_height = max(x$pxl_row_in_fullres),
             sampleID = nm)
    }))
  } else {
    # Validate dims
    if (!any(class(dims) %in% "tbl")) abort(glue("Invalid class of dims object '{class(dims)[1]}'"))
    if (!any(c("full_width", "full_height") %in% colnames(dims))) abort(glue("Couldn't find 'full_width' or 'full_height' in dims object."))
    dims <- dims |>
      select(full_width, full_height, sampleID) |>
      mutate(x_start = 0, y_start = 0) |>
      select(x_start, y_start, full_width, full_height, sampleID)
  }

  return(dims)

}

#' Get feature limits
#'
#' @param data a list of tibbles containing coordinates and feature values
#' @param coords_columns character vector specifying the column names of the coordinates
#' @param scale one of "free" or "shared". "free" will fix the value limits for each
#' feature separately and "shared" will fix the value limits for all plots to be identical
#'
#' @importFrom dplyr select summarize_all contains all_of
#'
#' @return a list of tibbles with feature value limits
#'
.get_feature_limits <- function (
  data,
  coords_columns,
  scale = c("shared", "free")
) {

  if (scale == "free") {
    feature_limits <- lapply(data, function(x) {
      x |>
        select(-barcode, -all_of(coords_columns), -contains("encoded_cols"), -sampleID) |>
        summarize_all(range)
    })
  } else if (scale == "shared") {
    feature_limits <- do.call(bind_rows, data) |>
      select(-barcode, -all_of(coords_columns), -contains("encoded_cols"), -sampleID) |>
      summarize_all(range)
    feature_limits <- setNames(lapply(names(data), function(nm) {
      feature_limits
    }), nm = names(data))
  }

  return(feature_limits)
}

#' Blend values
#'
#' @param data a list of tibbles containing coordinates and feature values
#' @param features a character vector of feature names
#' @param blend_order an integer vector specifying the order to blend values by.
#' This vector can only take values 1, 2 or 3 and at least two values needs to
#' be provided.
#' @param feature_limits a list of tibbles containing information about the
#' plot dimensions
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#' @importFrom scales rescale
#' @importFrom dplyr select mutate bind_cols
#' @importFrom farver encode_colour
#'
#' @return a list of tibbles similar to input data but in which the feature
#' columns have been replaced by a color vector with blended colors called
#' \code{encoded_cols}
#'
.color_blender <- function (
    data,
    features,
    blend_order,
    feature_limits
) {

  if (!length(features) == length(blend_order)) abort("'features' and 'blend_order' need to have the same length")
  if (!length(features) %in% 2:3) abort(glue("Feature blending only works with 2 or 3 features,",
                                             " but {length(features)} features were provided."))
  if (!all(blend_order %in% 1:3)) abort(glue("'blend_order' should be a vector with values 1, 2 or 3."))
  if (sum(duplicated(blend_order)) > 0) abort(glue("'blend_order' cannot have repeated values."))

  data <- setNames(lapply(names(data), function(nm) {
    x <- data[[nm]]
    feature_values <- x |>
      select(all_of(features)) |>
      mutate(across(everything(), ~ rescale(.x, from = c(feature_limits[[nm]][1, cur_column(), drop = TRUE],
                                                                 feature_limits[[nm]][2, cur_column(), drop = TRUE]),
                                                    to = c(0, 255))))
    mat <- matrix(0, ncol = 3, nrow = nrow(feature_values))
    mat[, blend_order[1:length(features)]] <- feature_values |> as.matrix()
    encoded_cols <- mat |>
      encode_colour(from = "rgb")
    x <- bind_cols(x, encoded_cols = encoded_cols)
    x <- x |> select(-all_of(features))
    return(x)
  }), nm = names(data))

  return(data)
}

#' Arrange list of plots
#'
#' @param plots list of ggplot objects
#' @param by.col should the plot be arranged in columns?
#'
#' @importFrom patchwork wrap_plots
#'
#' @return a patchwork object containing a grid of plots
#'
.arrange_plots <- function (
    plots,
    by.col = TRUE
) {

  if (!by.col) {
    wrapped_plots <- wrap_plots(plots, ncol = 1)
  } else {
    wrapped_plots <- wrap_plots(plots, nrow = 1)
  }

  return(wrapped_plots)
}


.gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
