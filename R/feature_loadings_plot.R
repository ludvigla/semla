#' Plot feature loadings for dimensional reduction data
#'
#' This function can be used to visualize the relative contribution of features (e.g. genes)
#' to dimensionality reduction vectors. The function provides three modes to draw a
#' bar plot, a dot plot or a heatmap.
#'
#' @section Select type:
#' For centered dimensionality reduction vectors, such as principal components,
#' it is best to select features with the highest and lowest loadings. For other types
#' of dimensionality reduction results, it is likely better to only select the features with
#' the highest loadings, in which case \code{type="positive"} is the appropriate choice.
#'
#' @section Select mode:
#' barplots or dotplots can be used to get detailed information about the loadings
#' for individual factors whereas heatmap is useful to summarize the loadings for
#' multiple factors. The heatmap option will override the \code{type} options and only
#' select the top features. With the heatmap mode, the values data will also be scaled
#' within each dimensionality reduction vectors to range between 0 and 1.
#'
#' @param object An object of class \code{Seurat}
#' @param dims An integer vector of dimensions to plot feature loadings for
#' @param reduction A character specifying the dimensionality reduction to use
#' @param nfeatures Number of features to show
#' @param mode Plot mode. One of "barplot", "dotplot" or "heatmap"
#' @param type Mode used to select features:
#' \itemize{
#'    \item{"positive" : select features with highest loadings}
#'    \item{"negative" : select features with lowest loadings}
#'    \item{"centered" : split selection of features to include top and bottom loadings}
#' }
#' @param fill Fill color for barplot and dotplot
#' @param color Border color for barplot and dotplot
#' @param bar_width With of barplot provided as a proportion
#' @param pt_size Size of points in dotplot
#' @param pt_stroke Width of border for points in dotplot
#' @param linetype Select a line type to be used for dotplot, e.g. "solid", "longdash" or "blank"
#' @param color_by_loadings Should the fill color of barplot or points reflect the feature loadings?
#' @param gradient_colors Colors to use for gradient if \code{color_by_loadings=TRUE}
#' @param ncols Number of columns used for final patchwork
#'
#' @import dplyr
#' @import rlang
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom patchwork wrap_plots
#'
#' @return An object of class \code{patchwork}
#'
#' @examples {
#' library(semla)
#'
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
#'
#' # Run PCA
#' se_mbrain <- se_mbrain |>
#'                 ScaleData() |>
#'                 RunPCA()
#'
#' # Plot feature loadings for PC_1 as a dotplot
#' PlotFeatureLoadings(se_mbrain, reduction = "pca", dims = 1,
#'                     mode = "dotplot", type = "centered")
#'
#' # Plot feature loadings for PC_1 and PC_2 as barplots
#' PlotFeatureLoadings(se_mbrain, reduction = "pca", dims = 1:2,
#'                     mode = "barplot", type = "centered")
#'
#' # Plot feature loadings for PC_1 and color bars by loading
#' PlotFeatureLoadings(se_mbrain, reduction = "pca", dims = 1:2,
#'                     mode = "barplot", type = "centered",
#'                     color_by_loadings = TRUE,
#'                     gradient_colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev())
#' }
#'
#' @export
PlotFeatureLoadings <- function (
  object,
  dims = 1,
  reduction = "pca",
  nfeatures = 30,
  mode = c("dotplot", "barplot", "heatmap"),
  type = c("positive", "negative", "centered"),
  fill = "lightgrey",
  color = "black",
  bar_width = 0.9,
  pt_size = 4,
  pt_stroke = 0.5,
  linetype = "dashed",
  color_by_loadings = FALSE,
  gradient_colors = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
  ncols = NULL
) {

  # Set global variables to NULL
  name <- value <- value_scaled <- NULL

  # Run checks
  mode <- match.arg(mode, choices = c("dotplot", "barplot", "heatmap"))
  type <- match.arg(type, choices = c("positive", "negative", "centered"))
  stopifnot(inherits(dims, what = c("numeric", "integer")),
            length(dims) > 0,
            inherits(nfeatures, what = c("numeric", "integer")),
            nfeatures > 1,
            inherits(fill, what = "character"),
            length(fill) == 1,
            inherits(color, what = "character"),
            length(color) == 1,
            inherits(bar_width, what = c("numeric", "integer")),
            length(bar_width) == 1,
            inherits(pt_size, what = c("numeric", "integer")),
            length(pt_size) == 1,
            inherits(pt_stroke, what = c("numeric", "integer")),
            length(pt_stroke) == 1,
            inherits(color_by_loadings, what = c("logical")),
            length(color_by_loadings) == 1,
            inherits(gradient_colors, what = "character"),
            length(gradient_colors) > 1)

  # Fetch data from Seurat object
  dimred_data_raw <- object[[reduction]]@feature.loadings[, dims, drop = FALSE]
  dimred_names <- colnames(dimred_data_raw)

  dimred_data <- dimred_data_raw |>
    as.data.frame() |>
    rownames_to_column(var = "name") |>
    as_tibble() |>
    pivot_longer(all_of(dimred_names), names_to = "dim", values_to = "value")

  if (mode != "heatmap") {

    # Select features based on type
    if (type == "centered") {
      top_n <- ceiling(nfeatures/2)
      bottom_n <- nfeatures - top_n
      dimred_data_top <- dimred_data |>
        group_by(dim) |>
        arrange(-value) |>
        slice_head(n = top_n)
      dimred_data_bottom <- dimred_data |>
        group_by(dim) |>
        arrange(value) |>
        slice_head(n = bottom_n)
      dimred_data <- bind_rows(dimred_data_top, dimred_data_bottom)
    } else {
      if (type == "positive") {
        dimred_data <- dimred_data |>
          arrange(-value)
      }
      if (type == "negative") {
        gradient_colors <- gradient_colors |> rev()
        dimred_data <- dimred_data |>
          arrange(value)
      }
      dimred_data <- dimred_data |>
        group_by(dim) |>
        slice_head(n = nfeatures)
    }

    # Rearrange values
    dimred_data_split <- dimred_data |>
      arrange(value) |>
      group_by(dim) |>
      group_split() |>
      setNames(nm = dimred_names)

    plots <- lapply(dimred_names, function(nm) {
      cur_data <- dimred_data_split[[nm]] |>
        mutate(name = factor(name, levels = name))
      p <- ggplot()
      if (mode == "barplot") {
        p <- p +
          {
            if (color_by_loadings) {
              geom_col(data = cur_data,
                       aes(name, value, fill = value), color = color, width = bar_width)
            } else {
              geom_col(data = cur_data,
                       aes(name, value), fill = fill, color = color, width = bar_width)
            }
          }
      }
      if (mode == "dotplot") {
        p <- p  +
          geom_segment(data = cur_data,
                       aes(x = name, xend = name, y = 0, yend = value),
                       linetype = linetype) +
          {
            if (color_by_loadings) {
              geom_point(data = cur_data,
                         aes(name, value, fill = value),
                         shape = 21, size = pt_size, color = color, stroke = pt_stroke)
            } else {
              geom_point(data = cur_data,
                         aes(name, value),
                         shape = 21, fill = fill, size = pt_size, color = color, stroke = pt_stroke)
            }
          }
      }
      p <- p +
        coord_flip() +
        theme_minimal() +
        labs(y = "Feature loading", x = "Feature name", title = nm) +
        theme(axis.text = element_text(colour = "black"))
      return(p)
    })

    ncols <- ncols %||% ceiling(sqrt(length(plots)))
    p <- wrap_plots(plots, ncol = ncols)

    # Add color scale if data should be colored by loadings
    if (color_by_loadings) {
      p <- p &
        scale_fill_gradientn(colours = gradient_colors) &
        labs(fill = "Feature\nloading", y = "")
    }
  } else {
    stopifnot(length(dims) > 1)
    
    # Convert dim to factor
    dimred_data$dim <- factor(dimred_data$dim, levels = dimred_data$dim |> unique())
    
    top_features <- dimred_data |>
      group_by(dim) |>
      arrange(-value) |>
      slice_head(n = nfeatures) |>
      pull(name) |>
      unique() |>
      rev()
    
    # Subset dimred data
    dimred_data_subset <- dimred_data |>
      group_by(dim) |>
      mutate(value_scaled = scales::rescale(value)) |>
      filter(name %in% top_features) |>
      mutate(name = factor(name, levels = top_features))

    p <- ggplot(dimred_data_subset, aes(dim, name, fill = value_scaled)) +
      geom_tile() +
      scale_fill_gradientn(colours = gradient_colors) +
      theme(panel.background = element_blank(), 
            axis.text.y = element_text(colour = "black"), 
            axis.text.x = element_text(colour = "black", angle = 60, hjust = 1)) +
      labs(x = "", y = "Feature name", fill = "scaled\nFeature\nLoading")
  }

  return(p)
}
