#' @include generics.R
#' @include checks.R
#' @include spatial_utils.R
#'
NULL


#' @section default method:
#' Takes a \code{tibble} and set of spot IDs and returns a named character vector with new labels.
#' The names of this vector corresponds to the input spot IDs.
#'
#' @section Seurat:
#' A categorical variable is selected from the \code{Seurat} object meta data
#' slot using \code{column_name}. From this column, one can specify what groups
#' to disconnect with \code{selected_groups}. If \code{selected_groups} isn't specified,
#' all groups in \code{selected_groups} will be disconnected separately.
#' The function returns a \code{Seurat} object with additional columns in the meta data
#' slot, one for each group in \code{selected_groups}. The suffix to these columns is
#' "_split", so a group in \code{selected_groups} called "tissue" will get a column
#' called "tissue_split" with new labels for each spatially disconnected region.
#'
#' @param spots A character vector with spot IDs present \code{object}
#'
#' @import dplyr
#' @import cli
#' @importFrom rlang abort
#' @importFrom glue glue
#'
#' @rdname disconnect-regions
#' @family spatial-methods
#'
#' @examples
#'
#' library(semla)
#' library(dplyr)
#' library(ggplot2)
#' library(patchwork)
#'
#' galt_spots_file <- system.file("extdata/mousecolon",
#'                                "galt_spots.csv",
#'                                 package = "semla")
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
#' # Find disconnected regions in GALT spots
#' disconnected_spot_labels <- DisconnectRegions(coords, spots)
#'
#' # Add information to coords and plot
#' gg <- coords |>
#'   mutate(galt = NA, galt_disconnected = NA)
#' gg$galt[match(spots, gg$barcode)] <- "galt"
#' gg$galt_disconnected[match(names(disconnected_spot_labels), gg$barcode)] <- disconnected_spot_labels
#'
#' p1 <- ggplot(gg, aes(x, y, color = galt))
#' p2 <- ggplot(gg, aes(x, y, color = galt_disconnected))
#' p <- p1 + p2 &
#'   geom_point() &
#'   theme_void() &
#'   coord_fixed()
#' p
#'
#' @export
#'
DisconnectRegions.default <- function (
  object,
  spots,
  verbose = TRUE,
  ...
) {

  # Set global variables to NULL
  barcode <- x <- y <- sampleID <- from <- to <- NULL

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
  if (verbose) cli_alert_info("Detecting disconnected regions for {length(spots)} spots")

  # Find border
  spatnet <- GetSpatialNetwork(object, ...)

  # Subset spatnet by selected spots
  spatnet <- lapply(spatnet, function(x) {
    x |>
      filter(from %in% spots, to %in% spots)
  })

  # load tidygraph
  if (!requireNamespace("tidygraph"))
    abort(glue("Package {cli::col_br_magenta('tidygraph')} is required. Please install it with: \n",
               "install.packages('tidygraph')"))

  # Check if graph is disconnected
  if (!.is_disconnected(spatnet)) {
    abort("Found no disconnected components.")
  }

  # Convert spatnet to tidygraph
  tidygraphs <- do.call(bind_rows, lapply(names(spatnet), function(nm) {
    x <- spatnet[[nm]]
    diconnected_graphs <- tidygraph::as_tbl_graph(x) |>
      tidygraph::to_components() |>
      lapply(as_tibble)
    sizes <- order(sapply(diconnected_graphs, nrow), decreasing = TRUE)
    diconnected_graphs <- diconnected_graphs[sizes]
    diconnected_graphs <- do.call(bind_rows, lapply(seq_along(diconnected_graphs), function(i) {
      diconnected_graphs[[i]] |>
        select(-contains("id")) |>
        bind_cols(id = paste0("S", nm, "_region", i))
    }))
    return(diconnected_graphs)
  }))
  if (verbose) cli_alert_info("Found {length(unique(tidygraphs$id))} disconnected graph(s) in data")
  if (verbose) cli_alert_info("Sorting disconnected regions by decreasing size")

  # return results as a character vector
  labeler <- setNames(tidygraphs$id, nm = tidygraphs$name)
  singleton_spots <- setdiff(spots, names(labeler))
  if (verbose)
    cli_alert_info("Found {length(spots) - length(labeler)} singletons in data")
  if (verbose & (length(singleton_spots) > 0))
    cli_alert("  These will be labeled as 'singletons'")
  labeler <- c(labeler, setNames(rep("singleton", length(singleton_spots)), nm = singleton_spots))
  labels <- labeler[spots]

  return(labels)

}


#' @param column_name A character specifying the name of a column in the meta data slot of a `Seurat` 
#' object that contains categorical data, e.g. clusters or manual selections
#' @param selected_groups A character vector to select specific groups in \code{column_name} with.
#' All groups are selected by default, but the common use case is to select a region of interest.
#' @param verbose Print messages
#'
#' @import dplyr
#' @importFrom rlang abort
#' @importFrom glue glue
#' @import cli
#' @importFrom tibble column_to_rownames
#'
#' @rdname disconnect-regions
#'
#' @author Ludvig Larsson
#'
#' @examples
#'
#' library(semla)
#' library(ggplot2)
#' library(patchwork)
#'
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
#'
#' # Plot selected variable
#' MapLabels(se_mcolon, column_name = "selection",
#'           pt_size = 3, override_plot_dims = TRUE)
#'
#' # Disconnect regions
#' se_mcolon <- DisconnectRegions(se_mcolon, column_name = "selection", selected_groups = "GALT")
#'
#' # Plot split regions
#' MapLabels(se_mcolon, column_name = "GALT_split",
#'           pt_size = 3, override_plot_dims = TRUE)
#'
#' # Note that if multiple sections are present, each section will be given
#' # it's own prefix in the disconnected groups.
#' se_merged <- MergeSTData(se_mcolon, se_mcolon)
#'
#' # Plot selected variable
#' MapLabels(se_merged, column_name = "selection",
#'           pt_size = 3, override_plot_dims = TRUE) +
#'   plot_layout(guides = "collect") &
#'   theme(legend.position = "top")
#'
#' # Disconnect regions
#' se_merged <- DisconnectRegions(se_merged, column_name = "selection", selected_groups = "GALT")
#'
#' # Plot split regions
#' MapLabels(se_merged, column_name = "GALT_split",
#'           pt_size = 3, override_plot_dims = TRUE) +
#'   plot_layout(guides = "collect") &
#'   theme(legend.position = "top")
#'
#' @export
#'
DisconnectRegions.Seurat <- function (
  object,
  column_name,
  selected_groups = NULL,
  verbose = TRUE,
  ...
) {

  # Set global variables to NULL
  barcode <- pxl_col_in_fullres <- pxl_row_in_fullres <- sampleID <- NULL

  # validate input
  .check_seurat_object(object)
  .validate_column_name(object, column_name)
  selected_groups <- .validate_selected_labels(object, selected_groups, column_name)

  # Select spots
  spots_list <- .get_spots_list(object, selected_groups, column_name, split_by_sample = FALSE)

  # Get coordinates
  coords <- GetStaffli(object)@meta_data |>
    select(barcode, pxl_col_in_fullres, pxl_row_in_fullres, sampleID) |>
    setNames(nm = c("barcode", "x", "y", "sampleID"))

  # Get disconnected groups
  disconnected_groups <- setNames(lapply(names(spots_list), function(lbl) {
    if (verbose) cli_alert_info("Extracting disconnected components for group '{lbl}'")
    DisconnectRegions(coords, spots = spots_list[[lbl]], verbose = verbose, ...)
  }), nm = names(spots_list))

  # Convert results
  column_labels <- paste0(names(disconnected_groups), "_split")
  disconnected_groups <- lapply(seq_along(disconnected_groups), function(i) {
    x <- disconnected_groups[[i]]
    tibble(barcode = names(x), label = x) |>
      setNames(nm = c("barcode", column_labels[i]))
  }) |>
    setNames(nm = names(disconnected_groups))

  # Add new labels to Seurat object
  tmp_mData <- object@meta.data |>
    select(-contains(column_labels), -contains("barcode")) |>
    bind_cols(barcode = colnames(object))

  # Add new labels to Seurat meta data
  for (lbl in names(disconnected_groups)) {
    tmp_mData <- tmp_mData |>
      left_join(y = disconnected_groups[[lbl]], by = "barcode")
  }
  tmp_mData <- tmp_mData |>
    column_to_rownames(var = "barcode")

  # Place new meta data in Seurat object
  object@meta.data <- tmp_mData

  return(object)

}


#' Check if spatial network is disconnected
#'
#' @param spatnet A list of tibbles with spot IDs for neighboring spots
#' generate with \code{\link{GetSpatialNetwork}}.
#'
#' @noRd
.is_disconnected <- function (
    spatnet
) {
  checks <- sapply(names(spatnet), function(nm) {
    gr <- tidygraph::as_tbl_graph(x = spatnet[[nm]])
    is_connected <- tidygraph::with_graph(graph = gr, expr = tidygraph::graph_is_connected())
    return(!is_connected)
  })
  return(any(checks))
}
