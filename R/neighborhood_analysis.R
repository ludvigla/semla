#' @include generics.R
#'
NULL


#' Autodetect region neighbors
#'
#' This function allows you to automatically identify neighbors of a selected region.
#'
#' @section Seurat:
#' If a Seurat object is provided, the \code{RegionNeighbors} takes a meta
#' data column (chosen with \code{column_name}) with categorical labels,
#' finds the nearest neighbors of spots for a selected group in this columns
#' (chosen with \code{column_labels}) and returns new meta data column with
#' labels for the nearest neighbors of the selected group. If no \code{column_labels}
#' are specified, the method will return a column for each separate category in
#' the \code{column_name} vector.
#'
#' Note that the prefix to the returned column names will be selected based on the \code{mode}.
#' You can overwrite this behavior by manually setting \code{column_key}.
#'
#' Below is some additional information about the behavior of different \code{mode}s:
#'
#' \itemize{
#'    \item{return outer border (default): \code{mode="outer"}}
#'    \item{return inner border: \code{mode="inner"}}
#'    \item{return inner and outer borders: \code{mode="inner_outer"}}
#'    \item{return all selected spots and outer border: \code{mode="all_inner_outer"}}
#' }
#'
#' @section default method:
#' The default method takes a list of spatial networks generated with
#' \code{\link{GetSpatialNetwork}} together with a vector of spot IDs
#' and returns the spot IDs for border spots. The behavior for border
#' spot selection is determined by the \code{mode}.
#'
#' @param spots A character vector with spot IDs present in `spatnet`
#' @param mode Select mode (see details)
#' @param verbose Print messages
#'
#' @rdname region-neighbors
#' @family spatial-methods
#'
#' @author Ludvig Larsson
#'
#' @importFrom rlang warn
#' @import cli
#' @importFrom dplyr filter mutate
#'
#' @export
#'
RegionNeighbors.default <- function (
    object,
    spots,
    mode = c("outer", "inner", "inner_outer", "all_inner_outer"),
    verbose = FALSE,
    ...
) {

  # Set global variables to NULL
  from <- to <- NULL

  # Get mode
  mode <- match.arg(mode, choices = c("outer", "inner", "inner_outer", "all_inner_outer"))

  # Check object
  stopifnot(
    inherits(object, what = "list"),
    length(object) > 0,
    inherits(spots, what = "character"),
    length(spots) > 0,
    inherits(mode, what = "character"),
    length(mode) == 1
  )

  # Combine spatial networks into one tibble
  spatnet_combined <- do.call(rbind, lapply(names(object), function(nm) {
    spnet <- object[[nm]]
    spnet$sampleID <- nm
    return(spnet)
  }))

  # Check if spots are in spatnet_combined
  if (!is.character(spots)) abort(glue("Invalid class {class(spots)} for 'spots', expected 'character'"))
  if (length(spots) == 0) abort(glue("'spots' character vector is empty"))
  if (!all(spots %in% spatnet_combined$from)) {
    warn(glue("{sum(!spots %in% spatnet_combined$from)} spots had 0 neighbors."))
    spots <- intersect(spots, spatnet_combined$from)
  }

  # Filter spatial network to contain from spots in selected group
  spatnet_combined <- spatnet_combined |>
    filter(from %in% spots)

  # If outer_border = TRUE, return neighboring spots
  # otherwise, return inner border
  if (mode == "outer") {
    if (verbose) cli_alert("  Excluding neighbors from the same group")
    spatnet_combined <- spatnet_combined |>
      filter(!to %in% spots)
    if (nrow(spatnet_combined) == 0) abort("0 neighbors found after filtering")
    if (verbose) cli_alert("  {nrow(spatnet_combined)} neighbors left")
    spots_keep <- unique(spatnet_combined$to)
  }
  if (mode == "inner") {
    spatnet_combined <- spatnet_combined |>
      filter(!((from %in% spots) & (to %in% spots)))
    spots_keep <- unique(spatnet_combined$from)
  }
  if (mode == "inner_outer") {
    spatnet_combined <- spatnet_combined |>
      filter(!((from %in% spots) & (to %in% spots)))
    spots_keep <- unique(c(spatnet_combined$from, spatnet_combined$to))
  }
  if (mode == "all_inner_outer") {
    spots_keep <- unique(c(spatnet_combined$from, spatnet_combined$to))
  }

  if (verbose) cli_alert("  Returning neighbors")

  return(spots_keep)
}

#' @param column_name string specifying a column name in your meta data
#' with labels, e.g. clusters or manual selections
#' @param column_labels character vector with labels to find nearest neighbors for.
#' These labels need to be present in the meta data columns specified by \code{column_name}
#' @param column_key prefix to columns returned in the Seurat object
#'
#' @rdname region-neighbors
#'
#' @import cli
#' @importFrom glue glue
#' @importFrom rlang abort %||%
#' @import dplyr
#' @importFrom tibble tibble
#'
#' @examples
#'
#' library(semla)
#' library(dplyr)
#'
#' se_mbrain <-
#'   readRDS(system.file("extdata",
#'   "/mousebrain/se_mbrain",
#'   package = "semla"))
#'
#' # Create Seurat object
#' se_mbrain <- se_mbrain |>
#'   ScaleData(verbose = FALSE) |>
#'   RunPCA(verbose = FALSE) |>
#'   FindNeighbors(verbose = FALSE) |>
#'   FindClusters(verbose = FALSE)
#'
#' # Find neighbors to cluster 10
#' se_mbrain <- RegionNeighbors(se_mbrain,
#'                              column_name = "seurat_clusters",
#'                              column_labels = "10")
#'
#' # Plot cluster 10 and its neighbors
#' se_mbrain$selected_clusters <- se_mbrain[[]] |>
#'   mutate(across(where(is.factor), as.character)) |>
#'   mutate(cl = case_when(seurat_clusters %in% "10" ~ seurat_clusters,
#'                         TRUE ~ NA_character_)) |>
#'   pull(cl)
#'
#' MapLabels(se_mbrain, column_name = "selected_clusters") |
#'   MapLabels(se_mbrain, column_name = "nb_to_10")
#'
#' \donttest{
#' # Find neighbors to clusters 8 and 10
#' se_mbrain$selected_clusters <- se_mbrain[[]] |>
#'   mutate(across(where(is.factor), as.character)) |>
#'   mutate(cl = case_when(seurat_clusters %in% c("8", "10") ~ seurat_clusters,
#'                         TRUE ~ NA_character_)) |>
#'   pull(cl)
#' se_mbrain <- RegionNeighbors(se_mbrain,
#'                              column_name = "seurat_clusters",
#'                              column_labels = c("8", "10"))
#'
#' # Plot cluster 8, 10 and its neighbors
#' library(patchwork)
#' MapLabels(se_mbrain, column_name = "selected_clusters") +
#'   MapLabels(se_mbrain, column_name = "nb_to_8") +
#'   MapLabels(se_mbrain, column_name = "nb_to_10") +
#'   plot_layout(design = c(area(1, 1, 1, 1),
#'                          area(1, 2, 1, 2),
#'                          area(1, 3, 1, 3)))
#'
#' # it is also possible to pass additional parameters to GetSpatialNetwork
#' # to make it find more neighbors at a larger distances
#' se_mbrain <- RegionNeighbors(se_mbrain,
#'                              column_name = "seurat_clusters",
#'                              column_labels = "10",
#'                              nNeighbors = 40,
#'                              maxDist = Inf)
#' MapLabels(se_mbrain, column_name = "nb_to_10")
#' }
#' 
#' @export
#'
RegionNeighbors.Seurat <- function (
    object,
    column_name,
    column_labels = NULL,
    mode = c("outer", "inner", "inner_outer", "all_inner_outer"),
    column_key = NULL,
    verbose = TRUE,
    ...
) {

  # Set global variables to NULL
  barcode <- var2 <- NULL

  # Check mode
  mode <- match.arg(mode, choices = c("outer", "inner", "inner_outer", "all_inner_outer"))

  # Define column key
  if (is.null(column_key)) {
    if (mode %in% c("outer", "inner_outer", "all_inner_outer")) {
      column_key <- "nb_to_"
    } else {
      column_key <- "inner_border_"
    }
  } else {
    if ((!inherits(column_key, what = "character")) | (length(column_key) != 1))
      abort(glue("Invalid column_key, expected a 'character' of length 1"))
  }

  # validate Seurat object
  .check_seurat_object(object)

  # Check if column_name is available in meta data
  if (!inherits(column_name, what = "character"))
    abort(glue("Invalid class '{class(column_name)}' for 'column_name', expected a 'character' of length 1"))
  if (length(column_name) != 1)
    abort(glue("Invalid length {length(column_name)} for 'column_name', expected a vector of length 1"))
  if (!column_name %in% colnames(object[[]]))
    abort("'column_name' is not present on the Seurat object meta data")

  # Check column label
  column_labels <- column_labels %||% {
    cli_alert_info("No 'column_labels' provided. Using all groups in '{column_name}' column")
    column_labels <- unique(object[[]] |> pull(all_of(column_name))) |> as.character()
  }
  if (!is.character(column_labels))
    abort(glue("Invalid class '{class(column_labels)}', expected a 'character'"))
  if (!all(column_labels %in% (object[[]] |> pull(all_of(column_name)))))
    abort(glue("Some 'column_labels' are not present in '{column_name}'"))
  if (any(grepl(" ", column_labels)))
    abort(glue("Some 'column_labels' in '{column_name}' contain blank spaces which is not permitted. Remove/replace blankspaces before proceeding."))
  
  # Select spots
  spots_list <- lapply(column_labels, function(lbl) {
    spots <- GetStaffli(object)@meta_data |>
      bind_cols(object[[]] |> select(all_of(column_name))) |>
      mutate_if(is.factor, as.character) |>
      filter(!! sym(column_name) == lbl) |>
      pull(barcode)
  }) |>
    setNames(column_labels)

  # get spatial networks
  spatnet <- GetSpatialNetwork(object, ...)

  # Find neighbors
  nbs <- setNames(lapply(names(spots_list), function(nm) {
    if (verbose) cli_alert_info("Finding neighboring spots for '{nm}'")
    spots <- spots_list[[nm]]
    to_spots <- RegionNeighbors(spatnet, spots, mode, verbose)
    to_spots <- tibble(to_spots, nm) |>
      setNames(nm = c("barcode", column_name)) |>
      distinct()
    return(to_spots)
  }), nm = column_labels)

  # Add data to Seurat object
  nbs_rearranged <- do.call(bind_cols, lapply(column_labels, function(lbl) {
    left_join(x = GetStaffli(object)@meta_data |>
                select(barcode) |>
                bind_cols(object[[]] |> select(all_of(column_name))),
              y = nbs[[lbl]],
              by = "barcode") |>
      select(-barcode) |>
      setNames(nm = c("var1", "var2")) |>
      mutate_if(is.factor, as.character) |>
      mutate(var2 = case_when((!is.na(var2)) & (var1 != var2) ~ paste0(column_key, lbl),
                              TRUE ~ var2)) |>
      select(var2) |>
      setNames(nm = paste0(column_key, lbl))
  }))

  # Remove columns from meta data and add results
  object@meta.data <- object@meta.data |>
    select(-contains(colnames(nbs_rearranged))) |>
    bind_cols(nbs_rearranged)

  return(object)
}
