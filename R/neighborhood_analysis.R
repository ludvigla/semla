#' @include generics.R
#'
NULL

#' Autodetect region neighbors
#'
#' This function allows you to automatically identify neighbors of a selected region.
#'
#' One way of using method this is to find spots surrounding a certain cluster. First, you need to make sure
#' the identity of the Seurat object is set to the meta.data column that you want to use, so for example
#' `se <- SetIdent(se, value = "seurat_clusters")` if you want to use the default seurat clusters.
#' Then you select the label that defined the region of interest using the `id` parameter, so for example
#' `ìd = "1"` will use cluster 1 as the region. If you set the `keep.idents` parameter to TRUE, the cluster ids
#' of the neighboring spots will be kept in the result, otherwise they will be returned as one single goup.
#' You can also activate the `keep.within.id` parameter to include all spots of the selected region in the output,
#' otherwise only the spots along the region border will be kept.
#'
#' @param spatnet A list of spatial networks generated with \code{\link{GetSpatialNetwork}}
#' @param spots A character vector with spot IDs present in `spatnet`
#' @param keep.within.id If set to TRUE, all id spots are kept, otherwise only the spots with outside neighbors are kept
#' @param verbose Print messages
#'
#' @rdname region-neighbors
#'
#' @importFrom rlang inform
#' @importFrom dplyr filter mutate
#'
#' @export
#'
#' @examples
#' \dontrun{
#' }
#'
RegionNeighbors.default <- function (
    spatnet,
    spots,
    keep.within.id = FALSE,
    verbose = FALSE
) {

  # Combine spatial networks into one tibble
  spatnet_combined <- do.call(rbind, lapply(seq_along(spatnet), function(i) {
    spnet <- spatnet[[i]]
    spnet$sampleID <- i
    return(spnet)
  }))

  # Check if spots are in spatnet_combined
  if (!is.character(spots)) abort(glue("Invalid class {class(spots)} for 'spots', expected 'character'"))
  if (length(spots) == 0) abort(glue("'spots' character vector is empty"))
  if (!all(spots %in% spatnet_combined$from)) {
    warn(glue("{sum(!spots %in% spatnet_combined$from)} spots had 0 neighbors."))
    spots <- intersect(spots, spatnet_combined$from)
  }

  # Obtain group labels
  spatnet_combined <- spatnet_combined |>
    filter(from %in% spots)
  if (nrow(spatnet_combined) == 0) abort("0 neighbors found")
  if (verbose) inform(c(">" = glue("Found {nrow(spatnet_combined)} neighbors for selected spots")))

  if (!keep.within.id) {
    if (verbose) inform(c(">" = "Excluding neighbors from the same group"))
    spatnet_combined <- spatnet_combined |>
      filter(!to %in% spots)
    if (nrow(spatnet_combined) == 0) abort("0 neighbors found after filtering")
    if (verbose) inform(c(">" = glue("{nrow(spatnet_combined)} neighbors left")))
  }

  if (verbose) inform(c(">" = "Returning neighbors"))
  return(spatnet_combined$to)
}

#' @param object an object
#' @param column_name string specifying a column name in your meta data
#' with labels, e.g. clusters or manual selections
#' @param column_labels character vector with labels to find nearest neighbors for.
#' These labels need to be present in the meta data columns specified by \code{column_name}
#' @param column_key prefix to columns returned in the Seurat object
#' @param keep.within.id logical specifying whether neighboring spots of the same label
#' should be kept
#' @param verbose print messages
#'
#' @rdname region-neighbors
#'
#' @importFrom glue glue
#' @importFrom rlang inform abort
#' @importFrom dplyr bind_cols select left_join filter bind_cols pull distinct sym all_of
#' @importFrom tibble tibble
#'
#' @examples
#'
#' library(STUtility2)
#'
#' se_mbrain <- readRDS(Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"), "/mousebrain/se_mbrain")))
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
#'                               column_name = "seurat_clusters",
#'                               column_labels = "10")
#'
#' # Plot cluster 10 and its neighbors
#' se_mbrain$selected_clusters <- se_mbrain[[]] |>
#'   mutate(across(where(is.factor), as.character)) |>
#'   mutate(cl = case_when(seurat_clusters %in% "10" ~ seurat_clusters,
#'                         TRUE ~ NA_character_)) |>
#'   pull(cl)
#'
#' MapLabels(se_mbrain, column_name = "selected_clusters") | MapLabels(se_mbrain, column_name = "nb_to_10")
#'
#' # Find neighbors to clusters 10 and 13
#' se_mbrain$selected_clusters <- se_mbrain[[]] |>
#'   mutate(across(where(is.factor), as.character)) |>
#'   mutate(cl = case_when(seurat_clusters %in% c("10", "13") ~ seurat_clusters,
#'                         TRUE ~ NA_character_)) |>
#'   pull(cl)
#' se_mbrain <- RegionNeighbors(se_mbrain,
#'                               column_name = "seurat_clusters",
#'                               column_labels = c("10", "13"))
#'
#' # Plot cluster 10, 13 and its neighbors
#' library(patchwork)
#' MapLabels(se_mbrain, column_name = "selected_clusters") +
#'   MapLabels(se_mbrain, column_name = "nb_to_10") +
#'   MapLabels(se_mbrain, column_name = "nb_to_13") +
#'   plot_layout(design = c(area(1, 1, 1, 1),
#'                          area(1, 2, 1, 2),
#'                          area(1, 3, 1, 3)))
#'
#' # it is also possible to pass additional parameters to GetSpatialNetwork
#' # to make it find more neighbors at a larger distances
#' se_mbrain <- RegionNeighbors(se_mbrain,
#'                               column_name = "seurat_clusters",
#'                               column_labels = "10",
#'                               nNeighbors = 40,
#'                               maxDist = Inf)
#' MapLabels(se_mbrain, column_name = "nb_to_10")
#'
#' @export
#'
RegionNeighbors.Seurat <- function (
    object,
    column_name,
    column_labels = NULL,
    column_key = "nb_to_",
    keep.within.id = FALSE,
    verbose = TRUE,
    ...
) {

  # validate Seurat object
  .check_seurat_object(object)

  # Check if column_name is available in meta data
  if (!is.character(column_name)) abort(glue("Invalid class '{class(column_name)}', expected a 'character'"))
  if (length(column_name) != 1) abort(glue("Invalid length {length(column_name)} for 'column_name', expected a character vector of length 1"))
  if (!column_name %in% colnames(object[[]])) abort("'column_name' is not present on the Seurat object meta data")

  # Check column label
  column_labels <- column_labels %||% {
    inform(glue("No 'column_labels' provided. Using all groups in '{column_name}' column"))
    column_labels <- unique(object[[]] |> pull(all_of(column_name)))
  }
  if (!is.character(column_labels)) abort(glue("Invalid class '{class(column_labels)}', expected a 'character'"))
  if (!all(column_labels %in% (object[[]] |> pull(all_of(column_name))))) abort(glue("Some 'column_labels' are not present in '{column_name}'"))

  # Select spots
  spots_list <- lapply(column_labels, function(lbl) {
    spots <- GetStaffli(object)@meta_data |>
      bind_cols(object[[]] |> select(all_of(column_name))) |>
      filter(!! sym(column_name) == lbl) |>
      pull(barcode)
  }) |>
    setNames(column_labels)

  # get spatial networks
  spatnet <- GetSpatialNetwork(object, ...)

  # Find neighbors
  nbs <- setNames(lapply(names(spots_list), function(nm) {
    if (verbose) inform(c("i" = glue("Finding neighboring spots for '{nm}'")))
    spots <- spots_list[[nm]]
    to_spots <- RegionNeighbors(spatnet, spots, keep.within.id, verbose)
    to_spots <- tibble(to_spots, nm) |>
      setNames(nm = c("barcode", column_name)) |>
      distinct()
    return(to_spots)
  }), nm = column_labels)

  # Add data to Seurat object
  nbs_rearranged <- do.call(bind_cols, lapply(column_labels, function(lbl) {
    left_join(x = GetStaffli(object)@meta_data |> select(barcode),
              y = nbs[[lbl]],
              by = "barcode") |>
      select(-barcode) |>
      setNames(nm = paste0(column_key, lbl))
  }))

  # Remove columns from meta data and add results
  object@meta.data <- object@meta.data |>
    select(-contains(colnames(nbs_rearranged))) |>
    bind_cols(nbs_rearranged)

  return(object)
}
