#' Autodetect region neighbors
#'
#' This function allows you to automatically identify neighbors of a selected region.
#'
#' One way of using method this is to find spots surrounding a certain cluster. First, you need to make sure
#' the identity of the Seurat object is set to the meta.data column that you want to use, so for example
#' `se <- SetIdent(se, value = "seurat_clusters")` if you want to use the default seurat clusters.
#' Then you select the label that defined the region of interest using the `id` parameter, so for example
#' `ìd = "1"` will use cluster 1 as the region. If you set the `keep.idents` parameter to TRUE, the cluster ids
#' of the neighbouring spots will be kept in the result, otherwise they will be returned as one single goup.
#' You can also activate the `keep.within.id` parameter to include all spots of the selected region in the output,
#' otherwise only the spots along the region border will be kept.
#'
#' @param spatnet A list of spatial networks generated with \code{\link{GetSpatialNetwork}}
#' @param spots A character vector with spot IDs present in `spatnet`
#' @param keep.within.id If set to TRUE, all id spots are kept, otherwise only the spots with outside neighbours are kept
#' @param verbose Print messages
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
RegionNeighbours.default <- function (
    spatnet,
    spots,
    keep.within.id = FALSE,
    verbose = FALSE
) {

  # Combine spatial networks into one tibble
  spatnet_combined <- do.call(rbind, lapply(seq_along(spatnet), function(i) {
    spnet <- spatnet[[i]]
    spnet$sample <- paste0(i)
    return(spnet)
  }))

  # Obtain group labels
  spatnet_combined <- spatnet_combined |>
    filter(from %in% spots)
  if (nrow(spatnet_combined) == 0) abort("0 neighbors found")
  if (verbose) inform(c(">" = glue("Found {nrow(spatnet_combined)} neighbours for selected spots")))

  if (!keep.within.id) {
    if (verbose) inform(c(">" = "Excluding neighbours from the same group"))
    spatnet_combined <- spatnet_combined |>
      filter(!to %in% spots)
    if (nrow(spatnet_combined) == 0) abort("0 neighbors found after filtering")
    if (verbose) inform(c(">" = glue("{nrow(spatnet_combined)} neighbours left")))
  }

  if (verbose) inform(c(">" = "Returning neighbors"))
  return(spatnet_combined$to)
}


# TESTS
# -----------------------

# samples <- Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"), "/*/filtered_feature_bc_matrix.h5"))
# imgs <- Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"), "/*/spatial/tissue_hires_image.png"))
# spotfiles <- Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"), "/*/spatial/tissue_positions_list.csv"))
# json <- Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"), "/*/spatial/scalefactors_json.json"))
#
# # Create a tibble/data.frame with file paths
# library(tibble)
# infoTable <- tibble(samples, imgs, spotfiles, json, sample_id = c("mousebrain", "mousecolon"))
#
# # Create Seurat object
# se <- ReadVisiumData(infoTable = infoTable)
#
# # Create a spatial network from a tibble with barcodes, (x, y) coordinates and sample IDs
# coordfiles <- c(system.file("extdata/mousebrain/spatial", "tissue_positions_list.csv", package = "STUtility2"),
#                 system.file("extdata/mousecolon/spatial", "tissue_positions_list.csv", package = "STUtility2"))
#
# # Load coordinate data into a tibble
# xys <- do.call(rbind, lapply(seq_along(coordfiles), function(i) {
#   coords <- setNames(read.csv(coordfiles[i], header = FALSE), nm = c("barcode", "selection", "grid_y", "grid_x", "y", "x"))
#   coords$sample <- paste0(i)
#   coords <- coords |>
#     dplyr::filter(selection == 1) |>
#     dplyr::select(barcode, x, y, sample) |>
#     tibble::as_tibble()
#   return(coords)
# }))
#
# # Create spatial networks
# spatnet <- GetSpatialNetwork(xys)
