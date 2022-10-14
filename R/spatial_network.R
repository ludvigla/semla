#' @include generics.R
#'
NULL

#' @importFrom dplyr group_by mutate ungroup filter left_join summarize n add_count
#' @importFrom tibble tibble
#'
#' @param object An object
#' @param nNeighbors Number of nearest neighbors to calculate for each spot. The default
#' number of neighbors is 6 given the hexagonal pattern of 10x Visium arrays.
#' @param maxDist Distance cut-off for nearest neighbors to consider. If set to NULL (default),
#' `maxDist` is estimated from the data by taking the minimum neighbor distance multiplied by
#' a factor of `1.2`.
#' @param minK Minimum nearest neighbors [default: 0]. Spots with fewer neighbors will be discarded.
#' Useful if you want to remove spots with few or no neighbors.
#'
#' @rdname get-network
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # Create a spatial network from a tibble with barcodes, (x, y) coordinates and sample IDs
#' coordfiles <- c(system.file("extdata/mousebrain/spatial", "tissue_positions_list.csv", package = "STUtility2"),
#'                 system.file("extdata/mousecolon/spatial", "tissue_positions_list.csv", package = "STUtility2"))
#'
#' # Load coordinate data into a tibble
#' xys <- do.call(rbind, lapply(seq_along(coordfiles), function(i) {
#'   coords <- setNames(read.csv(coordfiles[i], header = FALSE), nm = c("barcode", "selection", "grid_y", "grid_x", "y", "x"))
#'   coords$sampleID <- i
#'   coords <- coords |>
#'     dplyr::filter(selection == 1) |>
#'     dplyr::select(barcode, x, y, sampleID) |>
#'     tibble::as_tibble()
#'   return(coords)
#' }))
#'
#' # Create spatial networks from xys coordinates
#' spatnet <- GetSpatialNetwork(xys)
#'
#' # Plot network
#' p1 <- ggplot(spatnet[["1"]], aes(x = x_start, xend = x_end, y = y_start, yend = y_end)) +
#'   geom_segment() +
#'   scale_y_reverse()
#'
#' p2 <- ggplot(spatnet[["2"]], aes(x = x_start, xend = x_end, y = y_start, yend = y_end)) +
#'   geom_segment() +
#'   scale_y_reverse()
#'
#' p1 + p2
#' }
#'
#' @export
#'
GetSpatialNetwork.default <- function (
    object,
    nNeighbors = 6,
    maxDist = NULL,
    minK = 0
) {

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

  # install dbscan if not already installed
  if (!requireNamespace("dbscan")) install.packages("dbscan")

  # Set number of nearest neighbors if NULL
  nNeighbors <- nNeighbors %||% 6
  # Split coordinates by sample
  xys.list <- split(object, object$sampleID)[unique(object$sampleID)]

  # Compute network
  knn_long.list <- setNames(lapply(names(xys.list), function(sampleID) {

    xys_subset <- xys.list[[sampleID]]
    spotnames <- setNames(xys_subset$barcode, nm = c(1:nrow(xys_subset)) |> paste0())
    knn_spatial <- dbscan::kNN(x = xys_subset[, c("x", "y")] |> as.matrix(), k = nNeighbors)
    maxDist <- maxDist %||% (apply(knn_spatial$dist, 1, min) |> min())*1.2
    knn_long <- tibble::tibble(from = rep(1:nrow(knn_spatial$id), nNeighbors),
                                       to = as.vector(knn_spatial$id),
                                       distance = as.vector(knn_spatial$dist))

    knn_long$from <- spotnames[knn_long$from]
    knn_long$to <- spotnames[knn_long$to]
    knn_long <- knn_long |>
      add_count(from) |>
      filter(distance <= maxDist, n > minK) |>
      add_count(from, name = "nn") |>
      select(-n)

    # Merge with coordinates
    knn_long <-
      left_join(
        x = knn_long,
        y = setNames(xys_subset[, 1:3], nm = c("barcode", "x_start", "y_start")),
        by = c("from" = "barcode")
      )
    knn_long <-
      left_join(
        x = knn_long,
        y = setNames(xys_subset[, 1:3], nm = c("barcode", "x_end", "y_end")),
        by = c("to" = "barcode")
      )

    return(knn_long)
  }), nm = names(xys.list))

  return(knn_long.list)
}

#' @rdname get-network
#'
#' @importFrom dplyr select rename
#'
#' @export
#' @method GetSpatialNetwork Seurat
#'
GetSpatialNetwork.Seurat <- function (
    object,
    nNeighbors = 6,
    maxDist = NULL,
    minK = 0
) {

  # Check Seurat object
  .check_seurat_object(object)

  # Get coordinates
  xys <- GetStaffli(object)@meta_data |>
    select(barcode, pxl_col_in_fullres, pxl_row_in_fullres, sampleID) |>
    rename(x = pxl_col_in_fullres, y = pxl_row_in_fullres)

  # get spatial networks
  spatnet <- GetSpatialNetwork(xys, nNeighbors = nNeighbors, maxDist = maxDist, minK = minK)

  # Return spatial networks
  return(spatnet)
}
