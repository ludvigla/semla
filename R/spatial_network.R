#' Create Spatial Networks
#'
#' Create spatial networks from spatial coordinates. The function expects a `data.frame` or `tibble`
#' containing barcodes, coordinates and sample ids. One spatial network is created per sample and
#' the results are returned as a list.
#'
#' @param xys An object of class `tbl` or `data.frame` with four columns "barcode", "x", "y"
#' and "sample" holding the coordinates for a set of spots. The "barcode" column is a character
#' vector with spatial barcodes, "x", "y" hold numeric values representing the spot coordinates and
#' "sample" is a character vector with unique sample IDs.
#' @param nNeighbours Number of nearest neighbors to calculate for each spot. The default
#' number of neighbors is 6 for the 'Visium' platform and 4 for the '1k' and '2k' platforms.
#' @param maxDist Distance cut-off for nearest neighbors to consider. If set to NULL (default),
#' `maxDist` is estimated from the data by taking the minimum neighbor distance multiplied by
#' a factor of `1.2`.
#' @param minK Minimum nearest neighbors [default: 0]. Spots with fewer neighbors will be discarded.
#' Useful if you want to remove spots with few or no neighbors.
#'
#' @family network-methods
#'
#' @return A list of tibbles, each containing information about the nearest neighbors of each spot.
#' For one spot in the column "from", its nearest neighboring spots are provided in the "to" column.
#' Distances correspond to distances between "to" and "from", and usually correspond to H&E image
#' pixels. "numK" defines the number of nearest neighbors for "from" spots selected by `GetSpatialNetwork`.
#' "x_start", "y_start" are the spatial coordinates for "from" spots while "x_end", "y_end" are the
#' spatial coordinates for the neighboring "to" spots.
#'
#' @importFrom dplyr group_by mutate ungroup filter left_join
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # Create a spatial network from a data.frame (xys)
#' spatnet <- GetSpatialNetwork(xys)
#'
#' # Plot network
#' ggplot(spatnet[["1"]], aes(x = x_start, xend = x_end, y = y_start, yend = y_end)) +
#'    geom_segment()
#'
#' }
#'
#' @export
GetSpatialNetwork <- function (
    xys,
    nNeighbours = 6,
    maxDist = NULL,
    minK = 0
) {

  # Check xys object class
  if (!class(xys) %in% c("data.frame", "matrix", "tbl"))
    abort(glue("Invalid class '{class(xys)}'."))
  if (ncol(xys) != 4)
    abort(glue("Invalid number of columns '{ncol(xys)}'. Expected 4."))
  if (!all(
    xys |> summarize(
      check_barcode = is.character(barcode),
      check_x = is.numeric(x),
      check_y = is.numeric(y),
      check_sample = is.character(sample)
    ) |>
    unlist()
  )) {
    abort(glue("Invalid column class(es)."))
  }

  # install dbscan if not already installed
  if (!requireNamespace("dbscan")) install.packages("dbscan")

  # Set number of nearest neighbors if NULL
  nNeighbours <- nNeighbours %||% 6
  # Split coordinates by sample
  xys.list <- split(xys, xys$sample)[unique(xys$sample)]

  # Compute network
  knn_long.list <- setNames(lapply(names(xys.list), function(sampleID) {

    xys_subset <- xys.list[[sampleID]]
    spotnames <- setNames(xys_subset$barcode, nm = c(1:nrow(xys_subset)) |> paste0())
    knn_spatial <- dbscan::kNN(x = xys_subset[, c("x", "y")] |> as.matrix(), k = nNeighbours)
    maxDist <- maxDist %||% (apply(knn_spatial$dist, 1, min) |> min())*1.2
    knn_long <- tibble::tibble(from = rep(1:nrow(knn_spatial$id), nNeighbours),
                                       to = as.vector(knn_spatial$id),
                                       distance = as.vector(knn_spatial$dist))

    knn_long$from <- spotnames[knn_long$from]
    knn_long$to <- spotnames[knn_long$to]
    knn_long <- knn_long |>
      group_by(from) |>
      mutate(numK = n()) |>
      ungroup() |>
      filter(distance <= maxDist, numK > minK)

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
