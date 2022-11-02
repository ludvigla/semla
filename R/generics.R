#' Create Spatial Networks
#'
#' Create spatial networks from spatial coordinates. The spatial networks are provided in a long
#' format which holds information about spot neighbors, their center-to-center distances and positions.
#'
#' @details The default method expects an object of class `tbl` or `data.frame` with four columns
#' "barcode", "x", "y" and "sample" holding the coordinates for a set of spots. The "barcode" column
#' is a character vector with spatial barcodes, "x", "y" hold numeric values representing the spot
#' coordinates and "sample" is a character vector with unique sample IDs.
#'
#' @param object An object (see details)
#' @param ... Arguments passed to other methods
#'
#' @return A list of tibbles, each containing information about the nearest neighbors of each spot.
#' For one spot in the column "from", its nearest neighboring spots are provided in the "to" column.
#' Distances correspond to distances between "to" and "from", and usually correspond to H&E image
#' pixels. "numK" defines the number of nearest neighbors for "from" spots selected by `GetSpatialNetwork`.
#' "x_start", "y_start" are the spatial coordinates for "from" spots while "x_end", "y_end" are the
#' spatial coordinates for the neighboring "to" spots.
#'
#' @family network-methods
#' @rdname get-network
#'
#' @export GetSpatialNetwork
#'
GetSpatialNetwork <- function(object, ...) {
  UseMethod(generic = 'GetSpatialNetwork', object = object)
}


#' Find features with high spatial autocorrelation
#'
#' @param object An object (see details)
#' @param ... Arguments passed to other methods
#'
#' @family network-methods
#' @rdname cor-features
#'
#' @export
#'
#' @md
CorSpatialFeatures <- function(object, ...) {
  UseMethod(generic = 'CorSpatialFeatures', object = object)
}


#' Map numeric features or categorical labels in 2D space
#'
#' \code{MapFeatures} can be used to map numeric features as colors on spot
#' coordinates in 2D space. If multiple features and samples are provided, these
#' will be plotted individually and arranged into a grid of subplots.
#'
#' @details Note that you can only plot numeric features with \code{MapFeatures},
#' for example: gene expression, QC metrics or dimensionality reduction vectors.
#' If you want to plot categorical features, use \code{\link{MapLabels}} instead.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @family visualization methods
#' @rdname visualize-features
#'
#' @export
#'
MapFeatures <- function(object, ...) {
  UseMethod(generic = 'MapFeatures', object = object)
}


#' Map labels in 2D
#'
#' \code{MapLabels} colors spot coordinates in 2D space based on the values of
#' a categorical feature. Only 1 feature can be provided.
#'
#' @details Note that you can only plot categorical features with \code{MapLabels},
#' for example: spot annotations or clusters.
#' If you want to plot numerical features, use \code{\link{MapFeatures}} instead.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @family visualization methods
#' @rdname visualize-labels
#'
#' @export
#'
MapLabels <- function(object, ...) {
  UseMethod(generic = 'MapLabels', object = object)
}


#' Load images
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname load-images
#'
#' @export
#'
LoadImages <- function(object, ...) {
  UseMethod(generic = 'LoadImages', object = object)
}


#' Apply rigid transformations to images
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @family image transformations
#' @rdname transform-images
#'
#' @export
#'
RigidTransformImages <- function(object, ...) {
  UseMethod(generic = 'RigidTransformImages', object = object)
}

#' Mask images
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @family image transformations
#' @rdname mask-images
#'
#' @export
#'
MaskImages <- function(object, ...) {
  UseMethod(generic = 'MaskImages', object = object)
}


#' Find region neighbors
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname region-neighbors
#'
#' @export
#'
RegionNeighbors <- function(object, ...) {
  UseMethod(generic = 'RegionNeighbors', object = object)
}


#' @title Calculate radial distances from a region border
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname radial-distance
#'
#' @export
#'
RadialDistance <- function(object, ...) {
  UseMethod(generic = "RadialDistance", object = object)
}



#' @title Disconnect regions
#'
#' @description This function allows you split spatially disconnected regions
#' belonging to the same category. In cases where a certain tissue type create
#' isolated "islands" in the tissue section, these islands can be separated. A
#' common example is tertiary lymphoid structures (TLS) which are typically
#' dispersed across a tissue section.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname disconnect-regions
#'
#' @export
#'
DisconnectRegions <- function(object, ...) {
  UseMethod(generic = "DisconnectRegions", object = object)
}
