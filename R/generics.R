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


#' @title Calculate local G statistic
#'
#' @description This local spatial statistic measures the concentration of
#' high or low values for a given region. This can for example be used to
#' define spatial structures in a tissue section based on the values of
#' selected features. Furthermore, the local G statistic can be used for
#' high/low clustering of data.
#'
#' NB: This function only calculates G, not G star.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @family spatial methods
#' @rdname local-G
#'
#' @references
#' \itemize{
#'    \item{
#'       Ord, J. K. and Getis, A. 1995 Local spatial autocorrelation statistics:
#'       distributional issues and an application. Geographical Analysis, 27, 286–306
#'    }
#'    \item{
#'    Bivand RS, Wong DWS 2018 Comparing implementations of global and local indicators
#'    of spatial association. TEST, 27(3), 716–748 [DOI link](10.1007/s11749-018-0599-x)
#'    }
#' }
#'
#'
#' @seealso
#' [Emerging Hot Spot Analysis](https://sfdep.josiahparry.com/articles/understanding-emerging-hotspots.html)
#' [G and G star local spatial statistics](https://r-spatial.github.io/spdep/reference/localG.html)
#' [ESRI Getis-ord General G](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/h-how-high-low-clustering-getis-ord-general-g-spat.htm)
#'
#' @export
#'
#' @md
RunLocalG <- function(object, ...) {
  UseMethod(generic = "RunLocalG", object = object)
}


#' Manually apply rigid transformations to images
#'
#' Opens an interactive viewer where images are interactive objects
#' that can be rotated, moved, scaled and mirrored. It is often difficult
#' to automate image registration for H&E images, in particular if tissue
#' sections are damaged, distorted or folded. In such situations, manual
#' registration can be very useful.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @family image transformations
#' @rdname manual-transform-images
#'
#' @export
#'
RunAlignment <- function(object, ...) {
  UseMethod(generic = 'RunAlignment', object = object)
}


#' @title Cell type prediction with NNLS
#'
#' @description
#'
#' This function can be used to project cell type expression
#' profiles onto a 10x Visium gene expression matrix to obtain
#' cell type proportion estimates.
#'
#' @details This method is suitable for
#' paired 10x Visium and scRNA-seq data. In other words, the
#' scRNA-seq data should represent the cell types present in
#' the 10x Visium spatial transcriptomics data.
#' Redundant cell types (i.e. cells that are present in the
#' scRNA-seq data but not in the data 10x Visium) or missing
#' cell types (i.e. cells that are present in the 10x Visium
#' data but not in the data scRNA-seq data) might affect the
#' results.
#'
#' Preferably, the scRNA-seq data should be composed of the
#' same cell types as the tissue section(s) processed by 10x
#' Visium and be collected from the same tissue specimen.
#' This method is intended to be used as a fast cell type
#' mapping which is useful for data driven exploration.
#'
#' We encourage users to explore alternative methods such as
#' stereoscope, cell2location or RCTD.
#'
#' The method uses the NNLS implementation
#' from the `RcppML` R package developed by Zachary DeBruine.
#'
#' @references
#' [RcppML](https://doi.org/10.1101/2021.09.01.458620)
#' [RcppML GitHub page](https://github.com/zdebruine/RcppML)
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname celltype-prediction
#'
#' @export
#'
#' @md
RunNNLS <- function(object, ...) {
  UseMethod(generic = "RunNNLS", object = object)
}
