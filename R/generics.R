#' Load images
#' 
#' Load H&E images (or any other custom images) required for visualization methods in \code{semla}.
#' 
#' @param object An object
#' @param ... Arguments passed to other methods
#' 
#' @family pre-process
#' @rdname load-images
#' 
#' @return An object with images in \code{raster} format
#' 
#' @export
#' 
LoadImages <- function(object, ...) {
  UseMethod(generic = 'LoadImages', object = object)
}


#' Create Spatial Networks
#'
#' Create spatial networks from spatial coordinates. The spatial networks are provided in a long
#' format which holds information about spot neighbors, their center-to-center distances and positions.
#'
#' @details The default method expects an object of class \code{tbl} or \code{data.frame} with four columns
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
#' pixels. \code{nNeighbors} defines the number of nearest neighbors for "from" spots selected by \code{GetSpatialNetwork}.
#' "x_start", "y_start" are the spatial coordinates for "from" spots while "x_end", "y_end" are the
#' spatial coordinates for the neighboring "to" spots.
#'
#' @family network-methods
#' @rdname get-network
#'
#' @export
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


#' Map numeric features
#'
#' \code{MapFeatures} can be used to map numeric features to spots where the
#' values are encoded as colors. If multiple features and samples are provided, these
#' will be plotted individually and arranged into a grid of subplots.
#'
#' @details Note that you can only plot numeric features with \code{MapFeatures},
#' for example: gene expression, QC metrics or dimensionality reduction vectors.
#' If you want to plot categorical features, use \code{\link{MapLabels}} instead.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @family spatial-visualization-methods
#' @rdname visualize-features
#'
#' @export
#'
MapFeatures <- function(object, ...) {
  UseMethod(generic = 'MapFeatures', object = object)
}

#' Map numeric features or categorical labels
#'
#' \code{MapMultipleFeatures} can be used to map multiple numeric feature on the
#' same tissue section(s). \code{MapFeatures} provides an option to map 2 or 3
#' features using RGB color blending, whereas \code{MapMultipleFeatures} can handle
#' more than 3 features. See details below for more information.
#'
#' @details RGB color blending is only possible with 2 or 3 features. To visualize more
#' than 3 colors in the same plot, we can instead assign a specific color to each spot
#' depending on what feature has the highest value. This means that only the "dominant"
#' feature will be shown in each spot.
#'
#' Before using \code{MapMultipleFeatures}, you should be aware that this type of visualization
#' will bias what is being shown. Preferably, the selected features should be selecting in
#' a way that they are mutually exclusive, i.e. expressed in different spatial locations.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @family spatial-visualization-methods
#' @rdname visualize-multiple-features
#'
#' @export
#'
MapMultipleFeatures <- function(object, ...) {
  UseMethod(generic = 'MapMultipleFeatures', object = object)
}


#' Map categorical features
#'
#' \code{MapLabels} colors spots based on the values of a categorical feature.
#'
#' @details You can only plot categorical features with \code{MapLabels},
#' for example: spot annotations or clusters.
#' If you want to plot numerical features, use \code{\link{MapFeatures}} instead.
#' Only 1 categorical feature can be provided.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @family visualization-visualization-methods
#' @rdname visualize-labels
#'
#' @export
#'
MapLabels <- function(object, ...) {
  UseMethod(generic = 'MapLabels', object = object)
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
#' @returns An object with masked images
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
#' @returns An object with labeled region neighbors
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
#' @returns An object with radial distances
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
#' @returns An object with disconnected regions
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
#' @param ... Parameters passed to \code{\link{GetSpatialNework}}
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
#'    of spatial association. TEST, 27(3), 716–748 \doi{10.1007/s11749-018-0599-x}
#'    }
#' }
#' 
#' @returns An object with local G statistics
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
#' @author Ludvig Larsson
#' 
#' @returns An object with aligned images
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
#' from the \code{RcppML} R package developed by Zachary DeBruine.
#' 
#' @section Examples:
#' A tutorial can be found on our \href{https://ludvigla.github.io/semla/}{package website}.
#' Got to tutorials -> Cell type mapping
#'
#' @references
#' \doi{https://doi.org/10.1101/2021.09.01.458620}
#' [RcppML GitHub page](https://github.com/zdebruine/RcppML)
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#' 
#' @author Ludvig Larsson
#'
#' @rdname celltype-prediction
#' 
#' @return An object with cell type proportion estimates
#'
#' @export
#'
#' @md
RunNNLS <- function(object, ...) {
  UseMethod(generic = "RunNNLS", object = object)
}

#' @title Mixed cell type prediction with NNLS (experimental)
#' 
#' @description
#'
#' This function can be used to obtain proportion estimates for a mixture of cell types 
#' and unknown factors. The input spatial expression matrix is first deconvolved using
#' Non-negative Matrix Factorization (NMF) with the rank set to number of cell types + k 
#' additional factors. The estimated factor expression profiles are are then compared with 
#' the cell type expression profiles by computing pairwise correlation scores. 
#' Factor expression profiles with the lowest correlation scores are kept and combined
#' with the cell type expression profiles to predict proportions directly from the spatial
#' data using Non-negative Least Squares (NNLS).
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname celltype-prediction-mixed
#' 
#' @return An object with mixed cell type and factor proportion estimates
#'
#' @export
#'
#' @md
RunMixedNNLS <- function(object, ...) {
  UseMethod(generic = "RunMixedNNLS", object = object)
}
