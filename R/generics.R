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
#' This function can be used to find genes with high spatial autocorrelation in SRT data.
#' A more detailed description of the algorithm is outlined in the Details section below.
#'
#' @details The default method expects a matrix-like object with features in columns and spots in rows
#' and a list of spatial networks generated with \code{\link{GetSpatialNetwork}}.
#'
#' If \code{across_all} is set to \code{TRUE}, the spatial autocorrelation scores will be computed
#' across all samples. Otherwise, the scores will be calculated for each sample separately, and returns
#' a list.
#'
#' @section Spatial autocorrelation:
#' Spatial autocorrelation is the term used to describe the presence of systematic spatial
#' variation. Positive spatial autocorrelation of a feature is the tendency for regions that
#' are close together in space to have similar values for that feature.
#'
#' A simple example is when you have an anatomical structure or a tissue type that spans
#' across multiple neighboring spots in an SRT experiment, for example a gland, an immune
#' infiltrate or a region of the brain. Inside such structures, you might find that the
#' expression levels of certain genes (or other features) are highly similar and hence
#' these genes have a positive spatial autocorrelation.
#'
#' The method provided in `STUtility2` is relatively simple and fast. For each feature and spot,
#' the expression is averaged across all neighboring spots (typically the 6 closest neighbors)
#' to produce a lag expression vector. Since this vector represents the average of the surrounding
#' spots, we can use it to test if the expression in those spots is similar to the center spot.
#' One simple strategy is to calculate the pearson correlation between a genes' lag vector and
#' the original expression vector which typically captures the spatial autocorrelation well.
#'
#' @section Method steps:
#' \itemize{
#'    \item{Load a matrix with features in rows and spots in columns: \eqn{X_{expr}}}
#'    \item{Convert the corresponding spatial network to wide format and construct a nearest
#'    neighbor matrix \eqn{N_{neighbors}} in which neighboring spots have a value of 1
#'    and the remaining spots have a value of 0
#'    }
#'    \item{\eqn{N_{neighbors}} is then multiplied with the \eqn{X_{expr}} to
#'    calculate a lag vector for each feature: \cr \cr
#'    \eqn{X_{lagexpr} = (N_{neighbors}*X_{expr})/n_{neighbors}} \cr \cr
#'    where \eqn{n_{neighbors}} is the number of neighbors for each spot.
#'    }
#'    \item{The spatial autocorrelation score for a genes is the 'pearson' correlation of the
#'    lag vector and the initial expression vector: \cr \cr
#'    \eqn{spatcor_{feature} = cor(X_{lagexpr}[feature, ], X_{expr}[feature, ])}
#'    }
#' }
#'
#' @param object An object (see details)
#' @param ... Arguments passed to other methods
#'
#' @family network-methods
#' @rdname cor-features
#'
#' @return Either a list of tibbles or a tibble with feature names and correlation scores for
#' each feature in the input feature matrix
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
#' @section Seurat:
#' If a \code{Seurat} object is provided, the identified neighbors will be
#' stored as new columns in the meta data column with the names prefixed by
#' \code{column_key}
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
