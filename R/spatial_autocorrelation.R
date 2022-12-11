#' @include generics.R
#'
NULL

# TODO: dealing with NA values, test running with 0 cores
#' @description
#' This function can be used to find genes with high spatial autocorrelation in SRT data.
#' A more detailed description of the algorithm is outlined in the Details section below.
#'
#' @details The default method expects a matrix-like object with features in columns and spots in rows
#' and a list of spatial networks generated with \code{\link{GetSpatialNetwork}}.
#'
#' If \code{across_all} is set to \code{TRUE}, the spatial autocorrelation scores will be computed
#' across all samples. Otherwise, the scores will be calculated for each sample separately, and returns
#' a list with one `tibble` per sample.
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
#' The method provided in `STUtility2` works as follows. For each feature and spot,
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
#'    \item{The spatial autocorrelation score for a feature is the 'pearson' correlation of the
#'    lag vector and the initial expression vector: \cr \cr
#'    \eqn{spatcor_{feature} = cor(X_{lagexpr}[feature, ], X_{expr}[feature, ])}
#'    }
#' }
#' @param spatnet A list of spatial networks created with \code{\link{GetSpatialNetwork}}. The spots in these
#' networks should match the spots in the feature matrix.
#' @param across_all Should the autocorrelation scores be calculated across all samples?
#' @param nCores Number of cores to use for the spatial autocorrelation calculation
#' @param verbose Print messages
#'
#' @rdname cor-features
#'
#' @importFrom stats cor
#' @importFrom parallel mclapply detectCores
#' @importFrom tibble tibble
#' @importFrom Matrix rowSums
#' @importFrom tidyr pivot_wider
#' @import dplyr
#' @importFrom glue glue
#' @importFrom rlang abort
#' @importFrom methods as
#' @import cli
#'
#' @return Either a list of tibbles or a tibble with feature names and correlation scores
#'
#' @examples
#' \dontrun{
#' se_mbrain <- readRDS(system.file("extdata/mousebrain",
#'                                  "se_mbrain",
#'                                  package = "STUtility2"))
#' featureMat <- FetchData(se_mbrain, vars = VariableFeatures(se_mbrain)[1:100])
#'
#' coordfile <-
#'   system.file("extdata/mousebrain/spatial",
#'               "tissue_positions_list.csv",
#'               package = "STUtility2")
#'
#' # Load coordinate data into a tibble
#' xys <- setNames(read.csv(coordfile, header = FALSE),
#'                 nm = c("barcode", "selection", "grid_y", "grid_x", "y", "x"))
#' xys$sample <- paste0(1)
#' xys <- xys |>
#'   dplyr::mutate(barcode = paste0(barcode, "_", 1)) |>
#'   dplyr::filter(selection == 1) |>
#'   dplyr::select(barcode, x, y, sample) |>
#'   tibble::as_tibble()
#'
#' # Create spatial networks
#' spatnet <- GetSpatialNetwork(xys)
#' spatgenes <- CorSpatialFeatures(featureMat, spatnet, nCores = 1)
#'
#' # Check genes with highest spatial autocorrelation
#' head(spatgenes[[1]])
#' }
#'
#' @export
#'
#' @md
CorSpatialFeatures.default <- function (
    object,
    spatnet,
    across_all = FALSE,
    nCores = NULL,
    verbose = TRUE,
    ...
) {

  # Set global variables to NULL
  from <- to <- NULL

  if (verbose) cli_h2("Computing spatial autocorrelation")

  # Check objects
  if (!inherits(object, what = c("dgRMatrix", "dgCMatrix", "matrix", "data.frame"))) abort(glue("Invalid format of feature matrix: '{class(object)}'"))
  if (any(dim(object) == 0)) abort(glue("Invalid dimensions of feature matrix: '{paste(dim(object), collapse = 'x')}'"))
  if (!inherits(spatnet, what = "list")) abort(glue("Invalid format of spatnet: '{class(spatnet)}'"))

  # Convert object to sparse matrix
  if (!inherits(object, what = "dgCMatrix")) {
    object <- as(as.matrix(object), "dgCMatrix")
  }

  # get all spatial network spot barcode IDs
  all_spatial_network_spots <- Reduce(c, lapply(spatnet, function(x) x$from))

  # Check that objects match
  spots_in_spatnets <- !unique(all_spatial_network_spots) %in% rownames(object)
  if (any(spots_in_spatnets)) abort(glue("{sum(spots_in_spatnets)} spots in the spatial networks could not be found in the feature matrix.",
                                    "i" = "Make sure that the spatial networks share spot IDs with the feature matrix."))

  results <- lapply(seq_along(spatnet), function (i) {

    if (verbose) cli_alert_info("Sample {i}:")

    # pivot spatial network in long format to a wide format
    wide_spatial_network <- pivot_wider(spatnet[[i]] |> select(from, to) |> mutate(value = 1),
                                        names_from = "from", values_from = "value", values_fill = 0)

    # Convert wide spatial network to a matrix
    CN <- as.matrix(wide_spatial_network[, 2:ncol(wide_spatial_network)])
    rownames(CN) <- wide_spatial_network$to
    CN <- CN[colnames(CN), ]
    rm(wide_spatial_network)
    CN <- as(CN, "dgCMatrix")

    # Subset feature data to only include spots with neighbors
    x_subset <- object[colnames(CN), ]
    if (verbose) cli_alert("  Cleaned out spots without neighbors")

    # Calculate lag matrix
    lagMat <- (CN %*% x_subset) / rowSums(CN)
    if (verbose) cli_alert("  Computed feature lag expression")

    # return lagMat if the autocorrelation should be calculated across all samples
    if (across_all) {
      return(lagMat)
    }

    if (is.null(nCores)) {
      # Calculate spatial autocorrelation for each gene
      spatial_autocorrelation <- .colCors(x_subset, lagMat)
    } else {
      if (nCores > (detectCores() - 1)) {
        nCores <- detectCores() - 1
        cli_alert_info("Using {nCores} threads")
      }
      chunks <- ceiling((1:ncol(x_subset))/100)
      chunks <- split(1:ncol(x_subset), chunks)
      spatial_autocorrelation <- unlist(mclapply(seq_along(chunks), function(i) {
        inds <- chunks[[i]]
        .colCors(x_subset[, inds], lagMat[, inds])
      }, mc.cores = nCores))
    }

    if (verbose) cli_alert("  Computed feature spatial autocorrelation scores")

    # Summarize results
    results <- tibble(gene = names(spatial_autocorrelation), cor = spatial_autocorrelation) |>
      arrange(-cor)
  })

  # If across_all is set, calculate autocorrelations across all samples instead
  if (across_all) {
    if (verbose) cli_alert_info("Computing spatial autocorrelation scores across all samples")
    lagMat <- do.call(rbind, results)

    if (is.null(nCores)) {
      # Calculate spatial autocorrelation for each gene
      spatial_autocorrelation <- .colCors(object, lagMat)
    } else {
      if (nCores > (detectCores() - 1)) {
        nCores <- detectCores() - 1
        cli_alert_info("Using {nCores} threads")
      }
      chunks <- ceiling((1:ncol(object))/100)
      chunks <- split(1:ncol(object), chunks)
      spatial_autocorrelation <- unlist(mclapply(seq_along(chunks), function(i) {
        inds <- chunks[[i]]
        .colCors(object[, inds], lagMat[, inds])
      }, mc.cores = nCores))
    }

    if (verbose) cli_alert("  Computed feature spatial autocorrelation scores")
    # Summarize results
    results <- tibble(gene = colnames(object), cor = spatial_autocorrelation) |>
      arrange(-cor)
  }

  if (verbose) cli_alert_success("  Returning results")
  return(results)
}


#' @param features A character vector with features present in `Seurat` object. These
#' features need to be accessible with \code{\link{FetchData}}
#' @param assay_use Select assay to use for computation. If not specified, the default
#' assay will be used.
#' @param slot_use Select slot to use from assay object.
#'
#' @importFrom Seurat FetchData VariableFeatures GetAssayData
#' @importFrom rlang %||%
#'
#' @rdname cor-features
#'
#' @examples
#'
#' se_mbrain <- readRDS(system.file("extdata/mousebrain",
#'                                  "se_mbrain",
#'                                  package = "STUtility2"))
#' se_mbrain <- se_mbrain |>
#'   ScaleData() |>
#'   RunPCA()
#'
#' # Compute spatial autocorrelation for variable features
#' spatgenes <- CorSpatialFeatures(se_mbrain,
#'                                 features = VariableFeatures(se_mbrain),
#'                                 nCores = 1)
#'
#' # Check genes with highest spatial autocorrelation
#' head(spatgenes[[1]])
#'
#' # Note that the top variable genes are blood related (hemoglobin genes)
#' # These genes have lower spatial autocorrelation since blood vessels
#' # typically only cover a few spots and more randomly dispersed throughput the tissue
#' head(VariableFeatures(se_mbrain))
#'
#' # The same principle can be used to estimate spatial autocorrelation for other features,
#' # for example dimensionality reduction vectors
#' spatpcs <- CorSpatialFeatures(se_mbrain,
#'                              features = paste0("PC_", 1:10),
#'                              nCores = 1)
#'
#' # Calculate spatial autocorrelation scores for principal components
#' head(spatpcs[[1]])
#'
#' # Compute spatial autocorrelation scores for multiple datasets
#' se_mcolon <- readRDS(system.file("extdata/mousecolon",
#'                                  "se_mcolon",
#'                                  package = "STUtility2"))
#' se_merged <- MergeSTData(se_mbrain, se_mcolon) |>
#'   FindVariableFeatures()
#'
#' spatgenes <- CorSpatialFeatures(se_merged,
#'                                 features = VariableFeatures(se_merged),
#'                                 nCores = 1)
#'
#' # Check spatial autocorrelation scores mouse brain data
#' head(spatgenes[[1]])
#' # Check spatial autocorrelation scores mouse colon data
#' head(spatgenes[[2]])
#'
#' @export
#'
CorSpatialFeatures.Seurat <- function (
    object,
    features = NULL,
    assay_use = NULL,
    slot_use = "data",
    across_all = FALSE,
    nCores = NULL,
    verbose = TRUE,
    ...
) {

  # Validate Seurat object
  .check_seurat_object(object)

  # Get variable features if features=NULL
  features <- features %||% VariableFeatures(object)

  # Fetch features
  if (!is.null(assay_use)) {
    if (!requireNamespace("MatrixExtra", quietly = TRUE)) {
      install.packages("MatrixExtra")
    }
    stopifnot(is.character(assay_use),
              length(assay_use) == 1)
    featureMat <- GetAssayData(object, assay = assay_use, slot = slot_use)
    featureMat <- MatrixExtra::t(featureMat[features, ])
  } else {
    featureMat <- FetchData(object, vars = features)
  }

  # Obtain spatial networks
  spatnet <- GetSpatialNetwork(object)

  # Compute spatial autocorrelation
  spatfeatures <- CorSpatialFeatures(featureMat, spatnet, across_all, nCores, verbose, ...)

  return(spatfeatures)
}




#' Calculate pairwise correlation across two matrices
#'
#' @param x,y Numeric matrices with identical dimensions
#'
#' @importFrom Matrix colMeans
#'
#' @return A numeric vector with correlation scores
#'
#' @noRd
.colCors = function(x, y) {
  x = sweep(x, 2, colMeans(x))
  y = sweep(y, 2, colMeans(y))
  cor = colSums(x*y) / sqrt(colSums(x*x)*colSums(y*y))
  return(cor)
}
