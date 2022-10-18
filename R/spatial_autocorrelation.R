#' @include generics.R
#'
NULL


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
#' @importFrom dplyr mutate select
#' @importFrom glue glue
#' @importFrom rlang abort
#' @importFrom methods as
#'
#' @examples
#' \dontrun{
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "STUtility2"))
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
#' spatgenes <- CorSpatialFeatures(featureMat, spatnet)
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
#' dimMat <- se_mbrain |>
#' ScaleData() |>
#'   RunPCA() |>
#'   FetchData(vars = paste0("PC_", 1:10))
#'
#' # Calculate spatial autocorrelation scores for principal components
#' spatPCs <- CorSpatialFeatures(dimMat, spatnet)
#' head(spatPCs)
#' }
#'
#' @export
#'
#' @md
CorSpatialFeatures.default <- function (
    object,
    spatnet,
    across_all = FALSE,
    nCores = detectCores() - 1,
    verbose = TRUE,
    ...
) {

  # Set global variables to NULL
  from <- to <- NULL

  if (verbose) inform(c("i" = "Checking objects..."))

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
                                    "i" = "Make sure that the spatial networks were generated from the same dataset as feature matrix."))
  if (verbose) inform(c("v" = "Passed check!"))

  results <- lapply(seq_along(spatnet), function (i) {

    if (verbose) inform(c("", "i" = glue("Sample {i}:")))

    # pivot spatial network in long format to a wide format
    wide_spatial_network <- pivot_wider(spatnet[[i]] |> select(from, to) |> mutate(value = 1),
                                        names_from = "from", values_from = "value", values_fill = 0)
    if (verbose) inform(c("v" = "  Successfully reformatted spatial network from long to wide format"))

    # Convert wide spatial network to a matrix
    CN <- as.matrix(wide_spatial_network[, 2:ncol(wide_spatial_network)])
    rownames(CN) <- wide_spatial_network$to
    CN <- CN[colnames(CN), ]
    rm(wide_spatial_network)
    CN <- as(CN, "dgCMatrix")

    # Subset feature data to only include spots with neighbors
    x_subset <- object[colnames(CN), ]
    if (verbose) inform(c("v" = "  Cleaned out spots without neighbors"))

    # Calculate lag matrix
    lagMat <- (CN %*% x_subset) / rowSums(CN)
    if (verbose) inform(c("v" = "  Computed lag expression for all features"))

    # return lagMat if the autocorrelation should be calculated across all samples
    if (across_all) {
      return(lagMat)
    }

    # Calculate spatial autocorrelation for each gene
    spatial_autocorrelation <- unlist(mclapply(1:ncol(x_subset), function(i) {
      cor(as.numeric(x_subset[rownames(lagMat), i]), as.numeric(lagMat[, i]))
    }, mc.cores = nCores))
    if (verbose) inform(c("v" = "  Computed spatial autocorrelation scores for features"))

    # Summarize results
    res <- tibble(gene = colnames(x_subset), cor = spatial_autocorrelation) |>
      arrange(-spatial_autocorrelation)
    if (verbose) inform(c("v" = "  Returning results"))

    return(res)
  })

  # If across_all is set, calculate autocorrelations across all samples instead
  if (across_all) {
    if (verbose) inform(c("", "i" = "Computing spatial autocorrelation scores across all samples"))
    lagMat <- do.call(rbind, results)
    # Calculate spatial autocorrelation for each gene
    spatial_autocorrelation <- unlist(mclapply(1:ncol(object), function(i) {
      cor(as.numeric(object[rownames(lagMat), i]), as.numeric(lagMat[, i]))
    }, mc.cores = nCores))
    if (verbose) inform(c("v" = "Computation finished"))
    # Summarize results
    results <- tibble(gene = colnames(object), cor = spatial_autocorrelation) |>
      arrange(-spatial_autocorrelation)
  }

  return(results)
}
