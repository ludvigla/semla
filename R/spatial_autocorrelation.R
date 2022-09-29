#' Find features with high spatial autocorrelation
#'
#' This function can be used to find genes with high spatial autocorrelation in SRT data.
#' A more detailed description of the algorithm is outlined in the Details section below.
#'
#' Method steps:
#' \itemize{
#'    \item{Load a matrix with features in rows and spots in columns: \eqn{X_{expr}}}
#'    \item{Convert the corresponding spatial network to wide format and construct a nearest
#'    neighbor matrix \eqn{N_{neighbors}} in which neighboring spots have a value of 1
#'    and the remaining spots have a value of 0
#'    }
#'    \item{\eqn{N_{neighbors}} is then multiplied with the \eqn{X_{expr}} to
#'    calculate a lag vector for each gene: \cr \cr
#'    \eqn{X_{lagexpr} = (N_{neighbors}*X_{expr})/n_{neighbors}} \cr \cr
#'    where \eqn{n_{neighbors}} is the number of neighbors for each spot.
#'    }
#'    \item{The spatial autocorrelation score for a genes is the 'pearson' correlation of the
#'    lag vector and the initial expression vector: \cr \cr
#'    \eqn{spatcor_{gene} = cor(X_{lagexpr}[gene, ], X_{expr}[gene, ])}
#'    }
#' }
#'
#'
#' @param x A matrix-like object with features in columns and spots in rows
#' @param spatnet A list of spatial networks created with `GetSpatialNetwork`. The spots in these
#' networks should match the spots in `x`.
#' @param across_all Should the autocorrelation scores be calculated across all samples?
#' @param verbose Print messages
#'
#' @family network-methods
#'
#' @inheritParams GetSpatialNetwork
#'
#' @return Either a list of tibbles or a tibble with gene names and correlation scores for each gene in `x`
#'
#' @importFrom stats cor
#' @importFrom parallel mclapply detectCores
#' @importFrom tibble tibble
#' @importFrom Matrix rowSums
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr mutate select
#' @importFrom glue glue
#' @importFrom rlang abort
#'
#' @export
#'
#' @md
CorSpatialGenes <- function (
    x,
    spatnet,
    across_all = FALSE,
    nCores = detectCores() - 1,
    verbose = TRUE
) {

  if (verbose) inform(c("i" = "Checking objects..."))

  # Check objects
  if (!class(x) %in% c("dgRMatrix", "dgCMatrix", "matrix", "data.frame")) abort(glue("Invalid format of x: '{class(x)}'"))
  if (any(dim(x) == 0)) abort(glue("Invalid dimensions of x: '{paste(dim(x), collapse = 'x')}'"))
  if (!class(spatnet) == "list") abort(glue("Invalid format of spatnet: '{class(spatnet)}'"))

  # get all spatial network spot barcode IDs
  all_spatial_network_spots <- Reduce(c, lapply(spatnet, function(x) x$from))

  # Check that objects match
  spots_in_spatnets <- !unique(all_spatial_network_spots) %in% rownames(x)
  if (any(spots_in_spatnets)) abort(glue("{sum(spots_in_spatnets)} spots in the spatial networks could not be found in x.",
                                    "i" = "Make sure that the spatial networks were generated from the same dataset as x."))
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
    x_subset <- x[colnames(CN), ]
    if (verbose) inform(c("v" = "  Cleaned out spots without neighbors"))

    # Calculate lag matrix
    lagMat <- (CN %*% x_subset) / rowSums(CN)
    if (verbose) inform(c("v" = "  Computed lag expression for all features"))

    # return lagMat if the autocorrelation should be calculated across all samples
    if (across_all) {
      return(lagMat)
    }

    # Calculate spatial autocorrelation for each gene
    spatial_autocorrelation <- unlist(mclapply(1:ncol(x), function(i) {
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
    spatial_autocorrelation <- unlist(mclapply(1:ncol(x), function(i) {
      cor(as.numeric(x[rownames(lagMat), i]), as.numeric(lagMat[, i]))
    }, mc.cores = nCores))
    if (verbose) inform(c("v" = "Computation finished"))
    # Summarize results
    results <- tibble(gene = colnames(x), cor = spatial_autocorrelation) |>
      arrange(-spatial_autocorrelation)
  }

  return(results)
}
