% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/spatial_autocorrelation.R
\name{CorSpatialFeatures}
\alias{CorSpatialFeatures}
\alias{CorSpatialFeatures.default}
\alias{CorSpatialFeatures.Seurat}
\title{Find features with high spatial autocorrelation}
\usage{
CorSpatialFeatures(object, ...)

\method{CorSpatialFeatures}{default}(
  object,
  spatnet,
  across_all = FALSE,
  calculate_pvalue = FALSE,
  nCores = NULL,
  verbose = TRUE,
  ...
)

\method{CorSpatialFeatures}{Seurat}(
  object,
  features = NULL,
  assay_use = NULL,
  slot_use = "data",
  across_all = FALSE,
  calculate_pvalue = FALSE,
  nCores = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{An object (see details)}

\item{...}{Arguments passed to other methods}

\item{spatnet}{A list of spatial networks created with \code{\link{GetSpatialNetwork}}. The spots in these
networks should match the spots in the feature matrix.}

\item{across_all}{Boolean specifying if autocorrelation scores be calculated across all samples (default
FALSE).}

\item{calculate_pvalue}{Boolean specifying if p-values should be calculated for each correlation value
(default FALSE). Raw and BH-adjusted p-values provided.}

\item{nCores}{Number of cores to use for the spatial autocorrelation calculation}

\item{verbose}{Print messages}

\item{features}{A character vector with features present in \code{Seurat} object. These
features need to be accessible with \code{\link{FetchData}}}

\item{assay_use}{Select assay to use for computation. If not specified, the default
assay will be used.}

\item{slot_use}{Select slot to use from assay object.}
}
\value{
Either a list of tibbles or a tibble with feature names and correlation
scores. See section \strong{Mode} for details.
}
\description{
This function can be used to find genes with high spatial autocorrelation in SRT data.
A more detailed description of the algorithm is outlined below.
}
\section{Spatial autocorrelation}{

Spatial autocorrelation is the term used to describe the presence of systematic spatial
variation. Positive spatial autocorrelation of a feature is the tendency for regions that
are close together in space to have similar values for that feature.

A simple example is when you have an anatomical structure or a tissue type that spans
across multiple neighboring spots in an SRT experiment, for example a gland, an immune
infiltrate or a region of the brain. Inside such structures, you might find that the
expression levels of certain genes (or other features) are highly similar and hence
these genes have a positive spatial autocorrelation.

The method provided in \code{semla} works as follows. For each feature and spot,
the expression is averaged across all neighboring spots (typically the 6 closest neighbors)
to produce a lag expression vector. Since this vector represents the average of the surrounding
spots, we can use it to test if the expression in those spots is similar to the center spot.
One simple strategy is to calculate the pearson correlation between a genes' lag vector and
the original expression vector which typically captures the spatial autocorrelation well.
}

\section{Method steps}{

\itemize{
\item{Load a matrix with features in rows and spots in columns: \eqn{X_{expr}}}
\item{Convert the corresponding spatial network to wide format and construct a nearest
neighbor matrix \eqn{N_{neighbors}} in which neighboring spots have a value of 1
and the remaining spots have a value of 0
}
\item{\eqn{N_{neighbors}} is then multiplied with the \eqn{X_{expr}} to
calculate a lag vector for each feature: \cr \cr
\eqn{X_{lagexpr} = (N_{neighbors}*X_{expr})/n_{neighbors}} \cr \cr
where \eqn{n_{neighbors}} is the number of neighbors for each spot.
}
\item{The spatial autocorrelation score for a feature is the 'pearson' correlation of the
lag vector and the initial expression vector: \cr \cr
\eqn{spatcor_{feature} = cor(X_{lagexpr}[feature, ], X_{expr}[feature, ])}
}
}
}

\section{Default method}{

The default method expects a matrix-like object with features in columns and spots in rows
and a list of spatial networks generated with \code{\link{GetSpatialNetwork}}.
}

\section{Mode}{

If \code{across_all} is set to \code{TRUE}, the spatial autocorrelation scores will be computed
across all samples. Otherwise, the scores will be calculated for each sample separately, and returns
a list with one \code{tibble} per sample.
}

\examples{
\donttest{
se_mbrain <- readRDS(system.file("extdata/mousebrain",
                                 "se_mbrain",
                                 package = "semla"))
featureMat <- FetchData(se_mbrain, vars = VariableFeatures(se_mbrain)[1:100])

coordfile <-
  system.file("extdata/mousebrain/spatial",
              "tissue_positions_list.csv",
              package = "semla")

# Load coordinate data into a tibble
xys <- setNames(read.csv(coordfile, header = FALSE),
                nm = c("barcode", "selection", "grid_y", "grid_x", "y", "x"))
xys$sampleID <- 1L
xys <- xys |>
  dplyr::filter(selection == 1) |>
  dplyr::select(barcode, x, y, sampleID) |>
  tibble::as_tibble()

# Create spatial networks
spatnet <- GetSpatialNetwork(xys)
spatgenes <- CorSpatialFeatures(featureMat, spatnet, nCores = 1)

# Check genes with highest spatial autocorrelation
head(spatgenes[[1]])
}


se_mbrain <- readRDS(system.file("extdata/mousebrain",
                                 "se_mbrain",
                                 package = "semla"))
se_mbrain <- se_mbrain |>
  ScaleData() |>
  RunPCA()

# Compute spatial autocorrelation for variable features
spatgenes <- CorSpatialFeatures(se_mbrain,
                                features = VariableFeatures(se_mbrain),
                                nCores = 1)

# Check genes with highest spatial autocorrelation
head(spatgenes[[1]])

# Note that the top variable genes are blood related (hemoglobin genes)
# These genes have lower spatial autocorrelation since blood vessels
# typically only cover a few spots and more randomly dispersed throughput the tissue
head(VariableFeatures(se_mbrain))

\donttest{
# The same principle can be used to estimate spatial autocorrelation for other features,
# for example dimensionality reduction vectors
spatpcs <- CorSpatialFeatures(se_mbrain,
                             features = paste0("PC_", 1:10),
                             nCores = 1)

# Calculate spatial autocorrelation scores for principal components
head(spatpcs[[1]])

# Compute spatial autocorrelation scores for multiple datasets
se_mcolon <- readRDS(system.file("extdata/mousecolon",
                                 "se_mcolon",
                                 package = "semla"))
se_merged <- MergeSTData(se_mbrain, se_mcolon) |>
  FindVariableFeatures()

spatgenes <- CorSpatialFeatures(se_merged,
                                features = VariableFeatures(se_merged),
                                nCores = 1)

# Check spatial autocorrelation scores mouse brain data
head(spatgenes[[1]])
# Check spatial autocorrelation scores mouse colon data
head(spatgenes[[2]])
}

}
\seealso{
Other network-methods: 
\code{\link{GetSpatialNetwork}()}

Other spatial-methods: 
\code{\link{CutSpatialNetwork}()},
\code{\link{DisconnectRegions}()},
\code{\link{GetSpatialNetwork}()},
\code{\link{RadialDistance}()},
\code{\link{RegionNeighbors}()},
\code{\link{RunLabelAssortativityTest}()},
\code{\link{RunLocalG}()},
\code{\link{RunNeighborhoodEnrichmentTest}()}
}
\author{
Ludvig Larsson
}
\concept{network-methods}
\concept{spatial-methods}
