# semla 1.1.6

*2023-09-19*

## Bug fixes

-   Fixed bug in `RunLocalG` when `store_in_metadata = FALSE` and no alternative hypothesis is provided.
-   Fixed bug in `TileImage` when using multiple threads on linux.

## Changes

-   Modified behaviour of `mode = "heatmap"` in `PlotFeatureLoadings` to make sure that features are ordered by their highest dimensionality reduction value.
-   Updated NNMF vignette with a new data set.
-   Changed deprecated `ggplot2` function arguments in `AnglePlot` and `MapFeaturesSummary`.
-   Changed example `Seurat` object to contain basic array coordinates.

## Added

-   Option to set `launch.browser` in `FeatureViewer`. For instance, users can se `launch.browser = getOption("viewer")` in `FeatureViewer` to launch the app in the RStudio viewer window.
-   `CreateStaffliObject` now accepts basic array coordinates to be stored in `meta_data` slot. These coordinates are defined on the capture array grid of a 10x Genomics Visium slide and do not match the H&E images.

# semla 1.1.5

*2023-09-06*

## Bug fixes

-   Added check for non-finite values when visualizing numeric data in `FeatureViewer` to avoid crash
-   Fixed bug in `FeatureViewer` when a slot other than 'data' is selected

# semla 1.1.4

*2023-09-04*

## Bug fixes

-   Fixed bug in `LoadScaleFactors`
-   Fixed bug in `FeatureViewer` when handling data sets with spots located outside the H&E image

## Changes

-   Improved function documentation

## Added

-   Added support for conversion of `Seurat` objects with Slide-seq data in `UpdateSeuratForSemla`
-   Added option to color background and titles in `ImagePlot`
-   Added support for Visium + IF data in `LoadAndMergeMatrices` and `ReadVisiumData`
-   Added `renv.lock` file for installing semla and its dependencies with `renv`

## Website

-   Updated image alignment vignette
-   Updated Staffli objects article
-   Added cell type deconvolution benchmark article
-   Added Visium + IF data vignette
-   Added example to digital unrolling vignette

# semla 1.1.0

*2023-08-14*

## Added

-   Added helper function for `Staffli` objects:
    -   `GetScaleFactors` to fetch `tibble` with scale factor related information
-   `GetImages` to fetch images
-   `GetImageInfo` to fetch `tibble` with image related information
-   `ReplaceImagePaths` to validate and update image paths
-   Added `UpdateSeuratForSemla` to make `Seurat` objects with VisiumV1 data compatible with `semla`
-   Added `LoadAnnotationCSV` to load annotations exported from Loupe Browser

## Changes

-   Improved verbosity
-   Reduced size of example data
-   Updated function documentation
-   Updated function examples
-   Added option to load images from URL instead of local path (`LoadImages` and `ReplaceImagePaths`)
-   Added option to use multiple cores for `ExportDataForViewer`
-   Fixed bug in `RadialDistance` when no spots are available after removing singletons
-   Fixed bug in `MergeSTData` when merging two or more objects containing multiple tissue sections
-   Fixed bug in `MapMultipleFeatures` when NA values are present
-   Updated data loaders for increased compatibility with CytAssist data. Minor bug fixes for edge cases when spots are located outside of the H&E image
-   Fixed bug in `LoadAndMergeMatrices` when loading data from a data folder output by Space Ranger
-   Changed output of `AdjustTissueCoordinates` which now returns unadjusted y-axis values

# semla 1.0.0

*2023-03-24*

We are excited to announce the first release of our R package, `semla`.

The package is designed to process, analyze and visualize Spatially Resolved Transcriptomics (SRT) data.

Key features of this release include:

-   Interactive data visualization and annotation with the `FeatureViewer` web app
-   Spatially aware analysis methods: `RegionNeighbors`, `RadialDistance`, `CorSpatialFeatures`, `RunLocalG`, `RunLabelAssortativityTest`, `RunNeighborhoodEnrichmentTest`, `CutSpatialNetwork`
-   Visualization methods: `MapFeatures`, `MapLabels`
-   Image processing: `MaskImages`, `RigidTransformImages`, `RunAlignment`
-   Cell type mapping: `RunNNLS`

We would like to express our gratitude to [Javier Morlanes](https://github.com/jemorlanes) and [Marcos Machado](https://github.com/mabreumac) who contributed to this release. We appreciate your feedback, suggestions, and bug reports.

We invite users to download and install the development version of the package from our GitHub [repo](https://github.com/ludvigla/semla) and explore the features and functionality. Please feel free to provide feedback and suggestions for future improvements.

Thank you for your support and we look forward to continuing to enhance and improve this package.

*Sincerely,*

*Ludvig Larsson and Lovisa Franzén*
