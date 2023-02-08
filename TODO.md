# semla

## TODO

### External data

- [x] : An external data resource with mouse brain data. Contains spaceranger output files and a prepared Seurat object. These 
files can be accessed with `system.file("extdata/mousebrain", "path_to_file", package = "semla")` 
- [x] : An external data resource with mouse colon data. Contains spaceranger output files and a prepared Seurat object. These 
files can be accessed with `system.file("extdata/mousecolon", "path_to_file", package = "semla")`

### Data loading

- [x] `LoadAndMergeMatrices` : load and merge gene expression matrices
- [x] `LoadSpatialCoordinates` : load and merge spatial coordinates
- [x] `ReadVisiumData` : wrapper for `LoadAndMergeMatrices` and `LoadSpatialCoordinates` to load all data into a Seurat object. Dropped 
compatibility with "1k" and "2k" arrays. The function is similar to `InputFromTable`.
- [x] `LoadImages` : load H&E images into Seurat object

### Subset and Merge

- [x] `SubsetSTData` : for subsetting `Seurat` objects created with `semla`
- [x] `MergeSTData` : for merging `Seurat` objects created with `semla`

### On load

- [x] `.onAttach` : on load message when loading R package

### External data

- [x] _Visium mouse brain dataset_ : mouse brain test data 
- [x] _Visium mouse colon dataset_ : mouse colon test data

### Transformations

- [x] `CoordMirror` : mirror coordinates at predefined center
- [x] `CoordTransform` : apply rotations and translations to coordinates 
around predefined center
- [x] `ImageTranslate` : apply translations to image
- [x] `ImageTransform` : apply rotation and translation to image
- [x] `CoordAndImageTransform` : wrapper for `ImageTransform` and 
`CoordTransform` to apply transformations to paired image and coordinates. 
Can do rotations, translations and mirroring

### Image processing

- [x] `MaskImages` uses a simple method to identify the tissue section and mask the background implemented with 
`magick` 
- [x] `RigidTransformImages` takes a tibble with transformations and applies these transformations to the H&E 
images together with the spot coordinates.

### Tissue viewer

- [x] make OSD viewer from tile sourse
- [x] draw spots on top of tile source
- [x] lasso selection tool for manual annotation with paper JS
- [x] annotate spots and save selection
- [x] make tile function `TileImage` compatible with OSD
- [x] send data to OSD viewer from R or read from R on selection?

### Tissue aligner

- [x] make paper JS tool for image tranformations
- [x] communicate transformations to R
- [x] apply transformations to images in R object with `CoordAndImageTransform`
- [x] add side bar for image selection

### Digital unrolling

- [x] `CutSpatialNetwork` : Open a shiny/react application that can be used for "digital unrolling"

### Visualization tools
 
- [x] `FeatureOverlay` : this function has now been removed and instead we use `MapFeatures` which now takes an argument
`image_use` to inject images under the plot.
- [x] `DimOverlay` : this function has now been removed and instead we use `MapFeatures` which now takes an argument
`image_use` to inject images under the plot.
- [x] `MapFeatures` : similar to `ST.FeaturePlot` and `ST.DimPlot` but only works with numeric features such as
gene expression, PCs, QC metrics etc.
- [x] `MapLabels` : similar to `ST.FeaturePlot` but only works with categorical features such as
clusters, conditions etc.
- [x] `FactorGeneLoadingPlot` : Not sure if needed with the singlet R package available
- [x] `ImagePlot` : plots H&E images stored in a Seurat object
- [x] `MapFeaturesSummary` : create spatial maps with `MapFeatures` and add a summary plot next to it
- [x] `MapLabelsSummary` : create spatial maps with `MapLabels` and add a summary plot next to it
- [x] `AnglePlot` : draw a set of arrows on top of a selected region of spots. Mean to be used as a helper for selecting 
angle intervals in `RadialDistance` 

### Matrix factorization

- [x] `RunNMF` : no longer needed with the singlet R package available
- [x] `SeededNMF` : seed matrix factorization with a priori information fro single-cell data

### Spatial autocorelation

- [x] migrate `GetSpatNet` from STUtility
- [x] migrate `CorSpatialGenes` from STUtility and rename it to `CorSpatialFeatures`. This is a 
more general function that takes any type of feature as input, not just genes. This could for 
example be dimensionality reduction results instead.

### Spatial methods

- [x] migrate `RegionNeighbours` from STUtility
- [x] `RadialDistance` : calculates the radial distances from the borders of a region of interest
- [x] `DisconnectRegions` : separates spatially disconnected regions 
- [x] `NeighbourhoodEnrichmentTest` : Pairwise testing of cluster co-occurance. Start with checking nearest neighbors. Later, add possibility of measuring by distance(?).
- [x] `LabelAssortativityTest` : Inspired by Newman's Assortativity. Check whether selection of spots follow a dispersed or clustered spatial pattern.

### Staffli

- [x] migrate and modify STUtility S4 class object `Staffli`

### semla website

- [x] _getting started_
- [x] _visualization_
- [x] _spatial methods_
