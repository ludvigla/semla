# STUtility2

## TODO

### External data

- [x] : An external data resource with mouse brain data. Contains spaceranger output files and a prepared Seurat object. These 
files can be accessed with `system.file("extdata/mousebrain", "path_to_file", package = "STUtility2")` 
- [x] : An external data resource with mouse colon data. Contains spaceranger output files and a prepared Seurat object. These 
files can be accessed with `system.file("extdata/mousecolon", "path_to_file", package = "STUtility2")`

### Data loading

- [x] `LoadAndMergeMatrices` : load and merge gene expression matrices
- [x] `LoadSpatialCoordinates` : load and merge spatial coordinates
- [x] `ReadVisiumData` : wrapper for `LoadAndMergeMatrices` and `LoadSpatialCoordinates` to load all data into a Seurat object. Dropped 
compatibility with "1k" and "2k" arrays. The function is similar to `InputFromTable`.
- [x] `LoadImages` : load H&E images into Seurat object

### Subset and Merge

- [x] : Add working examples including image data
- [x] : `SubsetSTData` for subsetting `Seurat` objects created with `STUtility2`
- [x] : `MergeSTData` for merging `Seurat` objects created with `STUtility2`

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

- [x] option to mask images? 
- [ ] more examples to highlight flexibility

### Tissue viewer

- [x] make OSD viewer from tile sourse
- [x] draw spots on top of tile source
- [x] lasso selection tool for manual annotation with paper JS
- [ ] annotate spots and save selection
- [ ] make tile function `TileImage` compatible with OSD
- [ ] send data to OSD viewer from R or read from R on selection?
- [ ] add option to select samples in viewer (shiny or react?)
- [ ] make R function that calls react app

### Tissue aligner

- [x] make paper JS tool for image tranformations
- [x] communicate transformations to R
- [ ] apply transformations to images in R object with `CoordAndImageTransform`
- [ ] add side bar for image selection

### Visualization tools
 
- [x] `FeatureOverlay` : this function has now been removed and instead we use `MapFeatures` which now takes an argument
`image_use` to inject images under the plot.
- [x] `DimOverlay` : this function has now been removed and instead we use `MapFeatures` which now takes an argument
`image_use` to inject images under the plot.
- [x] `MapFeatures` : similar to `ST.FeaturePlot` and `ST.DimPlot` but only works with numeric features such as
gene expression, PCs, QC metrics etc.
- [x] `MapLabels` : similar to `ST.FeaturePlot` but only works with categorical features such as
clusters, conditions etc.
- [ ] `FactorGeneLoadingPlot` : migrate STUtility function. Name change?
- [x] `ImagePlot` : plots H&E images stored in a Seurat object
- [x] `MapFeaturesSummary` : create spatial maps with `MapFeatures` and add a summary plot next to it

### Matrix factorization

- [ ] `RunNMF` : modify STUtility function to work with RcppML
- [x] `SeededNMF` : seed matrix factorization with a priori information fro single-cell data
- [ ] `RunMixedNMF` : under development

### Spatial autocorelation

- [x] : migrate `GetSpatNet` from STUtility
- [x] : add working example to `GetSpatNet`
- [x] : migrate `CorSpatialGenes` from STUtility and rename it to `CorSpatialFeatures`. This is a 
more general function that takes any type of feature as input, not just genes. This could for 
example be dimensionality reduction results instead.
- [x] : add working example to `CorSpatialFeatures`. Here I have provided an example for how to 
calculate spatial autocorrelation of genes and for principal components.

### Spatial methods

- [x] : `RadialDistance` can be used to calculate distances from the border of a selected region. 
This can for example be useful when looking at the increase/decrease in expression as a function of 
distance from a predefined region. Negative values points toward the center of the selected region.
- [x] : `DisconnectComponents` can be used to disconnect spots belonging to a certain region that 
are spatially disconnected. This can for example be useful if you want to analyze disconnected a 
components of a particular tissue structure. 
- [ ] : Add examples in vignettes

### Neighborhood analysis

- [x] : migrate `RegionNeighbours` from STUtility
- [ ] : `NeighbourhoodEnrichmentTest` : under development. Pairwise testing of cluster co-occurance. Start with checking nearest neighbors. Later, add possibility of measuring by distance(?).
- [ ] : `GroupClusteringTest` : (function name??) under development. Inspired by Newman's Assortativity. Check whether selection of spots follow a dispersed or clustered spatial pattern.

### Staffli

- [x] : migrate and modify STUtility S4 class object `Staffli`

### STUtility2 website

- [x] : _getting started_
- [x] : _visualization_
- [ ] : _matrix factorization_

### Advanced features

- [ ] : crop data based on image rectangle
- [ ] : switch resolution. Perhaps not necessary with the Tissue Viewer?
