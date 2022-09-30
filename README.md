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
- [ ] `CoordAndImageTransform` : wrapper for `ImageTransform` and 
`CoordTransform` to apply transformations to paired image and coordinates. 
Need to fix add mirror as an option.

### Image processing

- [ ] option to mask images?

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

- [ ] `FeatureOverlay` : modify STUtility function, fix axes
- [ ] `DimOverlay` : modify STUtility function, fix axes
- [ ] `ST.FeaturePlot` : modify STUtility function, fix axes`. Change name?
- [ ] `ST.DimPlot` : modify STUtility function, fix axes. Change name?
- [ ] `FactorGeneLoadingPlot` : migrate STUtility function. Name change?

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

### Neighborhood analysis

- [ ] : migrate `RegionNeighbours` from STUtility

### Staffli

- [ ] : migrate and modify STUtility S4 class object `Staffli`

### STUtility2 website

- [ ] : _getting started_
- [ ] : _visualization_
- [ ] : _matrix factorization_

### Advanced features

- [ ] : crop data based on image rectangle
- [ ] : switch resolution. Perhaps not necessary with the Tissue Viewer?
