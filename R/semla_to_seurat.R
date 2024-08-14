#' @include checks.R
#'
NULL

#' Update a Seurat object created with \code{semla} for compatibility with \code{Seurat}'s spatial functions
#' 
#' @param object An object of class \code{Seurat} with Visium data, created with \code{semla}
#' @param image_use Strig specifying image type to use, either "raw" or "transformed"
#' @param verbose Print messages
#' 
#' @import dplyr
#' @import cli
#' @import rlang
#' @import glue
#' 
#' @return A \code{Seurat} object compatible with with \code{Seurat}'s spatial functions
#' 
#' @examples 
#' 
#' library(semla)
#' 
#' # Load example Visium data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
#' se_merged <- MergeSTData(se_mbrain, se_mcolon)
#' 
#' # Update object
#' se_merged <- UpdateSeuratFromSemla(se_merged)
#' 
#' # Use Seurat functions
#' SpatialFeaturePlot(se_merged, features = "Nrgn")
#' 
#' @export
UpdateSeuratFromSemla <- function (
    object,
    image_use = "raw",
    verbose = TRUE
) {
  
  if (verbose) cli_h3("Adding @images data from {col_br_magenta('Staffli')} object")
  if (verbose) cat_line()
  
  # Set global variables to NULL
  pxl_row_in_fullres <- pxl_col_in_fullres <- sampleID <- NULL
  
  # Check object
  if (!inherits(object, what = "Seurat")) abort(glue("Invalid class '{class(object)}'. Expected a 'Seurat' object."))
  .check_seurat_object(object)
  
  image_use <- match.arg(image_use, choices = c("raw", "transformed"))
  
  # Check if the transformed image exists
  if (!image_use %in% names(st_object@rasterlists))
    abort(glue("Transformed images are not available in this object."))
  
  # Get Staffli object
  st_object <- GetStaffli(object)
  
  # Add step for test data
  if (any(st_object@imgs %in% c("mousebrain", "mousecolon"))) {
    object <- LoadImages(object, verbose = FALSE)
    st_object <- GetStaffli(object)
  }
  
  # Check image paths
  for (impath in st_object@imgs) {
    if (!file.exists(impath)) {
      abort(glue("Invalid image path '{impath}'"))
    }
  }
  
  # Load images if using raw. If transformed, reverse from raster to 3D numeric array
  if (image_use == "raw") {
    if (verbose) cli_alert_info("Loading images")
    imgs <- lapply(st_object@imgs, function(im) {
      ar <- image_read(im)[[1]] |> as.integer()
      ar <- ar/max(ar)
      return(ar)
    })
  } else if (image_use == "transformed") {
    if (verbose) cli_alert_info("Converting hexadecimal colors to RGB")
    imgs <- lapply(object@tools[["Staffli"]]@rasterlists[["transformed"]], function(im) {
      # retrieve dimensions
      w <- ncol(im)
      h <- nrow(im)
      # flatten
      hex_flat <- as.vector(im)
      # hex to rgb 
      rgb_mat <- grDevices::col2rgb(hex_flat)
      # reshape and normalize
      ar <- array(rgb_mat, dim = c(3, w, h))
      ar <- aperm(ar, c(3,2,1))
      ar <- ar / 255
    })
  }
  
  # Which coordinates do we put into the Seurat object
  if (image_use == "raw") {
    rows <- "pxl_row_in_fullres"
    cols <- "pxl_col_in_fullres"
  } else {
    rows <- "pxl_row_in_fullres_transformed"
    cols <- "pxl_col_in_fullres_transformed"
  }
  
  # Convert coordinates
  if (verbose) cli_alert_info(paste0("Converting ", image_use, " coordinates"))
  coords <- GetCoordinates(object) |> 
    rename(imagerow = all_of(rows), imagecol = all_of(cols)) |> 
    group_by(sampleID) |> 
    group_split()
  coords <- setNames(coords, nm = paste0(1:length(coords)))
  coords <- lapply(coords, function(xy) {
    xy |> select(-sampleID) |> data.frame(row.names = 1)
  })
  
  # Convert scalefactors. If using transformed image, rescale scalefactors to this image
  if (verbose) cli_alert_info("Converting scale factors")
  if (image_use == "transformed") {
    scalefactors <- GetScaleFactors(object) |> 
      cbind(GetImageInfo(object) |> dplyr::select(!sampleID)) |>
      mutate(sampleID = as.integer(sampleID)) |> 
      group_by(sampleID) |> 
      mutate(f_max = max(full_width, full_height),
             h_max = max(width, height),
             t_max = sapply(imgs, function(x) max(dim(x))),
             tissue_hires_scalef = t_max / f_max, # if you are working with the transformed image, we only need one scalefactor, but since I dont know how seurat decides which one to use, I am providing both scalef with the same value
             tissue_lowres_scalef = t_max / f_max) |>
      select(all_of(colnames(GetScaleFactors(object)))) |>
      group_split()
  } else if (image_use == "raw") {
    scalefactors <- GetScaleFactors(object) |> 
      mutate(sampleID = as.integer(sampleID)) |> 
      group_by(sampleID) |> 
      group_split()
  }
  scalefactors <- setNames(scalefactors, nm = paste0(1:length(scalefactors)))
  
  spot_radius_list <- lapply(seq_along(scalefactors), function(i) {
    unnormalized.radius <- scalefactors[[i]]$fiducial_diameter_fullres * 
      scalefactors[[i]]$tissue_lowres_scalef
    spot.radius <- unnormalized.radius/max(dim(x = imgs[[i]]))
    return(spot.radius)
  })
  
  # Assemble all data
  if (verbose) cli_alert_info("Creating VisiumV1 slices")
  visium_data <- lapply(seq_along(scalefactors), function(i) {
    scale.factors <- scalefactors[[i]]
    slice <- new(Class = "VisiumV1", 
                 image = imgs[[i]], 
                 scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                              fiducial = scale.factors$fiducial_diameter_fullres, 
                                              hires = scale.factors$tissue_hires_scalef, 
                                              scale.factors$tissue_lowres_scalef), 
                 coordinates = coords[[i]], 
                 spot.radius = spot_radius_list[[i]])
    DefaultAssay(slice) <- DefaultAssay(object)
    return(slice)
  }) |> setNames(nm = paste0("slice", 1:length(coords)))
  
  # Place slices in Seurat object
  if (verbose) cli_alert_info("Storing VisiumV1 slices in @images slot")
  for (nm in names(visium_data)) {
    object[[nm]] <- visium_data[[nm]]
  }
  
  if (verbose) cat_line()
  if (verbose) cli_alert_success("Returning updated {col_br_magenta('Seurat')} object.")
  return(object)
}

