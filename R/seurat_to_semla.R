#' @include checks.R
#'
NULL

#' Update a Seurat object for compatibility with semla
#' 
#' In order to make a \code{Seurat} object with \code{semla}, 
#' you can use \code{AddSemlaToSeurat} to add the required 
#' data to the \code{Seurat} object. This assumes that 
#' the \code{Seurat} object contains one or more "VisiumV1" 
#' object(s) in the \code{images} slot. 
#' 
#' @details 
#' Note that you need to specify what H&E image was loaded, 
#' one of "tissue_lowres" or "tissue_hires". If this argument 
#' is incorrect, the tissue coordinates will be misplaced.
#' 
#' Visit the \href{https://ludvigla.github.io/semla/articles/getting_started.html}{getting started}
#' tutorial on our package website for more details on how to 
#' use this function.
#' 
#' @param object An object of class \code{Seurat} with Visium data
#' @param image_type One of "tissue_lowres" or "tissue_hires", specifying
#' what H&E image was loaded into the \code{Seurat} object
#' @param verbose Print messages
#' 
#' @return A \code{Seurat} object compatible with semla
#' 
#' @import rlang
#' @import glue
#' @import cli
#' @import dplyr
#' @import tibble
#' @importFrom magick image_read image_info
#' 
#' @examples 
#' 
#' \dontrun{
#' library(semla)
#' library(SeuratData)
#' 
#' InstallData("stxBrain")
#' 
#' # Load example Seurat object
#' brain <- LoadData("stxBrain", type = "anterior1")
#' 
#' # Make Seurat object compatible with semla
#' brain_semla <- AddSemlaToSeurat(brain)
#' }
#' 
#' @export
AddSemlaToSeurat <- function (
    object,
    image_type = c("tissue_lowres", "tissue_hires"),
    verbose = TRUE
) {
  # Check object
  if (!inherits(object, what = "Seurat")) abort(glue("Invalid class '{class(object)}'. Expected a 'Seurat' object."))
  if (length(object@images) == 0) abort(glue("No images available in 'Seurat' object."))
  
  image_type <- match.arg(arg = image_type, choices = c("tissue_lowres", "tissue_hires"))
  if (verbose) cli_alert_info("Expecting '{image_type}' format.")
  
  # Get images
  imgs <- lapply(object@images, function(slice) {
    slice@image |> as.raster()
  })
  
  # Export images
  temp_dir <- tempdir()
  image_paths <- c()
  img_info_tibble <- tibble()
  scalefactors_tibble <- tibble()
  coordinates_tibble <- tibble()
  for (i in seq_along(imgs)) {
    
    # Check if the slice represent Visium data
    if (!inherits(object@images$anterior1, what = "VisiumV1")) abort("Only 'VisiumV1' data is supported.")
    
    # Collect image paths
    if (verbose) cli_alert_info("Collecting data for sample {i}:")
    dims <- dim(imgs[[i]])
    image_path <- paste0(temp_dir, "/", i, ".png")
    if (verbose) cli_alert("  Exporting H&E image width dimensions {dims[1]}x{dims[2]} to {image_path}")
    png(filename = image_path, height = dims[1], width = dims[2])
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar)) 
    par(mar = c(0, 0, 0, 0))
    plot(imgs[[i]])
    dev.off()
    image_paths <- c(image_paths, image_path)
    
    # Collect scale factors
    sf <- object@images[[i]]@scale.factors
    sf <- tibble(spot_diameter_fullres = sf$spot,
                 tissue_hires_scalef = sf$hires,
                 fiducial_diameter_fullres = sf$fiducial,
                 tissue_lowres_scalef = sf$lowres,
                 sampleID = paste0(i))
    if (verbose) cli_alert("  Collected scale factors.")
    scalefactors_tibble <- bind_rows(scalefactors_tibble, sf)
    
    # Collect image information
    sf_image <- ifelse(image_type == "tissue_lowres", sf$tissue_lowres_scalef, sf$tissue_hires_scalef)
    iminfo <- image_read(image_path) |> 
      image_info() |> 
      mutate(full_width = ceiling(dims[1]/sf_image), full_height = ceiling(dims[2]/sf_image),
             sampleID = paste0(i), type = image_type) |> 
      select(format, width, height, full_width, full_height, colorspace, filesize, density, sampleID, type)
    if (verbose) cli_alert("  Collected image information.")
    img_info_tibble <- bind_rows(img_info_tibble, iminfo)
    
    # Collect spatial coordinates
    slice <- object@images[[i]]
    x <- slice@coordinates
    coordinates <- tibble(barcode = paste0(gsub(pattern = '[0-9]+', replacement = "", x = rownames(x)), i), 
           pxl_col_in_fullres = x[, "imagecol", drop = TRUE],
           pxl_row_in_fullres = x[, "imagerow", drop = TRUE],
           sampleID = i)
    coordinates_tibble <- bind_rows(coordinates_tibble, coordinates)
    if (verbose) cli_alert("  Collected coordinates")
  }
  
  # Create Staffli object
  if (verbose) cli_alert_warning(glue("Paths to raw images are unavailable. See '?ReplaceImagePaths()' ",
                                 "for more information on how to update the paths."))
  st_object <- CreateStaffliObject(imgs = image_paths,
                                   meta_data = coordinates_tibble, 
                                   image_info = img_info_tibble, 
                                   scalefactors = scalefactors_tibble)
  object@tools$Staffli <- st_object
  return(object)
}
