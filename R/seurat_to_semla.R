#' @include checks.R
#'
NULL

#' Update a Seurat object for compatibility with semla
#' 
#' In order to make a \code{Seurat} object with \code{semla}, 
#' you can use \code{AddSemlaToSeurat} to add the required 
#' data to the \code{Seurat} object. This assumes that 
#' the \code{Seurat} object contains one or more "VisiumV1" 
#' object(s) in the \code{images} slot. Alternatively, you
#' can convert a \code{Seurat} object with "SlideSeq" data.
#' 
#' @section VisiumV1:
#' Note that you need to specify what H&E image was loaded, 
#' one of "tissue_lowres" or "tissue_hires". If this argument 
#' is incorrect, the tissue coordinates will be misplaced.
#' 
#' Visit the \href{https://ludvigla.github.io/semla/articles/getting_started.html}{getting started}
#' tutorial on our package website for an example on how to convert a \code{Seurat} object with
#' VisiumV1 data.
#' 
#' @section VisiumV2:
#' Note that you need to specify what H&E image was loaded, 
#' one of "tissue_lowres" or "tissue_hires". If this argument 
#' is incorrect, the tissue coordinates will be misplaced.
#' 
#' The \href{https://ludvigla.github.io/semla/articles/getting_started.html}{getting started} 
#' tutorial applies also for VisiumV2 assays. If you are working with VisiumHD data however, 
#' we recommend loading the data using \code{semla}'s own functions, as detailed in the
#' \href{https://ludvigla.github.io/semla/articles/visiumHD.html}{VisiumHD} tutorial.
#'  
#' @section SlideSeq:
#' For SlideSeq data, there's no additional H&E image provided. If you convert
#' a Seurat object containing SlideSeq data, all image related functionality 
#' of \code{semla} will be inaccessible.
#' 
#' @param object An object of class \code{Seurat} with Visium data
#' @param image_type One of "tissue_lowres" or "tissue_hires", specifying
#' what H&E image was loaded into the \code{Seurat} object. Only used for 
#' "VisiumV1" data.
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
#' @importFrom grDevices dev.off png
#' 
#' @examples 
#' \dontrun{
#' library(semla)
#' library(SeuratData)
#' 
#' # Load example Seurat object with VisiumV1/VisiumV2 data (depends on your Seurat version)
#' InstallData("stxBrain")
#' brain <- LoadData("stxBrain", type = "anterior1")
#' 
#' # Make Seurat object compatible with semla
#' brain_semla <- UpdateSeuratForSemla(brain)
#' 
#' # Load example Seurat object with SlideSeq data
#' InstallData("ssHippo")
#' slide_seq <- LoadData("ssHippo")
#' 
#' # Make Seurat object compatible with semla
#' slide_seq <- UpdateSeuratForSemla(slide_seq)
#' }
#' 
#' @export
UpdateSeuratForSemla <- function (
    object,
    image_type = c("tissue_lowres", "tissue_hires"),
    verbose = TRUE
) {
  
  # Set global variables to NULL
  width <- height <- full_width <- full_height <- colorspace <- filesize <- density <- NULL
  type <- y <- sampleID <- NULL
  
  # Check object
  if (!inherits(object, what = "Seurat")) abort(glue("Invalid class '{class(object)}'. Expected a 'Seurat' object."))
  if (length(object@images) == 0) abort(glue("No images available in 'Seurat' object."))
  
  # Check images slot
  slice_types <- sapply(object@images, function(image) {
    image |> class() |> as.character()
  })
  if (length(unique(slice_types)) > 1) abort(glue("Only one spatial data type allowed. Got {paste(slice_types, collapse = ',')}"))
  if (!all(slice_types %in% c("SlideSeq", "VisiumV1", "VisiumV2"))) abort(glue("Only 'SlideSeq', 'VisiumV1' and 'VisiumV2' are currently supported. Got {slice_types}"))
  slice_type <- unique(slice_types)
  cli_alert_info("Found {slice_type} object(s).")
  
  if (slice_type %in%  c("VisiumV1", "VisiumV2")) {
    image_type <- match.arg(arg = image_type, choices = c("tissue_lowres", "tissue_hires"))
    if (verbose) cli_alert_info("Expecting '{image_type}' format.")
    
    # Get images
    imgs <- lapply(object@images, function(slice) {
      slice@image |> as.raster()
    }) 
  } else {
    imgs <- NULL
  }
  
  # Export images
  temp_dir <- tempdir()
  image_paths <- c()
  img_info_tibble <- tibble()
  scalefactors_tibble <- tibble()
  coordinates_tibble <- tibble()
  cli_h3("Collecting data from @images slot")
  for (i in seq_along(object@images)) {
    
    # Check if the slice represent Visium data
    if (!inherits(object@images[[i]], what = c("VisiumV1", "VisiumV2", "SlideSeq"))) abort("Only 'VisiumV1', 'VisiumV2' and 'SlideSeq' data is supported.")
    
    # Collect image paths
    if (verbose) cli_alert_info("Collecting data for sample {i}:")
    if (slice_type %in%  c("VisiumV1", "VisiumV2")) {
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
        mutate(full_width = ceiling(dims[2]/sf_image), 
               full_height = ceiling(dims[1]/sf_image),
               sampleID = paste0(i), type = image_type) |> 
        select(format, width, height, full_width, full_height, colorspace, filesize, density, sampleID, type)
      if (verbose) cli_alert("  Collected image information.")
      img_info_tibble <- bind_rows(img_info_tibble, iminfo)
    }
    
    # Collect spatial coordinates.
    ## VisiumV1 and Slide-seq assays access the coordinates the same way. This changed for VisiumV2 assays
    slice <- object@images[[i]]
    if (slice_type %in% c("VisiumV1", "SlideSeq")) {
      x <- slice@coordinates
    } else if (slice_type == "VisiumV2") {
      x <- slice@boundaries$centroids@coords
      rownames(x) <- slice@boundaries$centroids@cells
    }
    
    if (slice_type == "VisiumV1") {
      coordinates <- tibble(barcode = rownames(x), 
                            pxl_col_in_fullres = x[, "imagecol", drop = TRUE],
                            pxl_row_in_fullres = x[, "imagerow", drop = TRUE],
                            x = x[, "col", drop = TRUE],
                            y = y[, "row", drop = TRUE],
                            sampleID = i)
      coordinates_tibble <- bind_rows(coordinates_tibble, coordinates)
    } else  if (slice_type == "VisiumV2") {
      coordinates <- tibble(barcode = rownames(x), 
                            pxl_col_in_fullres = x[, "y", drop = TRUE],
                            pxl_row_in_fullres = x[, "x", drop = TRUE],
                            sampleID = as.integer(i))
      coordinates_tibble <- bind_rows(coordinates_tibble, coordinates)
      # Create array coordinates
      ## extract the x and y coordiantes of each spot
      if (!requireNamespace("stringr", quietly = TRUE)) {
        abort(glue("Package {cli::col_br_magenta('stringr')} is required. Please install it with: \n",
                   "install.packages('stringr')"))
      }
      array_coords <- stringr::str_split_fixed(coordinates_tibble$barcode, pattern = "[[:punct:]]", n = 5)
      ## add to coordinates tibble
      coordinates_tibble <- coordinates_tibble |> 
        mutate(x = as.integer(array_coords[,4]),
               y = as.integer(array_coords[,3]),
               pxl_col_in_fullres = as.integer(.data$pxl_col_in_fullres),
               pxl_row_in_fullres = as.integer(.data$pxl_row_in_fullres))
    } else if (slice_type == "SlideSeq") {
      coordinates <- tibble(barcode = rownames(x), 
                            pxl_col_in_fullres = x[, "x", drop = TRUE],
                            pxl_row_in_fullres = x[, "y", drop = TRUE],
                            sampleID = i)
      coordinates_tibble <- bind_rows(coordinates_tibble, coordinates)
      wh <- sapply(x |> select(x, y), range)
      iminfo <- tibble(format = NA_character_,
                       width = NA_character_,
                       height = NA_character_,
                       full_width = wh[2, 1] + wh[1, 1],
                       full_height = wh[2, 2] + wh[1, 2],
                       colorspace = NA_character_,
                       filesize = NA_integer_,
                       density = NA_character_,
                       sampleID = paste0(i),
                       type = NA_character_)
      img_info_tibble <- bind_rows(img_info_tibble, iminfo)
      sf <- tibble(spot_diameter_fullres = NA,
                   tissue_hires_scalef = NA,
                   fiducial_diameter_fullres = NA,
                   tissue_lowres_scalef = NA,
                   sampleID = paste0(i))
      scalefactors_tibble <- bind_rows(scalefactors_tibble, sf)
    }
    if (verbose) cli_alert("  Collected coordinates.")
  }
  
  # Create Staffli object
  if (slice_type %in%  c("VisiumV1", "VisiumV2")) {
    if (verbose) cli_alert_warning(glue("Paths to raw images are unavailable. See '?ReplaceImagePaths()' ",
                                        "for more information on how to update the paths."))
  }
  
  st_object <- CreateStaffliObject(imgs = image_paths,
                                   meta_data = coordinates_tibble, 
                                   image_info = img_info_tibble, 
                                   scalefactors = scalefactors_tibble)
  object@tools$Staffli <- st_object
  if (verbose) cat_line()
  if (verbose) cli_alert_success("Returning updated {col_br_magenta('Seurat')} object.")
  
  return(object)
}