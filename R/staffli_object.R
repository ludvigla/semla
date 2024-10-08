#' @include generics.R
#' @include checks.R
NULL

#' The Staffli Class
#'
#' The \code{Staffli} object is designed to hold information about the spatial data generated in a 10x
#' Visium SRT experiment. This includes paths to images, spot coordinates as well as loaded images
#' in \code{raster} format and additional information about these images.
#'
#' @slot imgs A character vector of paths to the raw HE images
#' @slot rasterlists A list of lists containing images in 'raster' format
#' @slot meta_data A tibble with spot coordinates and additional meta data
#' @slot image_height The height of the scaled images in pixels
#' @slot image_info A tibble with information about the raw images
#' @slot scalefactors A tibble with information about scalefactors used to convert coordinates
#' between the original image and down-scaled representations such as "tissue_lowres_image.jpg"
#' @slot pixels_per_um Numeric vector specifying the number of pixels in the raw images that
#' corresponds to 1 micrometer
#' @slot version Package version.
#'
#' @name Staffli-class
#' @rdname Staffli-class
#'
#' @author Ludvig Larsson
#'
#' @exportClass Staffli
#'
Staffli <- setClass (
  Class = 'Staffli',
  slots = c(
    imgs = 'ANY',
    rasterlists = 'list',
    meta_data = 'tbl',
    image_height = 'ANY',
    image_info = 'ANY',
    scalefactors = 'ANY',
    pixels_per_um = 'numeric',
    version = 'package_version'
  )
)


#' Create a Staffli object
#'
#' Create a Staffli object from a set of images and associated spot coordinates
#'
#' @param imgs Character vector specifying paths to images in JPG, PNG or TIF format
#' @param meta_data Spot-level metadata to add to the \code{Staffli} object. This should be a \code{tbl} with
#' required columns 'barcode' representing the spot IDs, 'pxl_col_in_fullres' and 'pxl_row_in_fullres'
#' which specifies the 10x Visium array coordinates and a 'sampleID' column with sample IDs
#' @param image_height Specifies the height of the scaled images in pixels [default: 400 pixels]
#' @param image_info a tibble with image information
#' @param scalefactors a tibble with scalefactors sued to transform coordinates from
#' the original image space to the downscaled images
#'
#' @importFrom utils packageVersion
#' @importFrom Matrix colSums
#' @importFrom stats setNames
#' @importFrom rlang abort
#' @importFrom methods new
#' 
#' @author Ludvig Larsson
#' 
#' @returns A \code{Staffli} object
#' 
#' @examples 
#' 
#' library(semla)
#' 
#' # Multiple samples
#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' # Create an object with multiple samples
#' he_imgs <- c(system.file("extdata/mousebrain", 
#'                          "spatial/tissue_lowres_image.jpg", 
#'                          package = "semla"),
#'              system.file("extdata/mousecolon", 
#'                          "spatial/tissue_lowres_image.jpg", 
#'                          package = "semla"))
#' spotfiles <- c(system.file("extdata/mousebrain", 
#'                            "spatial/tissue_positions_list.csv", 
#'                            package = "semla"),
#'                system.file("extdata/mousecolon", 
#'                            "spatial/tissue_positions_list.csv", 
#'                            package = "semla"))
#' jsonfiles <- c(system.file("extdata/mousebrain", 
#'                            "spatial/scalefactors_json.json", 
#'                            package = "semla"),
#'                system.file("extdata/mousecolon", 
#'                            "spatial/scalefactors_json.json", 
#'                            package = "semla"))
#' 
#' # Read coordinates and select relevant columns
#' coordinates <- LoadSpatialCoordinates(spotfiles) |> 
#'   select(all_of(c("barcode", "x", "y", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID")))
#' 
#' # Create image_info
#' image_info <- LoadImageInfo(he_imgs)
#' 
#' # Read scalefactors
#' scalefactors <- LoadScaleFactors(jsonfiles)
#' 
#' # Add additional columns to image_info using scalefactors
#' image_info <- UpdateImageInfo(image_info, scalefactors)
#' 
#' # Create Staffli object
#' staffli_object <- CreateStaffliObject(imgs = he_imgs, 
#'                                       meta_data = coordinates, 
#'                                       image_info = image_info, 
#'                                       scalefactors = scalefactors)
#' staffli_object
#'
#' @export
#'
CreateStaffliObject <- function (
    imgs = NULL,
    meta_data,
    image_height = 400,
    image_info,
    scalefactors
) {
  
  # Check meta_data
  if (!all(c("barcode", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID") %in% colnames(meta_data))) abort("Invalid meta_data.")
  if (all(c("x", "y") %in% colnames(meta_data))) {
    meta_data <- meta_data |> 
      select(all_of(c("barcode", "x", "y", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID")))
  } else {
    meta_data <- meta_data |> 
      select(all_of(c("barcode", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID")))
  }
  checks <- sapply(meta_data, class)
  if (!checks["barcode"] == "character")
    abort(glue("Invalid class '{checks['barcode']}' for barcode. Expected a 'character' vector."))
  if (!checks["pxl_col_in_fullres"] %in% c("numeric", "integer"))
    abort(glue("Invalid class '{checks['pxl_col_in_fullres']}' for pxl_col_in_fullres. Expected a 'numeric' or 'integer' vector."))
  if (!checks["pxl_col_in_fullres"] %in% c("numeric", "integer"))
    abort(glue("Invalid class '{checks['pxl_row_in_fullres']}' for pxl_row_in_fullres. Expected a 'numeric' or 'integer' vector."))
  if (!checks["sampleID"] == "integer")
    abort(glue("Invalid class '{checks['sampleID']}' for sampleID. Expected an 'integer' vector."))
  sampleIDs_meta_data <- meta_data$sampleID |> paste0() |> unique() |> sort()
  
  # Check image_info
  if (!all(c("width", "height", "full_width", "full_height") %in% colnames(image_info))) abort("Invalid image_info.")
  sampleIDs_image_info <- image_info$sampleID |> sort()
  
  # Check image_info
  sampleIDs_scalefactors <- image_info$sampleID |> sort()
  
  # Check sampleIDs
  if(!all(sampleIDs_meta_data == sampleIDs_image_info)) abort("Invalid sampleIDs")
  if(!all(sampleIDs_meta_data == sampleIDs_scalefactors)) abort("Invalid sampleIDs")
  if(!all(sampleIDs_image_info == sampleIDs_scalefactors)) abort("Invalid sampleIDs")
  
  # Check H&E images
  if (!is.null(imgs)) {
    if (!inherits(imgs, what = "character"))
      abort(glue("Invalid class '{class(imgs)}'. Expected a 'character' vector."))
    if (length(imgs) != length(sampleIDs_meta_data)) 
      abort(glue("Number of images ({length(imgs)}) does not match the number of sampleIDs ({length(sampleIDs_meta_data)})"))
  }

  object <- new (
    Class = 'Staffli',
    imgs = imgs,
    meta_data = meta_data,
    image_height = image_height,
    image_info = image_info,
    scalefactors = scalefactors,
    version = packageVersion(pkg = 'semla')
  )

  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Staffli methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Staffli Methods
#'
#' Methods for \code{\link{Staffli}} objects for generics defined in other
#' packages
#'
#' @param x A \code{\link{Staffli}} object
#' @param i Name or index of one or more metadata columns
#' @param ... Arguments passed to other methods
#'
#' @name Staffli-methods
#' @rdname Staffli-methods
#'
#' @concept staffli
#'
NULL


#' @describeIn Staffli-methods Auto completion for \code{$} access on a
#' \code{Staffli} object
#'
#' @inheritParams utils::.DollarNames
#'
#' @importFrom utils .DollarNames
#' @export
#' @method .DollarNames Staffli
#'
".DollarNames.Staffli" <- function(x, pattern = '') {
  meta_data <- as.list(x = colnames(x = x@meta_data))
  names(x = meta_data) <- unlist(x = meta_data)
  return(.DollarNames(x = meta_data, pattern = pattern))
}

#' @describeIn Staffli-methods Metadata access for \code{Staffli} objects
#'
#' @return A selected metadata column \code{i} for object \code{x}
#'
#' @export
#' @method $ Staffli
#' 
"$.Staffli" <- function(x, i, ...) {
  return(x@meta_data[, i, drop = TRUE])
}

#' Method to extract image info
#'
#' @param object A \code{Staffli} or \code{Seurat} object
#' 
#' @return A \code{tbl} with image information
#'
#' @export
#' @docType methods
#' @rdname GetImageInfo
#'
setGeneric("GetImageInfo", function(object) {
  standardGeneric("GetImageInfo")
})
#' @rdname GetImageInfo
#' 
#' @examples 
#' 
#' # Load example data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla")) |> 
#'    LoadImages()
#' 
#' # Fetch Staffli object
#' staffli_object <- GetStaffli(se_mbrain)
#' 
#' # Fetch images from a Staffli object
#' image_info <- GetImageInfo(staffli_object)
#' image_info
#'
#' @export
#'
setMethod (
  f = "GetImageInfo",
  signature = "Staffli",
  definition = function(object) {
    object@image_info
  }
)
#' @rdname GetImageInfo
#' 
#' @examples 
#'
#' # Fetch images from a Seurat object
#' image_info <- GetImageInfo(se_mbrain)
#' image_info
#'
#' @export
#'
setMethod (
  f = "GetImageInfo",
  signature = "Seurat",
  definition = function(object) {
    .check_seurat_object(object)
    st_object <- object@tools$Staffli
    st_object@image_info
  }
)


#' Method to extract image scale factors
#'
#' @param object A \code{Staffli} or \code{Seurat} object
#' 
#' @return A \code{tbl} with image scale factors
#'
#' @export
#' @docType methods
#' @rdname GetScaleFactors
#'
setGeneric("GetScaleFactors", function(object) {
  standardGeneric("GetScaleFactors")
})
#' @rdname GetScaleFactors
#' 
#' @examples 
#' 
#' # Load example data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla")) |> 
#'    LoadImages()
#' 
#' # Fetch Staffli object
#' staffli_object <- GetStaffli(se_mbrain)
#' 
#' # Fetch images from a Staffli object
#' scalefactors <- GetScaleFactors(staffli_object)
#' scalefactors
#'
#' @export
#'
setMethod (
  f = "GetScaleFactors",
  signature = "Staffli",
  definition = function(object) {
    object@scalefactors
  }
)
#' @rdname GetScaleFactors
#' 
#' @examples 
#' 
#' # Fetch images from a Seurat object
#' scalefactors <- GetScaleFactors(se_mbrain)
#' scalefactors
#'
#' @export
#'
setMethod (
  f = "GetScaleFactors",
  signature = "Seurat",
  definition = function(object) {
    .check_seurat_object(object)
    st_object <- object@tools$Staffli
    st_object@scalefactors
  }
)


#' Method to extract coordinates
#'
#' @param object A \code{Staffli} or \code{Seurat} object
#' 
#' @return A \code{tbl} with spot coordinates
#'
#' @export
#' @docType methods
#' @rdname GetCoordinates
#'
setGeneric("GetCoordinates", function(object) {
  standardGeneric("GetCoordinates")
})
#' @rdname GetCoordinates
#' 
#' @examples 
#' 
#' # Load example data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla")) |> 
#'    LoadImages()
#' 
#' # Fetch Staffli object
#' staffli_object <- GetStaffli(se_mbrain)
#' 
#' # Fetch coordinates from a Staffli object
#' coordinates <- GetCoordinates(staffli_object)
#' coordinates
#'
#' @export
#'
setMethod (
  f = "GetCoordinates",
  signature = "Staffli",
  definition = function(object) {
    object@meta_data
  }
)
#' @rdname GetCoordinates
#' 
#' @examples 
#' 
#' # Fetch images from a Seurat object
#' coordinates <- GetCoordinates(se_mbrain)
#' coordinates
#'
#' @export
#'
setMethod (
  f = "GetCoordinates",
  signature = "Seurat",
  definition = function(object) {
    .check_seurat_object(object)
    st_object <- object@tools$Staffli
    st_object@meta_data
  }
)


#' Method to extract images
#'
#' @param object A \code{Staffli} or \code{Seurat} object
#' @param image_use A string specifying the image type to get
#' 
#' @return A \code{list} of images in \code{raster} format
#'
#' @export
#' @docType methods
#' @rdname GetImages
#'
setGeneric("GetImages", function(object, image_use = c("raw", "transformed")) {
  standardGeneric("GetImages")
})
#' @rdname GetImages
#' 
#' @examples 
#' 
#' # Load example data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla")) |> 
#'    LoadImages()
#' 
#' # Fetch Staffli object
#' staffli_object <- GetStaffli(se_mbrain)
#' 
#' # Fetch images from a Staffli object
#' images <- GetImages(staffli_object)
#'
#' @export
#'
setMethod (
  f = "GetImages",
  signature = "Staffli",
  definition = function(object, image_use = c("raw", "transformed")) {
    image_use <- match.arg(image_use, choices = c("raw", "transformed"))
    if (!inherits(image_use, what = "character"))
      abort(glue("Invalid class '{class(image_use)}' of 'image_use'. Expected a 'character' vector."))
    if (!image_use %in% names(object@rasterlists)) 
      abort(glue("'{image_use}' images are not available."))
    object@rasterlists[[image_use]]
  }
)
#' @rdname GetImages
#' 
#' @examples 
#'
#' # Fetch images from a Seurat object
#' images <- GetImages(se_mbrain)
#'
#' @export
#'
setMethod (
  f = "GetImages",
  signature = "Seurat",
  definition = function(object, image_use) {
    image_use <- match.arg(image_use, choices = c("raw", "transformed"))
    .check_seurat_object(object)
    st_object <- object@tools$Staffli
    if (!inherits(image_use, what = "character"))
      abort(glue("Invalid class '{class(image_use)}' of 'image_use'. Expected a 'character' vector."))
    if (!image_use %in% names(st_object@rasterlists)) 
      abort(glue("'{image_use}' images are not available."))
    st_object@rasterlists[[image_use]]
  }
)


#' Method used to extract a \code{Staffli} object from the tools slot of a
#' \code{Seurat} object
#'
#' @param object A \code{Seurat} object
#' 
#' @return A \code{Staffli} object
#'
#' @export
#' @docType methods
#' @rdname GetStaffli
#'
setGeneric("GetStaffli", function(object) {
  standardGeneric("GetStaffli")
})
#' @rdname GetStaffli
#' 
#' @examples
#' 
#' # Load example data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla")) |> 
#'    LoadImages()
#' 
#' # Fetch Staffli object from a Seurat object
#' staffli_object <- GetStaffli(se_mbrain)
#' 
#' @export
setMethod (
  f = "GetStaffli",
  signature = "Seurat",
  definition = function(object) {
    .check_seurat_object(object)
    object@tools$Staffli
  }
)

#' Replace image paths
#' 
#' @param object A \code{Staffli} or \code{Seurat} object
#' @param paths A character vector with image paths
#' 
#' @import rlang
#' @import glue
#' @import dplyr
#' @importFrom magick image_read image_info
#' @importFrom tibble tibble
#' 
#' @return A \code{Staffli} or \code{Seurat} object with updated image paths
#' 
#' @export
#' @docType methods
#' @rdname ReplaceImagePaths
#' 
setGeneric("ReplaceImagePaths", function(object, paths) {
  standardGeneric("ReplaceImagePaths")
})
#' @rdname ReplaceImagePaths
#' 
#' @export
#' 
setMethod (
  f = "ReplaceImagePaths",
  signature = "Staffli",
  definition = function(object, paths) {
    
    # Set global variables to NULL
    width <- height <- full_width <- full_height <- colorspace <- filesize <- density <- sampleID <- type <- NULL
    
    if (!inherits(paths, what = "character")) 
      abort(glue("Invalid class '{class(paths)}'. Expected a 'character' vector."))
    if (!length(paths) == nrow(object@image_info)) 
      abort(glue("Expected {nrow(object@image_info)} images, got {length(paths)}"))
    image_info <- object@image_info
    image_info_new <- tibble()
    for (i in seq_along(paths)) {
      f <- paths[i]
      
      # Check file extension
      if (file.exists(f)) {
        if (!file_ext(f) %in% c("png", "jpg", "jpeg"))
          abort(glue("Only PNG and JPEG images are supported, got file extension .{file_ext(f)}"))
      }
      
      im <- try({f |> image_read()}, silent = TRUE)
      if (inherits(im, "try-error"))
        abort(glue("Invalid path/url {f}. Make sure to have ",
                   "valid image paths before running {col_br_magenta('LoadImages')}"))
      info <- im |>
        image_info()
      
      image_info_sample <- info |> 
        mutate(sampleID = paste0(i),
               type = case_when(basename(paths[i]) %in% paste0("tissue_hires_image.", c("jpg", "png")) ~ "tissue_hires",
                                basename(paths[i]) %in% paste0("tissue_lowres_image.", c("jpg", "png")) ~ "tissue_lowres",
                                TRUE ~ "unknown"))
      image_info_new <- bind_rows(image_info_new, image_info_sample)
    }
    image_info_new$full_width <- image_info$full_width
    image_info_new$full_height <- image_info$full_height
    object@imgs <- paths
    object@image_info <- image_info_new |> 
      select(format, width, height, full_width, full_height, colorspace, filesize, density, sampleID, type)
    return(object)
  }
)
#' @rdname ReplaceImagePaths
#' 
#' @export
#' 
setMethod (
  f = "ReplaceImagePaths",
  signature = "Seurat",
  definition = function(object, paths) {
    .check_seurat_object(object)
    st_object <- object@tools$Staffli
    st_object <- ReplaceImagePaths(st_object, paths)
    object@tools$Staffli <- st_object
    return(object)
  }
)


#' Show method for \code{Staffli} objects
#' 
#' @importFrom methods show
#'
#' @rdname show
#' 
#' @return No return value, plots the spot coordinates and alternatively also
#' the images found in a \code{Staffli} object
#'
#' @param object object to print pre-selected attributes for
#' 
#' @export
#'
setMethod (
  f = "show",
  signature = c(object = "Staffli"),
  definition = function(object) {
    cat("An object of class", class(x = object), "\n")
    cat(
      nrow(object@meta_data),
      'spots across',
      length(unique(object@meta_data[, "sampleID", drop = TRUE])),
      'samples. \n'
    )
    if (length(object@rasterlists) > 0) {
      cat(
        '\nAvailable image representations: \n  ',
        paste(names(object@rasterlists), collapse = ", "),
        '\n'
      )
    }
  }
)



#' Plot method for \code{Staffli} objects
#' 
#' @import rlang
#' @import glue
#' @import dplyr
#' @importFrom patchwork wrap_plots inset_element 
#' 
#' @return No return value, plots the content of a \code{Staffli} object
#'
#' @rdname plot
#'
#' @param x A \code{Staffli} object
#' @param image_use A string specifying the image to plot
#' @param coords_use A character vector of length 2 specifying the coordinates to use
#' @param ncol Integer specifying the number of columns in the plot grid
#' @param ... Additional parameters passed to \code{geom_point}
#' 
#' @export
#'
setMethod (
  f = "plot",
  signature = c(x = "Staffli", y = "missing"),
  definition = function(x, image_use = NULL, coords_use = c("raw", "transformed"), ncol = NULL, ...) {
    
    # Set global variables to NULL
    sampleID <- NULL
    
    # Check image_use
    if (!is.null(image_use)) {
      image_use <- match.arg(image_use, choices = c("raw", "transformed"))
    }
    coords_use <- match.arg(coords_use, choices = c("raw", "transformed"))
    
    if (!is.null(image_use)) {
      if (!inherits(image_use, what = "character"))
        abort(glue("Invalid class '{class(image_use)}' for 'image_use'. Expected a 'character'."))
      if (!image_use %in% names(x@rasterlists))
        abort(glue("'{image_use}' images are not available."))
      if (image_use == "raw") {
        coords_columns <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
      }
      if (image_use == "transformed") {
        coords_columns <- c("pxl_col_in_fullres_transformed", "pxl_row_in_fullres_transformed")
      }
      rstrs <- x@rasterlists[[image_use]]
    } else {
      if (coords_use == "raw") {
        coords_columns <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
      }
      if (coords_use == "transformed") {
        coords_columns <- c("pxl_col_in_fullres_transformed", "pxl_row_in_fullres_transformed")
      }
    }
    # Split coordinates
    coords_split <- x@meta_data |> 
      group_by(sampleID) |> 
      group_split()
    
    # Get dims
    dims <- x@image_info |> 
      mutate(sampleID = 1:n()) |> 
      group_by(sampleID) |> 
      group_split()
    
    # Create plot
    plots <- lapply(seq_along(coords_split), function(i) {
      xy <- coords_split[[i]]
      ggplot(xy, aes(x = !! sym(coords_columns[1]), y = !! sym(coords_columns[2]))) +
        geom_point(shape = 21, fill = NA, ...) +
        scale_x_continuous(limits = c(0, dims[[i]]$full_width), expand = c(0, 0)) +
        scale_y_reverse(limits = c(dims[[i]]$full_height, 0), expand = c(0, 0)) +
        theme_void() +
        coord_fixed() +
        ggtitle(label = paste0("Section ", i))
    })
    
    # Add images
    if (!is.null(image_use)) {
      plots <- lapply(seq_along(plots), function(i) {
        plots[[i]] +
          inset_element(p = rstrs[[i]], left = 0, bottom = 0, right = 1, top = 1, on_top = FALSE)
      })
    }
    
    # Wrap plots
    ncol <- ncol %||% ceiling(sqrt(length(plots)))
    wrap_plots(plots, ncol = ncol)
  }
)
