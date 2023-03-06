#' @include generics.R
#' @include checks.R
NULL

#' The Staffli Class
#'
#' The Staffli object is designed to hold information about the spatial data generated in a 10x
#' Visium SRT experiment. This includes paths to images, spot coordinates as well as loaded images
#' in `raster` format and additional information about these images.
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
#' @param meta_data Spot-level metadata to add to the `Staffli` object. This should be a `tbl` with
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
#' @export
#'
CreateStaffliObject <- function (
    imgs = NULL,
    meta_data,
    image_height = 400,
    image_info,
    scalefactors
) {
  if (!all(c("barcode", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID") %in% colnames(meta_data))) abort("Invalid meta_data columns.")

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

#' Method for extracting the raw H&E image information from a 'Staffli' object
#'
#' @param object A Staffli object
#' @export
#' @docType methods
#' @rdname ImageInfo
setGeneric("ImageInfo", function(object) {
  standardGeneric("ImageInfo")
})
#' @rdname ImageInfo
#' @aliases ImageInfo,Staffli,Staffli-method
#'
#' @export
#'
setMethod (
  f = "ImageInfo",
  signature = "Staffli",
  definition = function(object) {
    object@image_info
  }
)


#' Method to extract images from a Staffli object
#'
#' @param object A Staffli or Seurat object
#' @param type A string specifying the mage type to get
#'
#' @export
#' @docType methods
#' @rdname GetImages
#'
setGeneric("GetImages", function(object, type = "raw") {
  standardGeneric("GetImages")
})
#' @rdname GetImages
#' @aliases GetImages,Staffli,Staffli-method
#'
#' @export
#'
setMethod (
  f = "GetImages",
  signature = "Staffli",
  definition = function(object, type) {
    object@rasterlists[[type]]
  }
)
#' @rdname GetImages
#' @aliases GetImages,Seurat,Seurat-method
#'
#' @export
#'
setMethod (
  f = "GetImages",
  signature = "Seurat",
  definition = function(object, type) {
    .check_seurat_object(object)
    st.object <- object@tools$Staffli
    st.object@rasterlists[[type]]
  }
)


#' Method used to extract a Staffli object from the tools slot of a
#' Seurat object
#'
#' @param object A Seurat object
#'
#' @export
#' @docType methods
#' @rdname GetStaffli
#'
setGeneric("GetStaffli", function(object) {
  standardGeneric("GetStaffli")
})
#' @rdname GetStaffli
#' @aliases GetStaffli,Seurat,Seurat-method
setMethod (
  f = "GetStaffli",
  signature = "Seurat",
  definition = function(object) {
    .check_seurat_object(object)
    object@tools$Staffli
  }
)

#' Method used to get meta data from a Staffli object
#'
#' @rdname Staffli-get-methods
#' @aliases `[[`,Staffli,Staffli-method
#'
#' @param x object from which to extract element(s).
#' @param i row indices specifying elements to extract.
#' @param j column indices specifying elements to extract.
#' @param drop If TRUE the result is coerced to the lowest possible dimension.
#' This only works for extracting elements, not for the replacement.
#'
#' @export
#'
setMethod (
  f = "[[",
  signature = "Staffli",
  definition = function(x, i, j, drop = F) {
    x@meta_data[i, j, drop]
  }
)

#' Method used to set meta data in a Staffli object
#' @rdname Staffli-set-methods
#' @aliases `[[<-`,Staffli,Staffli-method
#'
#' @param x object in which to replace element(s).
#' @param i row indices specifying elements to replace.
#' @param j column indices specifying elements to replace.
#' @param value Data to add to meta data data.frame.
#' @param ... additional parameters
#'
#' @export
#'
setMethod (
  f = "[[<-",
  signature = "Staffli",
  definition = function(x, i, j, ..., value) {
    x@meta_data[i, j] <- value
    return(x)
  }
)

#' Show method for Staffli objects
#'
#' @rdname show
#' @aliases show,Staffli,Staffli-method
#'
#' @param object object to print preselected attributes for
#'
setMethod (
  f = "show",
  signature = "Staffli",
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
