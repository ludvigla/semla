#' Check if a Seurat object has been processed with \code{semla}
#'
#' @param object A Seurat object
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#'
#' @return an error message if conditions are not met
#'
#' @noRd
.check_seurat_object <- function (
    object
) {
  if (!inherits(object, what = "Seurat")) abort(glue("invalid class '{class(object)}'"))
  if (!"Staffli" %in% names(object@tools)) abort(c("This Seurat object does not appear to have been processed with semla.",
                                                 "x" = "'Staffli' object is missing from tools slot."))
}


#' Check if a Seurat object contains loaded images
#'
#' @param object A Seurat object
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#'
#' @return an error message if conditions are not met
#'
#' @noRd
.check_seurat_images <- function (
    object
) {
  st_object <- GetStaffli(object)
  if (!"raw" %in% names(st_object@rasterlists)) abort("Images have not been loaded yet. Did you run 'LoadImages()'?")
}
