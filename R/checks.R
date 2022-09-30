#' Check if a Seurat object has been processed with STUtility
#'
#' @param object A Seurat object
#' @param message An error message
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#'
.check_seurat_object <- function (
    object
) {
  if (!class(object) == "Seurat") abort(glue("invalid class '{class(object)}'"))
  if (!"Staffli" %in% names(object@tools)) abort(c("This Seurat object does not appear to have been processed with STUtility2.",
                                                 "x" = "'Staffli' object is missing from tools slot."))
}
