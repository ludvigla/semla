#' Create a merged expression matrix from multiple expression matrices stored on disk
#'
#' The merging makes sure that all genes detected are present in the merged output.
#' This means that if a gene is missing in a certain dataset, the spots in that
#' dataset will be assigned with 0 expression.
#'
#' @param samplefiles File paths to
#' @param verbose Print messages
#'
#' @return A sparse matrix of class 'dgCMatrix'
#'
#' @examples
#' \donttest{
#' samples <- list.files(...)
#' mergedMatrix <- LoadAndMergeMatrices(samples)
#' }
#'
#' @export
LoadAndMergeMatrices <- function (
    samplefiles,
    verbose = TRUE
) {

  # Run checks
  if (!is.character(samplefiles)) stop("'samplefiles' must be a character vector.")
  checks <- tibble::tibble(samplefiles) |>
    dplyr::mutate(is = dplyr::case_when(file.exists(samplefiles) ~ "file", dir.exists(samplefiles) ~ "dir"))
  if (any(is.na(checks$is))) stop(sprintf("Invalid path(s): \n\t%s", paste(checks$samplefiles[is.na(checks$is)], collapse = "\n\t")))

  # Load expression matrices
  if (verbose) cat("Loading expression matrices...\n")
  exprMats <- lapply(seq_along(samplefiles), function(i) {
    if (checks$is[i] == "dir") {
      # Assumes that the directory contains matrix, barcodes and genes
      exprMat <- Seurat::Read10X(data.dir = samplefiles[i])
    } else if (checks$is[i] == "file") {
      ext <- tools::file_ext(samplefiles[i])
      if (ext == "h5") {
        exprMat <- Seurat::Read10X_h5(samplefiles[i])
      } else if (ext %in% c("tsv", tsv.gz)) {
        exprMat <- data.frame(data.table::fread(samplefiles[i], sep = "\t", header = TRUE), row.names = 1)
        exprMat <- as(as.matrix(exprMat), "dgCMatrix")
      }
    }
    if (verbose) cat(sprintf("  Finished loading expression matrix %s.\n", i))
    return(exprMat)
  })

  # Merge matrices
  if (verbose) cat("\nMerging matrices...\n")
  mergedMat <- SeuratObject::RowMergeSparseMatrices(mat1 = exprMats[[1]], mat2 = exprMats[2:length(exprMats)])

  # Return merged matrix
  if (verbose) cat(sprintf("There are %s features and %s spots in the merged matrix.", nrow(mergedMat), ncol(mergedMat)))
  return(mergedMat)
}
