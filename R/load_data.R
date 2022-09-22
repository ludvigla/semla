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
  if (!is.character(samplefiles)) abort("'samplefiles' must be a character vector.")
  checks <- tibble::tibble(samplefiles) |>
    dplyr::mutate(is = dplyr::case_when(file.exists(samplefiles) ~ "file", dir.exists(samplefiles) ~ "dir"))
  if (any(is.na(checks$is))) abort(c("Invalid path(s):", glue::glue("{checks$samplefiles[is.na(checks$is)]}")))

  # Load expression matrices
  if (verbose) inform(c("i" = "Loading expression matrices:"))
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
    if (verbose) inform(c("v" = sprintf("  Finished loading expression matrix %s", i)))
    return(exprMat)
  })

  # Check gene overlap
  all.genes <- lapply(exprMats, function(exprMat) {rownames(exprMat)})
  intersecting.genes <- Reduce(intersect, all.genes)
  if (length(intersecting.genes) == 0) abort("No shared genes shared across datasets.")
  if (length(intersecting.genes) < 1000) inform(c("x" = glue::glue("There are only {length(intersecting.genes)} gene shared across all matrices:"),
                                                  "x" = "  Are you sure that the matrices share the same gene IDs?",
                                                  "x" = "  Are the datasets from the same species?"))
  # Merge matrices
  if (verbose) inform(c("i" = "Merging matrices:"))
  mergedMat <- SeuratObject::RowMergeSparseMatrices(mat1 = exprMats[[1]], mat2 = exprMats[2:length(exprMats)])

  # Return merged matrix
  if (verbose) inform(c("v" = glue::glue("  There are {nrow(mergedMat)} features and {ncol(mergedMat)} spots in the merged matrix.")))
  return(mergedMat)
}

samples <- c("~/targeted_vs_untargeted/data_curated/spaceranger_output/colon/V10S29-108_B1/filtered_feature_bc_matrix.h5",
             "~/targeted_vs_untargeted/data_curated/spaceranger_output/mousebone/V11B18-362_C1/filtered_feature_bc_matrix.h5")
myMat <- LoadAndMergeMatrices(samples)
myMat <- LoadAndMergeMatrices(c("hubba", "bubba"))
