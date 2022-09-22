#' Load and merge multiple gene expression matrices
#'
#' Gene expression matrices should have features in rows and spots in columns.
#'
#' @details The merging process makes sure that all genes detected are present in the merged output.
#' This means that if a gene is missing in a certain dataset, the spots in that dataset will
#' be assigned with 0 expression.
#'
#' @family pre-processing
#' @describeIn load-data
#'
#' @param samplefiles Character vector of file/directory paths. Paths should specify .h5 or
#' .tsv/.tsv.gz files. Alternatively, the paths could specify directories including barcodes.tsv,
#' features.tsv and matrix.mtx files.
#' @param verbose Print messages
#'
#' @return A sparse matrix of class 'dgCMatrix'
#'
#' @examples
#' \donttest{
#' # Load and merge two gene expression matrices
#' samples <- c(system.file("extdata/mousebrain", "filtered_feature_bc_matrix.h5", package = "STUtility2"),
#'              system.file("extdata/mousecolon", "filtered_feature_bc_matrix.h5", package = "STUtility2"))
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
    if (verbose) inform(c("v" = glue::glue("  Finished loading expression matrix {i}")))
    colnames(exprMat) <- gsub(pattern = "-\\d+", replacement = paste0("-", i), x = colnames(exprMat))
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
  if (length(exprMats) > 1) {
    if (verbose) inform(c("i" = "Merging matrices:"))
    mergedMat <- SeuratObject::RowMergeSparseMatrices(mat1 = exprMats[[1]], mat2 = exprMats[2:length(exprMats)])
    if (verbose) inform(c("v" = glue::glue("  There are {nrow(mergedMat)} features and {ncol(mergedMat)} spots in the merged matrix.")))
    return(mergedMat)
  } else {
    inform(c("i" = "only 1 expression matrix loaded."))
    if (verbose) inform(c("v" = glue::glue("  There are {nrow(exprMats[[1]])} features and {ncol(exprMats[[1]])} spots in the matrix.")))
    return(exprMats[[1]])
  }

}
