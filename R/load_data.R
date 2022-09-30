#' Load and merge multiple gene expression matrices
#'
#' Gene expression matrices should have features in rows and spots in columns.
#'
#' @details The merging process makes sure that all genes detected are present in the merged output.
#' This means that if a gene is missing in a certain dataset, the spots in that dataset will
#' be assigned with 0 expression.
#'
#' Spot IDs are renamed to be unique. Usually, the spots are named something similar to:
#' \code{"ACGCCTGACACGCGCT-1", "TACCGATCCAACACTT-1"}
#'
#' Since spot barcodes are shared across datasets, there is a risk that some of the spot IDs
#' will be duplicated after merging. To avoid this, the prefix (e.g. "-1") is replaced by
#' a unique prefix for each loaded matrix: "-1", "-2", "-3", ...
#'
#' @family pre-process
#'
#' @param samplefiles Character vector of file/directory paths. Paths should specify `.h5` or
#' `.tsv`/`.tsv.gz` files. Alternatively, the paths could specify directories including `barcodes.tsv`,
#' `features.tsv` and `matrix.mtx` files.
#' @param verbose Print messages
#'
#' @return A sparse matrix of class `dgCMatrix`
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


#' Load and merge multiple coordinate tables
#'
#' Load coordinates from \strong{'tissue_positions_list.csv'} files and merge them into a `tibble`.
#'
#' @details The merging process makes sure that all genes detected are present in the merged output.
#' This means that if a gene is missing in a certain dataset, the spots in that dataset will
#' be assigned with 0 expression.
#'
#' @family pre-process
#'
#' @param coordinatefiles Character vector of file paths. Paths should specify `.csv` files output by spaceranger
#' @param remove_spots_outside_tissue Should spots outside the tissue be removed?
#' @param verbose Print messages
#'
#' @return An object of class `tbl` (`tibble`)
#'
#' @examples
#' \donttest{
#' # Load and merge coordinates from two samples
#' coordinatefiles <- c(system.file("extdata/mousebrain/spatial", "tissue_positions_list.csv", package = "STUtility2"),
#'              system.file("extdata/mousecolon/spatial", "tissue_positions_list.csv", package = "STUtility2"))
#' coordinates <- LoadSpatialCoordinates(coordinatefiles)
#' }
#'
#' @export
LoadSpatialCoordinates <- function (
    coordinatefiles,
    remove_spots_outside_tissue = TRUE,
    verbose = TRUE
) {

  # Run checks
  if (!is.character(coordinatefiles)) abort("'coordinatefiles' must be a character vector.")
  checks <- tibble::tibble(coordinatefiles) |>
    dplyr::mutate(is = dplyr::case_when(file.exists(coordinatefiles) ~ "file"))
  if (any(is.na(checks$is))) abort(c("Invalid path(s):", glue::glue("{checks$coordinatefiles[is.na(checks$is)]}")))

  # Load coordinates
  if (verbose) inform(c("i" = "Loading coordinates:"))
  coordDF <- inject(rbind(!!!lapply(seq_along(coordinatefiles), function(i) {
    coords <- read.csv(file = coordinatefiles[i], header = FALSE) |>
      tibble::as_tibble() |>
      setNames(nm = c("barcode", "selected", "y", "x", "pxl_col_in_fullres", "pxl_row_in_fullres"))
    if (remove_spots_outside_tissue) {
      coords <- coords |>
        filter(selected == 1)
    }
    if (verbose) inform(c("v" = glue::glue("  Finished loading coordinates for sample {i}")))
    coords$barcode <- gsub(pattern = "-\\d+", replacement = paste0("-", i), x = coords$barcode)
    coords$sampleID <- i
    return(coords)
  })))

  # Return coordinates
  if (verbose) inform(c("v" = glue::glue("  Collected coordinates for {nrow(coordDF)} spots.")))
  return(coordDF)
}

#' Read image data
#'
#' Load image related data from \strong{'tissue_hires_image.png'} files.
#'
#' @family pre-process
#'
#' @param images An object of class `tibble` containing paths to images in PNG format, with
#' one row per sample. Paths should specify `.png` files output by spaceranger such as
#' `tissue_lowres_image.png` or `tissue_hires_image.png`. You do not have to load both
#' H&E images.
#' @param jsonfiles A character vector with file paths. Paths should specify `.json` files containing
#' scalefactors output by spaceranger.
#' @param verbose Print messages
#'
#' @return An object of class `DFrame``
#'
#' @examples
#' \donttest{
#' # Collect spaeranger output files
#' lowres.imagefiles <- c(system.file("extdata/mousebrain/spatial", "tissue_lowres_image.png", package = "STUtility2"),
#'              system.file("extdata/mousecolon/spatial", "tissue_lowres_image.png", package = "STUtility2"))
#' hires.imagefiles <- c(system.file("extdata/mousebrain/spatial", "tissue_hires_image.png", package = "STUtility2"),
#'              system.file("extdata/mousecolon/spatial", "tissue_hires_image.png", package = "STUtility2"))
#' jsonfiles <- c(system.file("extdata/mousebrain/spatial", "scalefactors_json.json", package = "STUtility2"),
#'              system.file("extdata/mousecolon/spatial", "scalefactors_json.json", package = "STUtility2"))
#'
#' # Load lowres/hires image data for two samples
#' imgFiles <- tibble::tibble(lowres.imagefiles, hires.imagefiles)
#' imgData <- LoadImageData(imgFiles, jsonfiles)
#'
#' # Load lower image data only
#' imgFiles <- tibble::tibble(lowres.imagefiles)
#' imgData <- LoadImageData(imgFiles, jsonfiles)
#' }
#' @export
LoadImageData <- function (
    images,
    jsonfiles,
    verbose = TRUE
) {

  # Check files
  check <- sapply(images, file.exists)
  if (any(!check)) abort(message = c("Invalid path(s):", glue::glue("File path '{unlist(images[!check])}' does not exist")))

  # read image data
  if (verbose) inform(c("i" = glue::glue("Reading image data for {nrow(images)} samples.")))
  imgData <- inject(rbind(!!!lapply(1:nrow(images), function(i) {
    png.files <- unlist(images[i, ])
    exts <- tools::file_ext(png.files)
    check <- exts %in% "png"
    if (!all(check)) abort(glue::glue("Invalid image format: '{exts[!check]}'"))
    DF <- readImgData(
      imageSources = unlist(images[i, ]),
      scaleFactors = jsonfiles[i],
      sample_id = paste0(i),
      load = FALSE
    )
    if (verbose)
      sapply(basename(unlist(images[i, ])), function(impath) {
        inform(c(
          "v" = glue::glue("  Finished reading '{impath}' image data for sample {i}")
        ))
      })
    return(DF)
  })))

  return(imgData)
}


#' Read spaceranger output files
#'
#' This function serves as a wrapper for \code{\link{LoadAndMergeMatrices}} and
#' \code{\link{LoadSpatialCoordinates}} to load spaceranger output files and create
#' a Seurat object. The spatial information, i.e. images and spot coordinates, are
#' stored inside the tools slot of the `Seurat` object in an object called `Staffli`.
#'
#' @details
#' \code{ReadVisiumData} takes a `data.frame` like table as input that should hold
#' certain spaceranger output file paths. The table should consist of four columns:
#' "samples", "imgs", "spotfiles" and "json".
#'
#' \itemize{
#'    \item{"samples" : file paths to expression matrices, e.g. `filtered_bc_matrix.h5`}
#'    \item{"imgs" : file paths to images, e.g. `tissue_hires_image.png`}
#'    \item{"spotfiles" : file paths to spot coordinate CSV files `tissue_positions_list.csv`}
#'    \item{"samples" : file paths to scalfactor JSOn files, e.g. `scalefactors_json.json`}
#' }
#'
#' @section Load data outside tissue:
#' Sometimes it can be useful to load data for all spots in a 10x Visium dataset, if you
#' need to explore transcripts captured outside of the tissue. In this case, you can
#' provide paths to the `raw_feature_bc_matrix.h5` files in the spaceranger output folders
#' and set `remove_spots_outside_tissue = FALSE`.
#'
#' @family pre-process
#'
#' @param infoTable A `data.frame` or `tbl` with paths to spaceranger output files
#' @param assay Assay name (default = "Spatial")
#' @param remove_spots_outside_tissue Should spots outside the tissue be removed?
#' @param verbose Print messages
#' @param ... Parameters passed to \code{\link{CreateSeuratObject}}
#'
#' @importFrom rlang abort inform
#' @importFrom glue glue
#' @importFrom Seurat CreateSeuratObject
#'
#' @examples
#' # Assemble spaceranger output files
#' samples <- Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"), "/*/filtered_feature_bc_matrix.h5"))
#' imgs <- Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"), "/*/spatial/tissue_hires_image.png"))
#' spotfiles <- Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"), "/*/spatial/tissue_positions_list.csv"))
#' json <- Sys.glob(paths = paste0(system.file("extdata", package = "STUtility2"), "/*/spatial/scalefactors_json.json"))
#'
#' # Create a tibble/data.frame with file paths
#' library(tibble)
#' infoTable <- tibble(samples, imgs, spotfiles, json, sample_id = c("mousebrain", "mousecolon"))
#'
#' # Create Seurat object
#' se <- ReadVisiumData(infoTable = infoTable)
#' se
#'
#' @export
#'
ReadVisiumData <- function (
    infoTable,
    assay = "Spatial",
    remove_spots_outside_tissue = TRUE,
    verbose = TRUE,
    ...
) {

  # Check infoTable
  if (!any(class(infoTable) %in% c("tbl", "data.frame"))) abort(glue("Invalid class '{class(infoTable)}' of 'infoTable'."))
  if (nrow(infoTable) == 0) abort(glue("'infoTable' is empty."))
  if (!all(sapply(infoTable, class) %in% "character")) abort(glue("Invalid column classes in 'infoTable'. Expecting 'character vectors.'"))
  missing_files <- infoTable |>
    select(samples, imgs, spotfiles, json) |>
    mutate(across(samples:json, ~ file.exists(.x))) |>
    summarize(across(samples:json, ~ any(!.x))) |>
    unlist()
  if (any(missing_files)) {
    checks <- missing_files[missing_files]
    sapply(names(checks), function(n) {
      abort(glue("Missing file(s) in the '{n}' column."))
    })
  }

  # Read expression matrices
  mergedMat <- LoadAndMergeMatrices(samplefiles = infoTable$samples,
                                    verbose = verbose)

  # Read spot coordinates
  coordinates <- LoadSpatialCoordinates(coordinatefiles = infoTable$spotfiles,
                                        remove_spots_outside_tissue = remove_spots_outside_tissue,
                                        verbose = verbose)

  # Make sure that coordinates and expression matrix are compatible
  if (!all(colnames(mergedMat) %in% coordinates$barcode)) {
    abort(glue("{sum(!colnames(mergedMat) %in% coordinates$barcode)} spots found in expression data but not in coordinate files."))
  }
  if (!all(coordinates$barcode %in% colnames(mergedMat))) {
    abort(glue("{sum(!coordinates$barcode %in% colnames(mergedMat))} spots found in coordinate files but not in expression data."))
  }
  if (verbose) inform(c("v" = "Expression matrices and coordinates are compatible"))

  # Sort coordinates
  coordinates <- coordinates[match(colnames(mergedMat), coordinates$barcode), ]

  # Create a Seurat object from expression matrix
  object <- CreateSeuratObject(counts = mergedMat,
                               assay = assay,
                               ...)
  if (verbose) inform(c(">" = "  Created the 'Seurat' object"))

  # Create a Stafli object
  staffli_object <- CreateStaffliObject(imgs = infoTable$imgs,
                                        meta_data = coordinates |>
                                          select(barcode, pxl_col_in_fullres, pxl_row_in_fullres, sampleID))
  if (verbose) inform(c(">" = "  Created the 'Staffli' object"))

  # Place Staffli object inside the tools slot of the Seurat object
  object@tools$Staffli <- staffli_object
  if (verbose) inform(c(">" = "  Loaded spatial data into the 'Seurat' object"))

  if (verbose) inform(c("i" = glue("Returning a 'Seurat' object with {nrow(object)} features and {ncol(object)} spots")))
  return(object)
}
