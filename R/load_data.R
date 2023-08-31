#' @include extdata.R
#'
NULL

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
#' @param samplefiles Character vector of file/directory paths. Paths should specify \code{.h5} or
#' \code{.tsv}/\code{.tsv.gz} files. Alternatively, the paths could specify directories including \code{barcodes.tsv},
#' \code{features.tsv} and \code{matrix.mtx} files.
#' @param verbose Print messages
#' 
#' @section IF data:
#' If the provided h5 files store antibody capture data, \code{LoadAndMergeMatrices} will
#' return a list of matrices. If multiple samples are loaded, the RNA expression matrices and 
#' antibody capture matrices will be merged and returned as separate elements of the list.
#' Note that if one or more samples only have RNA expression data, the function will add empty
#' values for those samples in the merged antibody capture matrix.
#'
#' @importFrom methods as
#' @importFrom Seurat Read10X Read10X_h5
#' @importFrom SeuratObject RowMergeSparseMatrices
#' @importFrom glue glue
#' @import cli
#' @import rlang
#' @importFrom dplyr case_when mutate
#' @importFrom tibble tibble
#' @importFrom tools file_ext
#'
#' @return A sparse matrix of class \code{dgCMatrix} or a list of sparse matrices of class \code{dgCMatrix}
#'
#' @examples
#'
#' # Load and merge two gene expression matrices
#' samples <-
#'   c(
#'     system.file(
#'       "extdata/mousebrain",
#'       "filtered_feature_bc_matrix.h5",
#'       package = "semla"
#'     ),
#'     system.file(
#'       "extdata/mousecolon",
#'       "filtered_feature_bc_matrix.h5",
#'       package = "semla"
#'     )
#'   )
#' mergedMatrix <- LoadAndMergeMatrices(samples)
#'
#' @export
LoadAndMergeMatrices <- function (
    samplefiles,
    verbose = TRUE
) {

  # Run checks
  if (!is.character(samplefiles)) abort("'samplefiles' must be a character vector.")
  checks <- tibble(samplefiles) |>
    mutate(is = case_when(dir.exists(samplefiles) ~ "dir", file.exists(samplefiles) ~ "file"))
  if (any(is.na(checks$is))) abort(c("Invalid path(s):", glue::glue("{checks$samplefiles[is.na(checks$is)]}")))

  # Load expression matrices
  if (verbose) cli_alert_info("Loading matrices:")
  exprMats <- list()
  abMats <- list()
  ab_exists <- FALSE
  antibodyMatrix <- NULL
  for (i in seq_along(samplefiles)) {
    if (checks$is[i] == "dir") {
      # Assumes that the directory contains matrix, barcodes and genes
      exprMat <- Read10X(data.dir = samplefiles[i])
    } else if (checks$is[i] == "file") {
      ext <- file_ext(samplefiles[i])
      if (ext == "h5") {
        exprMat <- suppressMessages({Read10X_h5(samplefiles[i])})
        if (verbose) cli_alert("  Finished loading expression matrix {i}")
      } else if (ext %in% c("tsv", "tsv.gz")) {
        if (!requireNamespace("data.table")) 
          abort(glue("Package {cli::col_br_magenta('data.table')} is required. Please install it with: \n",
                     "install.packages('data.table')"))
        exprMat <- data.frame(data.table::fread(samplefiles[i], sep = "\t", header = TRUE), row.names = 1)
        exprMat <- as(as.matrix(exprMat), "dgCMatrix")
      } else {
        abort(glue("Invalid file format '{ext}'"))
      }
    }
    # Check for additional matrices
    if (inherits(exprMat, what = "list")) {
      if (!all(names(exprMat) %in% c("Gene Expression", "Antibody Capture")))
        abort(glue("Data sets not supported: {paste(names(exprMat), collapse=', ')}"))
      antibodyMatrix <- exprMat[["Antibody Capture"]]
      exprMat <- exprMat[["Gene Expression"]]
      colnames(exprMat) <- gsub(pattern = "-\\d+", replacement = paste0("-", i), x = colnames(exprMat))
      colnames(antibodyMatrix) <- gsub(pattern = "-\\d+", replacement = paste0("-", i), x = colnames(antibodyMatrix))
      cli_alert("  Finished loading antibody capture matrix {i}")
      ab_exists <- TRUE
    } else {
      colnames(exprMat) <- gsub(pattern = "-\\d+", replacement = paste0("-", i), x = colnames(exprMat))
      if (ab_exists) cli_alert_warning("  No antibody capture matrix found for sample {i}. Returning empty matrix.")
      if (!is.null(antibodyMatrix)) {
        capture_ids <- rownames(antibodyMatrix)
        antibodyMatrix <- Matrix::rsparsematrix(0, nrow = nrow(antibodyMatrix), ncol = ncol(exprMat))
        rownames(antibodyMatrix) <- capture_ids
        colnames(antibodyMatrix) <- colnames(exprMat)
      }
    }
    if (ab_exists) {
      abMats <- c(abMats, list(antibodyMatrix))
    }
    exprMats <- c(exprMats, list(exprMat))
  }

  # Check gene overlap
  all.genes <- lapply(exprMats, function(exprMat) rownames(exprMat))
  intersecting.genes <- Reduce(intersect, all.genes)
  if (length(intersecting.genes) == 0) abort("No shared genes shared across datasets.")
  if (length(intersecting.genes) < 1000) {
    cli_alert_warning(col_br_red("There are only {length(intersecting.genes)} gene shared across all matrices:"))
    cli_alert(col_br_red("  Are you sure that the matrices share the same gene IDs?"))
    cli_alert(col_br_red("  Are the datasets from the same species?"))
  }
  
  # Merge matrices
  if (length(exprMats) > 1) {
    if (verbose) {
      cli_text("")
      cli_alert_info("Merging expression matrices:")
    }
    mergedMat <- RowMergeSparseMatrices(mat1 = exprMats[[1]], mat2 = exprMats[2:length(exprMats)])
    if (verbose) cli_alert_success(glue("There are {cli::col_br_blue(nrow(mergedMat))} features and {cli::col_br_magenta(ncol(mergedMat))} ",
                                        "spots in the merged expression matrix."))
  } else {
    if (verbose) cli_alert_success(glue("  There are {cli::col_br_blue(nrow(exprMats[[1]]))} features and",
                                        " {cli::col_br_magenta(ncol(exprMats[[1]]))} spots in the expression matrix."))
    mergedMat <- exprMats[[1]]
  }
  
  # Check antibody name overlaps
  if (ab_exists) {
    all.targets <- lapply(abMats, function(antibodyMatrix) rownames(antibodyMatrix))
    intersecting.targets <- Reduce(intersect, all.targets)
    
    # Merge antibody matrices
    if (length(abMats) > 1) {
      if (verbose) {
        cli_text("")
        cli_alert_info("Merging antibody capture matrices:")
      }
      mergedAbMat <- RowMergeSparseMatrices(mat1 = abMats[[1]], mat2 = abMats[2:length(abMats)])
      if (verbose) cli_alert_success(glue("There are {cli::col_br_cyan(nrow(mergedAbMat))} targets and {cli::col_br_magenta(ncol(mergedAbMat))} ",
                                          "spots in the merged antibody capture matrix."))
    } else {
      if (verbose) cli_alert_success(glue("  There are {cli::col_br_blue(nrow(abMats[[1]]))} targets and",
                                          " {cli::col_br_magenta(ncol(abMats[[1]]))} spots in the antobody capture matrix."))
      mergedAbMat <- abMats[[1]]
    }
    # Check if dimensions of the two matrices are the same
    if (ncol(mergedAbMat) != ncol(mergedMat))
      cli_alert_danger("Detected different numbers of spots in the RNA expression matrix and the antobody capture matrix")
    if (!all(colnames(mergedAbMat) == colnames(mergedMat))) 
      cli_alert_danger("The columns names of the RNA expression matrix and the antobody capture matrix differ.")
    return(list(exprMat = mergedMat, antibodyMatrix = mergedAbMat))
  } else {
    return(mergedMat)
  }

}


#' Load and merge multiple coordinate tables
#'
#' Load coordinates from \strong{'tissue_positions_list.csv'} files and merge them into a \code{tibble}.
#'
#' @details The merging process makes sure that all genes detected are present in the merged output.
#' This means that if a gene is missing in a certain dataset, the spots in that dataset will
#' be assigned with 0 expression.
#'
#' @family pre-process
#'
#' @param coordinatefiles Character vector of file paths. Paths should specify \code{.csv} files output by spaceranger
#' @param remove_spots_outside_tissue Should spots outside the tissue be removed?
#' @param verbose Print messages
#'
#' @import rlang
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom utils read.csv
#'
#' @return An object of class \code{tbl} containing spot coordinates
#'
#' @examples
#'
#' library(semla)
#'
#' # Load and merge coordinates from two samples
#' coordinatefiles <-
#'   c(system.file("extdata/mousebrain/spatial",
#'                 "tissue_positions_list.csv",
#'                 package = "semla"),
#'     system.file("extdata/mousecolon/spatial",
#'                 "tissue_positions_list.csv",
#'                 package = "semla"))
#' coordinates <- LoadSpatialCoordinates(coordinatefiles)
#'
#' @export
LoadSpatialCoordinates <- function (
    coordinatefiles,
    remove_spots_outside_tissue = TRUE,
    verbose = TRUE
) {

  # Set global variables to NULL
  selected <- pxl_col_in_fullres <- pxl_row_in_fullres <- NULL

  # Run checks
  if (!is.character(coordinatefiles)) abort("'coordinatefiles' must be a character vector.")
  checks <- tibble::tibble(coordinatefiles) |>
    dplyr::mutate(type = dplyr::case_when(file.exists(coordinatefiles) ~ "file"))
  if (any(is.na(checks$type))) abort(c("Invalid path(s):", glue::glue("{checks$coordinatefiles[is.na(checks$type)]}")))

  # Load coordinates
  if (verbose) cli_alert_info("Loading coordinates:")
  coordDF <- do.call(bind_rows, lapply(seq_along(coordinatefiles), function(i) {
    #bn <- gsub(basename(coordinatefiles[i]), pattern = ".csv", replacement = "")
    # New format since Space Ranger 2.0.0
    #use_header <- ifelse(bn == "tissue_positions", TRUE, FALSE)
    coords <- read.csv(file = coordinatefiles[i], header = TRUE) 
    if (!c("barcode") %in% colnames(coords)) {
      coords <- read.csv(file = coordinatefiles[i], header = FALSE) 
    }
    coords <- coords |>
      as_tibble() |>
      setNames(nm = c("barcode", "selected", "y", "x", "pxl_row_in_fullres", "pxl_col_in_fullres")) |>
      mutate(across(pxl_col_in_fullres:pxl_row_in_fullres, ~as.integer(.x)))
    if (remove_spots_outside_tissue) {
      coords <- coords |>
        filter(selected == 1)
    }
    if (verbose) cli_alert("  Finished loading coordinates for sample {i}")
    coords$barcode <- gsub(pattern = "-\\d+", replacement = paste0("-", i), x = coords$barcode)
    coords$sampleID <- i
    return(coords)
  }))

  # Return coordinates
  if (verbose) cli_alert_info("Collected coordinates for {cli::col_br_magenta(nrow(coordDF))} spots.")
  return(coordDF)
}

#' Load scale factors
#' 
#' Load coordinates from \strong{'scalefactors_json.json'} files and merge them into a \code{tibble}.
#' 
#' @param scalefactorfiles A character vector with file paths. Paths should specify JSON files containing
#' scale factors created with spaceranger.
#' 
#' @family pre-process
#' 
#' @import dplyr
#' @import rlang
#' @import glue
#' @importFrom jsonlite read_json
#' @importFrom tools file_ext
#' 
#' @return A \code{tbl} with scale factors
#' 
#' @examples
#' library(semla)
#' 
#' jsonfiles <-
#'   Sys.glob(paths = paste0(system.file("extdata", package = "semla"),
#'                           "/*/spatial/scalefactors_json.json"))
#' scalefactors <- LoadScaleFactors(jsonfiles)
#' scalefactors
#' 
#' @export
LoadScaleFactors <- function (
  scalefactorfiles
) {
  # Check files
  for (f in scalefactorfiles) {
    if (!file.exists(f))
      abort(glue("Invalid path '{scalefactorfiles}'. File doesn't exist."))
    if (file_ext(f) != "json")
      abort(glue("Invalid file extension '{file_ext(f)}'. Expected a JSON file."))
  }
  
  scalefactors <- do.call(bind_rows, lapply(seq_along(scalefactorfiles), function(i) {
    json_data <- data.frame(read_json(scalefactorfiles[i])) |> mutate(sampleID = paste0(i))
    return(json_data)
  })) |> as_tibble()
  
  return(scalefactors)
}


#' Load image information
#' 
#' @param imgfiles A character vector with file paths. Paths should specify PNG or JPEG files containing
#' images created with spaceranger. Alternatively, the character vector can contain URLs to images. 
#' 
#' @family pre-process
#' 
#' @import rlang
#' @import glue
#' @import dplyr
#' @import cli
#' @importFrom magick image_info
#' @importFrom tibble as_tibble
#' 
#' @return A \code{tbl} object with image info
#' 
#' @examples
#' library(semla)
#' 
#' imgs <-
#'   Sys.glob(paths = paste0(system.file("extdata", package = "semla"),
#'                           "/*/spatial/tissue_lowres_image.jpg"))
#' image_info <- LoadImageInfo(imgs)
#' image_info
#' 
#' @export
LoadImageInfo <- function (
  imgfiles
) {
  image_info <- do.call(rbind, lapply(seq_along(imgfiles), function(i) {
    f <- imgfiles[i]
    im <- try({image_read(f)}, silent = TRUE)
    if (inherits(im, "try-error"))
      abort(glue("Invalid path/url found in 'imgs'. Make sure to have ",
                 "valid image paths before running {col_br_magenta('ReadVisiumData')}"))
    image_data <- im |>
      image_info() |>
      mutate(sampleID = paste0(i),
             type = case_when(basename(f) %in% paste0("tissue_hires_image.", c("jpg", "png")) ~ "tissue_hires",
                              basename(f) %in% paste0("tissue_lowres_image.", c("jpg", "png")) ~ "tissue_lowres",
                              TRUE ~ "unknown"))
    return(image_data)
  })) |> as_tibble()
  
  return(image_info)
}


#' Update an \code{image_info} tibble
#' 
#' Takes two \code{tbl} objects as input, one with image information produced with 
#' \code{\link{LoadImageInfo}} and one with scale factors produced with 
#' \code{\link{LoadScaleFactors}}. The scale factors are used to calculate the 
#' width and height of the original H&E image which will be added to the output
#' \code{tbl}.
#' 
#' @family pre-process
#' 
#' @param image_info A \code{tbl} object with image info
#' @param scalefactors A \code{tbl} object with scale factors
#' 
#' @import dplyr
#' 
#' @return A \code{tbl} object with image info
#' 
#' @export
UpdateImageInfo <- function (
  image_info,
  scalefactors
) {
  
  # Add full res image dimensions to image data
  image_info <- image_info |> left_join(y = scalefactors, by = "sampleID") |>
    mutate(full_width = case_when(type == "tissue_lowres" ~ width/tissue_lowres_scalef,
                                  type == "tissue_hires" ~ width/tissue_hires_scalef),
           full_height = case_when(type == "tissue_lowres" ~ height/tissue_lowres_scalef,
                                   type == "tissue_hires" ~ height/tissue_hires_scalef)) |>
    select(all_of(c("format", "width", "height", "full_width", "full_height", "colorspace", "filesize", "density", "sampleID", "type")))
  
  return(image_info)
}


#' Load annotations from a CSV file
#' 
#' Takes a character vector with paths to CSV files and 
#' returns a \code{data.frame}. The resulting \code{data.frame}
#' can be added to a \code{Seurat} object using the \code{AddMetaData}
#' function. 
#' 
#' NB: The files need to be loaded in the correct order. For example, 
#' if you have an object with 5 tissue sections, you should provide
#' 5 \code{paths} in the correct order.
#' 
#' @family pre-process
#' 
#' @param paths A character vector with paths to CSV files
#' 
#' @import dplyr
#' @import rlang
#' @import glue
#' @importFrom tibble rownames_to_column
#' @importFrom tools file_ext
#' 
#' @return A \code{data.frame} object with barcode IDs and annotations
#' 
#' @examples 
#' 
#' library(semla)
#' 
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
#' se_merged <- MergeSTData(se_mcolon, y = list(se_mcolon, se_mcolon))
#' 
#' # Get annotation file(s)
#' annotation_files <- system.file("extdata/mousecolon", 
#'                                 "galt_spots.csv", 
#'                                 package = "semla") |> rep(3)
#' 
#' # Load annotations
#' annotations <- LoadAnnotationCSV(annotation_files)
#' head(annotations)
#' 
#' # Edit data.frame if needed
#' annotations <- annotations |> 
#'   rename(test = selection)
#' 
#' # Add annotations to meta data
#' se_merged <- AddMetaData(se_merged, metadata = annotations)
#' 
#' # Plot new categorical variable
#' MapLabels(se_merged, column_name = "test")
#' 
#' @export
LoadAnnotationCSV <- function (
  paths
) {
  
  # Set global variables to NULL
  barcode <- NULL
  
  # Create empty tibble
  ann <- data.frame()
  
  # Read csv files
  for (i in seq_along(paths)) {
    if (is.na(paths[i]))
      next
    if (!file.exists(paths[i]))
      abort(glue("Invalid path '{paths[i]}'. File doesn't exist."))
    if (file_ext(paths[i]) != "csv")
      abort("Invalid file format. Expected a CSV file.")
    sample_ann <- read.csv(paths[i], check.names = FALSE)
    if (!colnames(sample_ann)[1] == "barcode") {
      sample_ann <- sample_ann |> 
        rename(barcode = 1)
    }
    sample_ann <- sample_ann |> 
      mutate(barcode = gsub(pattern = "-\\d+", replacement = paste0("-", i), x = barcode))
    ann <- bind_rows(ann, sample_ann)
  }
  ann <- ann |> data.frame(row.names = 1, check.names = FALSE)
  return(ann)
}


#' Check if coordinates are located inside tissue section
#' 
#' @param coordinates A \code{tbl} with spatial coordinates
#' @param mergedMat A sparse matrix with gene expression data
#' @param image_info A \code{tbl} with image info
#' @param metaData A \code{tbl} with meta data
#' @param remove_spots_outside_HE Should spots outside the H&E be removed?
#' 
#' @import dplyr
#' @import cli
#' 
#' @return A list with \code{coordinates}, \code{mergedMat}, \code{image_info}
#' and \code{metaData}
#' 
#' @noRd
.check_coordinates <- function (
  coordinates,
  mergedMat,
  image_info,
  metaData,
  remove_spots_outside_HE
) {
  
  # Set global variables to NULL
  full_width <- full_height <- sampleID <- pxl_col_in_fullres <- pxl_row_in_fullres <- barcode <- NULL
  check_x <- check_y <- min_x <- max_x <- min_y <- max_y <- pad_before_x <- pad_after_x <- pad_after_y <- pad_before_y <- pad <- NULL
  
  # Remove coordinates that are located outside of the tissue
  check_coordinates <- coordinates |>
    group_by(sampleID) |>
    mutate(check_x = case_when(between(x = !! sym("pxl_col_in_fullres"),
                                       left = 0,
                                       right = image_info[cur_group_id(), "full_width", drop = TRUE]) ~ "inside",
                               TRUE ~ "outside"),
           check_y = case_when(between(x = !! sym("pxl_row_in_fullres"),
                                       left = 0,
                                       right = image_info[cur_group_id(), "full_height", drop = TRUE]) ~ "inside",
                               TRUE ~ "outside"))
  remove_spots <- c()
  if (any(check_coordinates$check_x == "outside")) {
    check_coordinates_x <- check_coordinates |> summarize(nOutside = sum(check_x == "outside"))
    for (ID in check_coordinates_x$sampleID) {
      cli_alert_warning("Found {check_coordinates_x |> filter(sampleID == ID) |> pull(nOutside)} spot(s) with x coordinates outside of the H&E image for sample {ID}")
    }
    remove_spots <- c(remove_spots, check_coordinates |> filter(check_x == "outside") |> pull(barcode))
  }
  if (any(check_coordinates$check_y == "outside")) {
    check_coordinates_y <- check_coordinates |> summarize(nOutside = sum(check_y == "outside"))
    for (ID in check_coordinates_y$sampleID) {
      cli_alert_warning("Found {check_coordinates_y |> filter(sampleID == ID) |> pull(nOutside)} spot(s) with y coordinates outside of the H&E image for sample {ID}")
    }
    remove_spots <- c(remove_spots, check_coordinates |> filter(check_y == "outside") |> pull(barcode))
  }
  if (length(remove_spots) > 0) {
    if (remove_spots_outside_HE) {
      cli_alert_warning("Removing {length(remove_spots)} spots from the dataset")
      coordinates <- coordinates |> filter(!barcode %in% remove_spots)
      mergedMat <- mergedMat[, colnames(mergedMat) %in% coordinates$barcode]
      metaData <- metaData[colnames(mergedMat), , drop = FALSE]
    } else {
      # Adjust width and height to match missing data
      spot_ranges <- coordinates |>
        group_by(sampleID) |>
        summarize(min_x = min(pxl_col_in_fullres),
                  max_x = max(pxl_col_in_fullres),
                  min_y = min(pxl_row_in_fullres),
                  max_y = max(pxl_row_in_fullres)) |>
        mutate(sampleID = paste0(sampleID)) |>
        left_join(y = image_info |>
                    select(sampleID, full_width, full_height) |>
                    mutate(sampleID = paste0(sampleID)), by = "sampleID")
      image_pad <- spot_ranges |>
        mutate(pad_before_x = ifelse(min_x < 0, abs(min_x), 0),
               pad_after_x = ifelse(max_x > full_width, max_x - full_width, 0),
               pad_before_y = ifelse(min_y < 0, abs(min_y), 0),
               pad_after_y = ifelse(max_y > full_height, max_y - full_height, 0)) |>
        mutate(across(pad_before_x:pad_after_y, ~ceiling(.x))) |>
        mutate(full_width_new = full_width + pad_before_x + pad_after_x,
               full_height_new = full_height + pad_before_y + pad_after_y)
      
      image_info$pad <- image_pad |>
        tidyr::unite("pad", pad_before_x:pad_after_y, sep = "x") |>
        pull(pad)
      coordinates <- coordinates |> 
        group_by(sampleID) |> 
        mutate(pxl_col_in_fullres = pxl_col_in_fullres + image_pad[cur_group_id(), ]$pad_before_x,
               pxl_row_in_fullres = pxl_row_in_fullres + image_pad[cur_group_id(), ]$pad_before_y)
      sfx <- image_pad$full_width_new/image_info$full_width
      sfy <- image_pad$full_height_new/image_info$full_height
      image_info$width <- ceiling(sfx*image_info$width)
      image_info$height <- ceiling(sfy*image_info$height)
      image_info$full_width <- ceiling(image_pad$full_width_new)
      image_info$full_height <- ceiling(image_pad$full_height_new)
      coordinates <- coordinates[match(colnames(mergedMat), coordinates$barcode), ]
      metaData <- metaData[coordinates$barcode, , drop = FALSE]
    }
  }
  return(list("coordinates" = coordinates, "mergedMat" = mergedMat, 
              "image_info" = image_info, "metaData" = metaData))
}


#' Read spaceranger output files
#'
#' This function serves as a wrapper for \code{\link{LoadAndMergeMatrices}} and
#' \code{\link{LoadSpatialCoordinates}} to load spaceranger output files and create
#' a Seurat object. The spatial information, i.e. images and spot coordinates, are
#' stored inside the tools slot of the \code{Seurat} object in \code{Staffli} object.
#'
#' @details
#' \code{ReadVisiumData} takes a \code{data.frame} like table as input that should hold
#' certain spaceranger output file paths. The table should consist of four columns:
#' "samples", "imgs", "spotfiles" and "json".
#'
#' \itemize{
#'    \item{"samples" : file paths to expression matrices, e.g. \code{filtered_bc_matrix.h5}}
#'    \item{"imgs" : file paths to images, e.g. \code{tissue_hires_image.jpg}}
#'    \item{"spotfiles" : file paths to spot coordinate CSV files \code{tissue_positions_list.csv}}
#'    \item{"json" : file paths to scale factor JSON files, e.g. \code{scalefactors_json.json}. 
#'    It is also possible to provide custom scale factors (more info below).}
#' }
#' 
#' @section Custom scale factors:
#' You can provide an additional column named 'scalefactor' with custom scale factor values.
#' The values should be between 0 to 1, where a 1 would correspond to the original H&E image
#' used as input for spaceranger. For instance, if your image was scaled from a height of 30,000
#' pixels to a height of 3,000 pixels (with the same aspect ratio), the scalefactor would be 0.1.
#'
#' @section Load data outside tissue section:
#' Sometimes it can be useful to load data for all spots in a 10x Visium dataset, if you
#' need to explore transcripts captured outside of the tissue. In this case, you can
#' provide paths to the \code{raw_feature_bc_matrix.h5} files in the spaceranger output folders
#' and set \code{remove_spots_outside_tissue = FALSE}.
#'
#' @section Filter data:
#' If you want to filter out spots and features, you can pass the \code{min.cells} and
#' \code{min.features} parameters (see \code{\link{CreateSeuratObject}} for more details);
#' however, it is recommended to use the \code{\link{SubsetSTData}} function for filtering
#' after the object has been created.
#'
#' @section Additional annotations:
#' You can add a column named 'annotation_files' with paths to CSV files containing 
#' annotations. These files should have the same structure as the CSV the files exported 
#' by Loupe Browser. One column with the barcode IDs and one column with the annotation labels.
#' The annotations will be added to the \code{meta.data} slot of the returned \code{\link{Seurat}} object.
#'
#' @family pre-process
#'
#' @param infoTable A \code{data.frame} or \code{tbl} with paths to spaceranger output files
#' @param assay Assay name (default = "Spatial")
#' @param remove_spots_outside_HE Should spots outside the H&E be removed? This option
#' can be useful for CytAssist data when the H&E image only cover a smaller part of the
#' entire tissue section.
#' @param remove_spots_outside_tissue Should spots outside the tissue be removed?
#' @param verbose Print messages
#' @param ... Parameters passed to \code{\link{CreateSeuratObject}}
#'
#' @import cli
#' @import rlang
#' @import dplyr
#' @importFrom zeallot %<-%
#' @importFrom glue glue
#' @importFrom Seurat CreateSeuratObject
#' @importFrom tidyr uncount
#'
#' @return A \code{\link{Seurat}} object with additional spatial information stored in
#' a \code{Staffli} object
#'
#' @examples
#' # Assemble spaceranger output files
#' samples <-
#'   Sys.glob(paths = paste0(system.file("extdata", package = "semla"),
#'                           "/*/filtered_feature_bc_matrix.h5"))
#' imgs <-
#'   Sys.glob(paths = paste0(system.file("extdata", package = "semla"),
#'                           "/*/spatial/tissue_lowres_image.jpg"))
#' spotfiles <-
#'   Sys.glob(paths = paste0(system.file("extdata", package = "semla"),
#'                           "/*/spatial/tissue_positions_list.csv"))
#' json <-
#'   Sys.glob(paths = paste0(system.file("extdata", package = "semla"),
#'                           "/*/spatial/scalefactors_json.json"))
#'
#' # Create a tibble/data.frame with file paths
#' library(tibble)
#' infoTable <- tibble(samples, imgs, spotfiles, json, 
#'                     sample_id = c("mousebrain", "mousecolon"))
#'
#' # Create Seurat object
#' se <- ReadVisiumData(infoTable = infoTable)
#' se
#' 
#' # Add additional annotations from CSV files
#' annotation_file <- 
#'   Sys.glob(paths = paste0(system.file("extdata", package = "semla"),
#'                           "/*/galt_spots.csv"))
#' annotation_files <- c(NA_character_, annotation_file)
#' 
#' # Create a tibble/data.frame with file paths
#' infoTable <- tibble(samples, imgs, spotfiles, json, 
#'                     sample_id = c("mousebrain", "mousecolon"), annotation_files)
#' 
#' se <- ReadVisiumData(infoTable)
#'
#' @export
#'
ReadVisiumData <- function (
  infoTable,
  assay = "Spatial",
  remove_spots_outside_HE = FALSE,
  remove_spots_outside_tissue = TRUE,
  verbose = TRUE,
  ...
) {

  # Set global variables to NULL
  barcode <- width <- height <- scalefactor <- custom_scalef <- NULL

  if (verbose) cli_h2("Reading 10x Visium data")

  # Check infoTable
  if (!all(c("samples", "imgs", "spotfiles") %in% colnames(infoTable)))
    abort("One or several of 'samples', 'imgs' and 'spotfiles' are missing from infoTable.")
  if (!any(c("json", "scalefactor") %in% colnames(infoTable)))
    abort("One of 'json' or 'scalefactor' columns needs to be provided")
  if (!any(class(infoTable) %in% c("tbl", "data.frame"))) abort(glue("Invalid class '{class(infoTable)}' of 'infoTable'."))
  if (nrow(infoTable) == 0) abort(glue("'infoTable' is empty."))
  if (!all(sapply(infoTable[, c("samples", "imgs", "spotfiles")], class) %in% "character")) abort(glue("Invalid column classes in 'infoTable'. Expected 'character' vectors."))
  if ("scalefactor" %in% colnames(infoTable)) {
    if (!inherits(infoTable[, "scalefactor", drop = TRUE], what = c("numeric", "integer")))
      abort(glue("Invalid column class for 'scalefactor'. Expected a 'numeric'."))
    if (!all(between(x = infoTable[, "scalefactor", drop = TRUE], left = 0.001, right = 1)))
      abort(glue("Scalefactors need to be in the range 0.001-1."))
  }
  missing_files <- infoTable |>
    select(all_of(c("samples", "spotfiles")), contains("json")) |>
    mutate(across(everything(), ~ file.exists(.x))) |>
    summarize(across(everything(), ~ any(!.x))) |>
    unlist()
  if (any(missing_files)) {
    checks <- missing_files[missing_files]
    sapply(names(checks), function(n) {
      abort(glue("Missing file(s) in the '{n}' column."))
    })
  }
  if ("json" %in% colnames(infoTable)) {
    if (!all(file.exists(infoTable |> pull(all_of("json"))))) {
      abort(glue("Missing json file(s)."))
    }
  }
  
  # Check if an annotation_files column was provided
  if ("annotation_files" %in% colnames(infoTable)) {
    annfiles <- infoTable |> pull(all_of("annotation_files"))
    if (verbose) cli_alert_info("Loading annotatons from CSV files")
    ann <- LoadAnnotationCSV(annfiles)
    infoTable <- infoTable |> 
      select(-all_of("annotation_files"))
    add_annotations <- TRUE
  } else {
    add_annotations <- FALSE
  }

  # Keep additional infoTable columns
  if (ncol(infoTable) > 4) {
    additionalMetaData <- infoTable |>
      select(-all_of(c("samples", "imgs", "spotfiles", "json")))
  } else {
    additionalMetaData <- NULL
  }

  # Read expression matrices
  mergedMat <- LoadAndMergeMatrices(samplefiles = infoTable$samples,
                                    verbose = verbose)
  
  # Check if multiple matrices are loaded
  if (inherits(mergedMat, what = "list")) {
    antibodyMatrix <- mergedMat[["antibodyMatrix"]]
    mergedMat <- mergedMat[["exprMat"]]
  } else {
    antibodyMatrix <- NULL
  }

  # Read spot coordinates
  coordinates <- LoadSpatialCoordinates(coordinatefiles = infoTable$spotfiles,
                                        remove_spots_outside_tissue = remove_spots_outside_tissue,
                                        verbose = verbose)

  # Create Seurat meta data
  if (!is.null(additionalMetaData)) {
    metaData <- coordinates |>
      select(barcode) |>
      bind_cols(
        additionalMetaData |>
          mutate(count = table(coordinates$sampleID) |> as.numeric()) |>
          uncount(count)
      ) |>
      data.frame(row.names = 1)
  } else {
    metaData <- NULL
  }

  # Make sure that coordinates and expression matrix are compatible
  if (!all(colnames(mergedMat) %in% coordinates$barcode)) {
    cli_alert_danger("{sum(!colnames(mergedMat) %in% coordinates$barcode)} spots found in the merged expression matrix that are missing in the coordinate files.")
    cli_alert_danger("Removing {sum(!colnames(mergedMat) %in% coordinates$barcode)} spots from the merged expression matrix")
    mergedMat <- mergedMat[, colnames(mergedMat) %in% coordinates$barcode]
  }
  if (!all(coordinates$barcode %in% colnames(mergedMat))) {
    cli_alert_danger("{sum(!coordinates$barcode %in% colnames(mergedMat))} spots found in coordinate files that are missing in the expression data.")
    cli_alert_danger("Removing {sum(!coordinates$barcode %in% colnames(mergedMat))} spots from the meta data")
    metaData <- metaData[colnames(mergedMat), , drop = FALSE]
  }
  if (verbose) cli_h3(text = "Creating `Seurat` object")
  if (verbose) cli_alert_success("Expression matrices and coordinates are compatible")

  # Read image info
  image_info <- LoadImageInfo(infoTable$imgs)

  # Read scale factors
  if ("scalefactor" %in% colnames(infoTable)) {
    if (verbose) cli_alert_info("Using custom scalefactor(s) from infoTable")
    scalefactors <- infoTable |> 
      select(scalefactor) |> 
      dplyr::rename(custom_scalef = scalefactor) |> 
      mutate(sampleID = paste0(1:n()))
    # Add full res image dimensions to image data
    image_info <- image_info |> left_join(y = scalefactors, by = "sampleID") |>
      mutate(full_width = width/custom_scalef,
             full_height = height/custom_scalef) |>
      select(all_of(c("format", "width", "height", "full_width", "full_height", "colorspace", "filesize", "density", "sampleID", "type")))
  } else {
    scalefactors <- LoadScaleFactors(infoTable$json)
    # Add full res image dimensions to image data
    image_info <- UpdateImageInfo(image_info, scalefactors)
  }

  # Sort coordinates
  coordinates <- coordinates[match(colnames(mergedMat), coordinates$barcode), ]
  
  # Check coordinates to make sure that they are located inside the H&E image
  c(coordinates, mergedMat, image_info, metaData) %<-% 
    .check_coordinates(coordinates, mergedMat, image_info, metaData, remove_spots_outside_HE)

  # Create a Seurat object from expression matrix
  object <- CreateSeuratObject(counts = mergedMat,
                               assay = assay,
                               meta.data = metaData,
                               ...)
  if (verbose) cli_alert_info("Created `Seurat` object")

  # Create a Staffli object
  staffli_object <- CreateStaffliObject(imgs = infoTable$imgs,
                                        meta_data = coordinates |>
                                          ungroup() |> 
                                          select(all_of(c("barcode", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID"))),
                                        image_info = image_info,
                                        scalefactors = scalefactors)
  if (verbose) cli_alert_info("Created `Staffli` object")

  # Place Staffli object inside the tools slot of the Seurat object
  object@tools$Staffli <- staffli_object
  
  # Add additional annotations if available
  if (add_annotations) {
    object <- AddMetaData(object, metadata = ann)
  }
  
  # Add additional antibody capture matrix if available
  if (!is.null(antibodyMatrix)) {
    if (verbose) cli_alert_info("Adding antibody capture data to assay 'AbCapture'")
    object[["AbCapture"]] <- suppressWarnings({CreateAssayObject(counts = antibodyMatrix[, coordinates$barcode])})
  }

  if (verbose) cli_alert_success(glue("Returning a `Seurat` object with {cli::col_br_blue(nrow(object))}",
                                   " features and {cli::col_br_magenta(ncol(object))} spots"))
  return(object)
}
