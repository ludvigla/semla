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
#' @param samplefiles Character vector of file/directory paths. Paths should specify `.h5` or
#' `.tsv`/`.tsv.gz` files. Alternatively, the paths could specify directories including `barcodes.tsv`,
#' `features.tsv` and `matrix.mtx` files.
#' @param verbose Print messages
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
#' @return A sparse matrix of class `dgCMatrix`
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
    mutate(is = case_when(file.exists(samplefiles) ~ "file", dir.exists(samplefiles) ~ "dir"))
  if (any(is.na(checks$is))) abort(c("Invalid path(s):", glue::glue("{checks$samplefiles[is.na(checks$is)]}")))

  # Load expression matrices
  if (verbose) cli_alert_info("Loading matrices:")
  exprMats <- lapply(seq_along(samplefiles), function(i) {
    if (checks$is[i] == "dir") {
      # Assumes that the directory contains matrix, barcodes and genes
      exprMat <- Read10X(data.dir = samplefiles[i])
    } else if (checks$is[i] == "file") {
      ext <- file_ext(samplefiles[i])
      if (ext == "h5") {
        exprMat <- Read10X_h5(samplefiles[i])
      } else if (ext %in% c("tsv", "tsv.gz")) {
        if (!requireNamespace("data.table")) install.packages("data.table")
        exprMat <- data.frame(data.table::fread(samplefiles[i], sep = "\t", header = TRUE), row.names = 1)
        exprMat <- as(as.matrix(exprMat), "dgCMatrix")
      }
    }
    if (verbose) cli_alert("  Finished loading expression matrix {i}")
    colnames(exprMat) <- gsub(pattern = "-\\d+", replacement = paste0("-", i), x = colnames(exprMat))
    return(exprMat)
  })

  # Check gene overlap
  all.genes <- lapply(exprMats, function(exprMat) {rownames(exprMat)})
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
      cli_alert_info("Merging matrices:")
    }
    mergedMat <- RowMergeSparseMatrices(mat1 = exprMats[[1]], mat2 = exprMats[2:length(exprMats)])
    if (verbose) cli_alert_success(glue("There are {cli::col_br_blue(nrow(mergedMat))} features and {cli::col_br_magenta(ncol(mergedMat))} ",
                                        "spots in the merged matrix."))
    return(mergedMat)
  } else {
    if (verbose) cli_alert_info("only 1 expression matrix loaded.")
    if (verbose) cli_alert_success(glue("  There are {cli::col_br_blue(nrow(exprMats[[1]]))} features and",
                                     " {cli::col_br_magenta(ncol(exprMats[[1]]))} spots in the matrix."))
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
#' @import rlang
#' @import dplyr
#' @importFrom utils read.csv
#'
#' @return An object of class `tbl` (`tibble`)
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
    bn <- gsub(basename(coordinatefiles[i]), pattern = ".csv", replacement = "")
    # New format since Space Ranger 2.0.0
    use_header <- ifelse(bn == "tissue_positions", TRUE, FALSE)
    coords <- read.csv(file = coordinatefiles[i], header = use_header) |>
      tibble::as_tibble() |>
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

# TODO: not yet implemented, only to be used for SpatialExperiment class
#' Read image data
#'
#' Load image related data from \strong{'tissue_hires_image.jpg'} files.
#'
#' @family pre-process
#'
#' @param images An object of class `tibble` containing paths to images in PNG format, with
#' one row per sample. Paths should specify `.png` files output by spaceranger such as
#' `tissue_lowres_image.jpg` or `tissue_hires_image.jpg`. You do not have to load both
#' H&E images.
#' @param jsonfiles A character vector with file paths. Paths should specify `.json` files containing
#' scalefactors output by spaceranger.
#' @param verbose Print messages
#'
#' @import glue
#' @import cli
#' @import rlang
#' @import dplyr
#' @importFrom tools file_ext
#'
#' @return An object of class `DFrame`
#'
#' @noRd
LoadImageData <- function (
    images,
    jsonfiles,
    verbose = TRUE
) {

  # Check files
  check <- sapply(images, file.exists)
  if (any(!check)) abort(message = c("Invalid path(s):", glue::glue("File path '{unlist(images[!check])}' does not exist")))

  # read image data
  if (verbose) cli_alert_info(glue::glue("Reading image data for {nrow(images)} samples."))
  imgData <- do.call(rbind, lapply(1:nrow(images), function(i) {
    png.files <- unlist(images[i, ])
    exts <- file_ext(png.files)
    check <- exts %in% "png"
    if (!all(check)) abort(glue::glue("Invalid image format: '{exts[!check]}'"))
    if (!requireNamespace("BiocManager"))
      install.packages("BiocManager")
    if (!requireNamespace("SpatialExperiment"))
      install.packages("SpatialExperiment")
    DF <- SpatialExperiment::readImgData(
      imageSources = unlist(images[i, ]),
      scaleFactors = jsonfiles[i],
      sample_id = paste0(i),
      load = FALSE
    )
    if (verbose)
      sapply(basename(unlist(images[i, ])), function(impath) {
        cli_alert_info("Finished reading '{impath}' image data for sample {i}")
      })
    return(DF)
  }))

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
#'    \item{"imgs" : file paths to images, e.g. `tissue_hires_image.jpg`}
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
#' @section Filter data:
#' If you want to filter out spots and features, you can pass the `min.cells` and
#' `min.features` parameters (see \code{\link{CreateSeuratObject}} for more details);
#' however, it is recommended to use the \code{\link{SubsetSTData}} function for filtering
#' after the object has been created.
#'
#' @family pre-process
#'
#' @param infoTable A `data.frame` or `tbl` with paths to spaceranger output files
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
#' @importFrom glue glue
#' @importFrom Seurat CreateSeuratObject
#' @importFrom magick image_read image_info
#' @importFrom jsonlite read_json
#' @importFrom dplyr select bind_cols mutate case_when left_join
#' @importFrom tibble as_tibble
#' @importFrom tidyr uncount
#'
#' @return A \code{\link{Seurat}} object with additional spatial information
#'
#' @examples
#' # Assemble spaceranger output files
#' samples <-
#'   Sys.glob(paths = paste0(system.file("extdata", package = "semla"),
#'                           "/*/filtered_feature_bc_matrix.h5"))
#' imgs <-
#'   Sys.glob(paths = paste0(system.file("extdata", package = "semla"),
#'                           "/*/spatial/tissue_hires_image.jpg"))
#' spotfiles <-
#'   Sys.glob(paths = paste0(system.file("extdata", package = "semla"),
#'                           "/*/spatial/tissue_positions_list.csv"))
#' json <-
#'   Sys.glob(paths = paste0(system.file("extdata", package = "semla"),
#'                           "/*/spatial/scalefactors_json.json"))
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
    remove_spots_outside_HE = FALSE,
    remove_spots_outside_tissue = TRUE,
    verbose = TRUE,
    ...
) {

  # Set global variables to NULL
  samples <- imgs <- spotfiles <- json <- barcode <- width <- height <- full_width <- full_height <- NULL
  colorspace <- density <- sampleID <- type <- pxl_col_in_fullres <- pxl_row_in_fullres <- filesize <- NULL
  check_x <- check_y <- min_x <- max_x <- min_y <- max_y <- pad_before_x <- pad_after_x <- pad_after_y <- pad_before_y <- pad <- NULL

  if (verbose) cli_h2("Reading 10x Visium data")

  # Check infoTable
  if (!all(c("samples", "imgs", "spotfiles", "json") %in% colnames(infoTable)))
    abort("One or several of 'samples', 'imgs', 'spotfiles' and 'json' are missing from infoTable.")
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

  # Keep additional infoTable columns
  if (ncol(infoTable) > 4) {
    additionalMetaData <- infoTable |>
      select(-samples, -imgs, -spotfiles, -json)
  } else {
    additionalMetaData <- NULL
  }

  # Read expression matrices
  mergedMat <- LoadAndMergeMatrices(samplefiles = infoTable$samples,
                                    verbose = verbose)

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
    cli_alert_danger("{sum(!colnames(mergedMat) %in% coordinates$barcode)} spots found in the merged expression matrix but not in coordinate files.")
    cli_alert_danger("Removing {sum(!colnames(mergedMat) %in% coordinates$barcode)} spots from the merged expression matrix")
    mergedMat <- mergedMat[, colnames(mergedMat) %in% coordinates$barcode]
  }
  if (!all(coordinates$barcode %in% colnames(mergedMat))) {
    cli_alert_danger("{sum(!coordinates$barcode %in% colnames(mergedMat))} spots found in coordinate files but not in expression data.")
    cli_alert_danger("Removing {sum(!coordinates$barcode %in% colnames(mergedMat))} spots from the meta data")
    metaData <- metaData[colnames(mergedMat), ]
  }
  if (verbose) cli_h3(text = "Creating `Seurat` object")
  if (verbose) cli_alert_success("Expression matrices and coordinates are compatible")

  # Read image info
  image_info <- do.call(rbind, lapply(seq_along(infoTable$imgs), function(i) {
    f <- infoTable$imgs[i]
    image_data <- image_read(f) |>
      image_info() |>
      mutate(sampleID = paste0(i),
             type = case_when(basename(f) %in% paste0("tissue_hires_image.", c("jpg", "png")) ~ "tissue_hires",
                              basename(f) %in% paste0("tissue_lowres_image", c("jpg", "png")) ~ "tissue_lowres",
                              TRUE ~ "unknown"))
    return(image_data)
  })) |> as_tibble()

  # Read scale factors
  scalefactors <- do.call(rbind, lapply(seq_along(infoTable$json), function(i) {
    json_data <- data.frame(read_json(infoTable$json[i])) |> mutate(sampleID = paste0(i))
    return(json_data)
  })) |> as_tibble()

  # Add full res image dimensions to image data
  image_info <- image_info |> left_join(y = scalefactors, by = "sampleID") |>
    mutate(full_width = case_when(type == "tissue_lowres" ~ width/tissue_lowres_scalef,
                                  type == "tissue_hires" ~ width/tissue_hires_scalef),
           full_height = case_when(type == "tissue_lowres" ~ height/tissue_lowres_scalef,
                                   type == "tissue_hires" ~ height/tissue_hires_scalef)) |>
    select(format, width, height, full_width, full_height, colorspace, filesize, density, sampleID, type)


  # Sort coordinates
  coordinates <- coordinates[match(colnames(mergedMat), coordinates$barcode), ]

  # Remove coordinates that are located outside of the tissue
  check_coordinates <- coordinates |>
    group_by(sampleID) |>
    mutate(check_x = case_when(between(x = pxl_col_in_fullres,
                                       left = 0,
                                       right = image_info[cur_group_id(), "full_width", drop = TRUE]) ~ "inside",
                               TRUE ~ "outside"),
           check_y = case_when(between(x = pxl_row_in_fullres,
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
    }
  }

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
                                          select(barcode, pxl_col_in_fullres, pxl_row_in_fullres, sampleID),
                                        image_info = image_info,
                                        scalefactors = scalefactors)
  if (verbose) cli_alert_info("Created `Staffli` object")

  # Place Staffli object inside the tools slot of the Seurat object
  object@tools$Staffli <- staffli_object

  if (verbose) cli_alert_success(glue("Returning a `Seurat` object with {cli::col_br_blue(nrow(object))}",
                                   " features and {cli::col_br_magenta(ncol(object))} spots"))
  return(object)
}
