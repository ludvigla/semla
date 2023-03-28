#' @include checks.R
#'
NULL

#' Export spatial coordinates to a JSON file
#'
#' Utility function to prepare data for \code{\link{FeatureViewer}}.
#' The exported JSON file should be exported to the same directory as
#' the H&E image tiles generated with \code{\link{TileImage}}
#' 
#' Coordinates are located outside the H&E images will not be shown in the viewer.
#'
#' @param object A \code{Seurat} object created with \code{semla}
#' @param sampleNumber An integer specifying a sample ID to export
#' spatial network for
#' @param outdir Name of a directory to export JSON file to
#' @param overwrite Should an existing coordinate file be overwritten?
#' @param verbose Print messages
#'
#' @import rlang
#' @import cli
#' @import glue
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom jsonlite write_json
#' 
#' @returns No return value, writes coordinates to a file
#'
#' @examples
#' library(semla)
#'
#' se_mbrain <- readRDS(system.file("extdata/mousebrain",
#'                                  "se_mbrain",
#'                                  package = "semla"))
#'
#' export_coordinates(se_mbrain, outdir = tempdir())
#'
#' @export
#'
export_coordinates <- function (
    object,
    sampleNumber = 1L,
    outdir,
    overwrite = FALSE,
    verbose = TRUE
) {

  # Set global variables to NULL
  barcode <- pxl_col_in_fullres <- pxl_row_in_fullres <- x <- y <- sampleID <- NULL

  # Check Seurat object
  .check_seurat_object(object)

  # Check input parameters
  stopifnot(
    inherits(sampleNumber, what = c("numeric", "integer")),
    length(sampleNumber) == 1,
    inherits(outdir, what = "character"),
    dir.exists(outdir)
  )

  # Create spatial network
  if (verbose) cli_alert_info("Fetching coordinates for sample {sampleNumber}")
  spatial_coords <- GetStaffli(object)@meta_data |>
    filter(sampleID == sampleNumber) |>
    select(barcode, pxl_col_in_fullres, pxl_row_in_fullres, sampleID)

  # Rescale coordinates
  image_info <- GetStaffli(object)@image_info |>
    filter(sampleID == sampleNumber)
  max_imwidth <- image_info[image_info$sampleID == sampleNumber, ]$full_width
  max_imheight <- image_info[image_info$sampleID == sampleNumber, ]$full_height
  
  # Check if image has been padded
  if ("pad" %in% colnames(image_info)) {
    if (verbose) cli_alert_warning("Spots located outside of the H&E image(s) will be removed")
    pad <- strsplit(image_info[image_info$sampleID == sampleNumber, ]$pad, "x") |> unlist() |> as.integer()
    # Remove spots outside image and subtract padding from coordinates
    # This is to make sure that the tiled images behave as expected
    # If spots are located outside of the H&E images, the dimensions 
    # of the viewer will no longer correspond to the dimensions of the 
    # H&E image
    # See reference: https://openseadragon.github.io/examples/viewport-coordinates/
    spatial_coords <- spatial_coords |> 
      mutate(pxl_col_in_fullres = pxl_col_in_fullres - pad[1]) |> 
      mutate(pxl_row_in_fullres = pxl_row_in_fullres - pad[3]) |> 
      filter(between(x = pxl_col_in_fullres, left = 0, right = max_imwidth - (pad[1] + pad[2]))) |> 
      filter(between(x = pxl_row_in_fullres, left = 0, right = max_imheight - (pad[3] + pad[4])))
    # Reduce the image dimensions with padding
    max_imwidth_reduced <- max_imwidth - (pad[1] + pad[2])
    max_imheight_reduced <- max_imheight - (pad[3] + pad[4])
    asp_ratio <- max_imheight_reduced/max_imwidth_reduced
    # Set ranges for coordinate transformation with reduced image dimensions
    from_x <- c(0, max_imwidth_reduced)
    from_y <- c(0, max_imheight_reduced)
    # The width of the viewer goes from 0 to 1 in openseadragon
    to_x <- c(0, 1)
    # The height of the viewer goes from 0 to the aspect ratio
    to_y <- c(0, asp_ratio)
  } else {
    # Set ranges for coordinate transformation with original image dimensions
    # if no padding is found
    asp_ratio <- max_imheight/max_imwidth
    from_x <- c(0, max_imwidth)
    from_y <- c(0, max_imheight)
    to_x <- c(0, 1)
    to_y <- c(0, asp_ratio)
  }
  
  # Transform spatial coordinates using ranges from and to
  spatial_coords <- spatial_coords |>
    mutate(x = pxl_col_in_fullres, y = pxl_row_in_fullres) |>
    select(-pxl_col_in_fullres, -pxl_row_in_fullres) |>
    mutate(across(x, ~ scales::rescale(x = .x, to = to_x,
                                       from = from_x))) |>
    mutate(across(y, ~ scales::rescale(x = .x, to = to_y,
                                       from = from_y))) |>
    mutate(index = 1:n())

  # Put nodes and edges into a list
  data <- list(nodes = spatial_coords)

  # Export
  outpath <- file.path(outdir, paste0('coords_Visium_', sampleNumber, '.json'))
  if (file.exists(outpath)) {
    if (overwrite) {
      cli_alert_warning("  Replacing file {outpath |> normalizePath(winslash = '/')}")
      unlink(x = outpath, recursive = TRUE)
    } else {
      abort(glue("File {outpath |> normalizePath(winslash = '/')} already exists. Overwrite it with {col_br_magenta('overwrite=TRUE')} or change outdir path"))
    }
  }
  if (verbose) cli_alert_info("Exporting Visium coordinates to <DATADIR>/{paste0('coords_Visium_', sampleNumber, '.json')}")
  data_json <- data |>
    write_json(auto_unbox = TRUE,
               path = outpath)
}
