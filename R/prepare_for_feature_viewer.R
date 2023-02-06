#' @include checks.R
#'
NULL

#' Export spatial coordinates to a JSON file
#'
#' Utility function to prepare data for \code{\link{FeatureViewer}}.
#' The exported JSON file should be exported to the same directory as
#' the H&E image tiles generated with \code{\link{TileImage}}
#'
#' @param object A `Seurat` object created with `STUtility2`
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
#' @examples
#' \dontrun{
#' libary(STUtility2)
#'
#' se_mbrain <- readRDS(system.file("extdata/mousebrain",
#'                                  "se_mbrain",
#'                                  package = "STUtility2"))
#'
#' export_coordinates(se_mbrain, outdir = "./")
#' }
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
  spatial_coords <- spatial_coords |>
    mutate(x = pxl_col_in_fullres, y = pxl_row_in_fullres) |>
    select(-pxl_col_in_fullres, -pxl_row_in_fullres) |>
    mutate(across(x, ~ scales::rescale(x = .x, to = c(0, 1),
                                       from = c(0, max_imwidth)))) |>
    mutate(across(y, ~ scales::rescale(x = .x, to = c(0, 1),
                                       from = c(0, max_imwidth)))) |>
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
