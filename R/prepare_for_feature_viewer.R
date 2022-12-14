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
#' @param sampleID An integer specifying a sample ID to export
#' spatial network for
#' @param outdir Name of a directory to export JSON file to
#' @param verbose Print messages
#'
#' @import rlang
#' @import cli
#' @import glue
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom jsonlite write_json
#'
#' @export
#'
export_coordinates <- function (
    object,
    sampleID = 1L,
    outdir,
    verbose = TRUE
) {

  # Set global variables to NULL

  # Check Seurat object
  .check_seurat_object(object)

  # Check input parameters
  stopifnot(
    inherits(sampleID, what = c("numeric", "integer")),
    length(sampleID) == 1,
    inherits(outdir, what = "character"),
    dir.exists(outdir)
  )

  # Create spatial network
  if (verbose) cli_alert_info("Fetching coordinates for sample {sampleID}")
  spatial_coords <- GetStaffli(object)@meta_data |>
    filter(sampleID == sampleID) |>
    select(barcode, pxl_col_in_fullres, pxl_row_in_fullres, sampleID)

  # Rescale coordinates
  image_info <- GetStaffli(object)@image_info |>
    filter(sampleID == sampleID)
  spatial_coords <- spatial_coords |>
    mutate(x = pxl_col_in_fullres, y = pxl_row_in_fullres) |>
    select(-pxl_col_in_fullres, -pxl_row_in_fullres) |>
    mutate(across(x, ~ scales::rescale(x = .x, to = c(0, 1),
                                       from = c(0, max(image_info$full_width, image_info$full_height))))) |>
    mutate(across(y, ~ scales::rescale(x = .x, to = c(0, 1),
                                       from = c(0, max(image_info$full_width, image_info$full_height))))) |>
    mutate(index = 1:n())

  # Put nodes and edges into a list
  data <- list(nodes = spatial_coords)

  # Export
  outpath <- file.path(outdir, paste0('coords_Visium_', sampleID, '.json'))
  if (verbose) cli_alert_info("Exporting Visium coordinates to <DATADIR>/{paste0('coords_Visium_', sampleID, '.json')}")
  data_json <- data |>
    write_json(auto_unbox = TRUE,
               path = outpath)
  if (verbose) cli_alert_success("Finished!")
}
