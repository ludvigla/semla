#' Update an STUtility v1 object to work with `semla`
#'
#' Some features will no longer be available when updating an STUtility v1 object.
#' STUtility v1 objects do not store information about scaling factors to convert
#' between pixels and real distances which will set limitations on some of the
#' visualization methods and spatial functions. To mitigate these issues, you can
#' either reload the data from the raw space ranger output files or manually add
#' the missing information stored in the "scalefactors_json.json" files.
#'
#' @param object A `Seurat` object created with STUtility v1
#' @param verbose Print messages
#'
#' @return A `Seurat` object compatible with semla
#'
#' @author Ludvig Larsson
#'
#' @import rlang
#' @import glue
#' @import dplyr
#' @importFrom tibble as_tibble
#'
#' @export
UpdateSTUtilityV1Object <- function (
  object,
  verbose = TRUE
) {

  # Set global variables to NULL
  barcode <- pixel_x <- pixel_y <- sample <- width <- height <- colorspace <- filesize <- density <- format <- NULL
  sampleID <- full_width <- full_height <- NULL

  # Check object
  if (!"Staffli" %in% names(object@tools)) {
    abort(glue("Missing {col_br_magenta('Staffli')} object."))
  }
  if (!all(object@tools$Staffl@platforms == "Visium")) {
    abort(glue("Only Visium data is supported."))
  }

  # get image paths
  imgs <- object@tools$Staffli@imgs
  if (verbose) cli_alert_info("Found {length(imgs)} datasets in {col_br_magenta('STUtility v1')} object")

  # get rasters
  if (length(object@tools$Staffl@rasterlists) > 0) {
    raw_rasters <- object@tools$Staffl@rasterlists[["raw"]]
    if (verbose) cli_alert_info("Found loaded raw H&E images")
    other_images <- setdiff(names(object@tools$Staffl@rasterlists), "raw")
    if (length(other_images) > 1 & verbose) {
      cli_alert_warning("Dropping processed images: {paste(other_images, collapse = ', ')}")
    }
  } else {
    raw_rasters <- list()
  }

  # Fetch meta data and convert it to a tibble with the correct column names
  if (verbose) cli_alert_info("Updating meta_data slot")
  meta_data <- object@tools$Staffl@meta.data |>
    as_tibble() |>
    select(barcode, pixel_x, pixel_y, sample) |>
    rename(pxl_col_in_fullres = pixel_x, pxl_row_in_fullres = pixel_y) |>
    mutate(sampleID = as.integer(sample)) |>
    select(-sample) |>
    mutate(barcode = gsub(pattern = "-\\d+", replacement = paste0(""), x = barcode) |> paste0("-", sampleID))

  # Update names of Seurat object
  if (verbose) cli_alert_info("Renaming spots in {col_br_magenta('Seurat')} object")
  object <- RenameCells(object, new.names = meta_data$barcode)

  # Create an image_info tibble
  if (verbose) cli_alert_info("Updating image_info slot")
  image_info <- do.call(bind_rows, lapply(seq_along(object@tools$Staffl@dims), function(i) {
    object@tools$Staffl@dims[[i]] |>
      as_tibble() |>
      select(width, height, colorspace, filesize, density, format) |>
      mutate(sampleID = paste0(i), full_width = width, full_height = height) |>
      select(format, width, height, full_width, full_height, colorspace, filesize, density, sampleID)
  }))

  # Create new Staffli object for semla
  if (verbose) cli_alert_info("Creating new Staffli object")
  staffli <- semla::CreateStaffliObject(imgs = imgs, meta_data = meta_data, image_info = image_info, scalefactors = tibble())
  staffli@rasterlists <- list(raw = raw_rasters)

  # Save staffli object to Seurat object
  object@tools$Staffli <- staffli

  # Try to unload STUtility v1
  unload_STUtilityv1 <- try({detach("package:STutility", unload = TRUE)}, silent = TRUE)
  if (inherits(unload_STUtilityv1, what = "try-error") & verbose)
    cli_alert_warning("Unloaded {col_br_magenta('STUtility v1')} package")

  # Throw warning
  if (verbose) {
    cli_alert_success("Successfully updated object")
  }

  return(object)
}
