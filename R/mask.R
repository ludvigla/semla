#' @include generics.R
#'
NULL


#' @description
#' Image masking can sometimes be useful when you want to remove background from
#' an H&E image. This is usually only relevant when creating figures. This function
#' provides a simple image masking technique based on blob extraction. Since image
#' masking is a challenging task, it is not guaranteed that it will work on all
#' H&E images. Certain artefacts can be particularly difficult to remove, for example
#' bubbles or if the tissue stain hasn't been removed properly.
#'
#' Even if some artefacts are detected by the algorithm, you can provide a set of
#' spot coordinates that will be used to try to filter out artefacts that are not
#' covered by spots.
#'
#' Note that this method is only useful on H&E images. Other stains or image types
#' (e.g. immunofluorescence images) will most likely not work.
#'
#' @param xy_coords Optional tibble with spot coordinates matching the input image
#' @param minPixels Minimum area for blobs used to remove small artefacts given in pixels.
#' @param verbose Print messages
#'
#' @importFrom magick image_convert image_channel image_blur image_normalize
#' image_quantize image_threshold image_threshold image_connect image_convert
#' image_read image_transparent image_composite image_negate
#' @importFrom rlang try_fetch abort
#' @import cli
#'
#' @rdname mask-images
#'
#' @section default method:
#' returns a masked 'magick-image' object
#'
#' @examples
#'
#' library(STUtility2)
#' library(magick)
#' library(dplyr)
#'
#' # Load image
#' lowresimagefile <- system.file("extdata/mousebrain/spatial",
#'                                "tissue_lowres_image.png",
#'                                package = "STUtility2")
#' im <- image_read(lowresimagefile)
#'
#' # Load coordinates
#' coordinatesfile <- system.file("extdata/mousebrain/spatial",
#'                                "tissue_positions_list.csv",
#'                                package = "STUtility2")
#' xy <- LoadSpatialCoordinates(coordinatefiles = coordinatesfile)
#'
#' # Load scalefactors
#' json <- system.file("extdata/mousebrain/spatial",
#'                     "scalefactors_json.json",
#'                     package = "STUtility2")
#' scalefactors <- jsonlite::read_json(path = json)
#' xy <- xy |>
#'   mutate(across(pxl_row_in_fullres:pxl_col_in_fullres,
#'                 ~ round(.x*scalefactors$tissue_lowres_scalef))) |>
#'   select(pxl_row_in_fullres:pxl_col_in_fullres)
#'
#' im_masked <- MaskImages(im, xy_coords = xy)
#'
#' par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))
#' im |> as.raster() |> plot()
#' im_masked |> as.raster() |> plot()
#'
#' @export
#'
MaskImages.default <- function (
  object,
  xy_coords = NULL,
  minPixels = 100,
  verbose = TRUE,
  ...
) {

  # Input image must be of class 'magick-image'
  if (!inherits(object, what = "magick-image"))
    abort(glue("Invalid class '{class(object)}', expected a 'magick-image' object"))

  # Validate minPixels
  if (!inherits(minPixels, what = "numeric"))
    abort(glue("Invalid class '{class(minPixels)}' for 'minPixels', expected a 'numeric'"))
  if (length(minPixels) != 1)
    abort(glue("'minPixels' should be a 'numeric' of length 1"))

  # Apply image conversions
  if (verbose) cli_alert_info("Segmenting image using blob extraction")

  im_blobs <- object |>
    image_convert(colorspace = "cmyk") |>
    image_channel(channel = "Magenta") |>
    image_blur(sigma = 2) |>
    image_normalize() |>
    image_quantize(max = 5) |>
    image_threshold(type = "black", threshold = "2%") |>
    image_threshold(type = "white", threshold = "2%") |>
    image_connect(connectivity = 1) |>
    image_convert(colorspace = "gray")

  # Convert image to bitmap array
  im_blobs_bm <- as.integer(im_blobs[[1]])

  # Find blob sizes
  blob_sizes <- table(im_blobs_bm)
  blob_names <- names(blob_sizes)
  blob_sizes <- blob_sizes |>
    as.numeric() |>
    setNames(blob_names) |>
    sort(decreasing = TRUE)

  # Remove blobs with fewer than minPixels pixels
  blob_sizes <- blob_sizes[blob_sizes > minPixels]

  # Remove background blob
  blob_sizes <- blob_sizes[setdiff(names(blob_sizes), "0")]

  # Remove small blobs
  if (verbose) cli_alert_info("Filtering out blobs with fewer than {minPixels} pixels")
  im_blobs_bm[!im_blobs_bm %in% as.numeric(names(blob_sizes))] <- 0

  # Find blobs that are overlapping with spatial coordinates
  # if coordinates are provided
  if (!is.null(xy_coords)) {
    if (verbose) cli_alert_info("Filtering out blobs that do not overlap with provided coordinates")
    if (!inherits(xy_coords, what = "tbl"))
      abort(glue("'xy_coords' should be a tibble, got '{class(xy_coords)}'"))
    if (nrow(xy_coords) == 0)
      abort("'xy_coords' is empty")
    if (!all(c("pxl_row_in_fullres", "pxl_col_in_fullres") %in% colnames(xy_coords)))
      abort("Spot coordinates are missing from 'xy_coords'")

    blobs_overlapping <- sapply(1:nrow(xy_coords), function(i) {
      error <- try_fetch({
        im_blobs_bm[xy_coords[i, ]$pxl_row_in_fullres,
                    xy_coords[i, ]$pxl_col_in_fullres, ]
      },
      error = function(cnd) {
        abort(glue("Spot coordinates are outside of the image. Make sure that ",
                   "the spot coordinates are scaled to fit the input image."), parent = cnd)
      })
    }) |>
      unique() |>
      setdiff(0)

    # Remove non-overlapping blobs
    im_blobs_bm[!im_blobs_bm %in% blobs_overlapping] <- 0
  }

  # read image from modified bitmap array
  im_filtered <- image_read(im_blobs_bm) |>
    image_negate() |>
    image_transparent(color = "#000000")

  # Combine original image with filtered mask
  im_final <- image_composite(image = object, composite_image = im_filtered)
  if (verbose) cli_alert_success("Composed masked image from selected blobs")

  return(im_final)

}


#' @param section_numbers An integer vector specifying samples to mask
#'
#' @importFrom rlang %||% abort
#' @import cli
#' @importFrom dplyr filter mutate group_by group_split across
#' @importFrom magick image_read
#'
#' @rdname mask-images
#'
#' @section Seurat:
#' Returns a Seurat object with masked images
#'
#' @examples
#'
#' library(STUtility2)
#'
#' # Load example Visium data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "STUtility2"))
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "STUtility2"))
#' se_merged <- MergeSTData(se_mbrain, se_mcolon) |> LoadImages()
#' se_merged <- se_merged |> MaskImages()
#'
#' # Plot masked images
#' ImagePlot(se_merged)
#'
#' @export
#'
MaskImages.Seurat <- function (
    object,
    section_numbers = NULL,
    minPixels = 100,
    verbose = TRUE,
    ...
) {

  # Set global variables to NulL
  sampleID <- pxl_col_in_fullres <- pxl_row_in_fullres <- NULL

  if (verbose) cli_h2("Masking image(s)")

  # Check Seurat object
  .check_seurat_object(object)
  .check_seurat_images(object)

  # Get Staffli object
  st_object <- GetStaffli(object)
  if (verbose) cli_alert_info("Found {nrow(st_object@image_info)} samples")

  # Check section numbers
  section_numbers <- section_numbers %||% st_object@image_info$sampleID |> as.integer()
  if (!all(section_numbers %in% st_object@image_info$sampleID))
    abort(glue("Invalid numbers for 'section_numbers'. ",
    "Sections available: {st_object@image_info$sampleID}"))

  # Get raw images
  raw_images <- masked_images <- st_object@rasterlists[["raw"]][section_numbers]
  if (verbose) cli_alert_info("Fetched images")

  # Get coordinates and scalefactors
  sfs <- st_object@image_height/st_object@image_info$full_height
  xy_coords_list <- st_object@meta_data |>
    filter(sampleID %in% section_numbers) |>
    group_by(sampleID) |>
    group_split()
  xy_coords_list <- lapply(seq_along(xy_coords_list), function(i) {
    xy_coords <- xy_coords_list[[i]] |>
      select(pxl_col_in_fullres:pxl_row_in_fullres) |>
      mutate(across(pxl_col_in_fullres:pxl_row_in_fullres, ~ round(.x*sfs[i])))
    return(xy_coords)
  })
  if (verbose) cli_alert_info("Fetched spot coordinates")

  # Mask images
  selected_masked_images <- lapply(seq_along(raw_images), function(i) {
    if (verbose) cli_alert_info("Processing sample {i}")
    im_masked <- MaskImages(raw_images[[i]] |> image_read(),
                            xy_coords = xy_coords_list[[i]],
                            minPixels = minPixels,
                            verbose = verbose) |>
      as.raster()
    return(im_masked)
  })

  # Place selected masked images in raster list
  masked_images[section_numbers] <- selected_masked_images

  # Place masked images in Staffli object
  st_object@rasterlists[["raw"]] <- masked_images

  # Place modified Staffli object in Seurat object
  object@tools$Staffli <- st_object
  if (verbose) cli_alert_success("Returning Seurat object with masked images")

  # Return Seurat object
  return(object)

}
