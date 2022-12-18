#' Tile an H&E image
#'
#' This function takes an image of class `magick-image` and create a tile
#' map. The size of each tile is 256x256 pixels and the number of zoom levels
#' are determined from the automatically.
#'
#' @param im An image of class `magick-image`
#' @param sampleID The section number to use. This number will be appended to the output files names
#' @param outpath A string specifying an output directory to save the tiled image in.
#' If this is not provided, a temporary directory will be created
#' @param maxZoomLevel Max zoom level
#' @param maxImgWidth Safety threshold to make sure that the zoom level doesn't get too deep.
#' @param nCores Number of cores to use for threading
#' @param verbose Print messages
#'
#' @importFrom parallel detectCores mclapply
#' @importFrom magick image_read image_info image_crop image_scale image_blank image_composite image_transparent image_write
#' @importFrom tibble tibble
#' @importFrom rlang %||%
#' @import cli
#'
#' @return Path to tiles
#'
#' @examples
#' \dontrun{
#' if (!requireNamespace("shiny"))
#'   install.packages("shiny")
#' if (!requireNamespace("leaflet"))
#'   install.packages("leaflet")
#' library(magick)
#' library(shiny)
#' library(leaflet)
#'
#' # Load H&E image with magick
#' imfile <-
#'   system.file("extdata/mousecolon/spatial",
#'               "tissue_hires_image.png",
#'               package = "STUtility2")
#' im <- image_read(imfile)
#'
#' # tile image and return path to tiles
#' tile_res <- TileImage(im, nCores = 2)
#'
#' # Create a simple viewer with leaflet
#' ui <- fluidPage(
#'   leafletOutput("map", height = 512, width = 512),
#' )
#'
#' server <- function(input, output, session) {
#'   addResourcePath("mytiles", "/Users/ludviglarsson/Downloads/tiles")
#'   output$map <- renderLeaflet({
#'     leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
#'       addTiles(urlTemplate = "/mytiles/{z}/{x}_{y}.jpg",
#'                options = tileOptions(continuousWorld = TRUE,
#'                                      tileSize = "256",
#'                                      minZoom = paste0(res$minZoomLevel),
#'                                      maxZoom = paste0(res$maxZoomLevel)))
#'   })
#' }
#'
#' shinyApp(ui, server)
#' }
#'
#' @export
TileImage <- function (
    im,
    sampleID = 1,
    outpath = NULL,
    maxZoomLevel = 4,
    maxImgWidth = 1e4,
    nCores = detectCores() - 1,
    verbose = TRUE
) {

  # Validate input
  stopifnot(
    inherits(im, what = "magick-image") & length(im) == 1,
    inherits(sampleID, what = c("numeric", "integer")) & length(sampleID) == 1,
    inherits(maxZoomLevel, what = c("numeric", "integer")) & length(maxZoomLevel) == 1,
    inherits(maxImgWidth, what = c("numeric", "integer")) & length(maxImgWidth) == 1
  )

  # Define max Zoom level
  info <- image_info(im)
  if (verbose) cli_alert_info("Got an image of dimensions {info$height}x{info$width} for sample {sampleID}")
  nLevels <- floor(max(log2(ceiling(info$width/256)), log2(ceiling(info$height/256)))):0
  scaled_dims <- tibble(scaled_widths = info$width/(2^nLevels),
                        scaled_heigths = info$width/(2^nLevels))
  lower_limit <- apply(scaled_dims, 1, max) > 256
  upper_limit <- apply(scaled_dims, 1, max) < maxImgWidth
  nLevels <- nLevels[lower_limit & upper_limit]
  if (verbose) cli_alert("  Tiling image into {length(nLevels)} zoom levels")

  # Create output path
  outpath <- outpath %||% tempdir()
  outpath_data <- paste0(outpath, "/osd_data")
  dir.create(outpath_data, showWarnings = FALSE)
  outpath_tiles <- paste0(outpath_data, paste0("/tiles", sampleID))
  dir.create(outpath_tiles, showWarnings = FALSE)

  # Zoom levels
  tiles <- list()
  if (verbose) cli_alert("  Creating tiles")
  for (n in seq_along(nLevels)) {

    dir.create(paste0(outpath_tiles, "/", n), showWarnings = FALSE)

    if (nLevels[n] > 0) {
      im_scaled <- image_scale(im, geometry = paste0(info$width/(2^nLevels[n])))
    } else {
      im_scaled <- im
    }
    info_scaled <- image_info(im_scaled)

    # Create crops
    nRows <- floor(info_scaled$height/256)
    nRowsMar <- info_scaled$height %% 256
    nCols <- floor(info_scaled$width/256)
    nColsMar <- info_scaled$width %% 256

    # Create crop windows
    if (nRowsMar > 0) {
      rowHeights <- c(rep(256, nRows), nRowsMar)
    } else {
      rowHeights <- rep(256, nRows)
    }
    rowOffsets <- c(0, cumsum(rowHeights)[1:length(rowHeights) - 1])
    if (nColsMar > 0) {
      colWidths <- c(rep(256, nCols), nColsMar)
    } else {
      colWidths <- rep(256, nCols)
    }
    colOffsets <- c(0, cumsum(colWidths)[1:(length(colWidths) - 1)])

    # Create tiles
    for (i in seq_along(colWidths)) {
      for (j in seq_along(rowHeights)) {
        crop.geom <- paste0(colWidths[i], "x", rowHeights[j], "+", colOffsets[i], "+", rowOffsets[j])
        im_cropped <- image_crop(im_scaled, geometry = crop.geom)
        if (image_info(im_cropped)$width < 256 | image_info(im_cropped)$height < 256) {
          im_blank <- image_blank(color = "white", width = 256, height = 256)
          im_cropped <- image_composite(im_blank, composite_image = im_cropped)
          im_cropped <- image_transparent(im_cropped, color = "white")
        }
        tiles <- c(tiles, setNames(list(im_cropped), nm = paste0(n, "/", i - 1, "_", j - 1)))
      }
    }
  }

  # Export tiles
  if (verbose) cli_alert("  Exporting tiles")
  results <- mclapply(names(tiles), function(tileName) {
    image_write(tiles[[tileName]], path = paste0(outpath_tiles, "/", tileName, ".jpg"))
  }, mc.cores = nCores)

  # Export data as JSON
  if (verbose) cli_alert("  Exporting meta data")
  d <- list(tilepath = outpath_tiles,
       minZoomLevel = ifelse(length(nLevels) > 1,
                             nLevels[length(nLevels)] + 1, 1),
       maxZoomLevel = nLevels[1] + 1,
       image_width = info$width,
       image_height = info$height,
       tilesize = 256)
  image_info_outpath <- paste0(outpath_data, paste0("/image_info_", sampleID, ".json"))
  jsonlite::write_json(x = d, path = image_info_outpath, auto_unbox = TRUE)

  return(list(datapath = outpath_data, tilepath = outpath_tiles, infopath = image_info_outpath))
}
