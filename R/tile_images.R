#' Tile an H&E image
#'
#' This function takes an image of class `magick-image` and create a tile
#' map. The size of each tile is 256x256 pixels and the number of zoom levels
#' are determined rom the
#'
#' @param im An image of class `magick-image`
#' @param outpath A string specifying an output directory to save the tiled image in.
#' If this is not provided, a temporary directory will be created
#' @param maxZoomLevel Max zoom level
#' @param maxImgWidth Safety threshold to make sure that the zoom level doesn't get too deep.
#' @param nCores Number of cores to use for threading
#'
#' @importFrom parallel detectCores mclapply
#' @importFrom magick image_read image_info image_crop image_scale image_blank image_composite image_transparent image_write
#' @importFrom tibble tibble
#' @importFrom rlang %||%
#'
#' @return Path to tiles
#'
#' @examples
#' \donttest{
#' library(magick)
#' library(shiny)
#' library(leaflet)
#'
#' # Load H&E image with magick
#' imfile <- system.file("extdata/mousecolon/spatial", "tissue_hires_image.png", package = "STUtility2")
#' im <- image_read(imfile)
#'
#' # tile image and return path to tiles
#' tile_res <- TileImage(im)
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
    outpath = NULL,
    maxZoomLevel = 4,
    maxImgWidth = 1e4,
    nCores = detectCores() - 1
) {

  # Define max Zoom level
  info <- image_info(im)
  nLevels <- floor(max(log2(ceiling(info$width/256)), log2(ceiling(info$height/256)))):0
  scaled_dims <- tibble(scaled_widths = info$width/(2^nLevels),
                        scaled_heigths = info$width/(2^nLevels))
  lower_limit <- apply(scaled_dims, 1, max) > 256
  upper_limit <- apply(scaled_dims, 1, max) < maxImgWidth
  nLevels <- nLevels[lower_limit & upper_limit]

  # Create output path
  outpath <- outpath %||% tempdir()
  outpath <- paste0(outpath, "/tiles")
  dir.create(outpath)

  # Zoom levels
  tiles <- list()
  for (n in seq_along(nLevels)) {

    dir.create(paste0(outpath, "/", n))

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
  results <- mclapply(names(tiles), function(tileName) {
    image_write(tiles[[tileName]], path = paste0(outpath, "/", tileName, ".jpg"))
  }, mc.cores = nCores)

  return(list(tilepath = outpath, minZoomLevel = ifelse(length(nLevels) > 1, nLevels[length(nLevels)] + 1, 1), maxZoomLevel = nLevels[1] + 1))
}
