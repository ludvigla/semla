#' @include generics.R
#' @include checks.R
#'
NULL


#' @param spots A character vector with spot IDs present 'object'
#' @param verbose Print messages
#'
#' @import dplyr
#' @importFrom rlang abort inform
#' @importFrom glue glue
#'
#' @rdname disconnect-regions
#'
#' @examples
#'
#' library(STUtility2)
#' library(dplyr)
#' library(ggplot2)
#' library(patchwork)
#'
#' galt_spots_file <- "~/STUtility2/repo/STUtility2/inst/extdata/mousecolon/galt_spots.csv"
#' galt_spots <- read.csv(galt_spots_file) |>
#'   as_tibble()
#'
#' # read coordinates
#' coordfile <- system.file("extdata/mousecolon/spatial",
#'                          "tissue_positions_list.csv",
#'                          package = "STUtility2")
#' coords <- read.csv(coordfile, header = FALSE) |>
#'   filter(V2 == 1) |>
#'   select(V1, V6, V5) |>
#'   setNames(nm = c("barcode", "x", "y")) |>
#'   bind_cols(sampleID = 1) |>
#'   as_tibble()
#'
#' # Select spots
#' spots <- galt_spots$barcode[galt_spots$selection == "GALT"]
#' head(spots)
#'
#' # Find disconnected regions in GALT spots
#' disconnected_spot_labels <- DisconnectRegions(coords, spots)
#'
#' # Add information to coords and plot
#' gg <- coords |>
#'   mutate(galt = NA, galt_disconnected = NA)
#' gg$galt[match(spots, gg$barcode)] <- "galt"
#' gg$galt_disconnected[match(names(disconnected_spot_labels), gg$barcode)] <- disconnected_spot_labels
#'
#' p1 <- ggplot(gg, aes(x, y, color = galt))
#' p2 <- ggplot(gg, aes(x, y, color = galt_disconnected))
#' p <- p1 + p2 &
#'   geom_point() &
#'   theme_void()
#' p
#'
#' @export
#'
DisconnectRegions.default <- function (
  object,
  spots,
  verbose = TRUE,
  ...
) {

  # Set global variables to NULL
  barcode <- x <- y <- sampleID <- NULL

  # Check object object class
  if (!any(class(object) %in% c("data.frame", "matrix", "tbl")))
    abort(glue("Invalid class '{class(object)}'."))
  if (ncol(object) != 4)
    abort(glue("Invalid number of columns '{ncol(object)}'. Expected 4."))
  if (!all(
    object |> summarize(
      check_barcode = is.character(barcode),
      check_x = is.numeric(x),
      check_y = is.numeric(y),
      check_sample = is.numeric(sampleID)
    ) |>
    unlist()
  )) {
    abort(glue("Invalid column class(es)."))
  }

  # Check that spots are in object
  stopifnot(
    inherits(spots, what = "character"),
    length(spots) > 0,
    all(spots %in% object$barcode)
  )
  if (verbose) inform(c("i" = glue("Detecting disconnected regions for {length(spots)} spots")))

  # Find border
  spatnet <- GetSpatialNetwork(object, ...)

  # Subset spatnet by selected spots
  spatnet <- lapply(spatnet, function(x) {
    x |>
      filter(from %in% spots, to %in% spots)
  })

  # load tidygraph
  if (!requireNamespace("tidygraph"))
    install.packages("tidygraph")

  # Convert spatnet to tidygraph
  tidygraphs <- do.call(bind_rows, lapply(names(spatnet), function(nm) {
    x <- spatnet[[nm]]
    diconnected_graphs <- as_tbl_graph(x) |>
      to_components() |>
      lapply(as_tibble)
    sizes <- order(sapply(diconnected_graphs, nrow), decreasing = TRUE)
    diconnected_graphs <- diconnected_graphs[sizes]
    diconnected_graphs <- do.call(bind_rows, lapply(seq_along(diconnected_graphs), function(i) {
      diconnected_graphs[[i]] |>
        select(-contains("id")) |>
        bind_cols(id = paste0("S", nm, "-region", i))
    }))
    return(diconnected_graphs)
  }))
  if (verbose) inform(c("i" = glue("Found {length(unique(tidygraphs$id))} disconnected graphs in data")))
  if (verbose) inform(c("i" = glue("Sorting disconnected regions by decreasing size")))

  # return results as a character vector
  labeler <- setNames(tidygraphs$id, nm = tidygraphs$name)
  singleton_spots <- setdiff(spots, names(labeler))
  if (verbose) inform(c("i" = glue("Found {length(spots) - length(labeler)} singletons in data")))
  if (verbose & (length(singleton_spots) > 0)) inform(c(">" = "   These will be labeled as 'singletons'"))
  labeler <- c(labeler, setNames(rep("singleton", length(singleton_spots)), nm = singleton_spots))
  labels <- labeler[spots]

  return(labels)

}

