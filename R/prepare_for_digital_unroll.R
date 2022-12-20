#' @include generics.R
#' @include checks.R
#'
NULL

#' Export a spatial graph to a JSON file
#'
#' Utility function to prepare data for \code{\link{CutSpatialNetwork}}.
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
#' @seealso CutSpatialNetwork
#'
#' @export
#'
export_graph <- function (
  object,
  sampleID = 1L,
  outdir,
  verbose = TRUE
) {

  # Set global variables to NULL
  from <- to <- barcode <- pxl_col_in_fullres <- pxl_row_in_fullres <- name <- NULL
  edges <- weight <- x <- y <- NULL

  # Check Seurat object
  .check_seurat_object(object)

  # Check input parameters
  stopifnot(
    inherits(sampleID, what = c("numeric", "integer")),
    length(sampleID) == 1,
    inherits(outdir, what = "character"),
    dir.exists(outdir)
  )

  # Import tidygraph
  if (!requireNamespace("tidygraph"))
    install.packages("tidygraph")

  # Create spatial network
  if (verbose) cli_alert_info("Creating spatial network for sample {sampleID}")
  spatnet <- GetSpatialNetwork(object)[[sampleID]]

  # Create an adjacency matrix
  if (verbose) cli_alert_info("Creating adjacency matrix from spatial network")
  wide_spatial_network <- pivot_wider(
    spatnet |>
      select(from, to) |>
      mutate(value = 1L) |>
      add_count(from) |>
      filter(n > 1) |>
      select(-n),
    names_from = "from",
    values_from = "value",
    values_fill = 0
  ) |>
    data.frame(row.names = 1, check.names = FALSE) |>
    as.matrix()

  # Sort columns of matrix
  if (verbose) cli_alert_info("Rearranging {nrow(wide_spatial_network)}x{ncol(wide_spatial_network)} adjacency matrix")
  wide_spatial_network <- wide_spatial_network[, rownames(wide_spatial_network)]
  if (!all(rownames(wide_spatial_network) == colnames(wide_spatial_network)))
    abort("Something went wrong...")

  # remove duplicated edges
  if (verbose) cli_alert_info("Removing duplicated edges")
  wide_spatial_network <- wide_spatial_network*upper.tri(wide_spatial_network)

  # Add attributes to nodes
  if (verbose) {
    cli_alert_info("Creating tidy graph object")
    cli_alert_info("Adding node attributes")
  }
  network <- suppressWarnings({tidygraph::as_tbl_graph(wide_spatial_network, directed = FALSE)})

  st_object <- GetStaffli(object)
  network <- network |>
    tidygraph::activate(nodes) |>
    left_join(y = st_object@meta_data |>
                filter(sampleID == sampleID) |>
                select(barcode, pxl_col_in_fullres, pxl_row_in_fullres) |>
                dplyr::rename(name = barcode), by = "name") |>
    mutate(x = scales::rescale(x = pxl_col_in_fullres, to = c(0, 1),
                               from = c(0, st_object@image_info$full_width[sampleID])),
           y = scales::rescale(x = pxl_row_in_fullres, to = c(0, 1),
                               from = c(0, st_object@image_info$full_width[sampleID])),
           id = name,
           `_row` = name)

  # Add attributes to edges
  if (verbose) cli_alert_info("Adding edge attributes")
  network <- network |>
    tidygraph::activate(edges) |>
    select(-weight) |>
    mutate(keep = TRUE) |>
    left_join(y = network |>
                tidygraph::activate(nodes) |>
                as_tibble() |>
                mutate(from = 1:n()) |>
                select(from, x, y),
              by = "from") |>
    left_join(y = network |>
                tidygraph::activate(nodes) |>
                as_tibble() |>
                mutate(to = 1:n()) |>
                select(to, x, y) |>
                dplyr::rename(x_end = x, y_end = y),
              by = "to") |>
    mutate(index = 1:n())

  # Extract nodes and edges
  nodes = network |>
    tidygraph::activate(nodes) |>
    as_tibble()
  links = network |>
    tidygraph::activate(edges) |>
    as_tibble()
  links$source <- nodes$name[links$from]
  links$target <- nodes$name[links$to]
  links <- links |> select(-from, -to)

  # Put nodes and edges into a list
  data <- list(nodes = nodes, links = links)

  # Export
  networkpath <- file.path(outdir, paste0("network_Visium_", sampleID, ".json"))
  if (verbose) cli_alert_info("Exporting spatial network to {networkpath}")
  data_json <- data |>
    write_json(auto_unbox = TRUE,
               path = networkpath)
  if (verbose) cli_alert_success("Finished!")
}

#' @noRd
export_tidygraph <- function (
  network,
  sampleID = 1,
  outdir
) {

  # Set global variables to NULL
  name <- edges <- from <- to <- NULL

  network <- network |>
    tidygraph::activate(nodes) |>
    mutate(id = name, `_row` = name)

  network <- network |>
    tidygraph::activate(edges) |>
    mutate(index = 1:n())

  # Extract nodes and edges
  nodes = network |>
    tidygraph::activate(nodes) |>
    as_tibble()
  links = network |>
    tidygraph::activate(edges) |>
    as_tibble()
  links$source <- nodes$name[links$from]
  links$target <- nodes$name[links$to]
  links <- links |> select(-from, -to)

  # Put nodes and edges into a list
  data <- list(nodes = nodes, links = links)

  # Export
  data_json <- data |>
    write_json(auto_unbox = TRUE,
               path = file.path(outdir, "network_Visium_", sampleID, ".json"))
}
