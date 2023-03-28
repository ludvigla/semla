#' Adjust tissue coordinates for digital unrolling
#'
#' This function takes a \code{tbl_graph} object as input generated
#' with \code{\link{CutSpatialNetwork}} and attempts to calculate
#' the unrolled tissue coordinates.
#'
#' @details
#' The algorithm is described briefly below:
#'
#' First, the end point of the graph are identified. The outermost
#' end point is assumed to be the starting point and the point at the
#' center of the roll is assumed to be the end point.
#'
#' Then, the algorithm tries to find the geodesic between the
#' two end points, forming a band of nodes which represents the shortest
#' path between the starting point and the end point.
#'
#' The order of these shortest represents the distances along the x axis
#' in the new coordinate system and the geodesics to all other nodes are used
#' as distances on the y axis of the in the new coordinate system.
#'
#' By using information about the location of nodes relative to the position
#' of the shortest path nodes, the algorithm also tries to adjust the y axis
#' to ensure non-negative values.
#'
#' NB: The graph has to be fully connected. If multiple sub graphs are found,
#' only the largest graph will be kept. It is crucial that edges has been cut
#' properly with \code{\link{CutSpatialNetwork}}, otherwise the results will be
#' inaccurate. See the package website for examples.
#'
#' @param full_graph A \code{tbl_graph} object generated with \code{\link{CutSpatialNetwork}}
#' @param verbose Print messages
#'
#' @import rlang
#' @import glue
#' @import dplyr
#' @importFrom tibble tibble as_tibble
#' @import cli
#'
#' @return A \code{tibble} with "unrolled" tissue coordinates
#'
#' @export
#'
AdjustTissueCoordinates <- function (
  full_graph,
  verbose = TRUE
) {

  # Set global variables to NULL
  edges <- keep <- x <- y <- name <- x_dist <- y_dist <-NULL

  # Import tidygraph and igraph
  if (!requireNamespace("tidygraph"))
    abort(glue("Package {cli::col_br_magenta('tidygraph')} is required. Please install it with: \n",
               "install.packages('tidygraph')"))
  if (!requireNamespace("igraph"))
    abort(glue("Package {cli::col_br_magenta('igraph')} is required. Please install it with: \n",
               "install.packages('igraph')"))
  if (!inherits(full_graph, what = "tbl_graph"))
    abort(glue("Invalid class '{class(full_graph)[1]}', expected a 'tbl_graph'"))

  # Check if edges has been removed
  edges_to_be_removed <- sum(!(full_graph |> tidygraph::activate(edges) |> pull(keep)))
  if (edges_to_be_removed > 0) {

    if (verbose) cli_alert_info("Removing {edges_to_be_removed} edges")

    # Filter graph
    small_graph <- full_graph |>
      tidygraph::activate(edges) |>
      filter(keep)
  } else {
    small_graph <- full_graph
  }

  # Check if multiple subgraphs are present
  if (verbose) cli_alert_info("Checking for disconnected graphs")
  small_graph_split <- tidygraph::to_components(small_graph)
  if (length(small_graph_split) > 1) {
    if (verbose) cli_alert_info(col_br_magenta("More than 1 subgraph identified. Keeping the largest subgraph"))
    small_graph <- small_graph_split[[1]]
  }

  # calculate pairwise geodesics between nodes
  if (verbose) cli_alert_info("Calculating pairwise geodesics between nodes in graph")
  distMat <- igraph::distances(small_graph)

  # identify end nodes
  if (verbose) cli_alert_info("Identifying end points")
  inds <- which(distMat == max(distMat), arr.ind = TRUE)
  if (nrow(inds) > 2) {
    inds <- inds[1:2, ]
  }

  # Get nodes from small graph
  nodes <- small_graph |>
    tidygraph::activate(nodes) |>
    as_tibble()

  # order end points
  outer_end <- nodes[match(rownames(inds), nodes$name), ] |>
    mutate(x = x - 0.5, y = y - 0.5) |>
    select(x, y) |>
    as.matrix() |>
    apply(1, function(x) {
      sqrt(sum(x[1]^2 + x[2]^2))
    }) |>
    which.max()

  # obtain distances from end point from distMat
  dists_from_end <- tibble(name = colnames(distMat), dist = distMat[inds[outer_end, 1], ])

  # Add distances from end
  nodes <- nodes |>
    left_join(y = dists_from_end, by = "name")

  # Find shortest path between end points
  if (verbose) cli_alert_info("Finding shortest path beetween end points")
  spath <- igraph::shortest_paths(
    graph = small_graph,
    from = inds[outer_end, 1],
    to = inds[outer_end, 2],
    mode = "all"
  )$vpath[[1]]

  # Add shortest path order to nodes
  nodes$short_path <- ifelse(nodes$name %in% names(spath), "short_path", NA)
  nodes$short_path_ord <- NA
  nodes$short_path_ord[match(names(spath), nodes$name)] <- 1:length(spath)

  # Subset distance matrix to include shortest path nodes in rows
  # and all other nodes in columns
  dist_subset <- distMat[names(spath), !colnames(distMat) %in% names(spath)]

  # Find shortest distance for each node not in shortest path
  node_names <- colnames(dist_subset)
  spath_to_nodes <- do.call(bind_rows, lapply(1:ncol(dist_subset), function(i) {
    closest_spath_node = which.min(dist_subset[, i])
    tibble(node_name = node_names[i], spath_node_name = names(closest_spath_node))
  }))

  # Get coodinates for shortest path nodes and convert to matrices
  coords_spath_nodes <- nodes[match(rownames(dist_subset), nodes$name), c("name", "x", "y")] |>
    mutate(ord = 1:n()) |>
    data.frame(row.names = 1) |>
    as.matrix()
  coords_nodes <- nodes[match(colnames(dist_subset), nodes$name), c("name", "x", "y")] |>
    data.frame(row.names = 1) |>
    as.matrix()

  # Calculate sign for angles
  if (verbose) cli_alert_info("Checking location of nodes relative to shortest path nodes")
  nodes_sign <- .get_nodes_sign(distMat = dist_subset,
                                coords_nodes = coords_nodes,
                                coords_spath_nodes = coords_spath_nodes,
                                spath_to_nodes = spath_to_nodes)

  # Add distances
  nodes <- nodes |>
    left_join(y = nodes_sign, by = "name")

  # Add x distances
  nodes$x_dist <- NA
  nodes$x_dist[match(rownames(coords_spath_nodes), nodes$name)] <- coords_spath_nodes[, "ord"]
  nodes$x_dist[match(spath_to_nodes$node_name, nodes$name)] <- coords_spath_nodes[spath_to_nodes$spath_node_name, "ord"]

  # Select required columns
  nodes <- nodes |>
    select(name, x, y, x_dist, y_dist)

  # rescale y distances
  if (verbose) cli_alert_info("Rescaling y distances to ensure non-negative values")
  nodes <- nodes |>
    mutate(y_dist = case_when(is.na(y_dist) ~ 0,
                              TRUE ~ y_dist)) |>
    na.omit() |>
    group_by(x_dist) |>
    mutate(y_dist = y_dist - min(y_dist))

  # Return results
  if (verbose) cli_alert_success("Finished!")
  return(nodes)
}


#' Calculate sign for nodes relative to base nodes
#'
#' @importFrom tibble tibble
#'
#' @noRd
.get_nodes_sign <- function(
  distMat,
  coords_nodes,
  coords_spath_nodes,
  spath_to_nodes
) {

  nodes_signs <- do.call(bind_rows, lapply(colnames(distMat), function(s) {
    P <- coords_nodes[s, 1:2]
    to_node <- spath_to_nodes[match(s, spath_to_nodes$node_name), ]
    xs <- expand.range(coords_spath_nodes[to_node$spath_node_name, "ord"],
                       exp.factor = 2, maxval = nrow(coords_spath_nodes))
    xs <- coords_spath_nodes[xs, ]
    ls <- adjust.coords(x = xs[1, 1:2], y = xs[2, 1:2], P = P)
    x <- ls[1:2]; y <- ls[3:4]
    y <- sign_angle(x = x, y = y, P = as.numeric(P))
    y_dist <- distMat[to_node$spath_node_name, s]
    y_dist <- ifelse(y > 0, y_dist, -y_dist)
    return(tibble(name = s, sign_angle = ifelse(y > 0, 1, -1), y_dist))
  }))

  return(nodes_signs)
}


#' Finds angle for point relative to line
#'
#' @noRd
Angle <- function (
  A,
  B,
  C
) {
  vector1 = c(A[1] - B[1], A[2] - B[2])
  vector2 = c(C[1] - B[1], C[2] - B[2])
  num = (vector1[1] * vector2[1] + vector1[2] * vector2[2])
  den = sqrt(vector1[1]^2 + vector1[2]^2) * sqrt(vector2[1]^2 +
                                                   vector2[2]^2)
  angle = acos(num/den)
  angle = (360 * angle)/(2 * pi)
  return(angle)
}

#' Adjusts coordinates for point
#'
#' @noRd
adjust.coords <- function(x, y, P){
  avg <- c(mean(c(x[1], y[1])), mean(c(x[2], y[2])))
  df <- (P - avg)
  x_new <- x - df
  y_new <- y - df
  return(c(x_new, y_new))
}

#' Checks what side of a line a point is located
#'
#' @noRd
sider <- function(x, y, P) {
  (P[1] - x[1])*(y[2] - x[2]) - (P[2] - x[2])*(y[1] - x[1])
}

#' Finds angle for point relative to line 2
#'
#' @noRd
angle <- function(x, y, P){
  Angle(A = P, B = c(mean(c(x[1], y[1])), mean(c(x[2], y[2]))), C = y)
}

#' Get sign for angle
#'
#' @noRd
sign_angle <- function(x, y, P) {
  angle(x, y, P)*sign(sider(x, y, P))
}

#' Expand range
#'
#' @noRd
expand.range <- function(x, exp.factor = 5, maxval = 1e4) {
  y <- c(max(1, x - exp.factor), min(x + exp.factor, maxval))
  return(y)
}
