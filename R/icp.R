#' Kabsch Algorithm
#'
#' Aligns two sets of points via rotations and translations.
#'
#' @details
#' Given two sets of points, with one specified as the reference set,
#' the other set will be rotated so that the RMSD between the two is minimized.
#' The format of the matrix is that there should be one row for each of
#' n observations, and the number of columns, d, specifies the dimensionality
#' of the points. The point sets must be of equal size and with the same
#' ordering, i.e. point one of the second matrix is mapped to point one of
#' the reference matrix, point two of the second matrix is mapped to point two
#' of the reference matrix, and so on.
#'
#' @section Author:
#' The original code was written by James Melville. See references for link to
#' the GitHub gist where the code was taken from.
#'
#' @param pm n x 2 matrix of points to align to to \code{qm}.
#' @param qm n x 2 matrix of reference points.
#' @return Matrix \code{pm} rotated and translated so that the ith point
#' is aligned to the ith point of \code{qm} in the least-squares sense.
#'
#' @references
#' \url{https://gist.github.com/jlmelville/9b4e5d076e719a7541881e8cbf58a895}
#' \url{https://en.wikipedia.org/wiki/Kabsch_algorithm}
#'
kabsch <- function(pm, qm) {
  pm_dims <- dim(pm)
  if (!all(dim(qm) == pm_dims)) {
    stop(call. = TRUE, "Point sets must have the same dimensions")
  }
  # The rotation matrix will have (ncol - 1) leading ones in the diagonal
  diag_ones <- rep(1, pm_dims[2] - 1)

  # center the points
  pm <- scale(pm, center = TRUE, scale = FALSE)
  qm <- scale(qm, center = TRUE, scale = FALSE)

  am <- crossprod(pm, qm)

  svd_res <- svd(am)
  # use the sign of the determinant to ensure a right-hand coordinate system
  d <- determinant(tcrossprod(svd_res$v, svd_res$u))$sign
  dm <- diag(c(diag_ones, d))

  # rotation matrix
  um <- svd_res$v %*% tcrossprod(dm, svd_res$u)

  return(list(um, qm, pm))
}


#' Iterative Closest Point algorithm ICP
#'
#' Aligns two sets of unpaired point sets by applying rotations and translations.
#' The point sets can be of unequal length.
#'
#' @param xy_ref m x 2 numeric matrix of reference points
#' @param xy_query n x 2 numeric matrix of query points
#' @param iterations Number of iterations to run before stopping the ICP
#'
#' @importFrom dbscan kNN
#' @importFrom zeallot %<-%
#'
#' @return A list with the following objects:
#'
#' \itemize{
#'    \item{y_transf: n x 2 matrix of aligned query points}
#'    \item{rot_mat: 2 x 2 rotation matrix}
#' }
#'
#' @examples
#'
#' library(STUtility2)
#' library(dplyr)
#' library(ggplot2)
#'
#' # Load example mouse brain data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain",
#'                                  "se_mbrain",
#'                                  package = "STUtility2"))
#'
#' # Get spatial network
#' spatnet <- GetSpatialNetwork(se_mbrain)
#'
#' # Keep tissue border points
#' n1 <- spatnet[[1]] |>
#'   filter(nn < 6)
#' ggplot(n1, aes(x_start, y_start)) +
#'   geom_point() +
#'   scale_y_reverse()
#'
#' # get spot coordinates points
#' xy <- GetStaffli(se_mbrain)@meta_data |>
#'   filter(barcode %in% unique(n1$from)) |>
#'   select(pxl_col_in_fullres, pxl_row_in_fullres) |>
#'   setNames(nm = c("x", "y")) |>
#'   bind_cols(type = "set1")
#' xy_diff <- CoordTransform(xy_coords = xy[, 1:2],
#'                           angle = 30,
#'                           xy_offset = c(500, 500)) |>
#'   slice_sample(n = nrow(xy) - 100) |>
#'   setNames(nm = c("x", "y")) |>
#'   bind_cols(type = "set2")
#' xy_orig <- bind_rows(xy, xy_diff)
#'
#' # Plot point sets
#' ggplot(xy_orig, aes(x, y, color = type)) +
#'   geom_point() +
#'   scale_y_reverse() +
#'   labs(title = "Original point sets")
#'
#' res <- icp(xy_ref = xy[, 1:2], xy_query = xy_diff[, 1:2])
#'
#' xy_aligned <- xy |>
#'   bind_rows(bind_cols(setNames(res$y_transf|> as_tibble(),
#'                                nm = c("x", "y")),
#'                       type = "set2_aligned"))
#'
#' # Plot aligned point sets
#' ggplot(xy_aligned, aes(x, y, color = type)) +
#'   geom_point() +
#'   scale_y_reverse() +
#'   labs(title = "Aligned point sets")
#'
#' # We can obtain the rotation angle in degrees from
#' # the results
#' atan2(res$rot_mat[2, 1], res$rot_mat[1, 1])*(180/pi)
#'
#' @export
#'
icp <- function (
  xy_ref,
  xy_query,
  iterations = 100
) {

  # Set global variables to NULL
  um <- qm <- pm <- NULL

  # Check input objects
  stopifnot(inherits(xy_ref, what = c("matrix", "data.frame")),
            inherits(xy_query, what = c("matrix", "data.frame")),
            ncol(xy_ref) == 2,
            ncol(xy_query) == 2)

  # Save transformations in tr_matrices
  tr_matrices <- list()
  for (i in 1:iterations) {
    # Find closest neighbors in xy_query
    closest_nn <- kNN(x = xy_ref, query = xy_query, k = 1)
    # Run kabsch algorithm of selecetd points
    c(um, qm, pm) %<-% kabsch(pm = xy_query, qm = xy_ref[closest_nn$id, ])
    # Save 2x2 transformation matrix for current iteration
    tr_matrices <- c(tr_matrices, list(um))
    # Overwrite xy_query with transformed points
    xy_query <- sweep(t(tcrossprod(um, pm)), 2, -attr(qm, "scaled:center"))
  }

  # Get complete 2x2 transformation matrix
  rot_mat <- Reduce("%*%", tr_matrices)
  return(list(y_transf = as.matrix(xy_query), rot_mat = rot_mat))
}




# library(shiny)
# library(ggplot2)
#
# se_mbrain <- readRDS(system.file("extdata/mousebrain",
#                                  "se_mbrain",
#                                  package = "STUtility2"))
#
# se_mbrain <- LoadImages(se_mbrain)
# im <- se_mbrain@tools$Staffli@rasterlists$raw[[1]]
#
# ui <- pageWithSidebar(
#   headerPanel("Example"),
#   sidebarPanel(
#     radioButtons("color", "Pick Color", c("Pink", "Green", "Blue")),
#     selectInput("shape", "Select Shape:", c("Circle", "Triangle"))
#   ),
#   mainPanel(
#     fluidRow(column(width = 6,
#                     h4("Click plot to add points"),
#                     actionButton("rem_point", "Remove Last Point"),
#                     plotOutput("plot1", click = "plot_click")),
#              column(width = 6,
#                     h4("Table of points on plot"),
#                     tableOutput("table")))
#   )
# )
#
# server = function(input, output){
#
#   ## 1. set up reactive dataframe ##
#   values <- reactiveValues()
#   values$DT <- data.frame(x = numeric(),
#                           y = numeric(),
#                           color = factor(),
#                           shape = factor())
#
#   ## 2. Create a plot ##
#   output$plot1 = renderPlot({
#     p <- ggplot(values$DT, aes(x = x, y = y)) +
#       geom_point(aes(color = color), size = 5) +
#       scale_x_continuous(limits = c(0, ncol(im))) +
#       scale_y_continuous(limits = c(0, nrow(im))) +
#       theme_void() +
#       # include so that colors don't change as more color/shape chosen
#       scale_color_discrete(drop = FALSE) +
#       scale_shape_discrete(drop = FALSE)
#     p <- p + inset_element(p = im, left = 0, bottom = 0, right = 1, top = 1, on_top = FALSE)
#     return(p)
#   })
#
#   ## 3. add new row to reactive dataframe upon clicking plot ##
#   observeEvent(input$plot_click, {
#     # each input is a factor so levels are consistent for plotting characteristics
#     add_row <- data.frame(x = input$plot_click$x,
#                           y = input$plot_click$y,
#                           color = factor(input$color, levels = c("Pink", "Green", "Blue")),
#                           shape = factor(input$shape, levels = c("Circle", "Triangle")))
#     # add row to the data.frame
#     values$DT <- rbind(values$DT, add_row)
#   })
#
#   ## 4. remove row on actionButton click ##
#   observeEvent(input$rem_point, {
#     rem_row <- values$DT[-nrow(values$DT), ]
#     values$DT <- rem_row
#   })
#
#   ## 5. render a table of the growing dataframe ##
#   output$table <- renderTable({
#     values$DT
#   })
# }
#
# shinyApp(ui, server)
#
