#' Create a scale bar to use for plots
#'
#' @param x Width of scale bar in microns. If the width is higher than 1,000 microns,
#' the units will be converted to millimeters in the scale bar title.
#' @param breaks Number of breaks to cut scale bar into. For example, 6 breaks will
#' create 6 vertical lines and 5 even intervals.
#' @param highlight_breaks Highlight specific breaks by increasing their height.
#' @param title_position One of "top" or "bottom" for title placement
#' @param flip_bar Should the scale bar be flipped to a vertical orientation?
#' @param text_height Height of scale bar title
#' @param ... Parameters passed to geom_segment
#'
#' @return A `ggplot` object with a scalebar
#'
#' @author Ludvig Larsson
#'
#' @import ggplot2
#' @importFrom tibble tibble
#'
#' @examples
#'
#' # Draw a scale bar for a 500 micron distance with 6 breaks where
#' # the ends are highlighted
#' scalebar()
#'
#' # Draw a scale bar for a 1mm mm distance with 1 breaks where the
#' # ends and the mid point are highlighted
#' scalebar(breaks = 11, highlight_breaks = c(1, 6, 11))
#'
#' @export
scalebar <- function (
  x = 500,
  breaks = 6,
  highlight_breaks = c(1, 6),
  title_position = c("top", "bottom"),
  flip_bar = FALSE,
  text_height = 2,
  ...
) {

  # load ggfittext
  if (!requireNamespace("ggfittext")) {
    install.packages("ggfittext")
  }

  # Set global variables to NULL
  ord <- min_y <- max_y <- NULL

  # Check input
  if (!inherits(x = x, what = c("numeric", "integer"))) {
    abort(glue("Invalid class '{class(x)}' for {col_br_magenta('x')}. Expected a numeric or integer of length 1"))
  }
  if (!between(x = x, left = 100, right = 1.1e4)) {
    abort(glue("Invalid value for {col_br_magenta('x')}. Expected a numeric or integer of length 1"))
  }
  if (!inherits(breaks, what = c("numeric", "integer"))) {
    abort(glue("Invalid class '{class(breaks)}' for {col_br_magenta('breaks')}. Expected a numeric or integer of length 1"))
  }
  if (length(breaks) > 1) {
    abort(glue("Invalid length '{length(breaks)}' for {col_br_magenta('breaks')}. Expected a numeric or integer of length 1"))
  }
  stopifnot(inherits(highlight_breaks, what = c("numeric", "integer")),
            inherits(title_position, what = "character"),
            inherits(flip_bar, what = "logical"),
            length(flip_bar) == 1)

  # Match title position args
  title_position <- match.arg(title_position, choices = c("top", "bottom"))

  # Create title
  title = ifelse(x >= 1e3, paste0(x/1e3, "mm"), paste0(x, "\u00B5m"))

  sb_breaks <- tibble(min_y = rep(-0.5, breaks),
                      max_y = rep(0.5, breaks)) |>
    mutate(ord = 1:n()) |>
    mutate(min_y = case_when(ord %in% highlight_breaks ~ min_y*1.5,
                             TRUE ~ min_y),
           max_y = case_when(ord %in% highlight_breaks ~ max_y*1.5,
                             TRUE ~ max_y))

  p <- ggplot() +
    geom_segment(aes(x = 1, xend = max(sb_breaks$ord), y = 0, yend = 0), ...) +
    geom_segment(data = sb_breaks, aes(x = ord, xend = ord, y = min_y, yend = max_y), ...) +
    theme_void()
  if (title_position == "top") {
    p <- p +
      ggfittext::geom_fit_text(aes(xmin = 1, xmax = breaks, ymin = max(sb_breaks$max_y), ymax = max(sb_breaks$max_y) + text_height, label = title),
                    grow = TRUE, angle = ifelse(flip_bar, 90, 0), min.size = 0.1, padding.x = grid::unit(0, "mm"), ...) +
      scale_y_continuous(limits = c(min(sb_breaks$min_y), max(sb_breaks$max_y) + text_height))
  } else {
    p <- p +
      ggfittext::geom_fit_text(aes(xmin = 1, xmax = breaks, ymin = min(sb_breaks$min_y) - text_height, ymax = min(sb_breaks$min_y), label = title),
                    grow = TRUE, angle = ifelse(flip_bar, 90, 0), min.size = 0.1, padding.x = grid::unit(0, "mm"), ...) +
      scale_y_continuous(limits = c(min(sb_breaks$min_y) - text_height, max(sb_breaks$max_y)))
  }

  # Flip bar?
  if (flip_bar) {
    p <- p + coord_flip()
  }

  p$labels$scalebar_width <- x

  return(p)
}
