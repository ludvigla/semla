#' Function title, e.g. "Scatter plot"
#'
#' Function description, e.g. "This function takes two numeric vectors
#' of equal length and returns a scatter plot"
#'
#' @details More details about function
#'
#' @section Section 1:
#' Some more details with its own header
#'
#' @section Section 2:
#' Some more details with its own header
#'
#' @param df data.frame
#' @param colx,coly column names specifying numeric vectors to plot
#' @param col Color of points
#' @param ... Optional parameters passed to \code{geom_point}
#'
#' @return Decription of what the function returns, e.g. "A scatter plot"
#'
#' @importFrom ggplot2 ggplot geom_point
#'
#' @examples
#' MyPlot(df = mtcars, colx = "cyl", coly = "hp", col = "red")
#'
#' @export
MyPlot <- function (
    df,
    colx,
    coly,
    col,
    ...
) {

  # Check column names
  check <- c(colx, coly) %in% colnames(df)
  if (!all(check)) stop(sprintf("'%s' not found in data.frame", paste(c(colx, coly)[!check], collapse = "', '")))

  # Check class
  check <- sapply(df[, c(colx, coly)], is.numeric)
  if (!all(check)) stop(sprintf("'%s' not integer vectors", paste(c(colx, coly)[!check], collapse = "', '")))

  # Create plot
  ggplot(data = df, aes_string(colx, coly)) +
    geom_point(color = col, ...) +
    labs(title="Hello there")

}
