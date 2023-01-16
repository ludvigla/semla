#' @param X List
#' @param FUN Function to apply to list elements
#' @param nCores Number of cores
#' @param ... Other parameters passed to parLapply
#'
#' @importFrom parallel makeCluster parLapply stopCluster
#'
#' @noRd
.winLapply <- function(X, FUN, nCores, ...) {
  cl <- makeCluster(nCores)
  res <- parLapply(cl, X, FUN, ...)
  stopCluster(cl)
  return(res)
}

#' @param X List
#' @param FUN Function to apply to list elements
#' @param nCores Number of cores
#' @param ... Other parameters passed to mclapply
#'
#' @importFrom parallel mclapply
#'
#' @noRd
.unixLapply <- function(X, FUN, nCores, ...) {
  res <- mclapply(X, FUN, mc.preschedule = TRUE, ...)
  return(res)
}
