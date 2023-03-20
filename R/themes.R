#' Clean theme
#' 
#' Removes titles, subtitles and legends from a \code{patchwork} object
#' 
#' @returns A \code{theme} object
#' 
#' @examples
#' library(semla)
#' 
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", 
#'                                  "se_mbrain", 
#'                                  package = "semla"))
#'                                  
#' # With MapFeatures
#' MapFeatures(se_mbrain, c("nFeature_Spatial", "nCount_Spatial")) &
#'    ThemeClean()
#' 
#' @export
ThemeClean <- function (
) {
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        legend.position = "none")
}

#' Move legend to right side
#' 
#' Move the legends in a \code{patchwork} object to the right hand side of the plots
#' 
#' @returns A \code{theme} object
#' 
#' @examples 
#' library(semla)
#' 
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", 
#'                                  "se_mbrain", 
#'                                  package = "semla"))
#'                                  
#' # With MapFeatures
#' MapFeatures(se_mbrain, c("nFeature_Spatial", "nCount_Spatial")) &
#'    ThemeLegendRight()
#' 
#' @export
ThemeLegendRight <- function (
) {
  theme(legend.position = "right", legend.text = element_text(angle = 0, hjust = 1))
}


#' Modify patchwork titles
#' 
#' Takes a \code{patchwork} object as input produced with \code{\link{MapFeatures}}
#' or \code{\link{MapLabels}} and moves the legends to the right side.
#' 
#' @param p A \code{patchwork} object
#' @param titles A character vector matching the number of patches in \code{p}
#' 
#' @returns A \code{patchwork} object
#' 
#' @examples 
#' library(semla)
#' 
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", 
#'                                  "se_mbrain", 
#'                                  package = "semla"))
#'                                  
#' # With MapFeatures
#' p <- MapFeatures(se_mbrain, c("nFeature_Spatial", "nCount_Spatial"))
#' ModifyPatchworkTitles(p, titles = c("First", "Second"))
#' 
#' # With MapLabels
#' p <- MapLabels(se_mbrain, c("sample_id"))
#' ModifyPatchworkTitles(p, titles = "My title")
#' 
#' @export
ModifyPatchworkTitles <- function (
  p,
  titles
) {
  
  # Check p & titles
  if (!inherits(p, what = "patchwork")) 
    abort(glue("Expected a 'patcwork' object, got an object of class '{class(p)}'."))
  if (!inherits(titles, what = "character")) 
    abort(glue("Expected titles to be 'character' vector, got an object of class '{class(titles)}'."))
  
  l <- p$patches$plots |> length()
  if (l == 0) {
    if (!length(titles) == 1)
      abort(glue("Invalid length of 'titles'. Expected a character vector of length 1."))
    p <- p + ggtitle(label = titles)
  }
  err <- try({
    if (l > 0) {
      if (!length(titles) == (l + 1))
        abort(glue("Invalid length of 'titles'. Expected a character vector of length {l + 1}."))
      p$labels$title <- titles[length(titles)]
      p$patches$plots <- lapply(seq_along(p$patches$plots), function(i) {
        pl <- p$patches$plots[[i]]
        pl$labels$title <- titles[i]
        return(pl)
      })
    }
  }, silent = TRUE)
  if (inherits(err, what = "try-error")) {
    if (l > 0) {
      if (!length(titles) == l)
        abort(glue("Invalid length of 'titles'. Expected a character vector of length {l + 1}."))
      p$patches$plots <- lapply(seq_along(p$patches$plots), function(i) {
        pl <- p$patches$plots[[i]]
        pl$patches$plots[[1]]$labels$title <- titles[i]
        return(pl)
      })
    }
  }
  return(p)
}
