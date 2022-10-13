#' @include checks.R
#'
NULL

#' Plot H&E images
#'
#' @param object a Seurat object
#' @param type image type, for example 'raw'
#' @param section_numbers an integer vector with section numbers to plot
#' @param ncol integer value specifying the number of columns in the plot
#' grid
#' @param mar margins around each plot. See \code{\link{par}} for details.
#'
#' @importFrom rlang abort
#' @importFrom graphics layout
#' @importFrom dplyr select bind_cols
#'
#' @return draws a plot of the H&E images
#'
#' @examples
#'
#' library(STUtility2)
#'
#' # Load example Visium data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "STUtility2"))
#' se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "STUtility2"))
#' se_merged <- MergeSTData(se_mbrain, se_mcolon)
#'
#' # ImagePlot will throw an error if no images are loaded
#' ImagePlot(se_merged)
#'
#' # Load images
#' se_merged <- LoadImages(se_merged)
#' ImagePlot(se_merged)
#'
#' # Plot only selected tissue sections
#' ImagePlot(se_merged, section_numbers = 1)
#'
#' # Change order of plot
#' ImagePlot(se_merged, section_numbers = 2:1)
#'
#' # Add a sample_id column and label plots
#' se_merged$sample_id <- ifelse(GetStaffli(se_merged)@meta_data$sampleID == 1, "brain", "colon")
#' ImagePlot(se_merged, label_by = "sample_id")
#'
#' @export
#'
ImagePlot <- function (
    object,
    label_by = NULL,
    type = "raw",
    section_numbers = NULL,
    ncol = NULL,
    mar = c(1, 1, 1, 1)
) {
  # obtain Staffli object
  .check_seurat_object(object)
  st_object <- GetStaffli(object)

  # Check label_by
  if (!is.null(label_by)) {
    if (!is.character(label_by)) abort("'label_by' should be a string specifying a meta data column in the Seurat object.")
    if (!label_by %in% colnames(object[[]])) abort(glue("'{label_by}' is not present in the Seurat object meta data."))
    label_by_vec <- object[[]] |>
      pull(all_of(label_by))
    if (!class(label_by_vec) %in% c("character", "factor")) abort(glue("Invalid class '{class(label_by_vec)}' for 'label_by' column. ",
                                                                       "Expects a 'character' of 'factor'."))
    label_by_vec  <- sapply(split(label_by_vec, st_object@meta_data$sampleID), function(x) {
      y <- unique(x)
      if (length(y) > 1) abort("Invalid 'label_by' meta data column. ",
                               "Cannot have multiple labels per tissue section.")
      return(y)
    })
  }

  # Check if images are available
  if (is.null(st_object@rasterlists[["raw"]])) abort("Images have not been loaded into the Seurat object. ",
                                                     "Did you run LoadImages() yet?")

  # fetch images
  images <- st_object@rasterlists[[type]]

  # Use all images if indices are not specified
  section_numbers <- section_numbers %||% {
    seq_along(images)
  }

  # Check if indices are OK
  if (any(!section_numbers %in% seq_along(images))) abort(glue("'section_numbers' is out of range. Section numbers ",
                                                               "{paste(seq_along(images), collapse = ', ')} are available."))

  # Select sections
  images <- images[section_numbers]
  if (!is.null(label_by)) {
    label_by_vec <- label_by_vec[section_numbers]
  }

  # Add sample ID
  ncols <- ncol %||% ceiling(sqrt(length(x = images)))
  nrows <- ceiling(length(x = images)/ncols)

   #Create a plot layout
  layout.matrix <- t(matrix(c(1:length(images), rep(0, nrows*ncols - length(images))), nrow = ncols, ncol = nrows))
  layout(mat = layout.matrix)

  # plot images
  for (i in seq_along(images)) {
    rst <- images[[i]]
    par(mar = mar)
    plot(rst)
    if (!is.null(label_by)) {
      title(main = label_by_vec[i])
    }
  }
}
