#' validate a selected column name
#'
#' @param object A Seurat object
#' @param column_name A character specifying a column name
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#'
#' @return nothing
#'
#' @noRd
.validate_column_name <- function (
    object,
    column_name
) {
  # Check if column_name is available in meta data
  if (!is.character(column_name)) abort(glue("Invalid class '{class(column_name)}', expected a 'character'"))
  if (length(column_name) != 1) abort(glue("Invalid length {length(column_name)} for 'column_name', expected a 'character' vector of length 1"))
  if (!column_name %in% colnames(object[[]])) abort("'column_name' is not present in the Seurat object meta data")
}


#' validate selected labels
#'
#' @param object A Seurat object
#' @param selected_groups A character vector with selected labels in \code{column_name}
#' @param column_name A character vector specifying the name of a column in \code{object}
#' meta data slot
#'
#' @import dplyr
#' @import cli
#' @importFrom glue glue
#' @importFrom rlang abort
#'
#' @return A character vector with selected labels in \code{column_name}
#'
#' @noRd
.validate_selected_labels <- function (
    object,
    selected_groups,
    column_name
) {
  selected_groups <- selected_groups %||% {
    cli_alert_info("No 'selected_groups' provided. Using all groups in '{column_name}' column")
    selected_groups <- unique(object[[]] |> pull(all_of(column_name)))
  }
  if (!inherits(selected_groups, what = c("factor", "character")))
    abort(glue("Invalid class '{class(selected_groups)}', expected a 'character'"))
  if (!all(selected_groups %in% (object[[]] |> pull(all_of(column_name)))))
    abort(glue("Some 'selected_groups' are not present in '{column_name}'"))

  return(selected_groups)
}


#' Get spot IDs for selected groups and samples
#'
#' @param split_by_sample A logical specifying if the spots should be split
#' by sample
#'
#' @inheritParams .validate_selected_labels
#'
#' @import dplyr
#'
#' @return A list of lists with spot IDs for selected groups split by sample
#'
#' @noRd
.get_spots_list <- function (
  object,
  selected_groups,
  column_name,
  split_by_sample = TRUE
) {

  # Set global variables to NULL
  sampleID <- barcode <- NULL

  spots_list <- lapply(selected_groups, function(lbl) {
    spots <- GetStaffli(object)@meta_data |>
      bind_cols(object[[]] |> select(all_of(column_name))) |>
      filter(!! sym(column_name) == lbl)
    if (split_by_sample) {
      spots <- spots |>
        group_by(sampleID) |>
        group_split() |>
        lapply(function(x) x |> pull(barcode)) |>
        setNames(unique(spots$sampleID))
    } else {
      spots <- spots |>
        pull(barcode)
    }
    return(spots)
  }) |>
    setNames(selected_groups)

  return(spots_list)

}


#' Get spot coordinates for selected groups and samples
#'
#' @inheritParams .validate_selected_labels
#'
#' @import dplyr
#'
#' @return A list of spot coordinates for each sample
#'
#' @noRd
.get_coords_list <- function (
    object
) {

  # Set global variables to NULL
  barcode <- pxl_col_in_fullres <- pxl_row_in_fullres <- sampleID <- NULL

  coords_list <- GetStaffli(object)@meta_data |>
    select(barcode, pxl_col_in_fullres, pxl_row_in_fullres, sampleID) |>
    setNames(nm = c("barcode", "x", "y", "sampleID")) |>
    group_by(sampleID) |>
    group_split() |>
    setNames(nm = paste0(1:nrow(GetStaffli(object)@image_info)))

  return(coords_list)

}
