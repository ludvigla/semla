#' @include generics.R
#' @include checks.R
#'
NULL

#' @param object A matrix-like object with 10x Visium data
#' @param singlecell_matrix A matrix-like object with scRNA-seq data
#' @param groups A character vector of length \code{ncol(singlecell_matrix)}
#' with cell type labels. For `Seurat` objects, one can omit groups in which case
#' the current identity is used instead. Alternatively, a string can be provided
#' specifying a column to use from the `singlecell_object` meta.data slot.
#' @param nCells_per_group Each cell type is down sampled to this number before
#' computing cell type expression profiles. This is mainly used to speed up
#' computation. Note that if a cell type has fewer than \code{nCells_per_group}
#' cells, the cells will be sampled with replacement to obtain \code{nCells_per_group}
#' cells.
#' @param minCells_per_celltype Minimum number of cells allowed per cell type
#' @param min_prop Minimum proportion allowed. Any proportion values lower than
#' this threshold will be removed and the remaining proportions will be rescaled.
#' @param L1 lasso penalty
#' @param seed An integer to set seed for reproducibility
#' @param return_expression_profiles Logical specifying if the expression profile matrix
#' should also be returned
#' @param verbose Print messages
#'
#' @import rlang
#' @import cli
#' @import glue
#'
#' @rdname celltype-prediction
#'
#' @section default method:
#' Input `object` is a matrix-like object with 10x Visium data. The function returns
#' a matrix with estimated cell type proportions of dimensions nCellTypes x nSpots,
#' where nCellTypes is the number of cell types and nSpots is the number of spots.
#' If `return_expression_profiles=TRUE`, the returned object will be a list with estimated
#' proportions `prop` and the cell type expression profile matrix `W`.
#'
#' @export
RunNNLS.default <- function (
  object,
  singlecell_matrix,
  groups,
  nCells_per_group = 50L,
  minCells_per_celltype = 10L,
  min_prop = 0.01,
  L1 = 0.01,
  seed = 1337L,
  return_expression_profiles = FALSE,
  verbose = TRUE,
  ...
) {

  # Check input objects
  if (!inherits(object, what = c("dgRMatrix", "dgCMatrix", "matrix", "data.frame")))
    abort(glue("Invalid format of spatial expression matrix: '{class(object)}'"))
  if (any(dim(object) == 0))
    abort(glue("Invalid dimensions of spatial expression matrix: '{paste(dim(object), collapse = 'x')}'"))
  if (!inherits(singlecell_matrix, what = c("dgRMatrix", "dgCMatrix", "matrix", "data.frame")))
    abort(glue("Invalid format of single-cell expression matrix: '{class(singlecell_matrix)}'"))
  if (any(dim(singlecell_matrix) == 0))
    abort(glue("Invalid dimensions of single-cell expression matrix matrix: '{paste(dim(singlecell_matrix), collapse = 'x')}'"))


  # Check that groups are valid
  if (!inherits(groups, what = "character"))
    abort(glue("Invalid class '{class(groups)}' for groups, expected a 'character' vector"))
  if (length(groups) == 0)
    abort("'groups' is empty")
  if (length(groups) != ncol(singlecell_matrix))
    abort(glue("'groups' does not match the single-cell expression matrix. ",
               "Expected {ncol(singlecell_matrix)}, got {length(groups)}"))


  # Make sure that dev version of RcppML is installed
  if (!requireNamespace("RcppML"))
    install.packages("RcppML") # compile dev version
  if ((packageVersion("RcppML") |> as.character()) != "0.3.7")
    warn(c("The NNLS function might break if using a dev version of RcppML is used. ",
           "If RcppML::project(...) fails, try installing CRAN version 0.3.7 of RcppML."))

  # Prepare data
  if (verbose) cli_alert_info("Preparing data for NNLS")
  W <- .get_w_matrix(object = object,
                     singlecell_matrix = singlecell_matrix,
                     nCells_per_group = nCells_per_group,
                     minCells_per_celltype = minCells_per_celltype,
                     groups = groups,
                     seed = seed,
                     verbose = verbose)

  # Run NNLS
  if (verbose) cli_alert_info("Predicting cell type proportions with NNLS for {ncol(W)} cell types")
  proj_expr <- RcppML::project(W, object, L1 = L1, ...)

  # Convert predicted values to proportions
  prop <- apply(proj_expr, 2, function(x) {prop.table(x)})
  prop[prop < min_prop] <- 0
  prop <- apply(prop, 2, function(x) {prop.table(x)})
  rownames(prop) <- colnames(W)
  colnames(prop) <- colnames(object)

  if (return_expression_profiles) {
    return(list(prop = prop, W = W))
  } else {
    return(prop)
  }
}


#' @param features Features to use for computation. If no features are
#' provided, the intersection of variable features between sc_object and st_object
#' will be used. Note that the intersect can be quite small, in which case you can
#' increase the number of variable features by rerunning \code{FindVariableFeatures}
#' and increase the number with the `nfeatures` argument.
#' @param singlecell_object A `Seurat` object with single-cell gene expression data
#' @param singlecell_assay Assay in single cell data object to use for deconvolution
#' @param spatial_assay Assay in Visium data object to use for deconvolution
#' @param slot Name of slot in `singlecell_assay` and `spatial_assay` to use for
#' deconvolution
#' @param assay_name Sets the name of the returned `Assay` object. Only active
#' if return.as.dimred = FALSE.
#' @param dimred_name Sets the name of  the returned `DimReduc` object. Only
#' active if return.as.dimred = TRUE.
#' @param return_as_dimred By default, the results are returned as an `Assay`
#' named "celltypeprops", where each feature is named after the cell type IDs.
#' Alternatively, you can return the results as a `DimReduc` object. This will
#' also provide the cell gene loadings for each cell type. However, the proportions
#' will instead be named "factor_1", "factor_2", ...
#'
#' @section Seurat method:
#' Input `object` is a `Seurat` object with 10x Visium data. The function returns
#' the `Seurat` object with either a new `Assay` or `DimReduc` object containing
#' estimated cell type proportions.
#'
#' @importFrom Seurat GetAssayData CreateAssayObject CreateDimReducObject
#' @import cli
#'
#' @rdname celltype-prediction
#'
#' @return An object with cell type proportion estimates
#'
#' @export
#'
#' @md
RunNNLS.Seurat <- function (
    object,
    singlecell_object,
    groups = NULL,
    features = NULL,
    singlecell_assay = "RNA",
    spatial_assay = "Spatial",
    slot = "data",
    min_prop = 0.01,
    nCells_per_group = 50L,
    minCells_per_celltype = 10L,
    assay_name = "celltypeprops",
    dimred_name = "nnls",
    return_as_dimred = FALSE,
    L1 = 0.01,
    seed = 1337L,
    verbose = TRUE,
    ...
) {

  # Validate input objects
  stopifnot(inherits(object, what = "Seurat"),
            inherits(singlecell_object, what = "Seurat"))

  if (verbose) cli_h2("Predicting cell type proportions")
  if (verbose) cli_alert_info("Fetching data from Seurat objects")

  # Get features
  features <- .prep_features(object, singlecell_object, features, verbose)

  # Get groups
  groups <- .prep_groups(singlecell_object, groups, verbose)

  # Get expression matrix for single cell data and Visium data
  x_singlecell <- GetAssayData(singlecell_object, slot = slot, assay = singlecell_assay)[features, ]
  x_spatial <- GetAssayData(object, slot = slot, assay = spatial_assay)[features, ]
  # TODO: check if RNA data slot is normalized if slot = "data"

  results <- RunNNLS(
    object = x_spatial,
    singlecell_matrix = x_singlecell,
    groups = groups,
    nCells_per_group = nCells_per_group,
    minCells_per_celltype = minCells_per_celltype,
    min_prop = min_prop,
    L1 = L1,
    seed = seed,
    verbose = verbose,
    return_expression_profiles = TRUE,
    ... = ...
  )

  # Return results as a DimReduc object or an Assay object
  object <- .return_as(
    st_object = object,
    prop = results$prop,
    W = results$W,
    return_as_dimred = return_as_dimred,
    st_assay = spatial_assay,
    assay_name = assay_name,
    dimred_name = dimred_name,
    verbose = verbose
  )

  if (verbose) cli_alert_success("Finished")

  return(object)
}


#' Utility function to format returned data
#'
#' @param st_object A `Seurat` object woth Visium data
#' @param prop A matrix with cell type proportion estimates
#' @param W A matrix with cell type expression profiles
#' @param return_as_dimred Logical specifying if the data should be returned
#' as a `DimReduc` object or an `Assay` object
#' @param st_assay Assay used for spatial data
#' @param assay_name Assay name for returned data. Only used if `return_as_dimred=FALSE`
#' @param dimred_name Name for `DimReduc` object. Only used if `return_as_dimred=TRUE`
#' @param verbose Print messages
#'
#' @import glue
#' @import cli
#' @importFrom Seurat CreateDimReducObject CreateAssayObject DefaultAssay `DefaultAssay<-`
#'
#' @noRd
.return_as <- function (
    st_object,
    prop,
    W,
    return_as_dimred,
    st_assay,
    assay_name,
    dimred_name,
    verbose
) {
  # Return results as a DimReduc object or an Assay object
  if (return_as_dimred) {
    if (verbose) cli_alert_info("Returning results as a 'DimReduc' object")
    prop <- t(prop)
    colnames(prop) <- paste0("cFactor_", 1:ncol(prop))
    feature_loadings <- W
    colnames(feature_loadings) <- paste0("cFactor_", 1:ncol(prop))
    props_dimreduc <- CreateDimReducObject(embeddings = prop,
                                           loadings = feature_loadings,
                                           key = "cFactor_",
                                           assay = st_assay)
    st_object[[dimred_name]] <- props_dimreduc
  } else {
    if (verbose) cli_alert_info("Returning results in a new 'Assay' named '{assay_name}'")
    props_assay <- CreateAssayObject(data = as(prop, "dgCMatrix"))
    st_object[[assay_name]] <- props_assay
    if (verbose) cli_alert_info("Setting default assay to '{assay_name}'")
    DefaultAssay(st_object) <- assay_name
  }

  return(st_object)
}


#' Prepare features to keep for NNLS
#'
#' @import cli
#' @import rlang
#' @import glue
#' @importFrom Seurat VariableFeatures
#'
#' @noRd
.prep_features <- function (
  st_object,
  sc_object,
  features,
  verbose
) {
  if (is.null(features)) {
    if (length(VariableFeatures(sc_object)) == 0) abort("No variable features found in single cell Seurat object")
    if (length(VariableFeatures(st_object)) == 0) abort("No variable features found in spatial Seurat object")
    features <- features %||% intersect(VariableFeatures(sc_object), VariableFeatures(st_object))
    if (length(features) < 1e3) warn(glue("Only {length(features)} shared features detected"))
    if (length(features) < 100) abort(glue("At least 100 shared features required, found {length(features)}"))
  }

  # Filter features
  if (verbose) cli_alert("  Filtering out features that are only present in one data set")
  keep_features <- intersect(intersect(rownames(sc_object), rownames(st_object)), features)
  if (verbose) cli_alert("  Kept {length(keep_features)} features for deconvolution")

  return(keep_features)
}


#' @param sc_object A `Seurat` object with single-cell data
#' @param groups A character vector with group labels
#' @param verbose Print messages
#'
#' @import rlang
#' @import glue
#' @import cli
#' @import dplyr
#' @importFrom Seurat Idents
#'
#' @noRd
.prep_groups <- function (
  sc_object,
  groups,
  verbose
) {
  # Select groups
  if (is.null(groups)) {
    if (verbose) cli_alert_info("Using current identity as groups.")
    groups <- as.character(Idents(sc_object))
  } else if (is.character(groups) & length(groups) == 1) {
    if (!groups %in% colnames(sc_object[[]])) abort(glue("'{groups}' is not a valid column in Seurat object meta.data slot"))
    groups <- sc_object[[]] |>  pull(all_of(groups))
  } else {
    abort("Invalid value for 'groups'")
  }
  return(groups |> as.character())
}


#' Calculates expression profiles for groups of cells
#'
#' @param object An expression matrix with Visium data
#' @param singlecell_matrix An expression matrix with single-cell data
#' @param nCells_per_group Number of cells per cell type to keep
#' @param minCells_per_celltype Minimum number of cells accepted for a cell type
#' @param groups A character vector with group (cell type) labels
#' @param seed A seed for reproducibility
#' @param verbose Print messages
#'
#' @import rlang
#' @import glue
#' @import cli
#' @importFrom tibble tibble
#'
#' @noRd
.get_w_matrix <- function (
    object,
    singlecell_matrix,
    nCells_per_group,
    minCells_per_celltype,
    groups,
    seed,
    verbose
) {

  # Set global variables to NULL
  group <- barcode <- NULL

  # Rescale data to ensure positive values
  if (any(object < 0)) {
    if (verbose) warn("Found negative values in input matrix")
    if (verbose) cli_alert_info("Rescaling data to ensure positive values")
    object <- t(apply(object, 1, function(x) (x - min(x))/(max(x) - min(x))))
  }

  # Sample barcodes
  if (verbose) cli_alert(glue("  Downsampling scRNA-seq data to include a maximum of ",
                                   "{nCells_per_group} cells per cell type"))
  set.seed(seed)
  barcodes <- tibble(barcode = colnames(singlecell_matrix), group = groups) |>
    group_by(group) |>
    slice(sample(min(nCells_per_group, n())))

  # Check that all celltypes have enough cells
  cells_n <- barcodes |>
    summarize(n = n())
  if (any(cells_n$n < minCells_per_celltype)) {
    too_small <- cells_n |>
      filter(n < minCells_per_celltype)
    if (verbose) cli_alert(cli::col_br_magenta("  Cell type(s) {paste(too_small$group, collapse = ',')} ",
                           "have too few cells (<{minCells_per_celltype}) and will be excluded"))
    barcodes <- barcodes |>
      filter(!group %in% too_small$group)
  }
  if (verbose) cli_alert("  Kept {length(unique(barcodes$group))} cell types after filtering")
  barcodes <- barcodes |>
    pull(barcode)
  new_groups <- setNames(groups, nm = colnames(singlecell_matrix))[barcodes]

  # Compute cell type expression profiles
  if (verbose) cli_alert("  Calculating cell type expression profiles")
  W <- .get_expression_profiles(x = singlecell_matrix[, barcodes], groups =  new_groups)

  return(W)
}


#' Get single cell data expression profiles
#'
#' This function is used to obtain cell type expression profiles represented by
#' enrichment scores calculated for each cell.
#'
#' @param x Expression matrix to calculate expression profiles from
#' @param groups Character vector with group (cell type) labels
#'
#' @importFrom Matrix rowMeans
#' @importFrom tibble tibble
#' @import dplyr
#'
#' @noRd
.get_expression_profiles <- function (
    x,
    groups
){

  # Calculate means
  row_means <- do.call(cbind, lapply(unique(groups), function(grp) {
    rowMeans(x[, groups == grp])
  }))
  colnames(row_means) <- unique(groups)

  # Calculate enrichment scores
  W <- do.call(bind_cols, lapply(colnames(row_means), function(grp) {
    x1 <- row_means[, grp]
    x2 <- row_means[, -which(colnames(row_means) == grp)]
    x2 <- rowMeans(x2) + 1
    y <- tibble(x1/x2) |> setNames(nm = grp)
    return(y)
  })) |>
    as.matrix()
  rownames(W) <- rownames(x)

  W <- apply(W, 2, function(x) {
    x/max(x)
  })

  return(W)
}
