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
#' @param k An integer specifying the number of factors to compute 
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
#' @rdname celltype-prediction-mixed
#'
#' @section default method:
#' Input `object` is a matrix-like object with 10x Visium data. The function returns
#' a matrix with estimated proportions of dimensions (nCellTypes + k) x nSpots,
#' where nCellTypes is the number of cell types, k is the number of factors and nSpots is the number of spots.
#' If `return_expression_profiles=TRUE`, the returned object will be a list with estimated
#' proportions `prop` and the cell type/factor expression profile matrix `W`.
#'
#' @export
RunMixedNNLS.default <- function (
    object,
    singlecell_matrix,
    groups,
    k = 10L,
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
  if (!inherits(k, what = "integer"))
    abort(glue("Invalid class '{class(k)}' for k. Expected an integer of length 1"))
  if (length(k) != 1)
    abort(glue("Invalid length of k. Expected an integer of length 1"))
  if (k < 2)
    abort(glue("Invalid value for k={{k}}. k needs to be a positive integer larger than 1"))
  
  
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
  if ((packageVersion("RcppML") |> as.character()) != "0.3.7") {
    cli_alert_warning("  The NNLS function might break if a dev version of RcppML is used. ")
    cli_alert_warning("  If RcppML::project(...) fails, try installing CRAN version 0.3.7 of RcppML.")
  }
  
  # Prepare data
  if (verbose) cli_alert_info("Preparing data for NNLS")
  W <- .get_w_matrix(object = object,
                     singlecell_matrix = singlecell_matrix,
                     nCells_per_group = nCells_per_group,
                     minCells_per_celltype = minCells_per_celltype,
                     groups = groups,
                     seed = seed,
                     verbose = verbose)
  
  # Run nmf
  if (verbose) cli_alert_info("Computing NMF with rank {ncol(W) + k} (k={k})")
  nmf_model <- RcppML::nmf(data = object |> as.matrix(), k = ncol(W) + k, L1 = L1, maxit = 1, verbose = FALSE)
  W_nmf <- apply(nmf_model$w, 2, function(x) {
    x/max(x)
  })
  
  # Check correlations
  if (verbose) cli_alert_info("Keeping {k} factors with lowest correlation scores")
  corMat <- try({RcppML::cosine(W, W_nmf)})
  if (inherits(corMat, what = "try-error")) {
    corMat <- cor(W, W_nmf, method = "pearson") |> as.matrix()
  }
  
  cors_ordered <- apply(corMat, 2, function(x) {
   max(x)
  }) |> order()
  W_combined <- cbind(W, W_nmf[, cors_ordered[1:k]])
  
  # Run NNLS
  if (verbose) cli_alert_info("Predicting proportions with NNLS for {ncol(W)} cell types and {k} factors")
  proj_expr <- try({RcppML::project(object, W_combined, L1 = L1, ...)}, silent = TRUE)
  if (inherits(proj_expr, what = "try-error")) {
    proj_expr <- RcppML::project(w = W_combined, data = object, L1 = L1, ...)
  }
  
  # Convert predicted values to proportions
  prop <- apply(proj_expr, 2, function(x) {prop.table(x)})
  prop[prop < min_prop] <- 0
  prop <- apply(prop, 2, function(x) {prop.table(x)})
  rownames(prop) <- c(colnames(W), paste0("NMF-", 1:k))
  colnames(prop) <- colnames(object)
  
  if (return_expression_profiles) {
    return(list(prop = prop, W = W_combined))
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
#' estimated proportions for cell types and additional factors.
#'
#' @importFrom Seurat GetAssayData CreateAssayObject CreateDimReducObject
#' @import cli
#'
#' @rdname celltype-prediction-mixed
#'
#' @return An object with cell type proportion estimates
#'
#' @export
#'
#' @md
RunMixedNNLS.Seurat <- function (
    object,
    singlecell_object,
    groups = NULL,
    k = 10L,
    features = NULL,
    singlecell_assay = "RNA",
    spatial_assay = "Spatial",
    slot = "data",
    min_prop = 0.01,
    nCells_per_group = 50L,
    minCells_per_celltype = 10L,
    assay_name = "celltypepropsmixed",
    dimred_name = "nnlsmixed",
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

  results <- RunMixedNNLS.default(
    object = x_spatial,
    singlecell_matrix = x_singlecell,
    groups = groups,
    k = k,
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
    dimred_prefix = "cmFactor_",
    verbose = verbose
  )
  
  if (verbose) cli_alert_success("Finished")
  
  return(object)
}
