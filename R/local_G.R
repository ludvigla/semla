#' @include checks.R
#' @include generics.R
#'
NULL

#' @param spatnet A list of tibbles containing spatial networks generated with
#' \code{\link{GetSpatialNetwork}}
#' @param alternative A string specifying the alternative hypothesis: "two.sided",
#' "greater" or "less". By default, only the local G scores are returned. If an
#' alternative test is specified, the function will return both local G scores and
#' p-values. Note that p-values are adjusted for multiple hypothesis testing within
#' each feature using "BH" correction. If you want to adjust p-values with a different
#' strategy, you can compute the p-values directly from the local G z-scores.
#' @param return_as_tibble Logical specifying wthether the results should be returned
#' as an object of class `tbl`
#' @param verbose Print messages
#'
#' @importFrom glue glue
#' @importFrom rlang abort inform
#' @importFrom tidyr pivot_wider
#' @import dplyr
#' @importFrom stats pnorm p.adjust
#' @importFrom tibble rownames_to_column
#'
#' @rdname local-G
#'
#' @section default method:
#' Takes a matrix-like object with one feature per column and a list of spatial
#' networks generated with \code{\link{GetSpatialNetwork}} and computes the local G
#' scores and optionally p-values. The G scores are prefixed with "Gi" and p-values
#' are prefixed with one of "Pr(z <!=", "Pr(z " or "Pr(z <" depending on the chosen
#' test.
#'
#' @examples
#' library(STUtility2)
#' library(tibble)
#' library(dplyr)
#' library(ggplot2)
#' library(viridis)
#'
#' # read coordinates
#' coordfile <- system.file("extdata/mousebrain/spatial",
#'                          "tissue_positions_list.csv",
#'                          package = "STUtility2")
#' coords <- read.csv(coordfile, header = FALSE) |>
#'   filter(V2 == 1) |>
#'   select(V1, V6, V5) |>
#'   setNames(nm = c("barcode", "x", "y")) |>
#'   bind_cols(sampleID = 1) |>
#'   as_tibble()
#'
#' # get spatial network
#' spatnet <- GetSpatialNetwork(coords)
#'
#' # Load expression data
#' feature_matrix <- system.file("extdata/mousebrain",
#'                               "filtered_feature_bc_matrix.h5",
#'                               package = "STUtility2") |>
#'   Seurat::Read10X_h5()
#' featureMat <- t(as.matrix(feature_matrix[c("Mgp", "Th", "Nrgn"), ]))
#' featureMat <- featureMat[coords$barcode, ]
#' head(featureMat)
#'
#' # Calculate G scores
#' g_scores <- RunLocalG(log1p(featureMat), spatnet)
#' head(g_scores)
#'
#' # Bind results with coords for plotting
#' gg <- coords |>
#'   bind_cols(log1p(featureMat)) |>
#'   left_join(y = g_scores, by = "barcode")
#'
#' # Plot some results
#' p1 <- ggplot(gg, aes(x, y, color = Th))
#' p2 <- ggplot(gg, aes(x, y, color = `Gi[Th]`))
#' p3 <- ggplot(gg, aes(x, y, color = -log10(`Pr(z != E(Gi[Th]))`)))
#'
#' p <- p1 + p2 + p3 &
#'   geom_point() &
#'   theme_void() &
#'   scale_y_reverse() &
#'   coord_fixed() &
#'   scale_colour_gradientn(colours = viridis(n = 9))
#' p
#'
#' @export
#'
RunLocalG.default <- function (
  object,
  spatnet,
  alternative = c("two.sided", "greater", "less"),
  return_as_tibble = TRUE,
  verbose = TRUE,
  ...
) {

  # Set global variables to NULL
  from <- to <- barcode <- NULL

  if (!is.null(alternative)) {
    alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
    if (verbose) inform(c("i" = glue("Setting alternative hypothesis to '{alternative}'")))
  }

  # Check objects
  if (!inherits(object, what = c("dgRMatrix", "dgCMatrix", "matrix", "data.frame")))
    abort(glue("Invalid format of feature matrix: '{class(object)}'"))
  if (inherits(object, what = "data.frame")) {
    object <- as.matrix(object)
  }
  object <- as(object, "dgCMatrix")
  if (any(dim(object) == 0)) abort(glue("Invalid dimensions of feature matrix: '{paste(dim(object), collapse = 'x')}'"))
  if (!inherits(spatnet, what = "list")) abort(glue("Invalid format of spatnet: '{class(spatnet)}'"))

  # Check that all spots in spatnet are present in object data
  check <- all(sapply(spatnet, function(x) {
    c(all(x$from %in% rownames(object)), all(x$to %in% rownames(object)))
  }))
  if (!check) {
    abort(glue("Spots in 'spatnet' are not available in 'object'. ",
          "Make sure that 'spatnet' is generated from the same ",
          "dataset as 'object'"))
  }

  # pivot spatial networks in long format to a wide format
  wide_spatial_networks <- lapply(names(spatnet), function(nm) {
    wide_spatial_network <- pivot_wider(spatnet[[nm]] |>
                                          select(from, to) |>
                                          mutate(value = 1),
                                        names_from = "from",
                                        values_from = "value",
                                        values_fill = 0)
    # Convert wide spatial network to a matrix
    CN <- as.matrix(wide_spatial_network[, 2:ncol(wide_spatial_network)])
    rownames(CN) <- wide_spatial_network$to
    CN <- CN[colnames(CN), ]
    rm(wide_spatial_network)
    CN <- as(CN, "dgCMatrix")
    return(CN)
  }) |>
    setNames(nm = names(spatnet))


  Gi_stats <- lapply(names(wide_spatial_networks), function(nm) {

    wide_spatial_network <- wide_spatial_networks[[nm]]
    data_subset <- object[match(rownames(wide_spatial_network), rownames(object)), drop = FALSE, ]

    # Number of observations
    n <- nrow(data_subset)

    if (verbose) inform(c(">" = glue("  Calculating local G scores for {n} spots in sample {nm}")))

    # Calculate lag matrix
    lagMat <- wide_spatial_network %*% data_subset
    xibarMat <- apply(data_subset, 2, function(val) {
      (rep(sum(val), n) - val) / (n - 1)
    })
    si2Mat <- do.call(cbind, lapply(1:ncol(data_subset), function(i) {
      ((rep(sum(data_subset[, i]^2), n) - data_subset[, i]^2) / (n - 1)) - xibarMat[, i]^2
    }))

    # For Visium data, the weights for neighbors are 1
    weights_i <- rowSums(wide_spatial_network)
    weighted_xibarMat <- weights_i*xibarMat
    numMat <- (lagMat - weighted_xibarMat)
    denMat <- sqrt(si2Mat*(((n - 1)*weights_i - weights_i^2)/(n - 2)))
    GiMat <- numMat / denMat

    return(GiMat)

  })

  # Combines results
  GiMat_tot <- do.call(rbind, Gi_stats)
  colnames(GiMat_tot) <- glue("Gi[{colnames(GiMat_tot)}]")

  if (!is.null(alternative)) {

    if (verbose) inform(c("i" = glue("Calculating p-values for local G scores, MH-adjusted within each feature")))

    # Get function to calculate p-values
    calc_pval <- switch(alternative,
                        "two.sided" = function(x) 2*pnorm(abs(x), lower.tail = FALSE),
                        "greater" = function(x) pnorm(x, lower.tail = FALSE),
                        "less" = function(x) pnorm(x))
    alternative_symbol <- switch(alternative, "two.sided" = "!=", "greater" = ">", "less" = "<")
    col_label <- glue("Pr(z {alternative_symbol} E(Gi[{colnames(object)}]))")

    if (verbose) inform(c("i" = glue("G scores will be named Gi[Ftr] and adjusted p-values will ",
                             "be named Pr(z {alternative_symbol} E(Gi[Ftr]))")))

    # Calculate p-values
    pvs <-
      apply(GiMat_tot, 2, function(x) {
        calc_pval(x) |> p.adjust(method = "BH")
      }) |>
        matrix(ncol = ncol(GiMat_tot), byrow = FALSE)

    # Set colnames
    colnames(pvs) <- col_label

    # Merge results if return_as_tibble = TRUE
    if (return_as_tibble) {
      fullMat <- cbind(GiMat_tot, pvs)[, c(rbind(1:ncol(GiMat_tot), (ncol(GiMat_tot) + 1):(2*ncol(GiMat_tot))))] |>
        as_tibble() |>
        mutate(barcode = rownames(GiMat_tot)) |>
        relocate(barcode)
      return(fullMat)
    } else {
      # Otherwise, return results in a list
      return(list(GiMat = GiMat_tot, pvals = pvs))
    }
  }

  # If the calculation of p-values is omitted, return Gi scores
  if (return_as_tibble) {
    GiMat_tot_tibble <- GiMat_tot |> as_tibble()
    GiMat_tot_tibble <- GiMat_tot_tibble |>
      mutate(barcode = rownames(GiMat_tot)) |>
      relocate(barcode)
    return(GiMat_tot_tibble)
  } else {
    return(GiMat_tot)
  }

}



#' @param features A character vector of feature names fetchable with
#' \code{\link{FetchData}}
#' @param store_in_metadata A logical specifying if the results should be
#' returned in the meta data slot. If set to FALSE, the results will instead
#' be returned as an assay named by \code{assay_name}. If a statistical test
#' is applied, the adjusted p-values will be returned to the meta data slot if
#' \code{store_in_metadata = TRUE} and in the \code{@misc} slot of the assay
#' if \code{store_in_metadata = FALSE}.
#' @param assay_name Name of the assay if \code{store_in_metadata=FALSE}
#'
#' @section Seurat:
#' Takes a Seurat object as input and returns the local G scores for a
#' selected number of features and optionally p-values. The G scores are
#' prefixed with "Gi" and p-values are prefixed with one of "Pr(z <!=",
#' "Pr(z " or "Pr(z <" depending on the chosen test.
#'
#' @import dplyr
#' @importFrom Seurat FetchData CreateAssayObject
#' @importFrom rlang inform abort
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @rdname local-G
#'
#' @examples
#'
#' library(STUtility2)
#' library(dplyr)
#'
#' # Load Seurat object with mouse bain data
#' se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "STUtility2"))
#' se_mbrain <- se_mbrain |>
#'   ScaleData(verbose = FALSE) |>
#'   FindVariableFeatures(verbose = FALSE) |>
#'   RunPCA(verbose = FALSE)
#'
#' # Calculate G scores
#' se_mbrain <- RunLocalG(se_mbrain, features = c("Th", "Mgp"), alternative = "greater")
#'
#' # Plot G scores
#' MapFeatures(se_mbrain, features = "Gi[Th]")
#'
#' # high/low clustering
#' se_mbrain$cluster <- se_mbrain[[]] |>
#'   mutate(cluster = case_when(
#'     `Pr(z > E(Gi[Th]))` > 0.05 ~ "Not Significant",
#'     `Pr(z > E(Gi[Th]))` <= 0.05 & `Gi[Th]` > 0 ~ "High"
#'   )) |> pull(cluster)
#' MapLabels(se_mbrain, column_name = "cluster")
#'
#' @export
#'
RunLocalG.Seurat <- function (
  object,
  features,
  alternative = NULL,
  store_in_metadata = TRUE,
  assay_name = "GiScores",
  verbose = TRUE,
  ...
) {

  if (verbose) cli_h2("Calculating local G")

  # Set global variables to NULL
  barcode <- NULL

  # Check objects
  .check_seurat_object(object)

  # Get a spatial network
  spatnet <- GetSpatialNetwork(object)

  # Get data
  if (!inherits(features, what = "character"))
    abort("Invalid class '{class(features)}' for 'features', expected a 'character' vector")
  if (length(features) == 0)
    abort("'features' is empty")
  data <- FetchData(object, vars = features)
  if (verbose) inform(c("i" = glue("Got {length(features)} feature(s)")))

  # Calculate local G
  Gi_res <- RunLocalG(data,
                      spatnet = spatnet,
                      alternative = alternative,
                      return_as_tibble = TRUE,
                      verbose = verbose,
                      ... = ...)

  # Return in meta data if store_in_metadata = TRUE
  if (store_in_metadata) {
    if (verbose)
      inform(c("i" = "Placing results in 'Seurat' object meta.data slot"))
    newMdata <- object@meta.data |>
      select(-contains(colnames(Gi_res))) |>
      rownames_to_column(var = "barcode") |>
      left_join(y = Gi_res, by = "barcode") |>
      column_to_rownames(var = "barcode")
    object@meta.data <- newMdata
  } else {
    if (verbose)
      inform(c("i" = glue("Placing results in 'Seurat' object as an 'Assay' object named {assay_name}")))
    data <- Gi_res |>
      select(-starts_with("Pr(z")) |>
      column_to_rownames(var = "barcode") |>
      as.matrix()
    # Add empty rows if missing
    fillMat <- matrix(data = 0, ncol = ncol(object), nrow = length(features),
                      dimnames = list(colnames(data), colnames(object)))
    fillMat[, rownames(data)] <- t(data)
    gi_assay <- CreateAssayObject(data = fillMat)
    pvals <-  Gi_res |>
      select(barcode, starts_with("Pr(z")) |>
      column_to_rownames(var = "barcode") |>
      as.matrix()
    # Add empty rows to pvals if missing
    pvalMat <- matrix(data = 0, ncol = ncol(object), nrow = length(features),
                      dimnames = list(colnames(pvals), colnames(object)))
    pvalMat[, rownames(pvals)] <- pvals
    gi_assay@misc$pvals <- pvalMat
    object[[assay_name]] <- gi_assay
  }

  if (verbose) inform(c("v" = "Returning results"))
  return(object)

}

