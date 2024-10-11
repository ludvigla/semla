#' @include checks.R
#'
NULL

#' Map points to nearest reference spots within a distance
#' 
#' This function takes two dataframes containing point coordinates 
#' (pixel column and pixel row) in the first two columns and identifies the nearest 
#' "spot" in the second dataframe for each point in the first dataframe, considering 
#' a maximum distance and a maximum number of neighbors. The coordinate ID needs to
#' be in the rownames of the input dataframes.
#' 
#' @param ref_coords A dataframe containing the reference point coordinates 
#' with the first two columns corresponding to x- and y- pixel coordinates and
#' coordinate IDs as rownames.
#' @param map_coords A dataframe containing the coordinates of point to map
#' to the reference, with the first two columns corresponding to x- and y- pixel 
#' coordinates and coordinate IDs as rownames.
#' @param distance_max The maximum distance threshold for considering 
#' a spot as a neighbor (defaults to half the minimum distance between reference points).
#' @param n_neighbors The maximum number of neighbors to consider for 
#' each reference point (defaults to 100).
#' 
#' @return A dataframe with "map_coord_IDs" indicating the nearest map spot ID 
#' for each point in the reference dataframe.
#' 
#' @import tidyr
#' 
#' @export
MapPointsToReference <- function(
    ref_coords, 
    map_coords, 
    distance_max = NULL, 
    n_neighbors = 100) {
  
  if (n_neighbors<1) {
    abort(glue("Invalid {col_br_magenta('n_neighbors')}. Expected a numeric of at least 1"))
  }
  
  # Calculate minimum center-to-center distance between reference points
  min_cc_dist <- dbscan::kNN(ref_coords[, c(1, 2)], k = 1)$dist[, 1] |> min()
  half_cc_dist <- min_cc_dist / 2
  if (is.null(distance_max) || distance_max > half_cc_dist) {
    distance_max <- half_cc_dist
  }
  
  # Prepare dataframes with origin and ID information
  ref_df <- ref_coords %>%
    select(1:2) %>%
    mutate(origin = "ref", coord_ID = rownames(.))
  map_df <- map_coords %>%
    select(1:2) %>%
    mutate(origin = "map", coord_ID = rownames(.))
  rownames(map_df) <- paste0("map_", rownames(map_df))
  
  # Combine dataframes and calculate k nearest neighbors
  full_df <- rbind(ref_df, map_df)
  full_df$index <- 1:nrow(full_df)
  knn_results <- dbscan::kNN(as.matrix(full_df[, c(1, 2)]), k = n_neighbors, sort = TRUE)
  
  # Extract distances and IDs for reference points
  ref_res_dist <- knn_results$dist[ref_df$coord_ID, ]
  ref_res_id <- knn_results$id[ref_df$coord_ID, ]
  
  # Filter based on distance and remove reference points themselves by setting them as NA
  ref_idx <- full_df[rownames(ref_df), "index"]
  ref_res_id_filt <- ref_res_id
  ref_res_id_filt[ref_res_dist > distance_max] <- NA
  ref_res_id_filt[ref_res_id %in% ref_idx] <- NA
  
  # Replace neighbor IDs with corresponding map spot IDs
  map_df$index <- full_df[rownames(map_df), "index"]
  ref_res_id_filt <- apply(ref_res_id_filt, 2, function(x) map_df$coord_ID[match(x, map_df$index)])
  rownames(ref_res_id_filt) <- rownames(ref_res_id)
  
  # Find matching reference IDs for each map spot ID
  map_ref_df <- ref_res_id_filt |> 
    as.data.frame() |> 
    mutate(ref_id = rownames(ref_res_id)) |> 
    pivot_longer(cols = c(-ref_id), values_to = "map_id") |> 
    drop_na(map_id) |> 
    select(map_id, ref_id)
  
  map_ids_order <- map_df$coord_ID[map_df$coord_ID %in% ref_res_id_filt]
  map_ref_df <- map_ref_df[match(map_ids_order, map_ref_df$map_id),] |> as_tibble()
  
  return(map_ref_df)
}

#' Create a multimodal object in a shared spot coordinate framework
#'
#' This function takes two sets of spatial coordinates (the reference and the one that
#' should be mapped to the reference coordinate system), with matching pixel coordinates, 
#' along with a feature count matrix, and generates a new count matrix 
#' containing the aggregated feature values for each reference coordinate based 
#' on the mapping relationship between the two coordinate systems.
#'
#' @details This function assumes the reference coordinates represent fixed 
#' locations (e.g., Visium spots) and the mapped coordinates represent measurements 
#' from another modality (e.g., Mass Spec Imaging) that have already been aligned to 
#' the same pixel coordinate system. Unless something else is specified, the coordinates
#' from the second modality is mapped to the reference spots by including all coordinates
#' falling within half of the center-to-center distance between reference spots.
#' Count values are then aggregated for each feature using the chosen aggregation function 
#' ("mean" or "sum") into the reference coordinate system. 
#' The function can optionally utilize multi-threading for 
#' faster processing with large datasets.
#'
#'
#' @param object_ref Data frame containing reference coordinates. Required columns 
#' include "barcode", "x" (x-coordinate), "y" (y-coordinate), and "sampleID".
#' @param object_map Data frame containing mapped coordinates. Required columns 
#' include "barcode", "x" (x-coordinate), "y" (y-coordinate), and "sampleID".
#' @param map_counts Feature count matrix where rows represent features and columns 
#' represent mapped coordinates identified by their "barcode" values.
#' @param agg_func Character vector specifying either "mean" or "sum", corresponding to
#' the aggregation function to use for summarizing feature values mapped to each 
#' reference coordinate. Default set to "mean".
#' @param n_neighbors Integer specifying the number of nearest neighbors to consider 
#' for mapping each point in the mapped dataset to the reference dataset. Defaults 
#' to 100.
#' @param distance_max Numeric value specifying the maximum allowed distance 
#' between a mapped coordinate and its nearest reference neighbor(s). Points 
#' exceeding this distance will be excluded from the mapping. Defaults to NULL 
#' (maximum distance set to half the center-to-center distance between reference spots).
#' @param nCores Integer specifying the number of cores to use for parallel 
#' processing during data aggregation. Defaults to (detectCores() - 1), which uses 
#' all available cores minus one. Set to NULL to disable multi-threading.
#' @param verbose Logical indicating whether to print informative messages during 
#' execution. Defaults to TRUE.
#'
#' @return A list containing two elements:
#'   - `data`: A matrix representing the aggregated feature values for each 
#'   reference coordinate. Rows represent features, and columns represent reference 
#'   coordinates.
#'   - `map`: A data frame containing detailed information about the mapping 
#'   between reference and mapped coordinates, including the reference ID(s) for 
#'   each mapped point and the chosen aggregation function value(s) for each 
#'   feature and reference coordinate combination.
#' 
#' @importFrom dbscan kNN
#' @importFrom glue glue
#' @importFrom parallel detectCores
#' @importFrom Matrix Matrix
#' @import cli
#' @import rlang
#' @import dplyr
#'
#' @rdname multimodal-object
#' @family pre-process
#'
#'
#' @export
CreateMultiModalObject.default <- function (
    object_ref,
    object_map,
    map_counts,
    agg_func = "mean",
    n_neighbors = 100,
    distance_max = NULL,
    nCores = (parallel::detectCores() - 1),
    verbose = TRUE,
    ...
) {
  
  # Run checks
  # Check input objects
  if (!inherits(object_ref, what = c("dgRMatrix", "dgCMatrix", "matrix", "data.frame")))
    abort(glue("Invalid format of reference coordinate matrix: '{class(object_ref)}'"))
  if (any(dim(object_ref) == 0))
    abort(glue("Invalid dimensions of reference coordinate matrix: '{paste(dim(object_ref), collapse = 'x')}'"))
  if (!inherits(object_map, what = c("dgRMatrix", "dgCMatrix", "matrix", "data.frame")))
    abort(glue("Invalid format of mapping coordinate matrix: '{class(object_map)}'"))
  if (any(dim(object_map) == 0))
    abort(glue("Invalid dimensions of mapping coordinate matrix: '{paste(dim(object_map), collapse = 'x')}'"))
  if (!inherits(map_counts, what = c("dgRMatrix", "dgCMatrix", "matrix", "data.frame")))
    abort(glue("Invalid format of mapping count matrix: '{class(object_map)}'"))
  if (any(dim(map_counts) == 0))
    abort(glue("Invalid dimensions of mapping count matrix: '{paste(dim(object_map), collapse = 'x')}'"))
  
  # Check that arg input are valid
  agg_func <- match.arg(agg_func, choices = c("mean", "sum"), several.ok = FALSE)
  if (!inherits(n_neighbors, what = "numeric"))
    abort(glue("Invalid class '{class(n_neighbors)}' for number of neighbors, expected a positive 'numeric' value"))
  if (length(n_neighbors) == 0)
    abort("'n_neighbors' is empty")
  if (!n_neighbors > 0)
    abort("'n_neighbors' should be larger than 0")
  if (!is.null(distance_max))
    if (!inherits(distance_max, what = "numeric"))
      abort(glue("Invalid class '{class(distance_max)}' for maximum distance, expected a positive 'numeric' value"))
  
  
  # Split coordinates by sample
  xys_ref.list <- object_ref |>
    group_by(sampleID) |>
    group_split() |>
    setNames(paste0(unique(object_ref$sampleID)))
  
  xys_map.list <- object_map |>
    group_by(sampleID) |>
    group_split() |>
    setNames(paste0(unique(object_map$sampleID)))
  
  # Map coordinates per sample
  mapped_coords.list <- setNames(lapply(names(xys_ref.list), function(sampleID) {
    # Subset coordinates based on sampleID
    xys_ref_subset <- xys_ref.list[[sampleID]] |> as.data.frame()
    rownames(xys_ref_subset) <- xys_ref_subset$barcode
    xys_map_subset <- xys_map.list[[sampleID]] |> as.data.frame()
    rownames(xys_map_subset) <- xys_map_subset$barcode
    
    # Map points
    map_ref_df <- MapPointsToReference(ref_coords = xys_ref_subset[, c("x","y")], 
                                       map_coords = xys_map_subset[, c("x","y")], 
                                       distance_max = distance_max, 
                                       n_neighbors = n_neighbors)
    map_ref_df$sampleID <- sampleID
    return(map_ref_df)
  }), nm = names(xys_ref.list))
  mapped_coords.merged <- bind_rows(mapped_coords.list)
  
  
  # Aggregate map-data values based on ref-data mapping by;
  
  # 1. Filter map count data
  map_ids_keep <- mapped_coords.merged |> pull(map_id)
  map_counts_subset <- map_counts[, map_ids_keep]
  
  if (verbose) cli_alert_info(glue("Excluding {cli::col_br_cyan(ncol(map_counts)-ncol(map_counts_subset))}",
                                   " coordinates from the mapping dataset due to being outside of reference coordinate scope"))
  
  # 2. Compute sum and mean values for map coords per ref coord
  map_counts_subset_t <- map_counts_subset |> # Takes quite a long time (when using very big datasets)
    Matrix::t() |> 
    as.data.frame()
  
  
  # Aggregate data per coordinate (spot), with or without mutli-threading
  if (verbose) cli_alert_info(glue("Aggregating assay data in mapping coordinates per reference coordinate", 
                                   " (this step may take long if the data is large)"))
  start.time <- Sys.time()
  x_subset <- map_counts_subset_t
  mapping_df <- mapped_coords.merged
  
  if (is.null(nCores)) {
    agg_data_res <- .aggDataPerSpot(x_subset, mapping_df)
    
  } else {
    # Define multicore lapply function depending on OS
    if (!.Platform$OS.type %in% c("windows", "unix")) {
      # Skip threading and use lapply instead
      if (verbose) cli_alert_danger("Threading not supported. Using single thread.")
      parLapplier <-  function(X, FUN, nCores) {
        res <- lapply(X, FUN)
        return(res)
      }
    } else {
      parLapplier <- switch(.Platform$OS.type,
                            "windows" = .winLapply,
                            "unix" = .unixLapply)
    }
    if (nCores > (detectCores() - 1)) {
      nCores <- detectCores() - 1
    }
    if (verbose) cli_alert_info("Using {nCores} threads")
    
    # Define chunks based on reference coord IDs
    ref_id_unique <- mapping_df$ref_id |> unique()
    chunks <- ceiling((1:length(ref_id_unique))/100)
    chunks <- split(1:length(ref_id_unique), chunks)
    
    
    agg_data_res <- parLapplier(seq_along(chunks), function(i) {
      ref_id_ind <- chunks[[i]]
      map_ind <- mapping_df |> filter(ref_id %in% ref_id_unique[ref_id_ind])
      x_subset_ind <- x_subset[map_ind$map_id, ]
      # Aggregate data per chunk
      .aggDataPerSpot(x_subset_ind, map_ind)
    }, 
    nCores = nCores) |> 
      bind_rows()
  }
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units = 'mins') |> round(digits = 2)
  if (verbose) cli_alert_success("Data aggregation finished in {time.taken} minutes")
  
  # Prepare sum/mean mapping data for adding to assay
  if (verbose) cli_alert_info("Storing aggregated data summarized by {agg_func} values")
  agg_data_long <- agg_data_res |> select(all_of(c("ref_id", "feature", agg_func, "sampleID")))
  colnames(agg_data_long)[3] <- "value"
  
  agg_data_wide <- pivot_wider(agg_data_long |> select(-sampleID), values_from = "value", names_from = "feature") |> 
    as.data.frame() |> 
    column_to_rownames(var = "ref_id") |>
    as.matrix() |> 
    t()
  
  # Prepare metadata info about map-coordinate per ref-coordinate to object
  map_coords_per_ref_coord <- agg_data_res |>
    select(ref_id, coords_per_ref_coord) |> 
    group_by(ref_id) |> 
    slice(which.max(coords_per_ref_coord)) |> 
    column_to_rownames(var = "ref_id")
  
  res.list <- list(agg_data_wide, map_coords_per_ref_coord) |> setNames(nm = c("data", "map"))
  return(res.list)
}



#' @param object_ref Seurat/semla object containing reference coordinates.
#' @param object_map Seurat/semla object containing mapped coordinates.
#' @param new_assay_name Name to be assigned to the new Assay. If left empty, the same
#' name as the `object_map` Active Assay will be used (default).
#' @param agg_func Character vector specifying either "mean" or "sum", corresponding to
#' the aggregation function to use for summarizing feature values mapped to each 
#' reference coordinate. Default set to "mean".
#' @param n_neighbors Integer specifying the number of nearest neighbors to consider 
#' for mapping each point in the mapped dataset to the reference dataset. Defaults 
#' to 100.
#' @param distance_max Numeric value specifying the maximum allowed distance 
#' between a mapped coordinate and its nearest reference neighbor(s). Points 
#' exceeding this distance will be excluded from the mapping. Defaults to NULL 
#' (maximum distance set to half the center-to-center distance between reference spots).
#' @param assay Seurat object assay in `object_map` to fetch data from
#' @param slot  Seurat object slot in `object_map` to fetch data from (Seurat < v.5)
#' @param layer  Seurat object layer in `object_map` to fetch data from (Seurat ≥ v.5)
#' @param nCores Integer specifying the number of cores to use for parallel 
#' processing during data aggregation. Defaults to (detectCores() - 1), which uses 
#' all available cores minus one. Set to NULL to disable multi-threading.
#' @param verbose Logical indicating whether to print informative messages during 
#' execution. Defaults to TRUE.
#' 
#' @return A \code{\link{Seurat}/semla} multi-modal object with the added modality data
#' in a new Assay, all in a shared coordinate system.
#'
#' @importFrom dbscan kNN
#' @importFrom glue glue
#' @importFrom parallel detectCores
#' @importFrom Matrix Matrix
#' @import cli
#' @import rlang
#' @import dplyr
#'
#' @rdname multimodal-object
#' @family pre-process
#' 
#' @author Lovisa Franzén
#'
#' @examples
#'
#' se_mod1 <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
#' se_mod2 <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
#' 
#' # By default, the mean values of the mapping modality will be calculated per 
#' # reference coordinates
#' se_mmo <- CreateMultiModalObject(object_ref = se_mod1, object_map = se_mod2, 
#'                                  new_assay_name = "Modality2", 
#'                                  agg_func = "mean")
#'
#' # To change this and instead compute the sum, the `agg_func` argument can modified:
#' se_mmo <- CreateMultiModalObject(object_ref = se_mod1, object_map = se_mod2, 
#'                                  new_assay_name = "Modality2", 
#'                                  agg_func = "sum")
#'
#' @export
CreateMultiModalObject.Seurat <- function (
    object_ref,
    object_map,
    new_assay_name = NULL,
    agg_func = "mean",
    n_neighbors = 100,
    distance_max = NULL,
    assay = NULL, # if NULL, uses default assay
    slot = "counts", # Seurat < v.5
    layer = "counts", # Seurat ≥ v.5
    nCores = (parallel::detectCores() - 1),
    verbose = TRUE,
    ...
) {
  
  # Set global variable to NULL
  barcode <- pxl_col_in_fullres <- pxl_row_in_fullres <- sampleID <- NULL
  
  # Check that arg input are valid
  agg_func <- match.arg(agg_func, choices = c("mean", "sum"), several.ok = FALSE)
  if (!inherits(n_neighbors, what = "numeric"))
    abort(glue("Invalid class '{class(n_neighbors)}' for number of neighbors, expected a positive 'numeric' value"))
  if (length(n_neighbors) == 0)
    abort("'n_neighbors' is empty")
  if (!n_neighbors > 0)
    abort("'n_neighbors' should be larger than 0")
  if (!is.null(distance_max))
    if (!inherits(distance_max, what = "numeric"))
      abort(glue("Invalid class '{class(distance_max)}' for maximum distance, expected a positive 'numeric' value"))
  
  
  if (verbose) cli_h2("Joining second modality data to reference object")

  if (is.null(new_assay_name)) {
    new_assay_name <- object_map@active.assay
  }
  if (new_assay_name == object_ref@active.assay) {
    new_assay_name <- paste0(new_assay_name, "2")
  }
  
  # # Check Seurat object
  .check_seurat_object(object_ref)
  .check_seurat_object(object_map)
  
  map_counts <- suppressWarnings({GetAssayData(object_map, assay = assay, slot = slot, layer = layer)})
  
  # Get coordinates
  xys_ref <- GetStaffli(object_ref)@meta_data |>
    select(barcode, pxl_col_in_fullres, pxl_row_in_fullres, sampleID) |>
    rename(x = pxl_col_in_fullres, y = pxl_row_in_fullres)
  
  xys_map <- GetStaffli(object_map)@meta_data |>
    select(barcode, pxl_col_in_fullres, pxl_row_in_fullres, sampleID) |>
    rename(x = pxl_col_in_fullres, y = pxl_row_in_fullres)
  
  # Create the new Seurat assay
  new_assay.list <- CreateMultiModalObject.default(object_ref = xys_ref, 
                                                   object_map = xys_map, 
                                                   map_counts = map_counts,
                                                   agg_func = agg_func, 
                                                   n_neighbors = n_neighbors, 
                                                   distance_max = distance_max,
                                                   nCores = nCores,
                                                   verbose = verbose,
                                                   ...)
  
  agg_data_wide <- new_assay.list$data
  map_coords_per_ref_coord <- new_assay.list$map
  
  # Prepare and filter new multimodal object to only include intersecting coordinates 
  obj_combined <- object_ref
  obj_combined <- SubsetSTData(obj_combined, spots = colnames(agg_data_wide))
  
  # Create and add new assay
  if (verbose) cli_alert_info(glue("Adding second modality data to a new 'Assay' named '{new_assay_name}', containing", 
                                   " {cli::col_br_blue(nrow(agg_data_wide))} features"))
  obj_combined[[new_assay_name]] <- CreateAssayObject(counts = agg_data_wide)
  
  # Add metadata info about map-coordinate per ref-coordinate to object
  obj_combined <- AddMetaData(obj_combined, metadata = map_coords_per_ref_coord, col.name = "coords_per_ref")
  
  if (verbose) cli_alert_success(glue("Returning a `Seurat` object with {cli::col_br_blue(nrow(obj_combined))}",
                                      " features and {cli::col_br_magenta(ncol(obj_combined))} spots"))
  
  return(obj_combined)
}


#' Aggregate data from mapping matrix to reference coordinates
#' 
#' @param count_matrix A matrix with map-coordinate-IDs in rows and map-features in columns
#' @param mapping_df A \code{tbl} with columns \code{map_id}, \code{ref_id}, and \code{sampleID}
#' 
#' @import dplyr
#' @import cli
#' 
#' @return A \code{tbl} containing aggregated (sum and mean) map-feature data for each \code{ref_id}
#' 
#' @noRd
.aggDataPerSpot <- function (count_matrix, mapping_df) {
  feats <- colnames(count_matrix)
  count_matrix$map_id <- rownames(count_matrix)
  
  count_matrix_long <- pivot_longer(count_matrix,
                                    cols = all_of(feats),
                                    names_to = "feature", 
                                    values_to = "value" 
  )
  
  # Add reference coord IDs
  count_matrix_long$ref_id <- mapping_df$ref_id[match(count_matrix_long$map_id, mapping_df$map_id)]
  
  agg_count_data <- count_matrix_long |>  # Runs very very slow atm...
    group_by(ref_id, feature) |> 
    summarize(sum = sum(value),
              coords_per_ref_coord = n()) |> 
    mutate(mean = sum / coords_per_ref_coord) |> 
    select(ref_id, feature, sum, mean, coords_per_ref_coord) |> 
    ungroup()
  
  agg_count_data$sampleID <- mapping_df$sampleID[match(agg_count_data$ref_id, mapping_df$ref_id)] # Add sampleID info
  
  return(agg_count_data)
}
