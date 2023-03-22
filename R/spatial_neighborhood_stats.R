#' @include checks.R spatial_utils.R
#'
NULL

#' Neighborhood Enrichment Analysis
#'
#' @description
#' Performs Neighborhood Enrichment Analysis between spot labels, describing whether
#' spots of two categories lie next to each other spatially more often than expected
#' by chance.
#'
#' @details
#' This analysis calculates the enrichment score, z-score, based on how often spots
#' of different categorical labels (specified with \code{column_name}) lies
#' adjacent to each other. The observed number of edges between the labels is
#' then compared with the results from a set number of permutations (chosen
#' with \code{n_permutations}), allowing the calculation of a z-score.
#' The z-score will indicate if a label pair is overrepresented or underrepresented
#' as compared to what would be expected to see by chance.
#' The output of this function is a tibble table that for each label pair contains
#' information about the observed number of edges (\code{edges}), the mean of the permuted
#' results (\code{perm_mean}), the standard deviation of the permuted results (\code{perm_sd}),
#' and the z-score (\code{z_score}).
#'
#' @param object A Seurat object
#' @param column_name Column name in metadata corresponding to label ID of the spots.
#' @param column_labels Optional. Provide vector of label IDs to subset your analysis by.
#' Spots of all other labels will be excluded. Must be more than one. Default is all (NA).
#' @param n_permutations Integer specifying number of iterations the labels should be randomized
#' [default: 200]. Recommended to increase the number of permutations to >=200 for more robust
#' results. A lower number of permutations will result in high standard deviations and thus more
#' unreliable z-scores.
#' @param seed A seed for reproducibility [default: 123]
#' @param nCores Number of cores [default: parallel::detectCores() - 1]
#' @param verbose Print messages [default: TRUE]
#'
#' @return A tibble with scores for each label pair.
#'
#' @import dplyr
#' @import cli
#' @import tibble
#' @importFrom tidyr crossing replace_na
#' @importFrom rlang abort
#' @importFrom parallel detectCores
#'
#' @rdname neighborhood-enrichment
#' @family spatial-methods
#'
#' @author Lovisa Franzén
#'
#'
#' @examples
#' library(semla)
#'
#' # Read data
#' se <- readRDS(system.file("extdata/mousebrain",
#'                           "se_mbrain",
#'                           package = "semla"))
#'
#' # Generate clusters
#' se <- se |>
#'   NormalizeData() |>
#'   ScaleData() |>
#'   FindVariableFeatures() |>
#'   RunPCA() |>
#'   FindNeighbors(reduction = "pca", dims = 1:30) |>
#'   FindClusters()
#'
#' # Run Neigborhood Enrichment Analysis
#' res <- RunNeighborhoodEnrichmentTest(object = se,
#'                                      column_name = "seurat_clusters",
#'                                      n_permutations = 100,
#'                                      nCores = 2)
#'
#' res |> arrange(desc(abs(z_score)))
#'
#'
#' @export
RunNeighborhoodEnrichmentTest <- function(
  object,
  column_name,
  column_labels = NA,
  n_permutations = 200,
  nCores = parallel::detectCores() - 1,
  seed = 123,
  verbose = TRUE
){
  # Set global variables to NULL
  edges <- label_from <- label_1 <- label_2 <- label_label <- NULL

  # checks - import from checks.R and spatial_utils.R
  .check_seurat_object(object)
  .validate_column_name(object = object, column_name = column_name)

  # checks
  if (n_permutations != round(n_permutations)) rlang::abort("Invalid input for 'n_permutations', expected a numeric integer.")
  if (n_permutations <= 1) rlang::abort("Invalid input for 'n_permutations', needs to be more than 1.")

  # check column labels
  if (!is.na(column_labels[1])) {
    if (length(unique(column_labels))<=1) rlang::abort("Provided unique 'column_labels' needs to be more than 1.")
    if (sum(column_labels %in% unique(object[[column_name]][,1])) < length(column_labels)) rlang::abort("Some of the provided 'column_labels' can not be found among the 'column_name' labels. Please provide valid 'column_labels'.")
  }
  
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


  # Start analysis
  if (verbose) cli_h2("Running Neighborhood Enrichment Analysis")
  if (verbose) cli_alert_info("Generating neighborhood adjacency data from observed labels in column '{column_name}'")

  if (!is.na(column_labels[1])) {
    if (sum(duplicated(column_labels))>0) {
      cli_alert_warning("Provided 'column_labels' contains duplicated values. Continuing analysis with abundant labels removed.")
      column_labels <- unique(column_labels)
    }
    if (verbose) cli_alert_info("Analysis limited to study {length(column_labels)} unique labels")
  }

  # Filter object if custom column_labels are provided
  if (!is.na(column_labels[1])) {
    spots_keep <- rownames(object[[]])[object[[]][,column_name] %in% column_labels]
    object <- SubsetSTData(object, spots = spots_keep)
  }

  # Generate spatnet
  spatnet <- do.call(bind_rows, GetSpatialNetwork(object))
  spatnet$label_from <- object[[]][spatnet$from, column_name, drop = TRUE]
  spatnet$label_to <- object[[]][spatnet$to, column_name, drop = TRUE]


  # Create template with available label-label options
  if (is.na(column_labels[1])) {
    unique_labels <- paste0("Label_", levels(object[[column_name]][,1]))
  } else {
    unique_labels <- paste0("Label_", column_labels)
  }
  unique_labels <- factor(unique_labels, levels = unique_labels)
  if (length(unique_labels) <= 1) rlang::abort("Too few unique labels included in the analysis. Make sure there are at least two unique labels provided within 'column_name'.")

  label_label_template <- setNames(tidyr::crossing(unique_labels, unique_labels),
                                   nm = c("label_1", "label_2")) |>
    mutate(label_label = paste0(label_1, "-", label_2))


  # Generate weighted edge lists from region neighbors of observed labels
  edge_table <- spatnet |>
    tibble::add_column(
      label_1 = paste0("Label_", as.character(object[[]][spatnet$from, column_name, drop = TRUE])),
      label_2 = paste0("Label_", as.character(object[[]][spatnet$to, column_name, drop = TRUE]))
    ) |>
    mutate(label_label = paste0(label_1, "-", label_2))

  # Count number of edges per unique label-label pair
  edge_table_stats <- edge_table |>
    count(label_label) |>
    rename(edges = n)

  edge_table_stats <- label_label_template |>
    left_join(y = edge_table_stats |> select(label_label, edges), by = "label_label") |>
    mutate_at("edges", ~tidyr::replace_na(., 0))

  if (verbose) cli_alert_success("Observed label adjacency calculations complete")


  # Generate weighted edge lists from region neighbors of randomized labels
  if (verbose) cli_alert_info("Generating neighborhood adjacency data from randomized labels")
  if (verbose & n_permutations < 50) {
    cli_alert_warning("The number of selected permutations is low (<50). Consider increasing 'n_permutations' for more robust results")
  }
  if (nCores > (detectCores() - 1)) {
    nCores <- detectCores() - 1
    if (verbose) cli_alert_info("Using {nCores} threads")
  }
  # Create random labels
  perm_data <- .randomize_label_ids(object, column_name, n_permutations, seed)

  perm_edge_list <- parLapplier(1:n_permutations, function(i){

    perm_edge_table <- spatnet |>
      mutate(label_label = paste(paste0("Label_", perm_data[spatnet$from, i]),
                                 paste0("Label_", perm_data[spatnet$to, i]),
                                 sep = "-"))

    # Count number of edges per unique label-label pair
    perm_edge_table_stats <- perm_edge_table |>
      count(label_label) |>
      rename(edges = n)

    perm_edge_table_stats <- label_label_template |>
      left_join(y = perm_edge_table_stats |> select(label_label, edges), by = "label_label") |>
      mutate_at("edges", ~tidyr::replace_na(., 0))

    return(perm_edge_table_stats)
    },
    nCores = nCores
  )
  if (verbose) cli_alert_success("Randomized label adjacency calculations complete from {n_permutations} iterations")

  # Compute mean and sd from permutations
  perm_edge_df <- do.call(bind_rows, perm_edge_list)
  perm_res <- perm_edge_df |>
    group_by(label_label) |>
    summarise(.groups = "keep",
              perm_mean = mean(edges),
              perm_sd = stats::sd(edges)) |>
    ungroup()

  # Prepare results output and compute z-scores
  res_out <- full_join(x = edge_table_stats, y = perm_res, by = "label_label") |>
    filter(label_1 != label_2)
  res_out[res_out$perm_sd == 0, "perm_sd"] <- 0.0001  # add pseudocount if sd is 0
  res_out$z_score <- round((res_out$edges - res_out$perm_mean) / (res_out$perm_sd),
                           digits = 3)

  if (verbose) cli_alert_success("Scores calculated for each label pair and returned as output tibble")

  return(res_out)
}


#' Label Assortativity Analysis
#'
#' @description
#' Performs a network assortativity test for each label, describing whether the spots of that
#' label displays a clustered or dispersed spatial pattern.
#'
#' @details
#' This analysis is inspired by the Newman's Assortativity measure by looking at the
#' average degree of spots belonging to the same class or label (specified with \code{column_name}).
#' The hypothesis that we want to address is that spots of the same label will lie next
#' to each other in an aggregated fashion, as compared to being randomly dispered throughout
#' the tissue. As the average degree within a label may change based on the total number of
#' observations, we scale it towards the maximum average degree (k) possible in the network (for
#' Visium, the theoretical max avg k is 6) and the minimum avg k you would expect to see if spots
#' are randomly dispered. With the \code{n_permutations} argument, you can specify the number of
#' iterations for generating a mean of the minimum avg k for each label.
#' The output of this function is a tibble table with entries for each label corresponding to
#' the observed avg k (\code{avg_k}), the mean of the randomized avg k (\code{min_avg_k_mean}), the standard
#' deviation of the randomized avg k (\code{min_avg_k_sd}), and the scaled avg k (\code{avg_k_scaled}), which
#' has been scaled towards the network max avg k (\code{k_max}) and thus goes from 0 (completely randomly
#' dispersed) to 1 (fully connected).
#'
#'
#' @param object A Seurat object
#' @param column_name Column name in metadata corresponding to label ID of the spots.
#' @param n_permutations Integer specifying number of iterations the labels should be randomized
#' [default: 100]. Recommended to increase the number of permutations to >=100 for more robust
#' results. A lower number of permutations will result in high standard deviations and thus more
#' unreliable output
#' @param seed A seed to use for reproducibility [default: 123]
#' @param nCores Number of cores [default: parallel::detectCores() - 1]
#' @param verbose Print messages [default: TRUE]
#'
#' @return A tibble with scores for each label.
#'
#' @import dplyr
#' @import cli
#' @import tibble
#' @importFrom rlang abort
#' @importFrom parallel detectCores
#'
#' @rdname label-assortativity
#' @family spatial-methods
#'
#' @author Lovisa Franzén
#'
#'
#' @examples
#' library(semla)
#'
#' # Read data
#' se <- readRDS(system.file("extdata/mousebrain",
#'                           "se_mbrain",
#'                           package = "semla"))
#'
#' # Generate clusters
#' se <- se |>
#'   NormalizeData() |>
#'   ScaleData() |>
#'   FindVariableFeatures() |>
#'   RunPCA() |>
#'   FindNeighbors(reduction = "pca", dims = 1:30) |>
#'   FindClusters()
#'
#' # Run Label Assortativity Analysis
#' res <- RunLabelAssortativityTest(object = se,
#'                                  column_name = "seurat_clusters",
#'                                  n_permutations = 100,
#'                                  nCores = 2)
#'
#' res |> arrange(desc(avg_k_scaled))
#'
#'
#' @export
RunLabelAssortativityTest <- function(
    object,
    column_name,
    n_permutations = 100,
    nCores = parallel::detectCores() - 1,
    seed = 123,
    verbose = TRUE
){
  # Set global variables to NULL
  avg_k <- from <- k <- label <- label_from <- label_to <- random_label_from <- random_label_to <- NULL

  # checks - import from checks.R and spatial_utils.R
  .check_seurat_object(object)
  .validate_column_name(object = object, column_name = column_name)

  # checks
  if (n_permutations != round(n_permutations)) abort("Invalid input for 'n_permutations', expected a numeric integer.")
  if (n_permutations <= 1) abort("Invalid input for 'n_permutations', needs to be more than 1.")
  
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

  # Start analysis
  if (verbose) cli_h2("Running Label Assortativity Test")
  if (verbose) cli_alert_info("Generating neighborhood adjacency data from observed labels in column '{column_name}'")

  # Generate spatnet
  spatnet <- do.call(bind_rows, GetSpatialNetwork(object))
  spatnet$label_from <- object[[]][spatnet$from, column_name, drop = TRUE] |> as.character()
  spatnet$label_to <- object[[]][spatnet$to, column_name, drop = TRUE] |> as.character()

  # Unique labels to use
  unique_labels <- levels(object[[column_name]][,1])

  # Observed k for each label
  obs_label_k_list <- lapply(unique_labels, function(l){
    spatnet_subset <- subset(spatnet, label_from == l & label_to == l)
    spatnet_subset <- spatnet_subset |>
      group_by(from) |>
      summarize(.groups = "keep", k = n()) |>
      tibble::add_column(label = l)
  })
  obs_label_k_df <- do.call(bind_rows, obs_label_k_list)

  # Calculate average degree per label
  obs_k_stats <- obs_label_k_df |>
    group_by(label) |>
    summarise(.groups = "keep", avg_k = mean(k)) |>
    ungroup()

  obs_k_stats$label <- factor(obs_k_stats$label, levels = unique_labels)
  obs_k_stats <- obs_k_stats |> arrange(label)

  if (verbose) cli_alert_success("Observed label average degree calculations complete")

  # Randomise labels and count k per label
  if (verbose) cli_alert_info("Generating neighborhood adjacency data from randomized labels")
  if (verbose & n_permutations < 25) {
    cli_alert_warning("The number of selected permutations is low (<25). Consider increasing 'n_permutations' for more robust results")
  }
  if (nCores > (detectCores() - 1)) {
    nCores <- detectCores() - 1
    if (verbose) cli_alert_info("Using {nCores} threads")
  }
  # Create random labels
  perm_data <- .randomize_label_ids(object, column_name, n_permutations, seed)

  # Iterate randomized labels
  perm_k_stats_list <- parLapplier(1:n_permutations, function(i){

    spatnet$random_label_from <- perm_data[spatnet$from, i]
    spatnet$random_label_to <- perm_data[spatnet$to, i]

    # Permutation generated k for each label
    perm_label_k_df <- spatnet |>
      mutate(identical = random_label_from == random_label_to) |>
      filter(identical) |>
      mutate(label = random_label_from) |>
      group_by(from, label) |>
      summarize(.groups = "keep", k = n())

    # Calculate average degree per label
    perm_k_stats <- perm_label_k_df |>
      group_by(label) |>
      summarise(.groups = "keep", avg_k = mean(k)) |>
      ungroup() |>
      mutate(label = factor(label, levels = unique_labels)) |>
      arrange(label) |>
      tibble::add_column(iteration = i)
    },
    nCores = nCores
  )
  if (verbose) cli_alert_success("Randomized label adjacency calculations complete from {n_permutations} iterations")

  # Summarize permuted results
  perm_k_stats_df <- do.call(bind_rows, perm_k_stats_list)
  perm_k_res <- perm_k_stats_df |>
    group_by(label) |>
    summarise(.groups = "keep",
              min_avg_k_mean = mean(avg_k),
              min_avg_k_sd = stats::sd(avg_k)) |>
    ungroup()

  if (verbose & max(perm_k_res$min_avg_k_sd) > 0.1) {
    cli_alert_warning("The standard deviation of some of the permuted results is relatively high (>0.1) which may indicate unreliable results.\nA possible solution could be to increase the number of permutations ('n_permutations').")
  }

  # Prepare results output and compute scaled avg k
  network_k_max <- mean(spatnet$nn)
  res_out <- full_join(x = obs_k_stats, y = perm_k_res, by = "label")
  res_out$k_max <- network_k_max
  res_out$avg_k_scaled <- (res_out$avg_k-res_out$min_avg_k_mean)/(network_k_max-res_out$min_avg_k_mean)

  if (verbose) cli_alert_success("Scores calculated for each label and returned as output tibble")
  return(res_out)
}

#' Randomize Label IDs within each sample
#'
#' @param object Seurat object containing label and sample identities for each spot in the metadata.
#' @param column_name Column name in metadata corresponding to label ID of the spots.
#' @param n_permutations Number of times labels should be randomized.
#' @param seed A seed to use for reproducibility.
#'
#' @return A matrix of dimensions NxP with random labels, where N is the number of spots in the
#' Seurat object and P is the number of permutations
#'
#' @noRd
.randomize_label_ids <- function (
    object,
    column_name,
    n_permutations,
    seed
) {

  if (!inherits(seed, what = c("numeric", "integer")) | (length(seed) != 1)) abort("seed should be an integer of length 1")
  set.seed(seed)

  # Set global variables to NULL
  barcode <- sampleID <- NULL

  # Get barcodes and labels
  new_md <- GetStaffli(object)@meta_data |>
    select(barcode, sampleID) |>
    bind_cols(object[[]] |> select(all_of(column_name))) |>
    setNames(nm = c("barcode", "sampleID", "original_labels"))

  # Create random labels
  new_md_groups <- new_md |>
    group_by(sampleID) |>
    group_split()
  perm_data <- do.call(rbind, lapply(new_md_groups, function(x) {
    sapply(1:n_permutations, function(i) {
      x$original_labels[sample.int(nrow(x))]
    })
  }))
  rownames(perm_data) <- colnames(object)

  return(perm_data)
}
