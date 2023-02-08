#' @include generics.R
#'
NULL

#' Merge 10x Visium data
#'
#' Merges two or more Seurat objects containing SRT data while making sure that the
#' spatial data (images and spot coordinates) are handled correctly.
#'
#' NB: If you use the generic \code{\link{merge}} function you will not be able
#' to use any of the `semla` visualization methods on the output object as
#' the `Staffli` object will be broken.
#'
#' @param x A `Seurat` object
#' @param y A `Seurat` object or a list of `Seurat` objects
#' @param merge_data Merge the data slots instead of just merging the counts
#' (which requires re-normalization); this is recommended if the same normalization
#' approach was applied to all objects. See \code{\link{merge}} for details.
#' @param merge_dr Merge specified DimReducs that are present in all objects; will
#' only merge the embeddings slots for the first N dimensions that are shared across
#' all objects. See \code{\link{merge}} for details.
#' @param project \code{\link{Project}} name for the `Seurat` object
#'
#' @importFrom dplyr select mutate group_split ungroup cur_group_id group_by bind_rows
#' @importFrom tibble as_tibble
#' @importFrom tidyr separate unite drop_na
#' @importFrom Seurat RenameCells
#' @importFrom rlang warn
#'
#' @return A merged `Seurat` object
#'
#' @family subset-and-merge
#' @rdname merge
#'
#' @author Ludvig Larsson
#'
#' @examples
#' se_mbrain <- readRDS(system.file("extdata",
#'                                  "mousebrain/se_mbrain",
#'                                  package = "semla"))
#' se_mcolon <- readRDS(system.file("extdata",
#'                                  "mousecolon/se_mcolon",
#'                                  package = "semla"))
#'
#' se_mbrain
#' se_mcolon
#'
#' # Merge a mousebrain dataset with two mousecolon datasets
#' se_merged <- MergeSTData(x = se_mbrain, y = se_mcolon)
#' se_merged
#'
#' # Plot H&E images
#' ImagePlot(se_merged |> LoadImages(verbose = FALSE))
#'
#' @export
#'
MergeSTData <- function (
    x,
    y,
    merge_data = TRUE,
    merge_dr = NULL,
    project = "SeuratProject"
) {

  # Set global variables to NULL
  uniqueID <- barcode <- new_sampleID <- old_names <- sampleID <- NULL

  # Check Seurat object
  .check_seurat_object(x)
  if (is.list(y)) {
    for (obj in y) {
      .check_seurat_object(obj)
    }
  } else {
    .check_seurat_object(y)
    y <- list(y)
  }

  # Obtain Staffli objects
  st_x <- GetStaffli(x)
  if (is.list(y)) {
    st_y <- lapply(y, GetStaffli)
  } else {
    st_y <- list(y)
  }

  # Create new names
  se_objects <- c(x, y); rm(x); rm(y)
  st_objects <- c(st_x, st_y); rm(st_x); rm(st_y)
  mergedMetaData <- do.call(rbind, lapply(seq_along(se_objects), function(i) {
    obj <- se_objects[[i]]
    st_obj <- st_objects[[i]]
    mergedSampleMetaData <- cbind(obj@meta.data, st_obj@meta_data) |>
      as_tibble() |>
      select(barcode, sampleID) |>
      mutate(sample = i, uniqueID = paste0(i, "_", sampleID))
    return(mergedSampleMetaData)
  })) |>
    group_by(uniqueID) |>
    mutate(new_sampleID = cur_group_id(), old_names = barcode) |>
    separate(barcode, sep = "-", into = c("barcode", NA)) |>
    unite(col = "barcode", barcode, new_sampleID, sep = "-", remove = FALSE) |>
    ungroup() |>
    select(sample, sampleID, barcode, new_sampleID, old_names)

  mergedMetaData_split <- mergedMetaData |>
    group_by(sample) |>
    group_split()

  # Rename spots in Seurat object
  se_objects <- lapply(seq_along(se_objects), function(i) {
    obj <- se_objects[[i]]
    obj <- RenameCells(obj, new.names = mergedMetaData_split[[i]]$barcode, for.merge = TRUE)
    return(obj)
  })

  # Merge Seurat data
  object_merged <-
    merge(
      x = se_objects[[1]],
      y = se_objects[2:length(se_objects)],
      merge.data = merge_data,
      merge.dr = merge_dr,
      project = project
    )

  # Check that version matches
  versions <- sapply(st_objects, function(x) paste0(x@version))
  if (length(unique(versions)) > 1) warn(glue("Different versions; {paste(versions, collapse = ', '')} have been used to process the data"))

  # and check that the same heights have been used
  heights_check <- sapply(st_objects, function(x) x@image_height)
  if (length(unique(heights_check)) > 1) {
    warn(glue("Different heights have been used for the different objects; {paste(heights_check, collapse = ", ")}.",
              "Any loaded images will be removed and a default value of 400 pixels in height will be set."))
    image_height <- 400
  } else {
    image_height  <- unique(heights_check)
  }

  # Merge meta data
  spot_coordinates <- do.call(bind_rows, lapply(st_objects, function(x) x@meta_data)) |>
    drop_na()
  spot_coordinates$barcode <- mergedMetaData$barcode
  spot_coordinates$sampleID <- mergedMetaData$new_sampleID

  # Check that all objects have images
  imgs_class <- sapply(st_objects, function(x) class(x@imgs))
  if (all(imgs_class == "character")) {
    imgs <- sapply(st_objects, function(x) x@imgs)
  } else {
    imgs <- NULL
  }

  # Merge images
  available_images <- Reduce(intersect, sapply(st_objects, function(st_obj) {
    st_obj@rasterlists |> names()
  }))
  if (!is.null(available_images)) {
    rasterlists <- setNames(lapply(available_images, function(type) {
      Reduce(c, lapply(st_objects, function(st_obj) {
        st_obj@rasterlists[[type]]
      }))
    }), nm = available_images)
  }

  # Merge image_info
  image_info <- do.call(bind_rows, lapply(st_objects, function(st_obj) {
    st_obj@image_info
  })) |>
    mutate(sampleID = paste0(unique(mergedMetaData$new_sampleID)))

  # Merge scalefactors
  scalefactors <- do.call(bind_rows, lapply(st_objects, function(st_obj) {
    st_obj@scalefactors
  })) |>
    mutate(sampleID = paste0(unique(mergedMetaData$new_sampleID)))

  # Create Staffli object
  st_object_merged <-
    CreateStaffliObject(
      imgs = imgs,
      meta_data = spot_coordinates,
      image_height = image_height,
      image_info = image_info,
      scalefactors = scalefactors
    )

  # Add images
  if (!is.null(available_images)) {
    st_object_merged@rasterlists <- rasterlists
  }

  # Return Seurat object
  object_merged@tools$Staffli <- st_object_merged
  return(object_merged)
}


#' Subset 10x Visium data
#'
#' Subset a `Seurat` object while making sure that the spatial data
#' (images and spot coordinates) are handled correctly.
#'
#' If you use the default \code{\link{subset}} function you will most likely not be able
#' to use any of the `semla` visualization methods on the output object as
#' the `Staffli` object will be broken. The exception is when filtering is only
#' done at the feature level.
#'
#' @param object A Seurat object
#' @param expression Logical expression indicating features/variables to keep
#' @param spots A vector of spots to keep
#' @param features A vector of features to keep
#' @param idents A vector of identity classes to keep
#'
#' @importFrom dplyr mutate group_by ungroup filter pull arrange
#' @importFrom tidyr unite separate
#' @importFrom Seurat WhichCells RenameCells
#' @importFrom rlang enquo
#'
#' @return A filtered Seurat object
#'
#' @family subset-and-merge
#' @rdname subset
#'
#' @author Ludvig Larsson
#'
#' @examples
#'
#' se_mbrain <- readRDS(system.file("extdata", "mousebrain/se_mbrain", package = "semla"))
#'
#' # Create Seurat object
#' se_mbrain <- se_mbrain |> LoadImages()
#' se_mbrain
#'
#' # Subset by spot IDs (first 100)
#' se_fewspots <- SubsetSTData(se_mbrain, spots = colnames(se_mbrain)[1:1000])
#' se_fewspots
#'
#' # Subset by feature IDs (first 50)
#' se_fewgenes <- SubsetSTData(se_mbrain, features = rownames(se_mbrain)[1:50])
#' se_fewgenes
#'
#' # Subset using an expression
#' se_filtered <- SubsetSTData(se_mbrain, expression = nFeature_Spatial > 20)
#' se_filtered
#'
#' @export
#'
SubsetSTData <- function (
    object,
    expression,
    spots = NULL,
    features = NULL,
    idents = NULL
) {

  # Set global variables to NULL
  barcode <- sampleID <- NULL

  # Check that a Staffli object is present
  .check_seurat_object(object)

  # Obtain Staffli object
  st_object <- GetStaffli(object)
  all_sampleIDs <- unique(st_object@meta_data$sampleID)

  # Subset Staffli object and rearrange rows
  if (!is.null(spots)) {
    st_meta_data <- st_object@meta_data |>
      filter(barcode %in% spots) |>
      arrange(sampleID)
    st_object@meta_data <- st_meta_data
    spots <- st_meta_data$barcode
  }

  if (!missing(x = expression)) {
    expression <- enquo(expression)
    spots <- WhichCells(object, cells = spots, idents = idents, expression = expression)
  }
  object <- subset(object, features = features, cells = spots, idents = idents)

  # Subset ST meta data
  st_meta_data <- st_object@meta_data |>
    filter(barcode %in% colnames(object))

  # Get a list of remaining samples
  remaining_samples <- st_meta_data |> pull(sampleID) |> unique()
  new_sampleIDs <- setNames(seq_along(remaining_samples), nm = remaining_samples)

  # rename samples if some samples were dropped
  if (length(all_sampleIDs) > length(remaining_samples)) {
    st_meta_data <- st_meta_data |>
      group_by(sampleID) |>
      mutate(sampleID = cur_group_id()) |>
      separate(barcode, sep = "-", into = c("barcode", NA)) |>
      unite(col = "barcode", barcode, sampleID, sep = "-", remove = FALSE) |>
      ungroup()
    # rename spots in Seurat object
    object <- RenameCells(object, new.names = st_meta_data$barcode)
  }

  # Add st_meta_data to STaffli object
  st_object@meta_data <- st_meta_data

  # Subset each slot in Staffli object
  if (length(st_object@imgs) > 0) {
    st_object@imgs <- st_object@imgs[remaining_samples]
  }
  if (length(st_object@rasterlists) > 0) {
    rl <- st_object@rasterlists
    rl <- lapply(rl, function(ls) {
      setNames(ls[remaining_samples], nm = new_sampleIDs)
    })
    st_object@rasterlists <- rl
  }

  # Subset image data
  st_object@image_info <- st_object@image_info |>
    filter(sampleID %in% remaining_samples) |>
    mutate(sampleID = new_sampleIDs[as.character(sampleID)])

  # Subset scalefactors
  st_object@scalefactors <- st_object@scalefactors |>
    filter(sampleID %in% remaining_samples) |>
    mutate(sampleID = new_sampleIDs[sampleID])

  object@tools$Staffli <- st_object
  return(object)
}
