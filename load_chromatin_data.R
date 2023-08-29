#' @include extdata.R
#' @include load_data.R
#'
NULL

#' Load and merge ATAC-seq data from feature matrices
#'
#' Chromatin feature matrices should have features in rows and spots in columns.
#'
#' @details The merging process makes sure that all peaks detected are present in the merged output.
#' This means that if a feature is missing in a certain dataset, the spots in that dataset will
#' be assigned with 0 expression.
#'
#' Spot IDs are renamed to be unique. Usually, the spots are named something similar to:
#' \code{"ACGCCTGACACGCGCT-1", "TACCGATCCAACACTT-1"}
#'
#' Since spot barcodes are shared across datasets, there is a risk that some of the spot IDs
#' will be duplicated after merging. To avoid this, the prefix (e.g. "-1") is replaced by
#' a unique prefix for each loaded matrix: "-1", "-2", "-3", ...
#' Only barcodes in the feature matrices are changed. 
#' Barcodes in the fragmentfiles stay the same.
#'
#' @family pre-process
#'
#' @param samplefiles Character vector of file/directory paths. Paths should specify \code{.h5} or
#' \code{.tsv}/\code{.tsv.gz} files. Alternatively, the paths could specify directories including \code{barcodes.tsv},
#' \code{features.tsv} and \code{matrix.mtx} files.
#' @param fragmentfiles Character vector of file paths. Paths should specify \code{.tsv.gz}
#' @param cells Character vector of wanted cell barcodes for each sample
#' @param verbose Print messages
#'
#' @importFrom methods as
#' @importFrom Seurat Read10X Read10X_h5
#' @importFrom SeuratObject RowMergeSparseMatrices
#' @importFrom glue glue
#' @import cli
#' @import rlang
#' @importFrom dplyr case_when mutate
#' @importFrom tibble tibble
#' @importFrom tools file_ext
#'
#' @return A \code{ChromatinAssay} object 
#'
#' @examples
#'
#' # Load and merge two chromatin feature matrices
#' samples <-
#'   c(
#'     path/to/raw_peak_bc_matrix1.h5,
#'     path/to/raw_peak_bc_matrix2.h5
#'   )
#' 
#' frags <-
#'   c(
#'     path/to/fragments1.tsv.gz,
#'     path/to/fragments2.tsv.gz
#'   )
#' 
#' mergedChromatinData <- LoadChromatinFromMatrix(
#'                         samplefiles = samples, 
#'                         fragmentfiles = frags, 
#'                         sep = c(":", "-")
#'                         )
#'
#' @export

LoadChromatinFromMatrix <- function (
    samplefiles,
    fragmentfiles,
    cells = NULL,
    verbose = TRUE,
    sep = c(":", "-"),
    ...
) {
  
  # Run checks
  if (!is.character(samplefiles)) abort("'samplefiles' must be a character vector.")
  if (!is.character(fragmentfiles)) abort("'fragmentfiles' must be a character vector.")
  checks <- tibble(samplefiles = samplefiles,
                   fragmentfiles = fragmentfiles) |>
    dplyr::mutate(is_samplefiles = case_when(dir.exists(samplefiles) ~ "dir", 
                                      file.exists(samplefiles) ~ "file"),
           is_fragments = case_when(file.exists(fragmentfiles) ~ "fragments"),
           index = paste0(fragmentfiles, ".tbi"))
  if (any(is.na(checks$is_fragments))) abort(c("Invalid path(s):", glue::glue("{checks$fragmentfiles[is.na(checks$is_fragments)]}")))
  if (any(is.na(checks$is_samplefiles))) abort(c("Invalid path(s):", glue::glue("{checks$samplefiles[is.na(checks$is_samplefiles)]}")))
  if (any(is.na(checks$index))) abort(c("Index file(s) don't exist:", glue::glue("{checks$index[is.na(checks$is_fragments)]}")))
  
    # Load expression matrices
  if (verbose) cli_alert_info("Loading matrices:")
  chromatinMats <- lapply(seq_along(samplefiles), function(i) {
    if (checks$is_samplefiles[i] == "dir") {
      # Assumes that the directory contains matrix, barcodes and genes
      if (!file.exists(paste0(samplefiles[i], "/matrix.mtx"))) abort("matrix.mtx file is missing.")
      if (!file.exists(paste0(samplefiles[i], "/barcodes.tsv"))) abort("barcodes.tsv file is missing.")
      if (!file.exists(paste0(samplefiles[i], "/peaks.bed"))) abort("peaks.bed file is missing.")
      
      chromatinMat <- Matrix::readMM(paste0(samplefiles[i], "/matrix.mtx"))
      colnames(chromatinMat) <- read.table(paste0(samplefiles[1], "/barcodes.tsv"))[,1]
      bed <- read.table(paste0(samplefiles[i], "/peaks.bed"))
      rownames(chromatinMat) <- paste0(bed[,1], sep[1], bed[,2], sep[2], bed[,3])
    
      } else if (checks$is_samplefiles[i] == "file") {
      ext <- file_ext(samplefiles[i])
      if (ext == "h5") {
        chromatinMat <- Read10X_h5(samplefiles[i])
      } else if (ext %in% c("tsv", "tsv.gz")) {
        if (!requireNamespace("data.table")) 
          abort(glue("Package {cli::col_br_magenta('data.table')} is required. Please install it with: \n",
                     "install.packages('data.table')"))
        chromatinMat <- read.table(samplefiles[i], header = TRUE, check.names = FALSE) |>
          as.matrix() |>
          as("dgCMatrix")
      } else {
        abort(glue("Invalid file format '{ext}'"))
      }
    }
    if (verbose) cli_alert("Finished loading expression matrix {i}")
    return(chromatinMat)
  })
  
  # check if cells exist
  cells <- lapply(seq_along(along.with = chromatinMats), function(i) {
    if(!is.null(cells[[i]])){
      frag_cells <- cells[[i]]
    }else{
      frag_cells <- colnames(chromatinMats[[i]])
    }
    return(frag_cells)
  })
  
  chromatinMats <- lapply(seq_along(along.with = chromatinMats), function(i){
    chromatinMats[[i]] <- chromatinMats[[i]][, colnames(chromatinMats[[i]]) %in% cells[[i]]]
  })

  # read fragments
  fragmentObj <- lapply(seq_along(fragmentfiles), function(i) {
    fragmentObj <- CreateFragmentObject(
      path = fragmentfiles[i],
      cells = cells[[i]]
    )
    return(fragmentObj)
  })


  # Merge matrices
  if (length(chromatinMats) > 1) {
    if (verbose) {
      cli_text("")
      cli_alert_info("Merging matrices:")
    }
    
    # change cell names
    for (i in seq_along(along.with = chromatinMats)){
      new.names <- gsub(pattern = "-\\d+", replacement = paste0("-", i), x = colnames(chromatinMats[[i]]))
      colnames(chromatinMats[[i]]) <- new.names
      # add the information of the match between old and new names to the fragment objects
      frag_cells <- slot(object = fragmentObj[[i]], name = "cells")
      names(frag_cells) <- new.names
      slot(object = fragmentObj[[i]], name = "cells") <- frag_cells
    }
    
    # direct merging form matrix is not suggested unless the matrices share common features/peak sets.
    # no peak merging is done here.
    # warning message will be given if overlapping ranges detected by CreateChromatinAssay function
    mergedMat <- RowMergeSparseMatrices(mat1 = chromatinMats[[1]], mat2 = chromatinMats[2:length(chromatinMats)])
    if (verbose) cli_alert_success(glue("There are {cli::col_br_blue(nrow(mergedMat))} features and {cli::col_br_magenta(ncol(mergedMat))} ",
                                        "spots in the merged matrix."))
    # do not change cell names in the fragments before, so here unable the validation
    chromatinAssay <- CreateChromatinAssay(counts = mergedMat, 
                                           fragments = fragmentObj,
                                           validate.fragments = FALSE, 
                                           sep = sep,
                                           ...)
  } else {
    if (verbose) cli_alert_info("only 1 feature matrix loaded.")
    if (verbose) cli_alert_success(glue("  There are {cli::col_br_blue(nrow(chromatinMats[[1]]))} features and",
                                        " {cli::col_br_magenta(ncol(chromatinMats[[1]]))} spots in the matrix."))
    chromatinAssay <- CreateChromatinAssay(counts = chromatinMats[[1]], 
                                           fragments = fragmentObj[[1]],
                                           sep = sep,
                                           ...)
  }
  return(chromatinAssay)
}


#' Load and merge chromatin datasets from fragments
#' 
#' create common peak sets for all samples
#'
#' @details The merging process create common peak sets from provided bed files or called peaks from fragments,
#' then feature matrices are calculated using the common peak sets
#'
#' Spot IDs are renamed to be unique. Usually, the spots are named something similar to:
#' \code{"ACGCCTGACACGCGCT-1", "TACCGATCCAACACTT-1"}
#'
#' Since spot barcodes are shared across datasets, there is a risk that some of the spot IDs
#' will be duplicated after merging. To avoid this, the prefix (e.g. "-1") is replaced by
#' a unique prefix for each loaded matrix: "-1", "-2", "-3", ...
#'
#' @family pre-process
#'
#' @param fragmentfiles Character vector of file paths. Paths should specify \code{.tsv.gz}
#' @param cells Character vector of wanted cell barcodes for each sample
#' @param bedfiles Character vector of file paths. Paths should specify \code{.bed}
#' @param mergePeaks Methods to use for merging common peaks. could be chosen form "reduce" or "disjoin"
#' @param verbose Print messages
#'
#' @importFrom methods as
#' @importFrom Seurat Read10X Read10X_h5
#' @importFrom SeuratObject RowMergeSparseMatrices
#' @importFrom glue glue
#' @import cli
#' @import rlang
#' @importFrom dplyr case_when mutate
#' @importFrom tibble tibble
#' @importFrom tools file_ext
#'
#' @return A \code{ChromatinAssay} object 
#'
#'
#' @export
LoadChromatinFromFragments <- function (
    fragmentfiles,
    cells = cells,
    bedfiles = NULL,
    mergePeaks = "reduce",
    verbose = TRUE,
    ...
) {
  
  # Run checks
  if (!is.character(fragmentfiles)) abort("'fragmentfiles' must be a character vector.")
  checks <- tibble(fragmentfiles = fragmentfiles,
                   bedfiles = bedfiles) |>
    dplyr::mutate(is_fragments = case_when(file.exists(fragmentfiles) ~ "fragments"),
           index = paste0(fragmentfiles, ".tbi"))
  if (any(is.na(checks$is_fragments))) abort(c("Invalid path(s):", glue::glue("{checks$fragmentfiles[is.na(checks$is_fragments)]}")))
  if (any(is.na(checks$index))) abort(c("Index file(s) don't exist:", glue::glue("{checks$index[is.na(checks$is_fragments)]}")))
  if (any(is.na(checks$bedfiles))) message(c("Bed file(s) don't exist, calling peaks with MACS2:", glue::glue("{checks$index[is.na(checks$is_fragments)]}")))
  
  # Create fragment objects
  if (verbose) cli_alert_info("Creating fragment objects:")
  fragmentObj <- lapply(seq_along(fragmentfiles), function(i) {
    # load metadata
    # specific for cellranger output
    #meta <- read.table(
    #  file = metadata,
    #  stringsAsFactors = FALSE,
    #  sep = ",",
    #  header = TRUE,
    #  row.names = 1
    #)[-1, ]
    #cells <- rownames(meta)
    
    # create fragment object
    fragmentObj <- CreateFragmentObject(
      path = fragmentfiles[i],
      cells = cells[[i]]
    )
    return(fragmentObj)
  })

  # get all the peak sets
  if(verbose) cli_alert_info("Get peak information")
  gr <- lapply(seq_along(fragmentfiles), function(i){
    if(is.null(bedfiles)){
      if (verbose) cli_alert_info("Call peaks")
      gr <- CallPeaks(fragmentfiles[i])
    }else{
      if(!is.na(bedfiles[i])){
        # Read peak sets
        # with path to bedfiles provided:      
        gr <- read.table(
          file = bedfiles[i],
          col.names = c("chr", "start", "end")) %>% 
          makeGRangesFromDataFrame()
      }else{
        if (verbose) cli_alert_info("Call peaks")
        gr <- CallPeaks(fragmentfiles[i])
      }
    } 
    return(gr)
  })
  
  
  # Calculate feature matrices and merge
  if (length(fragmentObj) > 1) {
    if (verbose) {
      cli_text("")
      cli_alert_info("Merging matrices:")
    }

    # Create a unified set of peaks to quantify in each dataset
    # reduce/disjoin
    
    if(mergePeaks == "reduce"){
      combinedPeaks <- lapply(gr, GenomicRanges::reduce) %>%
        Reduce(GenomicRanges::union, .)
      
    } else if(mergePeaks == "disjoin"){
      combinedPeaks <- lapply(gr, GenomicRanges::disjoin) %>%
        do.call(c, .)
    }
        
    # calculate matrix
    chromatinMats <- lapply(seq_along(fragmentObj), function(i){
      chromatinMat <- FeatureMatrix(
        fragments = fragmentObj[[i]],
        features = combinedPeaks,
        cells = cells[[i]]
      )
      return(chromatinMat)
    })
    
    # change cell names
    for (i in seq_along(along.with = chromatinMats)){
      new.names <- gsub(pattern = "-\\d+", replacement = paste0("-", i), x = colnames(chromatinMats[[i]]))
      colnames(chromatinMats[[i]]) <- new.names
      # add the information of the match between old and new names to the fragment objects
      cells <- slot(object = fragmentObj[[i]], name = "cells")
      names(cells) <- new.names
      slot(object = fragmentObj[[i]], name = "cells") <- cells
    }
    
    # merge matrices
    mergedMat <- RowMergeSparseMatrices(mat1 = chromatinMats[[1]], mat2 = chromatinMats[2:length(chromatinMats)])
    if (verbose) cli_alert_success(glue("There are {cli::col_br_blue(nrow(mergedMat))} features and {cli::col_br_magenta(ncol(mergedMat))} ",
                                        "spots in the merged matrix."))
    chromatinAssay <- CreateChromatinAssay(counts = mergedMat, 
                                           fragments = fragmentObj,
                                           ranges = combinedPeaks,
                                           validate.fragments = FALSE, 
                                           ...)
  } else {
    if (verbose) cli_alert_info("Only 1 fragment file loaded.")
    chromatinMat <- FeatureMatrix(
      fragments = fragmentObj[[1]],
      features = gr[[1]],
      cells = cells[[1]])
    if (verbose) cli_alert_success(glue("  There are {cli::col_br_blue(nrow(chromatinMat))} features and",
                                        " {cli::col_br_magenta(ncol(chromatinMat))} spots in the matrix."))
    chromatinAssay <- CreateChromatinAssay(counts = chromatinMat, 
                                           fragments = fragmentObj[[1]],
                                           ranges = gr[[1]],
                                           validate.fragments = FALSE, 
                                           ...)
  }
  return(chromatinAssay)
}
  
# cells <- read.table(
# file = "~/Dropbox/Postdoc/sequencing/sci_mouse_multiome/data/ALL_1dpi_1/outs/singlecell.csv",
# stringsAsFactors = FALSE,
# sep = ",",
# header = TRUE,
# row.names = 1
# )[-1, ] %>% rownames()
# for cellranger outputs cell names need to be filtered at least by passed_filters above 0
# or errors might appear








#' Load and merge multiple coordinate tables not in standard spaceranger format
#' Not used now, however could be integrated into the ReadChromatinData function for specific samples
#'
#' Load coordinates from \strong{'tissue.tsv'} files and merge them into a \code{tibble}.
#'
#' @details The output will have the same format as the function LoadSpatialCoordinates
#'
#' @family pre-process
#'
#' @param coordinatefiles Character vector of file paths. Paths should specify \code{.tsv} files containing spots positions
#' @param remove_spots_outside_tissue Should spots outside the tissue be removed?
#' @param verbose Print messages
#'
#' @import rlang
#' @import dplyr
#' @import cli
#' @importFrom tibble as_tibble
#'
#' @return An object of class \code{tbl} containing spot coordinates
#'
#'
#' @export
LoadAndProcessCoords <- function(
    spotfiles,
    image_info,
    remove_spots_outside_tissue = TRUE,
    verbose = TRUE
){
  # Run checks
  if (!is.character(spotfiles)) abort("'spotfiles' must be a character vector.")
  checks <- tibble::tibble(spotfiles) |>
    dplyr::mutate(type = dplyr::case_when(file.exists(spotfiles) ~ "file"))
  if (any(is.na(checks$type))) abort(c("Invalid path(s):", glue::glue("{checks$spotfiles[is.na(checks$type)]}")))
  
  # Load data from the tsv file
  if (verbose) cli_alert_info("Loading coordinates:")
  coordinates <- do.call(bind_rows, lapply(seq_along(spotfiles), function(i){
    spot_pos <- read.table(spotfiles[i], sep = "\t", header = T) 
    spot_pos$selected <- as.integer(as.logical(spot_pos$tissue))
    spot_pos$barcode <- paste0(spot_pos$barcode, "-", i)
    
    # Calculate the coordinates according to the scale factors
    coords <- spot_pos |> 
      dplyr::select(all_of(c("barcode", "selected"))) |> 
      dplyr::mutate(y = spot_pos$row, x = spot_pos$col)
    sf <- c(image_info$width[i]/max(coords$x), image_info$height[i]/max(coords$y))
    coords <- coords |> 
      dplyr::mutate(pxl_row_in_fullres = round(coords$y * sf[2]),
             pxl_col_in_fullres = round(coords$x * sf[1]),
             sampleID = i)
    
    # filter spots outside the tissue
    if (remove_spots_outside_tissue) {
      coords <- coords |>
        dplyr::filter(selected == 1)
    }
    
    if (verbose) cli_alert("Finished loading coordinates for sample {i}")
    return(coords)
    
  }))
  if (verbose) cli_alert_info("Collected coordinates for {cli::col_br_magenta(nrow(coordinates))} spots.")
  return(coordinates)
}




  #' Read spatial ATAC-seq data
  #'
  #' This function serves as a wrapper for \code{\link{LoadChromatinFromMatrix}} or 
  #' \code{\link{LoadChromatinFromFragments}}and \code{\link{LoadSpatialCoordinates}} to 
  #' load spaceranger output files and create a Seurat object. The spatial information, 
  #' i.e. images and spot coordinates, are stored inside the tools slot of the 
  #' \code{Seurat} object in \code{Staffli} object.
  #'
  #' @details
  #' \code{ReadChromatinData} takes a \code{data.frame} like table as input that should hold
  #' certain spaceranger output file paths. The table consists of columns "matrix", "fragments", 
  #' "bed",  "imgs", "spotfiles", "json", "scalefactor" and other additional information.
  #'
  #' \itemize{
  #'    \item{"matrix" : (optional) file paths to feature matrices, e.g. \code{peak_bc_matrix.h5}}
  #'    \item{"fragments" : file paths to fragments, e.g. \code{fragments.tsv.gz}}
  #'    \item{"bed" : (optional) file paths to bed like files containing peaks, e.g. \code{peaks.bed}}
  #'    \item{"imgs" : file paths to images, e.g. \code{tissue_hires_image.jpg}}
  #'    \item{"spotfiles" : file paths to spot coordinate CSV files \code{tissue_positions_list.csv}}
  #'    \item{"json" : file paths to scale factor JSON files, e.g. \code{scalefactors_json.json}. 
  #'    It is also possible to provide custom scale factors (more info below).}
  #' }
  #' 
  #' @section Custom scale factors:
  #' You can provide an additional column named 'scalefactor' with custom scale factor values.
  #' The values should be between 0 to 1, where a 1 would correspond to the original H&E image
  #' used as input for spaceranger. For instance, if your image was scaled from a height of 30,000
  #' pixels to a height of 3,000 pixels (with the same aspect ratio), the scalefactor would be 0.1.
  #'
  #' @section Load data outside tissue section:
  #' Sometimes it can be useful to load data for all spots in a 10x Visium dataset, if you
  #' need to explore transcripts captured outside of the tissue. In this case, you can
  #' provide paths to the \code{raw_feature_bc_matrix.h5} files in the spaceranger output folders
  #' and set \code{remove_spots_outside_tissue = FALSE}.
  #'
  #' @section Filter data:
  #' If you want to filter out spots and features, you can pass the \code{min.cells} and
  #' \code{min.features} parameters (see \code{\link{CreateSeuratObject}} for more details);
  #' however, it is recommended to use the \code{\link{SubsetSTData}} function for filtering
  #' after the object has been created.
  #'
  #' @section Additional annotations:
  #' You can add a column named 'annotation_files' with paths to CSV files containing 
  #' annotations. These files should have the same structure as the CSV the files exported 
  #' by Loupe Browser. One column with the barcode IDs and one column with the annotation labels.
  #' The annotations will be added to the \code{meta.data} slot of the returned \code{\link{Seurat}} object.
  #'
  #' @family pre-process
  #'
  #' @param infoTable A \code{data.frame} or \code{tbl} with paths to spaceranger output files
  #' @param assay Assay name (default = "ATAC")
  #' @param remove_spots_outside_tissue Should spots outside the tissue be removed or not
  #' @param verbose Print messages
  #' @param ... Parameters passed to \code{\link{CreateSeuratObject}}
  #'
  #' @import cli
  #' @import rlang
  #' @import dplyr
  #' @importFrom glue glue
  #' @importFrom Seurat CreateSeuratObject
  #' @importFrom tidyr uncount
  #'
  #' @return A \code{\link{Seurat}} object with additional spatial information stored in
  #' a \code{Staffli} object
  #'
  #' @export
  #'
  ReadChromatinData <- function (
    infoTable,
    assay = "ATAC",
    remove_spots_outside_tissue = TRUE,
    verbose = TRUE,
    ...
  ) {

    
    if (verbose) cli_h2("Reading spatial chromatin data")
    
    # Check infoTable
    if (!all(c("fragments", "imgs", "spotfiles") %in% colnames(infoTable)))
      abort("One or several of 'fragments', 'imgs' and 'spotfiles' are missing from infoTable.")
    if (!"matrix" %in% colnames(infoTable)) {
      if (!"bed" %in% colnames(infoTable)) cli_alert_info("A 'bed' file is missing from infoTable and is suggested to be provided to calculate the feature matrices. Calling peaks from the fragments.")
      else cli_alert_info(glue("Feature matrices not provided, calculating from fragment files."))
    }
    if (!any(c("json", "scalefactor") %in% colnames(infoTable)))
      abort("One of 'json' or 'scalefactor' columns needs to be provided")
    if (!any(class(infoTable) %in% c("tbl", "data.frame"))) abort(glue("Invalid class '{class(infoTable)}' of 'infoTable'."))
    if (nrow(infoTable) == 0) abort(glue("'infoTable' is empty."))
    if (!all(sapply(infoTable[, c("fragments", "imgs", "spotfiles")], class) %in% "character")) abort(glue("Invalid column classes in 'infoTable'. Expected 'character' vectors."))
    if ("scalefactor" %in% colnames(infoTable)) {
      if (!inherits(infoTable[, "scalefactor", drop = TRUE], what = c("numeric", "integer")))
        abort(glue("Invalid column class for 'scalefactor'. Expected a 'numeric'."))
      if (!all(between(x = infoTable[, "scalefactor", drop = TRUE], 0.001, 1)))
        abort(glue("Scalefactors need to be in the range 0.001-1."))
    }
    missing_files <- infoTable |>
      dplyr::select(all_of(c("fragments", "spotfiles")), contains("json")) |>
      dplyr::mutate(across(everything(), ~ file.exists(.x))) |>
      summarize(across(everything(), ~ any(!.x))) |>
      unlist()
    if (any(missing_files)) {
      checks <- missing_files[missing_files]
      sapply(names(checks), function(n) {
        abort(glue("Missing file(s) in the '{n}' column."))
      })
    }
    if ("json" %in% colnames(infoTable)) {
      if (!all(file.exists(infoTable |> pull(all_of("json"))))) {
        abort(glue("Missing json file(s)."))
      }
    }
    
    
    # adapted from AddAnnotationCSV file, remove the possible NO_BARCODE row
    # Check if an annotation_files column was provided
    if ("annotation_files" %in% colnames(infoTable)) {
      annfiles <- infoTable |> pull(all_of("annotation_files"))
      if (verbose) cli_alert_info("Loading annotatons from CSV files")
      ann <- data.frame()
      for(i in seq_along(along.with = annfiles)){
        if (is.na(annfiles[i]))
          next
        if (!file.exists(annfiles[i]))
          abort(glue("Invalid path '{annfiles[i]}'. File doesn't exist."))
        if (file_ext(annfiles[i]) != "csv")
          abort("Invalid file format. Expected a CSV file.")
        sample_ann <- read.csv(annfiles[i], check.names = FALSE)
        if (!colnames(sample_ann)[1] == "barcode") {
          sample_ann <- sample_ann |> 
            rename(barcode = 1)
        }
        sample_ann <- sample_ann[grep(pattern = "-\\d+", sample_ann$barcode),] |> 
          dplyr::mutate(barcode = gsub(pattern = "-\\d+", replacement = paste0("-", i), x = barcode))
        ann <- bind_rows(ann, sample_ann)
      }
      ann <- ann |> data.frame(row.names = 1, check.names = FALSE)
      infoTable <- infoTable |> 
        dplyr::select(-all_of("annotation_files"))
      add_annotations <- TRUE
    } else {
      add_annotations <- FALSE
    }
    
    # Keep additional infoTable columns
    additionalMetaData <- infoTable |> dplyr::select(-any_of(c("fragments", "bed", "matrix", "imgs", "spotfiles", "json", "scalefactor")))
    if(ncol(additionalMetaData) < 1){
      additionalMetaData <- NULL
    }
    
    # Read image info
    image_info <- LoadImageInfo(infoTable$imgs)
    
    # Read spot coordinates
    # The standard spaceranger output 
    coordinates <- LoadSpatialCoordinates(coordinatefiles = infoTable$spotfiles,
                                          remove_spots_outside_tissue = remove_spots_outside_tissue,
                                          verbose = verbose)
    
    
    # Create Seurat meta data
    if (!is.null(additionalMetaData)) {
      metaData <- coordinates |>
        dplyr::select(barcode) |>
        bind_cols(
          additionalMetaData |>
            dplyr::mutate(count = table(coordinates$sampleID) |> as.numeric()) |>
            uncount(count)
        ) |>
        data.frame(row.names = 1)
    } else {
      metaData <- NULL
    }
    
    metaData <- metaData[coordinates$barcode, , drop = FALSE]
    
    # The following check is not necessarily needed since loading chromatin data needs a cell list as input,
    # the loading functions do filtering with the list.
    # # Make sure that coordinates and expression matrix are compatible
    # if (!all(colnames(mergedMat) %in% coordinates$barcode)) {
    #   cli_alert_danger("{sum(!colnames(mergedMat) %in% coordinates$barcode)} spots found in the merged expression matrix that are missing in the coordinate files.")
    #   cli_alert_danger("Removing {sum(!colnames(mergedMat) %in% coordinates$barcode)} spots from the merged expression matrix")
    #   mergedMat <- mergedMat[, colnames(mergedMat) %in% coordinates$barcode]
    # }
    # if (!all(coordinates$barcode %in% colnames(mergedMat))) {
    #   cli_alert_danger("{sum(!coordinates$barcode %in% colnames(mergedMat))} spots found in coordinate files that are missing in the expression data.")
    #   cli_alert_danger("Removing {sum(!coordinates$barcode %in% colnames(mergedMat))} spots from the meta data")
    #   metaData <- metaData[colnames(mergedMat), , drop = FALSE]
    # }
    # if (verbose) cli_h3(text = "Creating `Seurat` object")
    # if (verbose) cli_alert_success("Expression matrices and coordinates are compatible")
    
    # Read scale factors
    if ("scalefactor" %in% colnames(infoTable)) {
      if (verbose) cli_alert_info("Using custom scalefactor(s) from infoTable")
      scalefactors <- infoTable |> 
        dplyr::select(scalefactor) |> 
        dplyr::rename(custom_scalef = scalefactor) |> 
        dplyr::mutate(sampleID = paste0(1:n()))
      # Add full res image dimensions to image data
      image_info <- image_info |> left_join(y = scalefactors, by = "sampleID") |>
        dplyr::mutate(full_width = width/custom_scalef,
               full_height = height/custom_scalef) |>
        dplyr::select(all_of(c("format", "width", "height", "full_width", "full_height", "colorspace", "filesize", "density", "sampleID", "type")))
    } else {
      scalefactors <- LoadScaleFactors(infoTable$json)
      # Add full res image dimensions to image data
      image_info <- UpdateImageInfo(image_info, scalefactors)
    }
    

    
    # check coordinates in case the they are outside the image and make errors when plotting
    check_coordinates <- coordinates |>
      group_by(sampleID) |>
      dplyr::mutate(check = case_when(between(x = !! sym("pxl_col_in_fullres"),
                                       0, image_info[cur_group_id(), "full_width", drop = TRUE]) &
                               between(x = !! sym("pxl_row_in_fullres"),
                                       0, image_info[cur_group_id(), "full_height", drop = TRUE])
                               ~ "inside",
                               TRUE ~ "outside"))
    
    if (any(check_coordinates$check == "outside")) {
      check_coordinates <- check_coordinates |> summarize(nOutside = sum(check == "outside"))
      for (ID in check_coordinates$sampleID) {
        cli_alert_danger("Found {check_coordinates |> filter(sampleID == ID) |> pull(nOutside)} spot(s) with coordinates outside of the image for sample {ID}. Coordinates and images are incompatible.")
      }
    }
    
    
    # load cells in tissue only
    cells <- coordinates |> 
      dplyr::mutate(raw_barcode = gsub(pattern = "-\\d+", replacement = "-1", x = coordinates$barcode)) |>
      dplyr::filter(selected == 1) |> dplyr::select(sampleID, raw_barcode)
    cells <- split(cells$raw_barcode, cells$sampleID)
    
    
    # Read chromatin data
    if("matrix" %in% colnames(infoTable)){
      chromatinAssay <- LoadChromatinFromMatrix(samplefiles = infoTable$matrix,
                                                fragmentfiles = infoTable$fragments,
                                                cells = cells,
                                                verbose = verbose,
                                                sep = c(":", "-"))
    } else {
      chromatinAssay <- LoadChromatinFromFragments(fragmentfiles = infoTable$fragments,
                                                bedfiles = infoTable$bed,
                                                cells = cells,
                                                verbose = verbose)
    }
    
     

    
    # Create a Seurat object from expression matrix
    object <- CreateSeuratObject(chromatinAssay,
                                 assay = assay,
                                 meta.data = metaData,
                                 ...)
    if (verbose) cli_alert_info("Created `Seurat` object")
    
    # Create a Staffli object
    staffli_object <- CreateStaffliObject(imgs = infoTable$imgs,
                                          meta_data = coordinates |>
                                            ungroup() |> 
                                            dplyr::select(all_of(c("barcode", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID"))),
                                          image_info = image_info,
                                          scalefactors = scalefactors)
    if (verbose) cli_alert_info("Created `Staffli` object")
    
    # Place Staffli object inside the tools slot of the Seurat object
    object@tools$Staffli <- staffli_object
    
    # Add additional annotations if available
    if (add_annotations) {
      object <- AddMetaData(object, metadata = ann)
    }
    
    if (verbose) cli_alert_success(glue("Returning a `Seurat` object with {cli::col_br_blue(nrow(object))}",
                                        " features and {cli::col_br_magenta(ncol(object))} spots"))
    return(object)
  }
  