#' Export data for FeatureViewer
#'
#' @param object A `Seurat` object created with STUtility2
#' @param sampleIDs An integer vector specifying
#' @param outdir A character vector specifying the path to an existing directory
#' @param overwrite Overwrite files if they already exists
#' @param verbose Print messages
#'
#' @examples
#' \dontrun{
#' libary(STUtility2)
#'
#' se_mbrain <- readRDS(system.file("extdata/mousebrain",
#'                                  "se_mbrain",
#'                                  package = "STUtility2"))
#'
#' # Export files to current working directory
#' ExportDataForViewer(se_mbrain, outdir = "./")
#' }
#'
#' @export
#'
ExportDataForViewer <- function (
  object,
  sampleIDs = NULL,
  outdir = "./",
  overwrite = FALSE,
  verbose = TRUE
) {

  # Check object
  .check_seurat_object(object)

  # Check directory
  stopifnot(inherits(x = outdir, what = "character"),
            length(outdir) == 1)
  if (!dir.exists(outdir)) abort(glue("oudir '{outdir}' does not exist"))
  if (!file.access(names = outdir, mode = 2) == 0) abort(glue("You do not have permisssion to write to '{outdir}' "))

  # Check sampleIDs
  sampleIDs <- sampleIDs %||% {
    GetStaffli(object)@image_info$sampleID |> as.integer()
  }
  stopifnot(inherits(sampleIDs, what = c("integer", "numeric")),
            length(sampleIDs) > 0)
  available_sampleIDs <- GetStaffli(object)@image_info$sampleID |> as.integer()
  if (!all(sampleIDs %in% available_sampleIDs)) {
    abort(glue("Invalid sampleIDs. Possible sampleIDs are: {paste(available_sampleIDs, collapse = ', ')}"))
  }

  if (verbose) cli_alert_info("Attempting to tile H&E image(s)")
  imgs <- GetStaffli(object)@imgs[sampleIDs]

  for (i in seq_along(imgs)) {
    if (!file.exists(imgs[i])) {
      abort(glue("{imgs[i]} is not a valid path. Add a valid path to @imgs slot in 'Staffli' object"))
    }
    dirs <- TileImage(im = image_read(imgs[i]), outpath = outdir, sampleID = sampleIDs[i], overwrite = overwrite, verbose = verbose)
    datapath <- dirs$datapath
    if (verbose) cli_alert("  Exporting Visium coordinates")
    export_coordinates(object = object, sampleNumber = sampleIDs[i], outdir = datapath, overwrite = overwrite, verbose = FALSE)
  }

  return(datapath)
}
