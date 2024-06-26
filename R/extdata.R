#' Visium mouse brain dataset
#'
#' A Visium dataset obtained from a coronal tissue section of a mouse brain.
#' \itemize{
#'    \item{\strong{mousebrain/filtered_feature_bc_matrix.h5}:\cr light weight gene expression matrix with 188 genes in hdf5 format filtered to include spots under the tissue}
#'    \item{\strong{mousebrain/spatial/tissue_lowres_image.jpg}: \cr H&E image (600x565) pixels}
#'    \item{\strong{mousebrain/spatial/tissue_positions_list.csv}: \cr CSV file with spot coordinates}
#'    \item{\strong{mousebrain/spatial/scalefactors_json.json}: \cr JSON file with scalefactors}
#'    \item{\strong{mousebrain/se_mbrain}: \cr Seurat object stored as an .Rds file with the top 100 most variable features}
#' }
#'
#' @name mbrain_dataset
#' @family datasets
NULL

#' Visium mouse colon dataset
#'
#' A Visium dataset obtained from a "swiss roll" of a mouse colon
#' \itemize{
#'    \item{\strong{mousecolon/filtered_feature_bc_matrix.h5}:\cr light weight gene expression matrix with 188 genes in hdf5 format filtered to include spots under the tissue}
#'    \item{\strong{mousecolon/spatial/tissue_lowres_image.jpg}: \cr H&E image (600x541) pixels}
#'    \item{\strong{mousecolon/spatial/tissue_positions_list.csv}: \cr CSV file with spot coordinates}
#'    \item{\strong{mousecolon/spatial/scalefactors_json.json}: \cr JSON file with scalefactors}
#'    \item{\strong{mousecolon/se_mcolon}: \cr Seurat object stored as an .Rds file with the top 100 most variable features}
#'    \item{\strong{mousecolon/galt_spots.csv}: \cr CSV file with GALT spot barcodes}
#' }
#' 
#' @references 
#' \doi{10.1038/s41467-022-28497-0}
#'
#' @name mcolon_dataset
#' @family datasets
NULL
