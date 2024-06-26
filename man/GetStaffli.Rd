% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/staffli_object.R
\docType{methods}
\name{GetStaffli}
\alias{GetStaffli}
\alias{GetStaffli,Seurat-method}
\title{Method used to extract a \code{Staffli} object from the tools slot of a
\code{Seurat} object}
\usage{
GetStaffli(object)

\S4method{GetStaffli}{Seurat}(object)
}
\arguments{
\item{object}{A \code{Seurat} object}
}
\value{
A \code{Staffli} object
}
\description{
Method used to extract a \code{Staffli} object from the tools slot of a
\code{Seurat} object
}
\examples{

# Load example data
se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla")) |> 
   LoadImages()

# Fetch Staffli object from a Seurat object
staffli_object <- GetStaffli(se_mbrain)

}
