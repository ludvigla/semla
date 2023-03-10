% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/staffli_object.R
\name{CreateStaffliObject}
\alias{CreateStaffliObject}
\title{Create a Staffli object}
\usage{
CreateStaffliObject(
  imgs = NULL,
  meta_data,
  image_height = 400,
  image_info,
  scalefactors
)
}
\arguments{
\item{imgs}{Character vector specifying paths to images in JPG, PNG or TIF format}

\item{meta_data}{Spot-level metadata to add to the `Staffli` object. This should be a `tbl` with
required columns 'barcode' representing the spot IDs, 'pxl_col_in_fullres' and 'pxl_row_in_fullres'
which specifies the 10x Visium array coordinates and a 'sampleID' column with sample IDs}

\item{image_height}{Specifies the height of the scaled images in pixels [default: 400 pixels]}

\item{image_info}{a tibble with image information}

\item{scalefactors}{a tibble with scalefactors sued to transform coordinates from
the original image space to the downscaled images}
}
\description{
Create a Staffli object from a set of images and associated spot coordinates
}
