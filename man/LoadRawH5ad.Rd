% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing_h5ad.R
\name{LoadRawH5ad}
\alias{LoadRawH5ad}
\title{Load RawH5ad data}
\usage{
LoadRawH5ad(h5ad_path, assay = "Spatial")
}
\arguments{
\item{h5ad_path}{h5ad path that contains counts matrix,
gene names, meta information, and spatial informations.}

\item{assay}{Name of assay to associate spatial data to}
}
\value{
A \code{\link{Seurat}} object
}
\description{
Load RawH5ad data
}
\concept{preprocessing}
