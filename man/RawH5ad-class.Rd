% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects_h5ad.R
\docType{class}
\name{RawH5ad-class}
\alias{RawH5ad-class}
\alias{RawH5ad}
\title{The raw h5ad class}
\description{
The raw h5ad class represents spatial information from the any plate platform from h5ad
}
\section{Slots}{

\describe{
\item{\code{coordinates}}{...}
}}

\section{Slots}{


\describe{
\item{\code{assay}}{Name of assay to associate image data with; will give this image
priority for visualization when the assay is set as the active/default assay
in a \code{Seurat} object}

\item{\code{key}}{A one-length character vector with the object's key; keys must
be one or more alphanumeric characters followed by an underscore
\dQuote{\code{_}} (regex pattern
\dQuote{\code{^[a-zA-Z][a-zA-Z0-9]*_$}})}
}
}

\concept{spatial}
