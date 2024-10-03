#' @importFrom Rcpp evalCpp
#' @importFrom Matrix colSums rowSums colMeans rowMeans
#' @importFrom methods setClass setOldClass setClassUnion slot
#' slot<- setMethod new signature slotNames is setAs setValidity .hasSlot
#' @importFrom SeuratObject Cells GetImage GetTissueCoordinates Radius RenameCells
#' @importClassesFrom Matrix dgCMatrix
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions for raw h5ad
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The raw h5ad class
#'
#' The raw h5ad class represents spatial information from the any plate platform from h5ad
#'
#' @inheritSection SeuratObject::SpatialImage Slots
#' @slot coordinates ...
#' @concept spatial
#'
RawH5ad <- setClass(
  Class = 'RawH5ad',
  contains = 'SpatialImage',
  slots = list(
    'coordinates' = 'data.frame'
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get Cell Names
#'
#' @inheritParams SeuratObject::Cells
#'
#' @rdname Cells
#' @concept objects
#' @concept spatial
#' @method Cells RawH5ad
#' @export
#'
#' @seealso \code{\link[SeuratObject:Cells]{SeuratObject::Cells}}
#'
Cells.RawH5ad <- function(x, ...) {
  return(rownames(x = GetTissueCoordinates(object = x)))
}

#' Get Image Data
#'
#' @inheritParams SeuratObject::GetImage
#'
#' @rdname GetImage
#' @method GetImage RawH5ad
#' @concept objects
#' @concept spatial
#' @export
#'
#' @seealso \code{\link[SeuratObject:GetImage]{SeuratObject::GetImage}}
#'
GetImage.RawH5ad <- function(
    object,
    mode = c('grob', 'raster', 'plotly', 'raw'),
    ...
) {
  mode <- match.arg(arg = mode)
  return(NullImage(mode = mode))
}

#' Get Tissue Coordinates
#'
#' @inheritParams SeuratObject::GetTissueCoordinates
#'
#' @rdname GetTissueCoordinates
#' @method GetTissueCoordinates RawH5ad
#' @concept objects
#' @concept spatial
#' @export
#'
#' @seealso \code{\link[SeuratObject:GetTissueCoordinates]{SeuratObject::GetTissueCoordinates}}
#'
GetTissueCoordinates.RawH5ad <- function(object, ...) {
  coords <- slot(object = object, name = 'coordinates')
  colnames(x = coords) <- c('x', 'y')
  # coords$y <- -rev(x = coords$y) + 1
  # coords$y <- FlipCoords(x = coords$y)
  coords$cells <- rownames(x = coords)
  return(coords)
}

#' Get Spot Radius
#'
#' @inheritParams SeuratObject::Radius
#'
#' @rdname Radius
#' @concept objects
#' @concept spatial
#' @method Radius RawH5ad
#' @export
#'
#' @seealso \code{\link[SeuratObject:Radius]{SeuratObject::Radius}}
#'
Radius.RawH5ad <- function(object, ...) {
  return(0.005)
}

#' Rename Cells in an Object
#'
#' @inheritParams SeuratObject::RenameCells
#'
#' @rdname RenameCells
#' @concept objects
#' @method RenameCells RawH5ad
#' @export
#'
#' @seealso \code{\link[SeuratObject:RenameCells]{SeuratObject::RenameCells}}
#'
RenameCells.RawH5ad <- function(object, new.names = NULL, ...) {
  return(RenameCells.VisiumV1(object = object, new.names = new.names))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @method [ RawH5ad
#' @concept objects
#' @export
#'
"[.RawH5ad" <- function(x, i, ...) {
  return(subset(x = x, cells = i, ...))
}


#' @method dim RawH5ad
#' @concept objects
#' @export
#'
dim.RawH5ad <- function(x) {
  #  coords <- GetTissueCoordinates(object = x)
  #  return(c(
  #    max(coords[, 1]) - min(coords[, 1]),
  #    max(coords[, 2]) - min(coords[, 2])
  #  ))
  # return(dim(x = GetImage(object = x, mode = 'raw')))
  return(c(599, 600))
}

#' @method subset RawH5ad
#' @concept objects
#' @export
#'
subset.RawH5ad <- function(x, cells, ...) {
  x <- subset.VisiumV1(x = x, cells = cells, ...)
  return(x)
}

#' @method subset VisiumV1
#' @concept objects
#' @export
#'
subset.VisiumV1 <- function(x, cells, ...) {
  coordinates <- GetTissueCoordinates(object = x, scale = NULL, cols = NULL)
  cells <- cells[cells %in% rownames(x = coordinates)]
  coordinates <- coordinates[cells, ]
  slot(object = x, name = 'coordinates') <- coordinates
  return(x)
}

#' @rdname RenameCells
#' @concept objects
#' @method RenameCells VisiumV1
#' @export
#'
RenameCells.VisiumV1 <- function(object, new.names = NULL, ...) {
  if (is.null(x = new.names)) {
    return(object)
  } else if (length(x = new.names) != length(x = Cells(x = object))) {
    stop("Wrong number of cell/spot names", call. = FALSE)
  }
  names(x = new.names) <- Cells(x = object)
  coordinates <- GetTissueCoordinates(object = object, scale = NULL, cols = NULL)
  rownames(x = coordinates) <- new.names[rownames(x = coordinates)]
  slot(object = object, name = 'coordinates') <- coordinates
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Return a null image
#
# @param mode Image representation to return
# see \code{\link{GetImage}} for more details
#
#' @importFrom grid nullGrob
#' @importFrom grDevices as.raster
#
NullImage <- function(mode) {
  image <- switch(
    EXPR = mode,
    'grob' = nullGrob(),
    'raster' = as.raster(x = new(Class = 'matrix')),
    'plotly' = list('visible' = FALSE),
    'raw' = NULL,
    stop("Unknown image mode: ", mode, call. = FALSE)
  )
  return(image)
}

