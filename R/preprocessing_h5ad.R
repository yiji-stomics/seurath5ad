#' @importFrom progressr progressor
#' @importFrom anndata read_h5ad
#' @importFrom Seurat CreateSeuratObject
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Load RawH5ad data
#'
#' @param h5ad_path h5ad path that contains counts matrix,
#' gene names, meta information, and spatial informations.
#' @param assay Name of assay to associate spatial data to
#'
#' @return A \code{\link{Seurat}} object
#'
#' @importFrom anndata read_h5ad
#'
#' @export
#' @concept preprocessing
#'
LoadRawH5ad <- function(h5ad_path, assay = "Spatial") {
  # check and find input files
  if (length(x = h5ad_path) > 1) {
    warning("'LoadRawH5ad' accepts only one 'data.dir'",
            immediate. = TRUE)
    h5ad_path <- h5ad_path[1]
  } else if (length(x = h5ad_path) == 0) {
    stop("No file matched the pattern '*.h5ad'", call. = FALSE)
  }
  
  # load h5ad and create seurat object
  h5ad_obj <- read_h5ad(h5ad_path)
  # extract matrix
  expression_matrix <- Matrix::t(h5ad_obj$X)
  # extract cells and genes information
  cell_metadata <- h5ad_obj$obs
  gene_metadata <- h5ad_obj$var
  # create seurat object
  seurat_object <- CreateSeuratObject(counts = expression_matrix, assay = assay, meta.data = cell_metadata)
  # extract spatial information
  spatial_data  <- h5ad_obj$obsm$spatial
  
  colnames(spatial_data) <- c("x", "y")
  rownames(spatial_data) <- rownames(seurat_object@meta.data)
  spatial_data <- as.data.frame(spatial_data)
  image <- new(Class = 'RawH5ad', assay = assay, coordinates = spatial_data)
  seurat_object[["Slice"]] <- image
  
  return(seurat_object)
}
