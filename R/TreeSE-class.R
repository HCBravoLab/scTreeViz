
#' TreeSummarizedExperiment class wrapper for SummarizedExperiment objects
#' @import SummarizedExperiment
setClass("TreeSummarizedExperiment",
         contains = "SummarizedExperiment")

# validity
setValidity("TreeSummarizedExperiment", function(object) {
  msg <- NULL
  if(!(is(rowData(object), "TreeIndex") || is(colData(object), "TreeIndex"))) {
     msg <- "neither rowData nor colData are TreeIndex objects, use SummarizedExperiment instead."
  }

  if(is.null(msg)) TRUE else msg
})

#' The TreeSummarizedExperiment class.
#'
#' SummarizedExperiment-like class for datasets that have hierarchies on either rowData or colData.
#' For microbiome data, rowData is a tree hierarchy
#' For single cell data, colData is a tree hierarchy
#' @param assays simple list of counts
#' @param rowData rowData
#' @param colData colData
#' @param ... other parameters for SummarizedExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new
#' @export
TreeSummarizedExperiment <- function(assays = SimpleList(),
                   rowData = NULL,
                   colData = NULL,
                   ...) {

  if (missing(colData) || is.null(colData)) {
    assay <- assays[[1]]
    colData <- DataFrame(x=seq_len(ncol(assay)), row.names=colnames(assay))[, FALSE]
  }

  sumExp <- SummarizedExperiment(assays = assays, rowData = rowData, colData = colData, ...)

  new("TreeSummarizedExperiment", sumExp)
}

#' Import `Seurat` and `Clustree` into `TreeSummarizedExperiment`
#'
#' @param seurat seurat object that contains rowData
#' @param cluster_names sluter names from seurat if different from standard
#' @param clustree clustree object in graph format
#' @import tidygraph
#' @export
ImportFromSeurat <- function(seurat, clustree, cluster_names = NULL) {
  # get clusters from seurat
  clusters <- seurat@meta.data
  if(is.null(cluster_names)) {
    cluster_names <- colnames(clusters)
    cluster_names <- cluster_names[grepl("res", cluster_names)]
  }

  clusters <- clusters[, cluster_names]

  # get cluster nodes from clustree
  clusnodes <- as.data.frame(with_graph(clustree, .N()))
  clusnames <- clusnodes$cluster
  names(clusnames) <- clusnodes$node

  # parse and create a dataframish structure
  clusters <- lapply(colnames(clusters), function(cn) {
    col <- clusters[, cn]
    col <- paste0(cn, "C", col)
    for (i in names(clusnames)) {
      col <- replace(col, col==i, paste0(i, "-cluster", clusnames[[i]]))
    }
    col
  })

  names(clusters) <- cluster_names
  clusters <- as.data.frame(clusters)
  clusters$samples <- colnames(seurat@data)
  cluster_names <- c(cluster_names, "samples")

  tree <- TreeIndex(clusters, cluster_names)
  rownames(tree) <- colnames(seurat@data)

  TreeSummarizedExperiment(SimpleList(counts = seurat@data), colData = tree)
}


