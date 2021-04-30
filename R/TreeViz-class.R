
#' TreeViz class wrapper for SummarizedExperiment objects
#' @import SummarizedExperiment
#' @exportClass TreeViz
setClass("TreeViz",
         contains = "SummarizedExperiment")

# validity
setValidity("TreeViz", function(object) {
  msg <- NULL
  if(!(is(rowData(object), "TreeIndex") || is(colData(object), "TreeIndex"))) {
     msg <- "neither rowData nor colData are TreeIndex objects, use SummarizedExperiment instead."
  }

  if(is.null(msg)) TRUE else msg
})

#' The TreeViz class.
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
#' @examples
#' library(metagenomeSeq)
#' data(mouseData)
#' counts <- MRcounts(mouseData)
#' hierarchy <- fData(mouseData)
#' tree <- TreeIndex(hierarchy)
#' mbiome <- TreeViz(SimpleList(counts=counts), rowData=tree)
#' @export
TreeViz <- function(assays = SimpleList(),
                   rowData = NULL,
                   colData = NULL,
                   ...) {

  if (missing(colData) || is.null(colData)) {
    assay <- assays[[1]]
    colData <- DataFrame(x=seq_len(ncol(assay)), row.names=colnames(assay))[, FALSE]
  }

  sumExp <- SummarizedExperiment(assays = assays, rowData = rowData, colData = colData, ...)

  new("TreeViz", sumExp)
}
