
#' TreeSE class wrapper for SummarizedExperiment objects
#' @import SummarizedExperiment
setClass("TreeSE",
         contains = "SummarizedExperiment")

#' The TreeSE class.
#'
#' SummarizedExperiment-like class for datasets that have hierarchies on either rowData or colData.
#' For microbiome data, rowData is a tree hierarchy
#' For single cell data, colData is a tree hierarchy
#' @param assays simple list of counts
#' @param rowData rowData
#' @param colData colData
#' @param ... other parameters for SummarizedExperiment
#' @export
TreeSE <- function(assays = SimpleList(),
                   rowData = NULL,
                   colData = NULL,
                   ...) {
  sumExp <- SummarizedExperiment(assays = assays, ...)

  if (!missing(rowData) && !is.null(rowData)) {
    rowData(sumExp) <- rowData
  }

  if (!missing(colData) && !is.null(colData)) {
    colData(sumExp) <- colData

  }

  new("TreeSE", sumExp)
}
