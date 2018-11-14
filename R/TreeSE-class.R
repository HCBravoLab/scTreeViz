setClass("TreeSE",
         contains = "SummarizedExperiment")

#' The TreeSE class.
#'
#' SummarizedExperiment-like class for datasets that have hierarchies on either rowData or colData.
#' For microbiome data, rowData is a tree hierarchy
#' For single cell data, colData is a tree hierarchy
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData rowData<- colData<-
#' @import S4Vectors
#' @export
TreeSE <- function(assays = SimpleList(),
                   rowData = TreeIndex(),
                   colData = TreeIndex(),
                   ...) {

  sumExp <- SummarizedExperiment(assays = assays, ...)

  if(!missing(rowData)) {
    rowData(sumExp) <- rowData
  }

  if(!missing(colData)) {
    colData(sumExp) <- colData
  }

  new("TreeSE", sumExp)
}
