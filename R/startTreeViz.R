.register_all_treeviz_things <- function(app) {

  app$server$register_action("registerChartTypes", function(request_data) {
    app$chart_mgr$.register_available_chart_types(request_data$data)
  })

  app$server$register_action("getHierarchy", function(request_data) {
    datasource = request_data$datasourceGroup
    nodeId = request_data$nodeId
    obj <- app$data_mgr$.find_datasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getHierarchy(nodeId)
  })

  app$server$register_action("propagateHierarchyChanges", function(request_data) {

    datasource = request_data$datasourceGroup
    selection = request_data$selection
    order = request_data$order
    selectedLevels <- request_data$selectedLevels

    obj <- app$data_mgr$.find_datasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }

    obj$propagateHierarchyChanges(selection, order, selectedLevels)
  })

  app$server$register_action("getRows", function(request_data) {
    datasource <- request_data$datasource
    seqName <- request_data$seqName
    start <- request_data$start
    end <- request_data$end
    metadata <- request_data$metadata

    obj <- app$data_mgr$.find_datasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getRows(seqName, start, end, metadata)
  })

  app$server$register_action("getValues", function(request_data) {
    datasource <- request_data$datasource
    measurement <- request_data$measurement
    seqName <- request_data$seqName
    start <- request_data$start
    end <- request_data$end

    obj <- app$data_mgr$.find_datasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }
    obj$getValues(measurement, seqName, start, end)
  })

  app$server$register_action("getCombined", function(request_data) {

    measurementsList <- request_data$measurements
    result <- lapply(names(measurementsList), function(m) {
      seqName <- request_data$seqName
      start <- request_data$start
      end <- request_data$end
      order <- request_data$order
      nodeSelection <- request_data$selection
      selectedLevels <- request_data$selectedLevels
      measurements <- measurementsList[[m]]

      if(is.null(start)) {
        start <- 1
      }
      else {
        start <- start - 1
      }

      if(is.null(end)) {
        end <- 100000
      }
      else {
        end <- end
      }

      obj <- app$data_mgr$.find_datasource(m)
      if (is.null(obj)) {
        stop("cannot find datasource", m)
      }
      res <- obj$getCombined(measurements, seqName, start, end, order, nodeSelection, selectedLevels)
      if (class(obj) == "EpivizMetagenomicsDataTimeSeries"){
        res$rows$metadata$splines <- "true"
      }
      res
    })
    names(result) <- names(measurementsList)
    result
  })

  app$server$register_action("splinesSettings", function(request_data) {
    updateAlpha <- as.double(request_data$settings$alpha)
    names_list <- ls(app$data_mgr$.ms_list)
    class_list <- lapply(names_list, function(i) {class(app$data_mgr$.get_ms_object(i))})
    obj <- app$data_mgr$.get_ms_object(names_list[which(class_list == "EpivizMetagenomicsDataTimeSeries")][1])
    obj$updateSplineAlpha(updateAlpha)
    return(list())
  })

  app$server$register_action("getSeqInfos", function(request_data) {
    return(list(
      "treevizr" = c(1, 100000)
    ))
  })

  app$server$register_action("partitions", function(request_data) {
    return(list(
      "treevizr" = c(1, 100000)
    ))
  })

  app$server$register_action("getPCA", function(request_data) {

    measurementsList <- request_data$measurements
    result <- lapply(names(measurementsList), function(m) {
      seqName <- request_data$seqName
      measurements <- measurementsList[[m]]

      obj <- app$data_mgr$.find_datasource(m)
      if (is.null(obj)) {
        stop("cannot find datasource", m)
      }
      obj$getPCA(measurements)
    })
    names(result) <- names(measurementsList)
    result

  })

  app$server$register_action("getDiversity", function(request_data) {

    measurementsList <- request_data$measurements
    result <- lapply(names(measurementsList), function(m) {
      seqName <- request_data$seqName
      measurements <- measurementsList[[m]]

      obj <- app$data_mgr$.find_datasource(m)
      if (is.null(obj)) {
        stop("cannot find datasource", m)
      }
      obj$getAlphaDiversity(measurements)
    })
    names(result) <- names(measurementsList)
    result

  })

  app$server$register_action("search", function(request_data) {
    query <- request_data$q
    max_results <- request_data$maxResults
    datasource = request_data$datasourceGroup

    obj <- app$data_mgr$.find_datasource(datasource)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource)
    }

    list(nodes = obj$searchTaxonomy(query, max_results))
  })
}

#' Start treeviz app and create \code{\link[treevizr]{TreevizApp}} object to manage connection.
#'
#' @param data TreeViz object to explore
#' @param host (character) host address to launch.
#' @param register_function (function) function used to register actions and charts on the treeviz app.
#' @param ... additional parameters passed to \code{\link[epivizr]{startEpiviz}}.
#'
#' @return An object of class \code{\link[treevizr]{TreevizApp}}
#'
#' @import epivizr
#' @import sys
#' 
#' @seealso \code{\link[treevizr]{TreevizApp}}
#' @examples
#' # see package vignette for example usage
#' app <- startTreeviz(non_interactive=TRUE, open_browser=FALSE)
#' app$stop_app()
#'
#' @export
startTreeviz <- function(data = NULL, host="http://metaviz.cbcb.umd.edu",
                         register_function = .register_all_treeviz_things,
                         ...) {
  chr="treevizr"
  start <- 1
  end <- 100000
  app <- startEpiviz(host = host, register_function = register_function,
                     chr=chr, start=start, end=end, ...)
  mApp <- TreeVizApp$new(.url_parms=app$.url_parms, .browser_fun=app$.browser_fun,
                         server=app$server, data_mgr=app$data_mgr, chart_mgr=app$chart_mgr)
  
  
  # if (is(data, "Seurat")) {
  #   treeViz <- createFromSeurat(data)
  # }
  # else if (is(data, "SingleCellExperiment")) {
  #   treeViz <- createFromSCE(data) 
  # }
  
  if (!is.null(data)) {
    # treeViz <- data 
    
    # find variable genes
    data <- find_top_variable_genes(data, 100)
    
    now <- Sys.time()
    while ((Sys.time() - now) < 10) {
      # just iterate
      # mApp$server$wait_to_clear_requests()
    }
    
    # add facetZoom
    facetZoom <- mApp$plot(data, datasource_name = "SCRNA", tree = "col")
    
    mApp$server$wait_to_clear_requests()
    
    mes <- mApp$get_ms_object(chart_id_or_object = facetZoom)
    # Get Measurements from the plot
    ms_list <- facetZoom$get_measurements()
    subset_ms_list <- Filter(function(ms) ms@id %in% metadata(data)$top_variable, ms_list)
    
    # add Heatmap
    mApp$chart_mgr$visualize(chart_type = "HeatmapPlot",  measurements = subset_ms_list)
    
    mApp$server$wait_to_clear_requests()
    
    # TSNE
    mApp$chart_mgr$revisualize(chart_type = "PCAScatterPlot", chart= facetZoom)
    
    mApp$server$wait_to_clear_requests() 
  }
  
  mApp
}


.viewer_option_browse_fun <- function(url) {
  viewer <- getOption("viewer")
  if (is.null(viewer)) {
    utils::browseURL(url)
  } else {
    viewer(url)
  }
}

#' Start treeviz app in standalone (locally) and create \code{\link[treeviz]{TreevizApp}} object to manage connection.
#'
#' @param register_function (function) function used to register actions and charts on the treeviz app.
#' @param use_viewer_option (function) run application in viewer defined by \code{getOption("viewer")}.
#'  This allows standalone app to run in Rstudio's viewer (FALSE by default)
#' @param ... additional parameters passed to \code{\link[epivizrStandalone]{startStandalone}}.
#'
#' @return An object of class \code{\link{TreevizApp}}
#'
#' @import epivizrStandalone
#' @examples
#'
#' #' # see package vignette for example usage
#' app <- startTreevizStandalone(non_interactive=TRUE)
#' app$stop_app()
#'
#'
#' @export
startTreevizStandalone <- function(register_function = .register_all_treeviz_things,
                                   use_viewer_option=FALSE, ...) {
  chr="treevizr"
  start=1
  end=100
  seq <- Seqinfo(seqnames=chr,
                 seqlengths=100000,
                 isCircular=FALSE,
                 genome="treevizr")

  path <- system.file("www", package="epivizrStandalone")



  app <- startStandalone(seqinfo=seq,
                         register_function=register_function,
                         use_viewer_option = use_viewer_option,
                         chr=chr, start=start, end=end, ...)

  mApp <- TreeVizApp$new(.url_parms=app$.url_parms, .browser_fun=app$.browser_fun,
                         server=app$server, data_mgr=app$data_mgr, chart_mgr=app$chart_mgr)
  mApp
}

#' set treeviz app standalone settings
#'
#' @param url (character) github url to use. defaults to (\url{"https://github.com/epiviz/epiviz.git"}).
#' @param branch (character) branch on the github repository. defaults to (master).
#' @param local_path (character) if you already have a local instance of treeviz and would like to run standalone use this.
#' @param non_interactive (logical) don't download repo, used for testing purposes.
#' @return path to the treeviz app git repository
#'
#' @import epivizrStandalone
#' @examples
#'
#' \dontrun{
#' #' # see package vignette for example usage
#' setTreevizStandalone()
#' }
#'
#' @export
setTreevizStandalone <- function(url="https://github.com/epiviz/epiviz.git", branch="metaviz-4.1", local_path=NULL, non_interactive=FALSE) {
  setStandalone(url = url, branch = branch, local_path = local_path, non_interactive = non_interactive)
}
