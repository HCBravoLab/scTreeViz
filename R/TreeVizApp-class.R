#' Class managing connection to metaviz application.
#' 
#' @importClassesFrom epivizr EpivizApp
#' @importClassesFrom epivizrServer EpivizServer
#' @importClassesFrom epivizrData EpivizDataMgr EpivizMeasurement EpivizData
#' @exportClass TreeVizApp
TreeVizApp <- setRefClass("TreeVizApp",
                        contains = "EpivizApp",
                        fields=list(
                          .url_parms="list",
                          .browser_fun="function",
                          server="EpivizServer",
                          data_mgr="EpivizDataMgr",
                          chart_mgr="EpivizChartMgr"
                        ),
                        methods=list(
                          initialize=function(.url_parms=.url_parms, .browser_fun=.browser_fun,
                                              server=server, data_mgr=data_mgr, chart_mgr=chart_mgr) {
                            callSuper(.url_parms=.url_parms, .browser_fun=.browser_fun,
                                      server=server, data_mgr=data_mgr, chart_mgr=chart_mgr)
                          }
                        )
)

# TreeVizApp$methods(
#   explore=function(data) {
#     'explore adds facetZoon, heatmap and TSNE to the app.'
#     # treeViz <- data 
#     
#     # find variable genes
#     data <- find_top_variable_genes(data, 100)
#     
#     now <- Sys.time()
#     while ((Sys.time() - now) < 10) {
#       # just iterate
#       # .self$server$wait_to_clear_requests()
#     }
#     
#     .self$server$wait_to_clear_requests()
#     
#     # add facetZoom
#     facetZoom <- .self$plot(data, datasource_name = "SCRNA", tree = "col")
#     
#     .self$server$wait_to_clear_requests()
#     
#     now <- Sys.time()
#     while ((Sys.time() - now) < 10) {
#       # just iterate
#       # .self$server$wait_to_clear_requests()
#     }
#     
#     mes <- .self$get_ms_object(chart_id_or_object = facetZoom)
#     # Get Measurements from the plot
#     ms_list <- facetZoom$get_measurements()
#     subset_ms_list <- Filter(function(ms) ms@id %in% metadata(data)$top_variable, ms_list)
#     
#     print("variable genes length")
#     print(length(subset_ms_list))
#     
#     # add Heatmap
#     .self$chart_mgr$visualize(chart_type = "HeatmapPlot",  measurements = subset_ms_list)
#     
#     .self$server$wait_to_clear_requests()
#     
#     now <- Sys.time()
#     while ((Sys.time() - now) < 10) {
#       # just iterate
#       # .self$server$wait_to_clear_requests()
#     }
#     
#     # TSNE
#     .self$chart_mgr$revisualize(chart_type = "PCAScatterPlot", chart= facetZoom)
#     
#     .self$server$wait_to_clear_requests() 
#   }
# )

# TreeVizApp$methods(
#   explore=function(data) {
#     'explore adds facetZoon, heatmap and TSNE to the app.'
#     # treeViz <- data 
#     
#     # find variable genes
#     data <- find_top_variable_genes(data, 100)
#     
#     now <- Sys.time()
#     while ((Sys.time() - now) < 10) {
#       # just iterate
#       # .self$server$wait_to_clear_requests()
#     }
#     
#     .self$server$wait_to_clear_requests()
#     
#     # add facetZoom
#     facetZoom <- .self$plot(data, datasource_name = "SCRNA", tree = "col")
#     
#     .self$server$wait_to_clear_requests()
#     
#     now <- Sys.time()
#     while ((Sys.time() - now) < 10) {
#       # just iterate
#       # .self$server$wait_to_clear_requests()
#     }
#     
#     mes <- .self$get_ms_object(chart_id_or_object = facetZoom)
#     # Get Measurements from the plot
#     ms_list <- facetZoom$get_measurements()
#     subset_ms_list <- Filter(function(ms) ms@id %in% metadata(data)$top_variable, ms_list)
#     
#     print("variable genes length")
#     print(length(subset_ms_list))
#     
#     # add Heatmap
#     .self$chart_mgr$visualize(chart_type = "HeatmapPlot",  measurements = subset_ms_list)
#     
#     .self$server$wait_to_clear_requests()
#     
#     now <- Sys.time()
#     while ((Sys.time() - now) < 10) {
#       # just iterate
#       # .self$server$wait_to_clear_requests()
#     }
#     
#     # TSNE
#     .self$chart_mgr$revisualize(chart_type = "PCAScatterPlot", chart= facetZoom)
#     
#     .self$server$wait_to_clear_requests() 
#   }
# )

# # navigation methods
# TreeVizApp$methods(
#   navigate=function(node) {
#     'Navigate to a given node on the metaviz app.'
#     if (.self$is_server_closed()) {
#       return(invisible())
#     }
#     
#     callback <- function(response_data) {
#       invisible()
#     }
#     request_data=list(action="navigate",
#                       range=epivizrServer::json_writer(
#                         list(seqName="metavizr",start=node$leafIndex(),end=(node$leafIndex() + node$nleaves()))))
#     .self$server$send_request(request_data, callback)
#   },
#   slideshow=function(nodes, n=length(nodes), .callback=NULL) {
#     'Navigate on metaviz app successively to given nodes.
#     
#     \\describe{
#     \\item{nodes}{An object of class indicating
#     set of nodes to navigate in epiviz app.}
#     \\item{n}{(integer) The number of regions in \\code{nodes} to navigate to.}
#     \\item{.callback}{(function) function to call after navigating to each region. Used for testing purposes.}
#     }'
# 
#     n <- min(n, length(nodes))
#     ind <- seq(len=n)
#     chr <- 'metavizr'
#     start <- nodes[ind]$leafIndex()
#     end <- (nodes[ind]$leafIndex() + nodes[ind]$nleaves())
#     for (i in ind) {
#       .self$navigate(chr=chr[i], start=start[i], end=end[i])
#       if (!is.null(.callback) && is.function(.callback)) {
#         .callback(chr[i], start[i], end[i])
#       } else {
#         if (.self$server$is_interactive()) {
#           cat("Node", i, "of", n, ". Press key to continue (ESC to stop)...\n")
#           readLines(n=1)
#         }
#       }
#       tryCatch(.self$server$service(), interrupt=function(int) invisible())
#     }
#     invisible()
#   }
# )