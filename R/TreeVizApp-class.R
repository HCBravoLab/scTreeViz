#' Class managing connection to metaviz application.
#'
#' @importClassesFrom epivizr EpivizApp
#' @importClassesFrom epivizrServer EpivizServer
#' @importClassesFrom epivizrData EpivizDataMgr EpivizMeasurement EpivizData
#' @exportClass TreeVizApp
TreeVizApp <- setRefClass(
  "TreeVizApp",
  contains = "EpivizApp",
  fields = list(
    .url_parms = "list",
    .browser_fun = "function",
    server = "EpivizServer",
    data_mgr = "EpivizDataMgr",
    chart_mgr = "EpivizChartMgr"
  ),
  methods = list(
    initialize = function(.url_parms = .url_parms,
                          .browser_fun = .browser_fun,
                          server = server,
                          data_mgr = data_mgr,
                          chart_mgr = chart_mgr) {
      callSuper(
        .url_parms = .url_parms,
        .browser_fun = .browser_fun,
        server = server,
        data_mgr = data_mgr,
        chart_mgr = chart_mgr
      )
    }
  )
)

TreeVizApp$methods(
  plotGene = function(gene = NULL,
                      datasource_name = "SCRNA_1") {
    if (is.null(gene)) {
      stop("gene symbol must be provided")
    }
    
    obj <- .self$data_mgr$.find_datasource(datasource_name)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource_name)
    }
    
    ms_list <- obj$get_measurements()
    subset_ms_list <- Filter(function(ms)
      ms@id == gene, ms_list)
    
    if (length(subset_ms_list) == 0) {
      stop("cannot find the gene in the dataset!")
    }
    
    .self$chart_mgr$visualize(chart_type = "DiversityScatterPlot",  measurements = subset_ms_list)
    .self$server$wait_to_clear_requests()
    .delay_requests(.self$server)
    
  },
  
  
  #' @import umap
  #' @import Rtsne
  extract_SCE = function(datasource_name = "SCRNA_1") {
    obj <- .self$data_mgr$.find_datasource(datasource_name)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource_name)
    }
    level <- obj$.levelSelected+1
    node <- obj$.nodeSelections
    sce<-obj$extract_SCE_epiviz(level, node)
    
    # aggtreeviz <-
    #   aggregateTree(obj$.object, level, node, by = "col", format = "TreeViz")
    # 
    
    # sce <-
    #   SingleCellExperiment(
    #     assays = list(counts = assays(aggtreeviz)$counts),
    #     colData = colData(aggtreeviz)@hierarchy_tree
    #   )
    # 
    # reduced_dims=list()
    # for( name in names(metadata(obj$.object)$reduced_dim)){
    #   if(name=="pca"){
    #     pca_data <- prcomp(t(assays(aggtreeviz)$counts), rank=2)
    #     reduced_dims[[name]]=pca_data$x
    #   }
    #   if(name=="tsne"){
    #     tsne_data <- Rtsne(t(assays(aggtreeviz)$counts),  perplexity = 1)
    #     reduced_dims[[name]]=tsne_data$Y
    #   }
    #   
    #   if(name=="umap"){
    #     custom.settings = umap.defaults
    #     custom.settings$n_neighbors = 2
    #     
    #     umap_data <- umap(t(assays(aggtreeviz)$counts),custom.settings)
    #     reduced_dims[[name]]=umap_data$layout
    #   }
    #   
    # }
    # sce <- SingleCellExperiment(assays = list(counts = assays(aggtreeviz)$counts
    # ),
    # colData = colData(aggtreeviz)@hierarchy_tree,
    # reducedDims =reduced_dims
    # )
    # 
  }
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