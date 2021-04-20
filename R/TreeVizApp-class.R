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

TreeVizApp$methods(
  plotGene = function(gene = NULL, datasource_name = "SCRNA_1") {
    " Plot a bar plot for a gene across cell types
    \\describe{
    \\item{gene}{gene to extract expression values}
    \\item{datasource_name}{object to extract from (automatically selected)}
    }
    "
    
    if (is.null(gene)) {
      stop("gene symbol must be provided")
    }
    
    obj <- .self$data_mgr$.find_datasource(datasource_name)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource_name)
    }
    
    ms_list <-obj$get_measurements()
    subset_ms_list <- Filter(function(ms) ms@id == gene, ms_list)
    
    if (length(subset_ms_list) == 0) {
      stop("cannot find the gene in the dataset!")
    }
    
    .self$chart_mgr$visualize(chart_type = "DiversityScatterPlot",  measurements = subset_ms_list)
    .self$server$wait_to_clear_requests()
    .delay_requests(.self$server)
    
  },
  
  extract_SCE = function(cluster_name="treeviz_clusters", datasource_name = "SCRNA_1") {
    obj <- .self$data_mgr$.find_datasource(datasource_name)
    if (is.null(obj)) {
      stop("cannot find datasource", datasource_name)
    }

    sce<-obj$extract_SCE_epiviz(cluster_name)
  }
)