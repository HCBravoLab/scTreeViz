#' ClusterHierarchy class to manage treeviz cluster data
setClass(
  "ClusterHierarchy",
  contains = c("DataFrame")
)


#' create a new ClusterHierarchy object. User can give either
#' col_regex or  columns option to filter the columns or specify 
#' the column order
#' @param hierarchy hierarchy as a dataFrame
#' @param col_regex Regular Expression for choosing columns
#' @param columns Vector containing list of columns to choose from with ordering
#' @return `ClusterHierarchy`` return an object of class ClusterHierarchy containing cluster information
#' @importFrom methods new
#' @importFrom S4Vectors DataFrame
#' @export
#' 
ClusterHierarchy <- function(hierarchy, col_regex=NULL, columns =NULL) {
  #will se later if this is needed
    # 
    # if (is.null(hierarchy)) {
    #   stop("No hirerarchy")
    #   
    # }
    # 
    # if (ncol(hierarchy) == 0) {
    #   return(
    #     new(
    #       "TreeCluster",
    #       DataFrame(hierarchy)
    #     )
    #   )
    # }
    # 

  #Filter key words
  if(!is.null(col_regex) && !is.null(columns)){
    message("Cannot use both")
  }
  if(!is.null(col_regex)){
      hierarchy <- hierarchy[, grep( col_regex, colnames(hierarchy))] 
  }
  
  # ordering
  if(!is.null(columns)){
    hierarchy <- hierarchy[columns]
  }
  for(cols in colnames(hierarchy)){
    hierarchy[[cols]]<- as.factor(hierarchy[[cols]])
  }
  
  hierarchy_dt <- as.data.table(hierarchy)
  hierarchy_dt$samples <- rownames(hierarchy_dt) <- rownames(hierarchy)
  cols <- colnames(hierarchy_dt)[1:length(colnames(hierarchy_dt))-1]
  order <- rep(1, length(hierarchy_dt)-1)
  hierarchy_dt <- setorderv(hierarchy_dt, cols = cols, order = order)
  hierarchy_df <- as.data.frame(hierarchy_dt)
  rownames(hierarchy_df) <- hierarchy_df$samples
  hierarchy <- hierarchy_df[,cols]
  
  # hierarchy<- as.data.frame(hierarchy)
  hierarchy <- rename_clusters(hierarchy)
  
  uniqeness <- check_unique_parent(hierarchy)
  
  #print(str(hierarchy))
  # Create clustree object
  hierarchy_graph <-
    clustree(
      hierarchy ,
      prefix = "cluster",
      prop_filter = 0,
      return = "graph"
    )
  
  if(uniqeness==FALSE){
    # prune the graph with only core edges (this makes it a ~tree)
    graph_df <- as_long_data_frame(hierarchy_graph)
    hierarchy <- prune_tree(graph_df, hierarchy)
    hierarchy_graph = hierarchy$Clustree_obj
    hierarchy <- DataFrame(hierarchy$Cluster_obj)
  }
  
  # collapses tree if the levels are the same at different resolutions
  collapsed_graph <- collapse_tree(hierarchy_graph)
  
  cluster_names <-
    unique(sapply(strsplit(collapsed_graph$node, "C"), '[', 1))
  
  #take only uncollapsed columns
  hierarchy <- hierarchy[, cluster_names]
  
  digits<- sub(pattern = "cluster", replacement = "", x=cluster_names)
  cluster_names <- paste0(digits, ".",cluster_names, sep="")
 
  #renaming nodes from numbers to cluster1C1 cluster1C2 so on..
  for (clusnames in names(hierarchy)) {
    hierarchy[[clusnames]] <-
      paste(clusnames, hierarchy[[clusnames]], sep = 'C')
  }
  
  names(hierarchy)<- cluster_names
 
  samples <- rownames(hierarchy)
  if(is.null(rownames(hierarchy))){
    samples<- 1:nrow(hierarchy)
  }
  
  hierarchy <- cbind(hierarchy, samples)
  hierarchy <- checkRoot(hierarchy)
  
  new(
    "ClusterHierarchy",
    hierarchy
  )
}