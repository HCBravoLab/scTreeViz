setClassUnion("CharacterOrNull", c("character", "NULL"))

#' TreeCluster class to manage treeviz cluster data
setClass(
  "TreeCluster",
  contains = c("DataFrame")
)


#' create a new Treecluster object
#' @param hierarchy hierarchy as a dataFrame
#' @importFrom methods new
#' @importFrom S4Vectors DataFrame
#' @export
TreeCluster <- function(hierarchy = NULL, filter=NULL, keyword="cluster", ordering =NULL) {
  #will se later if this is needed
    # 
    if (is.null(hierarchy)) {
      stop("No hirerarchy")
      
    }
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
  if(!is.null(filter)){
    if (is.null(keyword)) {
      stop("No Keyword giver for Filter")
      
    }
    if(filter=="prefix"){
      hierarchy<- hierarchy[, startsWith(colnames(hierarchy), keyword)] 
    }
    if(filter=="suffix"){
      hierarchy<- hierarchy[, endssWith(colnames(hierarchy), keyword)] 
    }
    if(filter=="match"){
      hierarchy<- hierarchy[, grep( keyword, colnames(hierarchy))] 
    }
    
  }
  
  #ordering
  #print(colnames(hierarchy))
  #column name, not indices
  if(!is.null(ordering)){
    hierarchy<-hierarchy[ordering]
    #print(hierarchy)
    #reorder(hierarchy,ordering)
  }
  for(cols in colnames(hierarchy)){
    hierarchy[[cols]]<- as.factor(hierarchy[[cols]])
  }
  hierarchy <- rename_clusters(hierarchy)
  
  uniqeness <-
    check_unique_parent(hierarchy)
  
  # Create clustree object
  hierarchy_graph <-
    clustree(
      hierarchy ,
      prefix = "cluster",
      prop_filter = 0,
      return = "graph"
    )
  
  
  
  if(uniqeness==FALSE){
    #print(hierarchy_graph)
    # prune the graph with only core edges (this makes it a ~tree)
    print("non-unique")
    graph_df <- as_long_data_frame(hierarchy_graph)
    hierarchy <- prune_tree(graph_df, hierarchy)
    #print(hierarchy)
    hierarchy_graph = hierarchy$Clustree_obj
    hierarchy <- DataFrame(hierarchy$Cluster_obj)
    
  }
  
  
  # collapses tree if the levels are the same at different resolutions
  #print(hierarchy_graph)
  collapsed_graph <- collapse_tree(hierarchy_graph)
  
  
  cluster_names <-
    unique(sapply(strsplit(collapsed_graph$node, "C"), '[', 1))
  
  
  hierarchy <- hierarchy[, cluster_names]
  
  for (clusnames in names(hierarchy)) {
    hierarchy[[clusnames]] <-
      paste(clusnames, hierarchy[[clusnames]], sep = 'C')
  }
  
  samples <- rownames(hierarchy)
  if(is.null(rownames(hierarchy))){
    samples<- 1:nrow(hierarchy)
  }
  print(samples)
  hierarchy <- cbind(hierarchy, samples)
  hierarchy <- checkRoot(hierarchy)
  
  
  new(
    "TreeCluster",
    hierarchy
  )
  
  # tree <- TreeIndex(hierarchy)
  # rownames(tree) <- rownames(hierarchy)
  # 
  # treeviz <- TreeViz(SimpleList(counts = counts), colData = tree)
  # # plot_tree(TreeSE_obj@colData@hierarchy_tree)
  # treeviz
  # 
}