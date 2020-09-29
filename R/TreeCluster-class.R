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
  #print(colnames(hierarchy))
  #column name, not indices
  if(!is.null(ordering)){
    hierarchy<-hierarchy[ordering]
    #print(hierarchy)
    #reorder(hierarchy,ordering)
  }
  uniqeness <-
    check_unique_parent(hierarchy)
  if(uniqeness==FALSE){
    clusters <- rename_clusters(hierarchy)
    
    # Create clustree object
    clustree_graph <-
      clustree(
        clusters ,
        prefix = "cluster",
        prop_filter = 0,
        return = "graph"
      )
    
    graph_df <- as_long_data_frame(clustree_graph)
    
    # prune the graph with only core edges (this makes it a ~tree)
    hierarchy <- prune_tree(graph_df, clusters)
  }
  
  hierarchy_df <- DataFrame(hierarchy$Cluster_obj)
  #hierarchy_df <- hierarchy$Cluster_obj
  
  
  new(
    "TreeCluster",
    hierarchy_df
  )
}