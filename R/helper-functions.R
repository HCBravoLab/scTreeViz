#' check if alternative edge exists
#' @param sub_df dataframe
#' @param all_df dataframe
#' @return Dataframe with assigned edges
#'
check_alternate <- function(sub_df, all_df)
{
  #check where clusters from not-core edges should actually belong
  assign_df <- data.frame(
    "Sink_node_res" = character(),
    "Sink_node_clust" = character(),
    "Assign_from_res" = character(),
    "Assign_from_clust" = character(),
    "Assign_to_res" = character(),
    "Assign_to_clust" = character(),
    stringsAsFactors = FALSE
  )
  for (i in seq(1, nrow(sub_df))) {
    for (j in seq(1, nrow(all_df))) {
      if (sub_df$to_node[i] == all_df$to_node[j]) {
        assign_df[nrow(assign_df) + 1, ] <-
          c(
            sub_df$to_cluster[i],
            sub_df$to_clust[i],
            sub_df$from_cluster[i],
            sub_df$from_clust[i],
            all_df$from_cluster[j],
            all_df$from_clust[j]
          )
      }
    }
  }
  
  assign_df
}


#' propagate cluster assignment changes from graph object to hierarchy
#'
#' @param graph graph from clustree
#' @param cluster_obj cluster hierarchies at different resolutions
#' @return hierarchy of cluster resolutions
#'
change_assignment <- function(graph, cluster_obj) {
  #for each row in the core edge false dataset, find out corresponding
  #source and sink in Seurat, and change accordingly
  for (i in seq(1:nrow(graph))) {
    sink_res = match(paste0('cluster', graph$`Sink_node_res`[i]),
                     colnames(cluster_obj))
    source_res = match(paste0('cluster', graph$`Assign_from_res`[i]),
                       colnames(cluster_obj))
    
    
    for (j in 1:nrow(cluster_obj)) {
      if (as.numeric(as.character(cluster_obj[[j, sink_res]])) == as.numeric(graph$Sink_node_clust[i])  &&
          as.numeric(as.character(cluster_obj[[j, source_res]])) == as.numeric(graph$Assign_from_clust[i])) {
        cluster_obj[[j, source_res]] <-
          as.factor(graph$`Assign_to_clust`[i])
      }
      
    }
  }
  return (cluster_obj)
}

#' check if graph is cyclic
#'
#' @param pruned_graph graph from clustree
#' @return graph object after removing cyclic edges
#'
check_cycle <- function(pruned_graph) {
  # Check if a diamond or circle still exists in tree
  delete_set_vertices <- vector('numeric')
  ver_list <- V(pruned_graph)
  
  for (nodes in ver_list) {
    adj_edge <- incident(pruned_graph, nodes, mode = "in")
    
    if (length(adj_edge) > 1) {
      # Not a Tree yet
      adj_ver <-
        adjacent_vertices(pruned_graph, nodes, mode = "in")
      adj_ver <- as_ids(adj_ver[[1]])
      str(adj_ver)
      remove_node <- sample(adj_ver, 1)
      
      remove_edge <-
        incident(pruned_graph, remove_node, mode = "all")
      
      pruned_graph <- delete_edges(pruned_graph, remove_edge)
      delete_set_vertices <-
        c(delete_set_vertices, remove_node)
      
    }
  }
  pruned_graph <- delete.vertices(pruned_graph, delete_set_vertices)
}

#' prune tree to remove non-core edges from clustree hierarchy
#'
#' @param graph graph from clustree
#' @param cluster_df Original cluster hierarchy as dataframe
#' @return Dataframe containing cluster information at different resolutions
#'
prune_tree <- function(graph, cluster_df) {
  # Prune the tree so only core edges remain
  repeat {
    # drop duplicated columns (Throws error otherwise)
    
    graph <- graph[!duplicated(names(graph))]
    
    # See number of core edges at each ites
    # print(nrow(graph_Df[graph_Df$is_core == FALSE,]))
    # No False edges acyclic tree
    if (nrow(graph[graph$is_core == FALSE,]) == 0)
      break
    
    # Apply function
    assign_df <-
      check_alternate(graph[graph$is_core == FALSE,], graph[graph$is_core == TRUE,])
    
    cluster_df <- change_assignment(assign_df, cluster_df)
    ngraph <-
      clustree(
        cluster_df ,
        prefix = "cluster",
        prop_filter = 0,
        return = "graph"
      )
    
    graph <- as_long_data_frame(ngraph)
  }
  
  ngraph <- check_cycle(ngraph)
  
  result <-
    list("Cluster_obj" = cluster_df, "Clustree_obj" = ngraph)
  result
}

#' collapse tree if resolutions does not change cluster assignments
#'
#' @param original_graph Dataframe containing cluster information at different resolutions
#' @return Dataframe containing cluster information at different resolutions
#'
collapse_tree <- function(original_graph) {
  node_dists <- distances(original_graph, to = 1)
  layered_dist <- unique(distances(original_graph, to = 1))
  
  delete_set_vertices <- vector('numeric')
  
  # Vertex and edge lists of the graph from which we will construct the collapsed graph
  ver_list <-
    igraph::as_data_frame(original_graph, what = "vertices")
  
  i <- length(layered_dist)
  while (i >= 2) {
    prev_layer <- which(node_dists == layered_dist[i - 1])
    current_layer <- which(node_dists == layered_dist[i])
    
    if (length(prev_layer) == length(current_layer)) {
      while (length(prev_layer) == length(current_layer)) {
        delete_set_vertices <- c(delete_set_vertices, prev_layer)
        i <- i - 1
        
        prev_layer <-
          which(node_dists == layered_dist[i - 1])
        
      }
    }
    i <- i - 1
  }
  
  ver_list <- ver_list[-delete_set_vertices,]
  ver_list
}

#' check if hierarchy has more than 1 root, if do create a root
#'
#' @param cluster_df Dataframe containing cluster information at different resolutions
#' @return Dataframe containing cluster information at different resolutions
#'
checkRoot <- function(cluster_df) {
  # Handle Forests
  cols <- colnames(cluster_df)
  if (length(unique(cluster_df[[1]])) > 1) {
    cluster_df$root <- "ClusterAllClusters"
    cols <- c("root", cols)
  }
  
  cluster_df[cols]
}

#' check if every node has a single parent (tree vs graph)
#'
#' @param clusterdata Dataframe containing cluster information at different resolutions
#' @return Dataframe containing cluster information at different resolutions
#'
check_unique_parent <- function(clusterdata) {
  #check if user provided data has unique parents at each level
  for (i in seq(2, ncol(clusterdata))) {
    childs <- unique(clusterdata[[i]])
    for (values in childs) {
      subsetted_list <- clusterdata[clusterdata[[colnames(clusterdata)[[i]]]] == values,]
      
      
      parent <- length(unique(subsetted_list[[i - 1]]))
      
      if (parent > 1) {
        stop("Not a tree, some nodes with multiple parents in level", i)
      }
      #cat(nrow(subsetted_list)," ",parent, "\n")
    }
  }
}


#' Renames cluster names for consistency in visualization
#'
#' @param clusdata Dataframe containing cluster information at different resolutions
#' @return Dataframe containing cluster information at different resolutions
#'
rename_clusters <- function(clusdata) {
  clusnames <- seq(length(colnames(clusdata)))
  
  clusnames <- paste0("cluster", clusnames)
  colnames(clusdata) <- clusnames
  clusdata
}


#' Creates the `TreeViz` object from cluster hierarchy and count matrix.
#' Also prunes and collapses the hierarchy to remove redundant levels.
#'
#' @param clusters Dataframe containing cluster information at different resolutions
#' @param count matrix Dense or sparse matrix containing the count matrix
#' @return TreeSummarizedExperiment Object that can be visualized with metavizr
#' @import clustree
#' @importFrom igraph as_long_data_frame
#'
preprocessAndCreateTreeViz <- function(clusters, counts) {
  clusters <- rename_clusters(clusters)
  
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
  modified_obj <- prune_tree(graph_df, clusters)
  
  # modified graph and seurat object
  modified_graph = modified_obj$Clustree_obj
  clusters_new <-  modified_obj$Cluster_obj
  
  # collapses tree if the levels are the same at different resolutions
  collapsed_graph <- collapse_tree(modified_graph)
  cluster_names <-
    unique(sapply(strsplit(collapsed_graph$node, "C"), '[', 1))
  clusters_new <- clusters_new[, cluster_names]
  
  for (clusnames in names(clusters_new)) {
    clusters_new[[clusnames]] <-
      paste(clusnames, clusters_new[[clusnames]], sep = 'C')
  }
  
  samples <- rownames(clusters_new)
  clusters_new <- cbind(clusters_new, samples)
  clusters_new <- checkRoot(clusters_new)
  
  tree <- TreeIndex(clusters_new)
  rownames(tree) <- rownames(clusters_new)
  
  treeviz <- TreeViz(SimpleList(counts = counts), colData = tree)
  # plot_tree(TreeSE_obj@colData@hierarchy_tree)
  treeviz
}


#' Creates a `TreeViz` object from `Seurat`
#'
#' @param object `Seurat` class containing cluster information at different resolutions
#' @return `TreeViz` Object
#' @export
#'
createFromSeurat <- function(object) {
  clusterdata <- object@meta.data
  clusterdata <- clusterdata[, grep("*snn*", colnames(clusterdata))]
  
  treeviz <-
    preprocessAndCreateTreeViz(clusterdata, GetAssayData(object))
  
  if ("tsne" %in% Reductions(object)) {
    reducdim <- Reductions(object, slot = "tsne")
    
    metadata(treeviz)$tsne <- reducdim@cell.embeddings
    rownames(metadata(treeviz)$tsne)<- colnames(object)
  }
  treeviz
}


#' Creates a `TreeViz`` object from `SingleCellExperiment`
#'
#' @param object `SingleCellExperiment` class containing cluster information at different resolutions
#' @return `TreeViz` Object
#' @export
#'
createFromSCE <- function(object) {
  clusterdata <- colData(object)
  clusterdata <-
    clusterdata[, grep("(cluster|sc3_)", colnames(clusterdata), ignore.case = TRUE)]
  
  count <- counts(object)
  rownames(count) <- rownames(counts(object))
  
  treeviz <- preprocessAndCreateTreeViz(as.data.frame(clusterdata), count)
  
  if ("TSNE" %in% reducedDimNames(object)) {
      metadata(treeviz)$tsne <- reducedDims(object)$"TSNE"
      rownames(metadata(treeviz)$tsne)<- colnames(object)
  }
  treeviz
}

#' Creates `TreeViz` object from hierarchy and count matrix
#'
#' This module returns the Tree Summarized Experiment if the dataframe is tree, throws error otherwise
#'
#' @param clusters Dataframe containing cluster information at different resolutions
#' @param counts matrix Dense or sparse matrix containing the count matrix
#' @return `TreeViz`` Object
#' @export
#'
createTreeViz <- function(clusters, counts) {
  clusters <- rename_clusters(clusters)
  
  for (clusnames in names(clusters)) {
    clusters[[clusnames]] <-
      paste(clusnames, clusters[[clusnames]], sep = 'C')
  }
  
  check_unique_parent(clusters)
  samples <- rownames(clusters)
  clusters <- cbind(clusters, samples)
  
  clusters <- checkRoot(clusters)
  
  tree <- TreeIndex(clusters)
  rownames(tree) <- rownames(clusters)
  
  treeviz <- TreeViz(SimpleList(counts = counts), colData = tree)
  # plot_tree(TreeSE_obj@colData@hierarchy_tree)
  treeviz
}


#' Finds top variable genes `TreeViz`
#'
#' Finds the n top variable genes within the genes present in `TreeViz` object
#'
#' @param treeVizobject TreeViz object
#' @param top number of top genes to be calculated
#' @return `TreeViz` Object with added top_variable_gene information in metadata
#' 
#' @import scran
#' @export
#'
find_top_variable_genes <- function(treeviz, top=100) {
  dec.treeviz <- modelGeneVar(assays(treeviz)$counts)
  top_n <- getTopHVGs(dec.treeviz, n = top)
  metadata(treeviz)$top_variable <- top_n
  
  treeviz
}
