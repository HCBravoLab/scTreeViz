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
#' @importFrom utils str
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
prune_tree <- function(graph, cluster_df) {
  # Prune the tree so only core edges remain
  
  repeat {
    # drop duplicated columns (Throws error otherwise)
    
    graph <- graph[!duplicated(names(graph))]
    
    # See number of core edges at each iter
    # No False edges acyclic tree
    if (nrow(graph[graph$is_core == FALSE, ]) == 0) {
      ngraph <-
        clustree(
          cluster_df ,
          prefix = "cluster",
          prop_filter = 0,
          return = "graph"
        )
      break
    }
    
    # Apply function
    assign_df <-
      check_alternate(graph[graph$is_core == FALSE, ], graph[graph$is_core == TRUE, ])
    
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
collapse_tree <- function(original_graph) {
  #find all roots
  root_list <- which(sapply(sapply(V(original_graph),
                                   function(x)
                                     neighbors(original_graph, x, mode = "in")), length) == 0)
  
  delete_set_vertices <- vector('numeric')
  for (roots in root_list) {
    node_dists <- distances(original_graph, to = roots)
    layered_dist <- unique(distances(original_graph, to = roots))
    
    layered_dist <- layered_dist[is.finite(layered_dist) == TRUE]
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
  }
  
  if (length(delete_set_vertices) != 0) {
    message("Collapsing the graph by removing redundant nodes")
    ver_list <- ver_list[-delete_set_vertices, ]
  }
  ver_list
}

#' check if hierarchy has more than 1 root, if do create a root
#'
#' @param cluster_df Dataframe containing cluster information at different resolutions
#' @return Dataframe containing cluster information at different resolutions
checkRoot <- function(cluster_df) {
  # Handle Forests
  cols <- colnames(cluster_df)
  if (length(unique(cluster_df[[1]])) > 1) {
    message("Given Tree is a Forest\n Adding extra root")
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
      subsetted_list <-
        clusterdata[clusterdata[[colnames(clusterdata)[[i]]]] == values, ]
      
      
      parent <- length(unique(subsetted_list[[i - 1]]))
      
      if (parent > 1) {
        message(
          "Not a tree, some nodes with multiple parents in level ",
          i,
          "\n Performing Cluster reassignment "
        )
        return (FALSE)
      }
      #cat(nrow(subsetted_list)," ",parent, "\n")
    }
  }
  return(TRUE)
}


#' Renames cluster names for consistency in visualization
#'
#' @param clusterdata Dataframe containing cluster information at different resolutions
#' @return Dataframe containing cluster information at different resolutions
#'
rename_clusters <- function(clusterdata) {
  message("Renaming cluster Levels...")
  message("Previous Level names ", paste(colnames(clusterdata), collapse = "\t"))

  clusnames <- seq(length(colnames(clusterdata)))
  clusnames <- paste0("cluster", clusnames)
  colnames(clusterdata) <- clusnames
  message("New Level names ", paste(colnames(clusterdata), collapse = "\t"))
  clusterdata
}


#' Creates the `TreeViz` object from cluster hierarchy and count matrix.
#' Also prunes and collapses the hierarchy to remove redundant levels.
#'
#' @param clusters Dataframe containing cluster information at different resolutions
#' @param counts matrix Dense or sparse matrix containing the count matrix
#' @return TreeSummarizedExperiment Object that can be visualized with metavizr
#' @import clustree
#' @importFrom data.table setorder
#' @importFrom igraph as_long_data_frame
#' @export
preprocessAndCreateTreeViz <- function(clusters, counts) {
  clusters <- rename_clusters(clusters)
  setorder(clusters)
  counts <- counts[, rownames(clusters)]
  
  # Create clustree object
  clustree_graph <-
    clustree(
      clusters,
      prefix = "cluster",
      prop_filter = 0,
      return = "graph"
    )
  
  graph_df <- as_long_data_frame(clustree_graph)
  
  # prune the graph with only core edges (this makes it a ~tree)
  modified_obj <- prune_tree(graph_df, clusters)
  
  # modified graph and seurat object
  modified_graph <- modified_obj$Clustree_obj
  clusters_new <-  modified_obj$Cluster_obj
  
  # collapses tree if the levels are the same at different resolutions
  collapsed_graph <- collapse_tree(modified_graph)
  cluster_names <- unique(sapply(strsplit(collapsed_graph$node, "C"), '[', 1))
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
  treeviz
}

#' Creates hierarchical clustering on `Single Cell Experiment` object via the walktrap algorithm.
#' walktrap returns clustering at highest modularity. The modularity value indicates quality of cluster division.
#' Intermediate cluster assignments are created based on monotonically increasing level of modularity
#' @param object `Single Cell Experiment` on which `WalkTrap` clustering will be computed
#' @param nsteps number of steps to use in `cluster_walktrap`
#' @return dataframe with hierarchy information
#' @export
#'
generate_walktrap_hierarchy <- function(object, nsteps = 7) {
  
  message("calculating walktrap clusters")
  SNN_Graph <- scran::buildSNNGraph(object)
  clusters <- igraph::cluster_walktrap(SNN_Graph, steps = nsteps)
  modularity <- c()
  
  for (i in 1:length(clusters)) {
    modularity[i] <-
      igraph::modularity(SNN_Graph, igraph::cut_at(clusters, n = i))
  }
  
  monotonic_index <- match(unique(cummax(modularity)), modularity)
  cluster_data =  list()
  for (i in 1:length(monotonic_index)) {
    cluster_data[[i]] =  list(igraph::cut_at(clusters, n = monotonic_index[i]))
  }
  
  cluster_data <- as.data.frame(cluster_data)
  colnames(cluster_data) <- paste0("cluster", monotonic_index)
  cluster_data$samples <- rownames(cluster_data) <- colnames(object)
  cluster_data
}

#' Creates a `TreeViz` object from `Seurat`
#'
#' @param object `Seurat` class containing cluster information at different resolutions
#' @param check_metadata whether to metaData of `Seurat` object for cluster information or not
#' @param col_regex common regular expression shared across all columns with cluster information
#' @param columns vector containing columns with cluster information
#' @param reduced_dim Vector of Dimensionality reduction information provided in `Seurat` object to be added in `TreeViz` (if exists)
#' @return `TreeViz` Object
#' @export
#' @import Seurat
createFromSeurat <- function(object,
   check_metadata = FALSE,
   col_regex = "*snn*",
   columns = NULL,
   reduced_dim = c("TSNE")) {
  
  if (check_metadata==FALSE) {
    message("No default clusters provided")
    object.sce <- as.SingleCellExperiment(object)
    clusterdata <- generate_walktrap_hierarchy(object.sce)
    clusterdata <- ClusterHierarchy(clusterdata)
  }
  else{
    clusterdata <- object@meta.data
    clusterdata$samples <- rownames(clusterdata) <- colnames(object)
    clusterdata <- ClusterHierarchy(clusterdata, col_regex, columns)
  }

  treeviz <- createTreeViz(clusterdata, GetAssayData(object))
  
  for (dim_names in reduced_dim) {
    if (dim_names %in% Reductions(object)) {
      reducdim <- Reductions(object, slot = dim_names)
      
      metadata(treeviz)$reduced_dim[[dim_names]] <- reducdim@cell.embeddings[, 1:2]
      rownames(metadata(treeviz)$reduced_dim[[dim_names]]) <- colnames(object)
    }
  }

  print(names(metadata(treeviz)$reduced_dim))
  treeviz
}

#' Creates a `TreeViz`` object from `SingleCellExperiment`. Generates
#' clusters based on Walktrap algorithm if no default is provided
#' @param object `SingleCellExperiment` object to be visualized
#' @param check_coldata whether to colData of `SingeCellExperiment` object for cluster information or not
#' @param col_regex common regular expression shared across all columns with cluster information
#' @param columns vector containing columns with cluster information
#' @param reduced_dim Vector of Dimensionality reduction information provided in `SingeCellExperiment` object to be added in `TreeViz` (if exists)
#' @return `TreeViz` Object
#' @export
#' @import SingleCellExperiment
createFromSCE <-
  function(object,
           check_coldata = FALSE,
           col_regex = NULL,
           columns = NULL,
           reduced_dim = c("TSNE")) {
    
    if (check_coldata == TRUE) {
      clusterdata <- colData(object)
      clusterdata$samples <- rownames(clusterdata) <- colnames(object)
      clusterdata <- ClusterHierarchy(clusterdata, col_regex, columns)
    }
    else{
      message("No default clusters provided")
      clusterdata <- generate_walktrap_hierarchy(object)
      clusterdata$samples <- rownames(clusterdata) <- colnames(object)
      clusterdata <- ClusterHierarchy(clusterdata)
    }

    count <- counts(object)
    rownames(count) <- rownames(counts(object))
    
    treeviz <- createTreeViz(clusterdata, count)
    for (dim_names in reduced_dim) {
      if (dim_names %in% reducedDimNames(object)) {
        metadata(treeviz)$reduced_dim[[dim_names]] <-
          reducedDims(object)[[dim_names]][, 1:2]
        
        rownames(metadata(treeviz)$reduced_dim[[dim_names]]) <- colnames(object)
      }
    }
    
  treeviz
}

#' Creates `TreeViz` object from hierarchy and count matrix
#' This module returns the Tree Summarized Experiment if the dataframe is tree, throws error otherwise
#' @param clusters `ClusterHierarchy` object containing cluster information at different resolutions
#' @param counts matrix Dense or sparse matrix containing the count matrix
#' @return `TreeViz`` Object
#' @export
#'
createTreeViz <- function(clusters, counts) {
  
  if (class(clusters) %in% "data.frame") {
    clusters$samples <- rownames(clusters) <- colnames(counts)
    clusters <- ClusterHierarchy(clusters)
  }
  
  
  if (!class(clusters) %in% "ClusterHierarchy") {
    stop("clusters is not a ClusterHierarchy class")
  }
  
  counts <- counts[,clusters$samples]
  column_names <- colnames(clusters)
  clusters <- as.data.frame(clusters)
  colnames(clusters) <- column_names
  
  tree <- TreeIndex(clusters)
  rownames(tree) <- rownames(clusters)
  treeviz <- TreeViz(SimpleList(counts = counts), colData = tree)
  treeviz
}


#' Finds top variable genes `TreeViz`
#' Finds the n top variable genes within the genes present in `TreeViz` object
#' @param treeviz TreeViz object
#' @param top number of top genes to be calculated
#' @return `TreeViz` Object with added top_variable_gene information in metadata
#' @importFrom scran getTopHVGs
#' @importFrom scran modelGeneVar
#' @export
#'
find_top_variable_genes <- function(treeviz, top = 100) {
  dec.treeviz <- modelGeneVar(assays(treeviz)$counts)
  top_n <- getTopHVGs(dec.treeviz, n = top)
  metadata(treeviz)$top_variable <- top_n
  
  treeviz
}

#' Sets gene list for visualization
#' @param treeviz TreeViz object
#' @param genes list of genes to use
#' @return TreeViz object set with gene list
set_gene_list <- function(treeviz, genes) {
  genes_in_assay  <- rownames(assays(treeviz)$counts)
  metadata(treeviz)$top_variable <- intersect(genes_in_assay, genes)
  
  treeviz
}

#' Calculates `tsne` Dimensionality Reduction on `TreeViz` object. 
#' The result is added to reduced_dim slot in metadata
#' @param treeviz TreeViz object
#' @return `TreeViz` Object with added 'TSNE'tnse`information in reduced_dim slot of metadata
#' @importFrom scater calculateTSNE
#' @export
calculate_tsne <- function(treeviz) {
  message("No defaults dimensionality reductions provided")
  message("Calculating TSNE")
  tsne <- calculateTSNE(assays(treeviz)$counts)
  metadata(treeviz)$reduced_dim[['TSNE']] <- tsne[, 1:2]
  rownames(metadata(treeviz)$reduced_dim[['TSNE']]) <- colnames(treeviz)
  message("adding tsne to reduced dim slots")
  treeviz
}