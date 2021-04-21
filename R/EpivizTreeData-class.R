#' Data container for MRexperiment objects
#'
#' Used to serve hierarchical data (used in e.g., icicle plots and heatmaps).
#' @importClassesFrom epivizrData EpivizData
#' @importClassesFrom epivizrData SparseEpivizMeasurement
#' @import data.table
#' @import digest
#' @import methods
#' @import httr
#' @import SummarizedExperiment
#' @import Rtsne
#' @exportClass EpivizTreeData
EpivizTreeData <- setRefClass("EpivizTreeData",
  contains = "EpivizData",
  fields = list(
    .levels = "ANY",
    .maxDepth = "numeric",
    .feature_order = "character",
    .minValue = "numeric",
    .maxValue = "numeric",
    .sampleAnnotation = "ANY",
    .nodeSelections = "ANY",
    .levelSelected = "ANY",
    .lastRootId = "character",
    .json_query = "ANY",
    .graph = "ANY",
    .treeIn = "character",
    .measurements = "ANY"
  ),
  
  methods=list(
    initialize=function(object, tree="row", columns=NULL, control=metavizControl(), feature_order=NULL, ...) {
      
      if(is.null(feature_order)) {
        if(tree == "row") {
          .self$.feature_order = colnames(rowData(object))
        } else {
          .self$.feature_order = colnames(colData(object))
        }
      } else {
        .self$.feature_order <- feature_order
      }
      
      .self$.levelSelected <- 3
      .self$.lastRootId <- "0-0"
      .self$.nodeSelections <- list()
      .self$.treeIn <- tree
      .self$.measurements <- NULL
      
      if(tree == "row") {
        .self$.graph <- rowData(object)
      } else {
        .self$.graph <- colData(object)
      }
      
      featureSelection = NULL
      
      if(!is.null(featureSelection)) {
        featureSelection <- featureSelection[which(names(featureSelection) != "NA")]
        featureSelection <- featureSelection[which(names(featureSelection) != "no_match")]
        
        if(tree == "row") {
          node_ids <- sapply(names(featureSelection), function(n) {
            as.character(rowData(object)$nodes_table[node_label==n,id])
          })
        } else {
          node_ids <- sapply(names(featureSelection), function(n) {
            as.character(colData(object)$nodes_table[node_label==n,id])
          })
        }
        
        temp_selections <- unname(featureSelection)
        names(temp_selections) <- node_ids
        .self$.nodeSelections = temp_selections
      }
      
      callSuper(object=object, ...)
    }
  )
)

# Epiviz Websockets Protocol
EpivizTreeData$methods(
  get_default_chart_type=function() {
    "epiviz.ui.charts.tree.Icicle"
  },
  
  get_measurements=function(top_genes=TRUE) {
    "Get all annotation info for all samples
    \\describe{
    \\item{chart_id_or_object}{An object of class \\code{EpivizChart} or an id for
    a chart loaded to the epiviz app.}
    }
    "
    
    if(top_genes) {
      
      if (!is.null(.self$.measurements)) {
        out <- .self$.measurements
      }
      else {
        genes <- metadata(.self$.object)$top_variable
        
        out <- lapply(genes, function(gene) {
          epivizrData:::SparseEpivizMeasurement(id=gene,
                                                datasourceId=.self$.id)
        })
        
        .self$.measurements <- out  
      }
      return(out)
    }
    
    if (.self$.treeIn == "row") {
      sample_table <- as.data.frame(colData(.self$.object))
      # .sampleAnnotation <- colData(.self$.object)
    } else {
      sample_table <- as.data.frame(rowData(.self$.object))
      # .sampleAnnotation <- rowData(.self$.object)
    }
    
    if (!is.null(.self$.measurements)) {
      out <- .self$.measurements
    }
    else {
      out <- lapply(rownames(sample_table), function(sample) {
        epivizrData:::SparseEpivizMeasurement(id=sample,
          datasourceId=.self$.id)
      })
      
      .self$.measurements <- out 
    }
    
    return(out)
  },
  
  row_to_dict=function(row){
    "Helper function to format each node entry for getHierarchy response
    \\describe{
    \\item{row}{Information for current node.}
    }
    "
    
    toRet = list()
    toRet['end'] = row['end']
    toRet['partition'] = "NA"
    toRet['leafIndex'] = row['leafIndex']
    toRet['nchildren'] = row['nchildren']
    toRet['label'] = row['label']
    toRet['name'] = row['label']
    toRet['start'] = row['start']
    toRet['depth'] = row['depth']
    toRet['globalDepth'] = row['depth']
    toRet['nleaves'] = row['nleaves']
    toRet['parentId'] = row['parentId']
    toRet['order'] = row['order']
    toRet['id'] = row['id']
    # if(toRet['id'] %in% names(.self$.nodeSelections)){
    #   toRet['selectionType'] = .self$.nodeSelections[[as.character(toRet['id'])]]
    # } else{
    #   toRet['selectionType'] = 1
    # }
    toRet['selectionType'] = 1
    toRet['taxonomy'] = row['taxonomy']
    toRet['size'] = 1
    toRet['children'] = NULL
    return(toRet)
  },
  
  df_to_tree=function(root, df){
    "Helper function to recursively build nested response for getHierarchy
    \\describe{
    \\item{root}{Root of subtree}
    \\item{df}{data.frame containing children to process}
    }
    "
    
    if(nrow(df) == 0) {
      root$children = NULL
      return(root)
    }
    
    children = df[which(df['parentId'] == as.character(unlist(root['id']))),]
    
    if(length(children) == 0){
      root$children = NULL
      return(root)
    }
    
    otherChildren = df[which(df['parentId'] != as.character(unlist(root['id']))),]
    
    children = children[order(children['order']),]
    
    if(nrow(children) > 0){
      for(row_index in seq_len(nrow(children))){
        childDict = row_to_dict(children[row_index,])
        subDict = df_to_tree(childDict, otherChildren)
        
        if(!is.null(subDict)){
          root$children[[row_index]] = subDict
        }
        else {
          root$children = NULL
        }
      }
    }
    return(root)
  },
  
  getHierarchy=function(nodeId = NULL) {
    "Retrieve feature hierarchy information for subtree with specified root
    \\describe{
    \\item{nodeId}{Feature identifier with level info}
    }
    "
    
    # getHierarchy can be called with NULL from App
    if(is.null(nodeId) || nodeId == ""){
      nodeId <- .self$.lastRootId
    }
    .self$.lastRootId <- nodeId
    root <- nodeId
    
    #Split the node id to get level and index
    split_res <- strsplit(nodeId, "-")[[1]]
    level <- as.integer(split_res[1])+1
    index <- which(.self$.graph@node_ids_table[,level, with=FALSE] == nodeId)
    
    graph_tree <- .self$.graph@hierarchy_tree[,-which(colnames(.self$.graph@hierarchy_tree) == "otu_index")]
    
    label <- as.character(unique(graph_tree[,level][index]))
    taxonomy <- colnames(graph_tree)[level]
    
    if(length(.self$.graph@feature_order) >= level+3){
      last_level_of_subtree <- level+3
    } else{
      last_level_of_subtree <- length(.self$.graph@feature_order)
    }
    
    hierarchy_slice <- unique(.self$.graph@node_ids_table[get(taxonomy)==nodeId, (level+1):last_level_of_subtree])
    
    nodes_of_subtree <- sapply(seq(1,length((level+1):last_level_of_subtree)), function(i) {
      unname(unlist(unique(hierarchy_slice[,i, with=FALSE])))
    })
    
    nodes_of_subtree <- unlist(nodes_of_subtree)
    
    if(level == 0 || nodeId == "0-0"){
      nodesToRet <- c(root, unlist(nodes_of_subtree))
    } else{
      parent_of_root_taxonomy <- colnames(graph_tree)[(level-1)]
      parent_of_root <- unique(.self$.graph@node_ids_table[get(taxonomy)==nodeId, get(parent_of_root_taxonomy)])
      nodesToRet <- c(parent_of_root,root, unlist(nodes_of_subtree))
    }
    
    num_rows <- length(nodesToRet)
    
    starts <- rep(1, num_rows)
    labels <- rep(1, num_rows)
    leafIndexes <- rep(1, num_rows)
    parentIds <- rep(1, num_rows)
    depths <- rep(0, num_rows)
    partitions <- rep(1, num_rows)
    ends <- rep(1, num_rows)
    ids <- rep(1, num_rows)
    nchildrens <- rep(1, num_rows)
    taxonomys <- rep(1, num_rows)
    nleaves <- rep(1, num_rows)
    orders <- rep(1, num_rows)
    
    leaf_ordering_table <- as.data.table(.self$.graph@hierarchy_tree[,c(.self$.graph@feature_order[length(.self$.graph@feature_order)], "otu_index")])
    setnames(leaf_ordering_table, c("leaf", "otu_index"))
    
    lineage_DF <- as.data.frame(.self$.graph@node_ids_table)
    lineage_table <- .self$.graph@node_ids_table
    lineage_DF[,.self$.graph@feature_order[1]] <- lineage_table[,get(.self$.graph@feature_order[1])]
    
    for(i in seq(2,length(.self$.graph@feature_order))){
      lineage_DF[,.self$.graph@feature_order[i]] <- paste(lineage_DF[,.self$.graph@feature_order[i-1]], lineage_table[,get(.self$.graph@feature_order[i])], sep=",")
    }
    lineage_DT <- as.data.table(lineage_DF)
    
    lineage_DT_long <- melt(lineage_DT, id.vars = c("otu_index"),
                            measure.vars = c(colnames(lineage_DT)[1:(length(colnames(lineage_DT))-1)]),
                            variable.name = "taxonomy", variable.factor = FALSE)
    
    setnames(lineage_DT_long, c("otu_index", "taxonomy", "lineage"))
    
    temp_nodes_table <- merge(.self$.graph@nodes_table, lineage_DT_long, by="lineage")
    temp_nodes_table <- temp_nodes_table[, otu_index:=as.integer(otu_index)]
    
    
    for(i in seq_len(num_rows)){
      if(as.integer(strsplit(nodesToRet[i], "-")[[1]][1]) == last_level_of_subtree){
        depths[i] = length(.self$.graph@feature_order)
        level = length(.self$.graph@feature_order)
        index <- which(.self$.graph@node_ids_table[,level,with=FALSE] == nodeId)
        
        label <- as.character(unique(graph_tree[,level][index]))
        labels[i] <- label
        
        partition <- "NA"
        partitions[i] <- partition
        
        otu_index_temp <- leaf_ordering_table[leaf == nodeId, otu_index]
        starts[i] <- otu_index_temp
        leafIndexes[i] <- otu_index_temp
        ends[i] <- otu_index_temp
        
        ids[i] <- nodeId
        
        taxonomy <- colnames(graph_tree)[level]
        taxonomys[i] <- taxonomy
        
        nchildrens[i] <- 0
        nleaves[i] <- 0
        
        orders[i] <- .self$.graph@nodes_table[get("id")==nodeId,get("order")][[1]]
      } else{
        nodeId <- nodesToRet[i]
        split_res <- strsplit(nodesToRet[i], "-")[[1]]
        depths[i] <- as.integer(split_res[1])
        level <- as.integer(split_res[1])+1
        
        index <- which(.self$.graph@node_ids_table[,level,with=FALSE] == nodeId)
        
        label <- as.character(unique(graph_tree[,level][index]))
        labels[i] <- label
        taxonomy <- colnames(graph_tree)[level]
        
        if(nodesToRet[i] != "0-0"){
          parentId_taxonomy <- colnames(graph_tree)[(level-1)]
          parentId <- unique(.self$.graph@node_ids_table[get(taxonomy)==nodesToRet[i], get(parentId_taxonomy)])[1]
          parentIds[i] <- parentId
        } else{
          parentIds[i] <- "NA"
        }
        
        partition <- "NA"
        partitions[i] <- partition
        
        if(nodesToRet[i] == "0-0"){
          leaf_indexes_temp <- temp_nodes_table[level == (length(.self$.graph@feature_order)-1), otu_index,]
        } else{
          ids_match <- .self$.graph@nodes_table[id == nodesToRet[i],]
          parents_match <- ids_match[parent == parentIds[i]]
          leaf_indexes_temp <- temp_nodes_table[lineage == parents_match[,lineage,], otu_index,]
        }
        #leaf_indexes_temp <- temp_nodes_table[lineage == .self$.graph@nodes_table[id == nodesToRet[i],lineage,], otu_index,]
        if(length(leaf_indexes_temp) > 0){
          start <- min(leaf_indexes_temp)
        }    else{
          start <- nodesToRet[i]
        }
        
        if(length(leaf_indexes_temp) > 0){
          end <- max(leaf_indexes_temp)
        }    else{
          end <- nodesToRet[i]
        }
        
        starts[i] <- start
        
        leafIndex <- start
        leafIndexes[i] <- leafIndex
        ends[i] <- end
        
        id <- nodesToRet[i]
        ids[i] <- id
        
        taxonomy <- colnames(graph_tree)[level]
        taxonomys[i] <- taxonomy
        
        nchildren <- length(unique(.self$.graph@node_ids_table[get(taxonomy)==nodesToRet[i],][[as.integer(level)+1]]))
        nchildrens[i] <- nchildren[1]
        
        nleaves_temp <- length(unname(unlist(unique(.self$.graph@leaf_of_table[node_label==label, leaf]))))
        nleaves[i] <- nleaves_temp[1]
        
        if(nodesToRet[i] != "0-0"){
          orders[i] <- .self$.graph@nodes_table[get("id")==nodesToRet[i],get("order")][[1]]
        } else {
          orders[i] <- 1
        }
      }
    }
    ret_data_frame <- data.frame(start = starts, label = labels, leafIndex = leafIndexes, parentId = parentIds,
                                 depth = depths, partition = partitions, end = ends, id = ids, nchildren = nchildrens,
                                 taxonomy = taxonomys, nleaves = nleaves, order = orders)
    
    if(length(ret_data_frame) > 0){
      # convert columns to int
      ret_data_frame['start'] = as.numeric(unlist(ret_data_frame['start']))
      ret_data_frame['end'] = as.numeric(unlist(ret_data_frame['end']))
      ret_data_frame['order'] = as.numeric(unlist(ret_data_frame['order']))
      ret_data_frame['leafIndex'] = as.numeric(unlist(ret_data_frame['leafIndex']))
      ret_data_frame['nchildren'] = as.numeric(unlist(ret_data_frame['nchildren']))
      ret_data_frame['nleaves'] = as.numeric(unlist(ret_data_frame['nleaves']))
      ret_data_frame['depth'] = as.numeric(unlist(ret_data_frame['depth']))
      ret_data_frame['id'] = as.character(unlist(ret_data_frame['id']))
      
      
      root = ret_data_frame[1,]
      rest = ret_data_frame[-1,]
      rootDict = .self$row_to_dict(root)
      result = .self$df_to_tree(rootDict, rest)
      
      result[["rootTaxonomies"]] = .self$.graph@feature_order
      lineage = .self$.graph@nodes_table[get("id")==nodesToRet[1],get("lineage")][[1]]
      
      lineageLabel <- sapply(strsplit(lineage, ",")[[1]], function(str_id) {
        .self$.graph@nodes_table[get("id") == str_id, get("node_label")][[1]]
      })
      
      result[["lineageLabel"]] = paste(lineageLabel, sep=", ")
      
      resultResp = list(nodeSelectionTypes = .self$.nodeSelections,
                        selectionLevel = .self$.levelSelected,
                        tree = result)
      
      return(resultResp)
    }
    
    return(ret_data_frame)
  },
  
  propagateHierarchyChanges=function(selection = NULL, order = NULL, selectedLevels = NULL, request_with_labels = FALSE) {
    "Update internal state for hierarchy
    \\describe{
    \\item{selection}{Node-id and selectionType pairs}
    \\item{order}{Ordering of features}
    \\item{selectedLevels}{Current aggregation level}
    \\item{request_with_labels}{For handling requests using fData entries from MRexperiment}
    }
    "
    
    if(request_with_labels && !is.null(selection)){
      selection_ids <- sapply(names(selection), function(i){
        .self$.graph@nodes_table[node_label==i,id]
      })
      names(selection) <- selection_ids
    }
    
    # update node selections types to metaviztree
    if(!is.null(selection)) {
      for(n in names(selection)){
        .self$.nodeSelections[[n]] = selection[[n]]
      }
    }
    
    if(!is.null(selectedLevels)) {
      .self$.levelSelected <- as.integer(names(selectedLevels)[1])
    }
    .self$.mgr$.clear_datasourceGroup_cache(.self)
  },
  
  getRows=function(measurements = NULL, start = 1, end = 1000, selectedLevels = 3, selections = NULL) {
    "Return the sample annotation and features within the specified range and level for a given sample and features
    \\describe{
    \\item{measurements}{Samples to retrieve for}
    \\item{start}{Start of feature range to query}
    \\item{end}{End of feature range to query}
    \\item{selections}{Node-id and selectionType pairs}
    \\item{selectedLevels}{Current aggregation level}
    }
    "
    
    nodes_at_level <- .self$.graph@nodes_table[level==selectedLevels, ]
    nodes_at_level_ids <- nodes_at_level[,id]
    
    if(!is.null(selections) && !(length(selections) == 0)){
      nodes_at_level_selections <- rep(2, length(nodes_at_level_ids))
      names(nodes_at_level_selections) <- nodes_at_level_ids
      selections <- c(selections, nodes_at_level_selections)
      
      expand_selections <- which(selections == 1)
      if(!is.null(expand_selections) && length(expand_selections) > 0){
        expand_selection_indices = which(names(selections) %in% names(expand_selections))
        selections <- selections[-expand_selection_indices]
      }
      
      expanded_children <- .self$.graph@nodes_table[parent %in% names(expand_selections),id]
      
      child_lineage <- .self$.graph@nodes_table[id %in% names(selections),]
      remove_selections <- which(selections == 0)
      if(length(remove_selections) > 0){
        kept_nodes <- child_lineage[!grepl(paste(paste(names(remove_selections), collapse=",|"), ",",sep=""), lineage),]
        kept_nodes <- kept_nodes[!(id %in% names(remove_selections)),]
      } else {
        kept_nodes <- child_lineage
      }
      
      agg_selections <- which(selections == 2)
      if(length(agg_selections) > 0){
        kept_nodes <- kept_nodes[!grepl(paste(paste(names(agg_selections), collapse=",|"), ",",sep=""), lineage),]
      }
      kept_nodes <- as.character(kept_nodes[,id])
      nodes_at_level <- .self$.graph@nodes_table[id %in% c(kept_nodes,expanded_children),]
    }
    
    leaf_ordering_table <- as.data.table(.self$.graph@hierarchy_tree[,c(.self$.feature_order[length(.self$.feature_order)], "otu_index")])
    setnames(leaf_ordering_table, c("leaf", "otu_index"))
    leaf_ordering_table <- leaf_ordering_table[,leaf:=as.character(leaf)]
    leaf_ordering_table <- leaf_ordering_table[otu_index >= start & otu_index <= end]
    
    first_join <- merge(leaf_ordering_table, merge(nodes_at_level, .self$.graph@leaf_of_table, by="id"), by="leaf")
    setorderv(first_join, "otu_index.x")
    
    nodes <- as.data.frame(unique(first_join[,c("node_label.x", "lineage.x")]))[,"node_label.x"]
    
    ends <- rep(0, length(nodes))
    starts <- rep(0, length(nodes))
    indexes <- rep(0, length(nodes))
    metadata <- list()
    metadata[['label']] <- list()
    metadata[['id']] <- list()
    metadata[['lineage']] <- list()
    
    for (i in seq_along(nodes)) {
      feature_label <- nodes[i]
      res <- first_join[node_label.x==feature_label,otu_index.x]
      if(length(res) > 0) {
        ends[i] <- max(res)
        starts[i] <- min(res)
        # if(i == length(nodes)) {
        #   starts[i] <- min(res) - 1
        # }
        indexes[i] <- min(res)
        metadata[['label']][i] <- feature_label
        metadata[['id']][i] <- unique(.self$.graph@nodes_table[node_label==nodes[i], id])[1]
        metadata[['lineage']][i] <- unique(.self$.graph@nodes_table[node_label==nodes[i], lineage])[1]
      }
    }
    
    data_rows = list(end=ends, start=starts, index =indexes, metadata=metadata)
    return(data_rows)
  },
  
  searchTaxonomy=function(query = NULL, max_results = 15) {
    "Return list of features matching a text-based query
    \\describe{
    \\item{query}{String of feature for which to search}
    \\item{max_results}{Maximum results to return}
    }
    "
    return(list())
  },
  
  getCombined=function(measurements = NULL,
                       seqName, start = 1, end = 1000,
                       order = NULL, nodeSelection = NULL, selectedLevels = NULL) {
    "Return the counts aggregated to selected nodes for the given samples
    \\describe{
    \\item{measurements}{Samples to get counts for}
    \\item{seqName}{name of datasource}
    \\item{start}{Start of feature range to query}
    \\item{end}{End of feature range to query}
    \\item{order}{Ordering of nodes}
    \\item{nodeSelection}{Node-id and selectionType pairs}
    \\item{selectedLevels}{Current aggregation level}
    }
    "
    
    # update node selections types to tree
    if(!is.null(nodeSelection)) {
      for(n in names(nodeSelection)){
        .self$.nodeSelections[[n]] = nodeSelection[[n]]
      }
    }
    
    if(is.null(selectedLevels)) {
      selectedLevels = .self$.levelSelected + 1
    }
    
    selections = .self$.nodeSelections
    measurements = unique(measurements)
    
    data_rows = .self$getRows(measurements = measurements, start = start, end = end, selectedLevels = selectedLevels, selections = selections)
    row_order = unlist(data_rows$metadata$label)
    aggcounts = aggregateTree(.self$.object, selectedLevel=selectedLevels, selectedNodes=selections, start = start, end = end, by=.self$.treeIn, format="counts")
    
    data_columns = list()
    
    if (.self$.treeIn == "row") {
      for(m in measurements) {
        if(m %in% colnames(aggcounts)){
          inner_result <- aggcounts[,m]
          data_columns[[m]] <- unname(unlist(inner_result))
        } else {
          inner_result <- rep(0.0, nrow(aggcounts))
          data_columns[[m]] <- inner_result
        }
      }
    } else if (.self$.treeIn == "col") {
      for(m in measurements) {
        if(m %in% rownames(aggcounts)){
          inner_result <- aggcounts[m,]
          data_columns[[m]] <- unname(unlist(inner_result))
        } else {
          inner_result <- rep(0.0, ncol(aggcounts))
          data_columns[[m]] <- inner_result
        }
      } 
    }
    
    result <- list(
      cols = data_columns,
      rows = data_rows,
      globalStartIndex = data_rows$start[[1]]
    )
    
    return(result)
  },
  
  getReducedDim=function(method = NULL, gene = NULL) {
    " Compute PCA over all features for given samples
    \\describe{
    \\item{method}{which dimension to access}
    \\item{gene}{send expression of a gene back with the dimensions}
    }
    "
    
    if (is.null(method)) {
      method <- names(metadata(.self$.object)$reduced_dim)[[1]]
    }
    
    data_rows = .self$getRows(measurements = NULL, start = 1, end = 100000,
                              selectedLevels = .self$.levelSelected + 1,
                              selections = .self$.nodeSelections)
    
    max_length <- ncol(.self$.object)
    # max(data_rows$end)
    
    cluster_names <- rep("removed", max_length)
    for (i in 1:length(data_rows$metadata$label)) {
      start <- data_rows$start[i]
      end <- data_rows$end[i]
      cluster_names[start:end] <- data_rows$metadata$label[i]
    }
    
    measurements <-  metadata(.self$.object)$reduced_dim[[method]]
    
    genes <- rep(100, length(measurements))
    if (!(is.null(gene)) && gene != "") {
      genes <- assays(.self$.object)$counts[gene, ]
    }
    
    data <- list()
    level <- .self$.levelSelected + 1
    i <- 1
    for (col in rownames(measurements)) {
      row_index = which(colData(.self$.object)$samples == col)
      
      # TODO: need to add sample attributes
      temp    <-
        list(
          sample_id = col,
          dim1 = unname(measurements[col, 1]),
          dim2 = unname(measurements[col, 2]),
          name = unlist(cluster_names[[row_index]]),
          gene = unname(genes[row_index])
          # name = name
        )
      data[[col]] <- temp
      i <- i+1
    }
    
    result <- list("data" = unname(data), "pca_variance_explained" = c(1,1),
                   "cluster_order" = data_rows$metadata$label,
                   "gene_min_max" = c(min(genes), max(genes)))
    return(result)
  },
  
  extract_SCE_epiviz = function(cluster_name="treeviz_clusters") {

    data_rows = .self$getRows(measurements = NULL, start = 1, end = 100000,
                              selectedLevels = .self$.levelSelected + 1,
                              selections = .self$.nodeSelections)
    
    max_length <- ncol(.self$.object)
    # max(data_rows$end)
    
    cluster_names <- rep("removed", max_length)
    for (i in 1:length(data_rows$metadata$label)) {
      start <- data_rows$start[i]
      end <- data_rows$end[i]
      cluster_names[start:end] <- data_rows$metadata$label[i]
    }
    
    # TODO: dumb but whatever, need a shorthand notation
    coldata <- data.frame(treeviz_clusters = unlist(cluster_names))
    colnames(coldata) <- c(cluster_name)
    
    sce <- SingleCellExperiment(assays = list(counts = assays(.self$.object)$counts),
                                colData = coldata,
                                rowData=rowData(.self$.object)
    )
  },
  
  getGeneExpr=function(measurements=NULL){
    
    if (!is.null(measurements)) {
      gene <- measurements[1]
    }
    
    data_rows = .self$getRows(measurements = NULL, start = 1, end = 100000,
                              selectedLevels = .self$.levelSelected + 1,
                              selections = .self$.nodeSelections)
    
    max_length <- ncol(.self$.object)
    # max(data_rows$end)
    
    cluster_names <- rep("removed", max_length)
    for (i in 1:length(data_rows$metadata$label)) {
      start <- data_rows$start[i]
      end <- data_rows$end[i]
      cluster_names[start:end] <- data_rows$metadata$label[i]
    }
    
    level<- .self$.levelSelected + 1
    selectedNodes <- .self$.nodeSelections
    
    counts <- assays(.self$.object)$counts
    
    gene_col <- counts[gene, ]
    data <- list()
    i <- 1
    for (col in names(gene_col)) {
      
      row_index = which(colData(.self$.object)$samples == col)
      
      data[[col]] <- list(
        alphaDiversity = gene_col[[col]],
        sample_id = col,
        name = unlist(cluster_names[[row_index]])
      )
      i <- i+1
    }
    
    result <- list(data = unname(data))
    return(result)
  }
)