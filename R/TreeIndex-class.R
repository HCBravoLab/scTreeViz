setClassUnion("CharacterOrNull", c("character", "NULL"))

#' TreeIndex class to manage and query hierarchical data
setClass(
  "TreeIndex",
  contains = c("DataFrame"),
  representation(
    feature_order = "CharacterOrNull",
    leaf_of_table = "data.table",
    hierarchy_tree = "data.frame",
    node_ids_table = "data.table",
    nodes_table = "data.table"
  )
)

#' create a new TreeIndex object
#' @param hierarchy hierarchy as a data.table
#' @param feature_order order of the tree if different from colnames
#' @importFrom methods new
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table setorderv
#' @importFrom data.table melt
#' @importFrom data.table setnames
#' @importFrom S4Vectors DataFrame
#' @import digest
#' @export
TreeIndex <- function(hierarchy = NULL,
                      feature_order = NULL) {
  if (is.null(hierarchy)) {
    return(
      new(
        "TreeIndex",
        DataFrame(),
        feature_order = feature_order,
        leaf_of_table = data.table(),
        hierarchy_tree = data.frame(),
        node_ids_table = data.table(),
        nodes_table = data.table()
      )
    )
  }

  if (ncol(hierarchy) == 0) {
    return(
      new(
        "TreeIndex",
        DataFrame(hierarchy),
        feature_order = feature_order,
        leaf_of_table = data.table(),
        hierarchy_tree = data.frame(),
        node_ids_table = data.table(),
        nodes_table = data.table()
      )
    )
  }

  if (is.null(feature_order)) {
    feature_order <- colnames(hierarchy)
  }

  hierarchy_tree <-
    .generate_hierarchy_tree(hierarchy, feature_order)
  node_ids_table <-
    .generate_node_ids(hierarchy_tree, feature_order)
  nodes_table <-
    .generate_nodes_table(hierarchy_tree, node_ids_table, feature_order)
  leaf_of_table <-
    .generate_leaf_of_table(hierarchy_tree, node_ids_table, nodes_table, feature_order)

  hierarchy_df <- DataFrame(hierarchy)

  new(
    "TreeIndex",
    hierarchy_df,
    feature_order = feature_order,
    leaf_of_table = leaf_of_table,
    hierarchy_tree = hierarchy_tree,
    node_ids_table = node_ids_table,
    nodes_table = nodes_table
  )
}

.generate_hierarchy_tree <- function(hierarchy, feature_order) {
  fd <- hierarchy
  for (i in seq(ncol(fd))) {
    fd[, i] = as.character(fd[, i])
  }
  hierarchy <- fd

  replacing_na_obj_fData <- hierarchy[, feature_order]

  nas_replaced <-
    .replaceNAFeatures(replacing_na_obj_fData, feature_order)

  obj_fData <- as.data.table(nas_replaced)
  cols <- feature_order[1:length(feature_order) - 1]
  order <- rep(1, length(feature_order) - 1)
  ordered_fData <- setorderv(obj_fData, cols = cols, order = order)

  otu_indexes <-
    seq(1:length(ordered_fData[, get(feature_order[length(feature_order)])]))
  ordered_fData <- ordered_fData[, otu_index := otu_indexes]
  ordered_fData_df <- as.data.frame(ordered_fData)

  if (length(unique(ordered_fData_df[, 1])) > 1) {
    allFeatures <- rep("AllFeatures", nrow(ordered_fData_df))
    ordered_fData_df <- cbind(allFeatures, ordered_fData_df)
    feature_order <- unlist(c("allFeatures", feature_order))
  }

  ordered_fData_df
}

.generate_node_ids <- function(hierarchy_tree, feature_order) {
  table_node_ids <- hierarchy_tree
  id_list <- sapply(feature_order, function(level) {
    depth <- which(feature_order == level)
    temp_level <-
      data.table(table_node_ids[, c(level, "otu_index")])
    temp_level_count <-
      temp_level[, .(leaf_index = .I[which.min(otu_index)], count = .N), by =
                   eval(level)]

    level_features <- as.character(table_node_ids[[level]])
    for (i in seq_len(nrow(temp_level_count))) {
      row <- temp_level_count[i, ]
      if (depth == 1 && i == 1) {
        id <- paste(depth - 1, 0, sep = "-")
      } else{
        id <-
          paste(depth - 1, paste(digest(row[, 1], algo = "crc32"), i, sep = ""), sep =
                  "-")
      }
      level_features <-
        replace(level_features, which(level_features == row[[level]]), id)
    }
    level_features
  })

  node_ids_dt <- as.data.table(id_list)
  node_ids_dt$otu_index <- as.character(table_node_ids$otu_index)

  node_ids_table <- node_ids_dt
}

.generate_nodes_table <-
  function(hierarchy_tree,
           node_ids_table,
           feature_order) {
    lineage_DF <- as.data.frame(node_ids_table)
    lineage_table <- node_ids_table
    lineage_DF[, feature_order[1]] <-
      lineage_table[, get(feature_order[1])]

    for (i in seq(2, length(feature_order))) {
      lineage_DF[, feature_order[i]] <-
        paste(lineage_DF[, feature_order[i - 1]], lineage_table[, get(feature_order[i])], sep =
                ",")
    }
    lineage_DT <- as.data.table(lineage_DF)

    root_parents <-
      rep("None", length(node_ids_table[, get(feature_order[1])]))
    nodes_tab <-
      data.frame(
        id = node_ids_table[, get(feature_order[1])],
        parent = root_parents,
        lineage = node_ids_table[, get(feature_order[1])],
        node_label = hierarchy_tree[, 1],
        level = rep(1, length(hierarchy_tree[, 1]))
      )

    for (i in seq(2, length(feature_order))) {
      temp_nodes_tab <-
        data.frame(
          id = node_ids_table[, get(feature_order[i])],
          parent = node_ids_table[, get(feature_order[i -
                                                        1])],
          lineage = lineage_DT[, get(feature_order[i])],
          node_label = hierarchy_tree[, i],
          level = rep(i, length(hierarchy_tree[, i]))
        )

      nodes_tab <-
        rbind(nodes_tab[rownames(unique(nodes_tab[, c("id", "parent")])), ], temp_nodes_tab[rownames(unique(temp_nodes_tab[, c("id", "parent")])), ])
    }

    ret_table <- as.data.table(nodes_tab)
    ret_table <- ret_table[, id := as.character(id)]
    ret_table <- ret_table[, parent := as.character(parent)]
    ret_table <- ret_table[, lineage := as.character(lineage)]
    ret_table <- ret_table[, node_label := as.character(node_label)]
    ret_table <- ret_table[, level := as.integer(level)]

    ret_table <- ret_table[order(parent)]
    parent_list <- ret_table[, parent]
    orders <- rep(1, length(parent_list))

    for (j in seq(2, length(parent_list))) {
      if (parent_list[j] == parent_list[j - 1]) {
        orders[j] = orders[j - 1] + 1
      }
    }
    ret_table[, order := orders]

    ret_table
  }

.generate_leaf_of_table <-
  function(hierarchy_tree,
           node_ids_table,
           nodes_table,
           feature_order) {
    temp_hiearchy_DT <- as.data.table(hierarchy_tree)
    num_features <- length(feature_order)
    hiearchy_cols <- colnames(hierarchy_tree)

    melt_res <-
      melt(
        temp_hiearchy_DT,
        id.vars = c(feature_order[num_features], "otu_index"),
        measure.vars = c(hiearchy_cols[1:(length(hiearchy_cols) -
                                            1)])
      )
    label_table <- melt_res[, c(1, 2, 4)]
    setnames(label_table, c("leaf", "otu_index", "node_label"))

    label_table <- label_table[, leaf := as.character(leaf)]
    label_table <-
      label_table[, otu_index := as.character(otu_index)]

    lineage_DF <- as.data.frame(node_ids_table)
    lineage_table <- node_ids_table
    lineage_DF[, feature_order[1]] <-
      lineage_table[, get(feature_order[1])]

    for (i in seq(2, length(feature_order))) {
      lineage_DF[, feature_order[i]] <-
        paste(lineage_DF[, feature_order[i - 1]], lineage_table[, get(feature_order[i])], sep =
                ",")
    }
    lineage_DT <- as.data.table(lineage_DF)

    melt_res_lineage <-
      melt(
        lineage_DT,
        id.vars = c(feature_order[num_features], "otu_index"),
        measure.vars = c(hiearchy_cols[1:(length(hiearchy_cols)) - 1])
      )

    lineage_leaf_of_table <- unique(melt_res_lineage[, c(2, 4)])
    setnames(lineage_leaf_of_table, c("otu_index", "lineage"))

    lineage_leaf_of_table <-
      lineage_leaf_of_table[, otu_index := as.character(otu_index)]

    lineage_df <- as.data.frame(lineage_leaf_of_table)
    leaf_node_label <-
      as.data.frame(label_table)[, c("leaf", "node_label")]

    ret_table <- as.data.table(cbind(lineage_df, leaf_node_label))

    leaf_of_table <- ret_table

    leaf_of_table <-
      merge(unique(nodes_table[, mget(c("lineage", "id"))]),
            unique(leaf_of_table) , by = "lineage")
    leaf_of_table[, id := as.character(id)]
  }

.replaceNAFeatures = function(replacing_na_obj_fData, feature_order) {
  for (i in seq(1, length(feature_order))) {
    na_indices <-
      which(is.na(replacing_na_obj_fData[, feature_order[i]]))
    for (j in seq(1, length(na_indices))) {
      if (i > 1) {
        replacing_na_obj_fData[, feature_order[i]][na_indices[j]] <-
          paste("Not_Annotated",
                feature_order[i],
                # replacing_na_obj_fData[, feature_order[1]][na_indices[j]],
                replacing_na_obj_fData[, feature_order[i-1]][na_indices[j]],
                sep = "_")
      } else {
        replacing_na_obj_fData[, feature_order[i]][na_indices[j]] <-
          paste("Not_Annotated", feature_order[i], sep = "_")
      }
    }
    na_indices <-
      which(replacing_na_obj_fData[, feature_order[i]] == "NA")
    for (j in seq(1, length(na_indices))) {
      if (i > 1) {
        replacing_na_obj_fData[, feature_order[i]][na_indices[j]] <-
          paste("Not_Annotated",
                feature_order[i],
                # replacing_na_obj_fData[, feature_order[1]][na_indices[j]],
                replacing_na_obj_fData[, feature_order[i-1]][na_indices[j]],
                sep = "_")
      } else{
        replacing_na_obj_fData[, feature_order[i]][na_indices[j]] <-
          paste("Not_Annotated", feature_order[i], sep = "_")
      }
    }
    null_indices <-
      which(replacing_na_obj_fData[, feature_order[i]] == "NULL")
    for (j in seq(1, length(null_indices))) {
      if (i > 1) {
        replacing_na_obj_fData[, feature_order[i]][null_indices[j]] <-
          paste("Not_Annotated",
                feature_order[i],
                # replacing_na_obj_fData[, feature_order[1]][na_indices[j]],
                replacing_na_obj_fData[, feature_order[i-1]][na_indices[j]],
                sep = "_")
      } else{
        replacing_na_obj_fData[, feature_order[i]][null_indices[j]] <-
          paste("Not_Annotated", feature_order[i], sep = "_")
      }
    }
  }

  replacing_na_obj_fData
}
