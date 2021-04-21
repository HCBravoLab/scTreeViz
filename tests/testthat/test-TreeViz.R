context("testing TreeViz Class")

n=32
df <- data.frame(cluster0=rep(seq(1:2),each=ceiling(n/(2)),len=n))

for(i in seq(1,2)){
  df[[paste0("cluster",2*i-1)]]<- rep(seq(1:(2**i)),each=ceiling(n/(2**i)),len=n)
  df[[paste0("cluster",2*i)]]<- rep(seq(1:(2**i)),each=ceiling(n/(2**i)),len=n)
}

df[[paste0("cluster",2*i+1)]]<- rep(seq(1:7),each=ceiling(n/7),len=n)
counts <- matrix(rpois(3200, lambda = 10), ncol=n, nrow=100)

for(cols in colnames(df)){
  df[[cols]] <- as.factor(df[[cols]])
}

colnames(counts) <- df$samples <- rownames(df) <- paste0("cell-", 1:ncol(counts))
rownames(counts) <- paste0("gene-", 1:nrow(counts))
hierarchy <- ClusterHierarchy(df)

treeviz <- createTreeViz(hierarchy, counts)

test_that("create TreeVizClass", {
  expect_is(treeviz, "TreeViz")

})

test_that("check tree_collapse", {
  expect_equal(6, ncol(colData(treeviz)@hierarchy_tree))
})

test_that("check single_root", {
  expect_equal("ClusterAllClusters", unique(colData(treeviz)@hierarchy_tree[[1]]))
})

test_that("check multiple_parent", {

  hierarchydf <- colData(treeviz)@hierarchy_tree
  hierarchydf <- hierarchydf[, !colnames(hierarchydf) %in% c("samples", "otu_index")]
  
  flag <- TRUE
  
  childs <- unique(hierarchydf[[ncol(hierarchydf)]])
  for (values in childs) {
    subsetted_list <- hierarchydf[hierarchydf[[ncol(hierarchydf)]] == values,]
    parent <- length(unique(subsetted_list[[ncol(hierarchydf) - 1]]))
    if (parent > 1) {
      flag <- FALSE
    }
  }
  
  expect_equal(flag, TRUE)
})

test_that("check aggregate_tree", {
  sums <- list()
  sums[[1]] <- list(rowSums(counts[,c(1:10)]))
  sums[[2]] <- list(rowSums(counts[,c(11:15)]))
  sums[[3]] <- list(rowSums(counts[,c(16:25)]))
  sums[[4]] <- list(rowSums(counts[,c((26:32))]))
  sum_df <- as.data.frame(sums)
  aggtree <- TreeViz:::aggregateTree(treeviz, by="col", selectedLevel=3)
  
  agg_df <- as.data.frame(assays(aggtree)$counts)
  # result <- all.equal(agg_df, sum_df, check.attributes= FALSE)

  expect_equal(dim(agg_df),dim(sum_df))
  result<- all.equal(agg_df, sum_df, check.attributes= FALSE, ignore.col.order=TRUE, ignore.row.order=TRUE)

  expect_equal(result, TRUE)
})

