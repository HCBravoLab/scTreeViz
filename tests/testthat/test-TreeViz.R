context("testing TreeViz Class")


n=32
df<- data.frame(cluster0=rep(seq(1:2),each=ceiling(n/(2)),len=n))
for(i in seq(1,2)){
  df[[paste0("cluster",2*i-1)]]<- rep(seq(1:(2**i)),each=ceiling(n/(2**i)),len=n)
  df[[paste0("cluster",2*i)]]<- rep(seq(1:(2**i)),each=ceiling(n/(2**i)),len=n)
}
df[[paste0("cluster",2*i+1)]]<- rep(seq(1:7),each=ceiling(n/7),len=n)
counts <- matrix(rpois(3200, lambda = 10), ncol=n, nrow=100)

for(cols in colnames(df)){
  df[[cols]]<- as.factor(df[[cols]])
}

treeviz <- preprocessAndCreateTreeViz(df, counts)

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
  expect_error(createTreeViz(df,counts),"Not a tree")
  
  hierarchydf <- colData(treeviz)@hierarchy_tree
  
  
  hierarchydf <- hierarchydf[, !colnames(hierarchydf) %in% c("samples", "otu_index")]
  
  flag=TRUE
  
  childs <- unique(hierarchydf[[ncol(hierarchydf)]])
  for (values in childs) {
    subsetted_list <- hierarchydf[hierarchydf[[ncol(hierarchydf)]] == values,]
    parent <- length(unique(subsetted_list[[ncol(hierarchydf) - 1]]))
    if (parent > 1) {
      flag=FALSE
    }
    expect_equal(flag, TRUE)
  }
  
  
  #expect_equal("ClusterAllClusters", unique(treeviz@colData@hierarchy_tree[[1]]))
})

test_that("check aggregate_tree", {
  sums<- list()
  sums[[1]]<- list(rowSums(counts[,c(1:10)]))
  sums[[2]]<- list(rowSums(counts[,c(11:15)]))
  sums[[3]]<- list(rowSums(counts[,c(16:25)]))
  sums[[4]]<- list(rowSums(counts[,c((26:32))]))
  sum_df<-as.data.frame(sums)
  aggtree<- TreeViz:::aggregateTree(treeviz, by="col", selectedLevel=3)
  
  agg_df<- as.data.frame(assays(aggtree)$counts)
  result<- all.equal(agg_df, sum_df, check.attributes= FALSE)
  expect_equal(result, TRUE)
})

# test_that("getHierarchy", {
#   res <- treeviz$getHierarchy(nodeId = NULL)
#   print(res)
#   expect_equal("cluster", as.character(res$tree$label))
#   expect_equal(2, res$tree$nchildren)
#   #Go to next level and make sure that the nchildren are correct, also feature names match
# })
# 
# test_that("getValues", {
#   sampleId<- "PM1:20080107"
#   res <- mObj$getValues(measurements = sampleId, start=0, end=10172, selectedLevels = 3)
#   
#   expected <- c(0.0, 0.0, 1812.977099, 0.0, 3.816794, 3.816794, 0.0, 1488.549618, 0.0, 3.816794, 72.519084, 187.022901, 11.450382, 76.335878, 0.0, 0.0, 0.0, 0.0, 68.702290)
#   
#   diff_result <- setdiff(round(unname(res[[sampleId]]), digits=3), round(expected,digits=3))
#   expect_equal(length(diff_result), 0)
# })
# 
# test_that("getRows", {
#   sampleId<- "PM1:20080107"
#   res <- mObj$getRows(measurement = sampleId, start=1, end=10172, selectedLevels = 3)
#   expected_label <- c("Actinomycetales","Coriobacteriales", "Bifidobacteriales","Lactobacillales",
#                       "Clostridiales","Erysipelotrichales","Rhizobiales","Campylobacterales",
#                       "Enterobacteriales","Pasteurellales","Not_Annotated_order_Bacteria")
#   
#   intersect_result <- intersect(res$metadata$label, expected_label)
#   expect_equal(length(intersect_result), length(expected_label))
# })
# 
# #getPCA
# test_that("getPCA", {
#   sampleIds<- c("PM9:20071217", "PM10:20080211")
#   res <- mObj$getPCA(measurement = sampleIds)
#   
#   expected_PM9_20071217 <- c(0.7675132, -0.6410331)
#   expect_equal(expected_PM9_20071217, unname(c(res$data[[1]]$PC1, res$data[[1]]$PC2)), tolerance = .0001)
#   
#   expected_PM10_20080211 <- c(0.6410331, 0.7675132)
#   expect_equal(expected_PM10_20080211, unname(c(res$data[[2]]$PC1, res$data[[2]]$PC2)), tolerance = .0001)
#   
#   expected_pca_variance_explained <- c(0.6640456, 0.4342690)
#   expect_equal(expected_pca_variance_explained, res$pca_variance_explained, tolerance = .0001)
# })
# 
# #getAlphaDiversity
# test_that("getAlphaDiversity", {
#   sampleIds<- c("PM9:20071217", "PM10:20080211")
#   res <- mObj$getAlphaDiversity(measurement = sampleIds)
#   
#   expected_PM9_20071217 <- 4.589864
#   expect_equal(expected_PM9_20071217, res$data[[1]]$alphaDiversity, tolerance = .0001)
#   
#   expected_PM10_20080211 <- 3.838823
#   expect_equal(expected_PM10_20080211, res$data[[2]]$alphaDiversity, tolerance = .0001)
# })
# 
# #searchTaxonomy
# test_that("searchTaxonomy", {
#   res <- mObj$searchTaxonomy(query = "bact", max_results = 10)
#   
#   expect_equal(10, length(res))
# })
