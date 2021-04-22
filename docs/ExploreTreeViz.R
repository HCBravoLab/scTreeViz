## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-packages, message=FALSE, warning=FALSE------------------------------
library(palmtree)
library(dplyr)
library(Seurat)
library(SC3)
library(scran)
library(scater)
library(clustree)
library(igraph)
library(scRNAseq)

## ---- eval=FALSE, results='hide', warning=FALSE, error=FALSE, message=FALSE----
#  # load dataset
#  sce<- LunSpikeInData('416b')
#  # Normalization
#  sce <- logNormCounts(sce)
#  # calculate umap and tsne
#  sce <- runUMAP(sce)
#  sce<- runTSNE(sce)
#  sce<- runPCA(sce)

## ---- eval=FALSE, warning=FALSE, error=FALSE, message=FALSE-------------------
#  treeViz <- createFromSCE(sce, reduced_dim = c("UMAP","PCA","TSNE"))
#  plot(treeViz)

## ---- eval=FALSE, warning=FALSE, error=FALSE, message=FALSE-------------------
#  # Forming clusters
#  for (i in  seq(10)) {
#    clust.kmeans <- kmeans(reducedDim(sce, "TSNE"), centers = i)
#    sce[[paste0("clust", i)]] <- factor(clust.kmeans$cluster)
#  }
#  
#  treeViz<- createFromSCE(sce, check_coldata = TRUE, col_regex = "clust")
#  plot(treeViz)

## ---- eval=FALSE, echo=TRUE, results='hide', warning=FALSE, error=FALSE, message=FALSE----
#  data(pbmc_small)
#  pbmc <- pbmc_small

## ---- eval=FALSE, echo=TRUE, results='hide', warning=FALSE, error=FALSE, message=FALSE----
#  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#  pbmc <- NormalizeData(pbmc)
#  all.genes <- rownames(pbmc)
#  pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
#  pbmc <- FindVariableFeatures(object = pbmc)
#  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#  pbmc <- FindNeighbors(pbmc, dims = 1:10)
#  pbmc <- FindClusters(pbmc, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), print.output = 0, save.SNN = TRUE)
#  pbmc

## ---- eval=FALSE, echo=TRUE, results='hide', warning=FALSE, error=FALSE, message=FALSE----
#  # pbmc<- RunTSNE(pbmc)
#  pbmc<- RunUMAP(pbmc, dims=1:3)
#  Reductions(pbmc)

## ---- eval=FALSE, echo=TRUE,  warning=FALSE, error=FALSE, message=FALSE-------
#  treeViz<- createFromSeurat(pbmc, check_metadata = TRUE, reduced_dim = c("umap","pca","tsne"))
#  plot(treeViz)

## ---- eval=FALSE, results='hide', warning=FALSE, error=FALSE, message=FALSE----
#  n=64
#  # create a hierarchy
#  df<- data.frame(cluster0=rep(1,n))
#  for(i in seq(1,5)){
#    df[[paste0("cluster",i)]]<- rep(seq(1:(2**i)),each=ceiling(n/(2**i)),len=n)
#  }
#  
#  # generate a count matrix
#  counts <- matrix(rpois(6400, lambda = 10), ncol=n, nrow=100)
#  
#  # create a `TreeViz` object
#  treeViz <- createTreeViz(df, counts)
#  plot(treeViz)

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  app <- startTreeviz(treeViz, top_genes = 500)

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  app$stop_app()

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  setTreevizStandalone()
#  app <- startTreevizStandalone(treeViz)

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  app$stop_app()

