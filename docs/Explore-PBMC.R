## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load-packages, message=FALSE, warning=FALSE------------------------------
library(TreeViz)

## -----------------------------------------------------------------------------
data(PBMC3k_TreeViz)

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  library(Seurat)
#  
#  pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#  pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
#  
#  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#  
#  pbmc <- NormalizeData(pbmc)
#  
#  all.genes <- rownames(pbmc)
#  pbmc <- ScaleData(pbmc, features = all.genes)
#  
#  pbmc <- FindVariableFeatures(pbmc)
#  pbmc <- FindNeighbors(pbmc, dims=1:10)
#  
#  pbmc <- FindClusters(pbmc, resolution=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), save.SNN=TRUE)
#  
#  pbmc <- RunTSNE(pbmc)
#  pbmc <- RunUMAP(pbmc, dims=1:5)
#  pbmc <- RunPCA(pbmc)
#  
#  PBMC3K_TreeViz <- createFromSeurat(pbmc, reduced_dim = c("pca", "tsne", "umap"))
#  PBMC3K_TreeViz

## ---- eval=FALSE--------------------------------------------------------------
#  app <- startTreeviz(PBMC3k_TreeViz, top_genes = 100)

## ---- eval=FALSE--------------------------------------------------------------
#  app$plotGene(gene="TYROBP")

## ---- eval=FALSE--------------------------------------------------------------
#  app$stop_app()

