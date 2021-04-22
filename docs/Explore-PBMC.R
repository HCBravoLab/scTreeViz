## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load-packages, message=FALSE, warning=FALSE------------------------------
library(palmtree)
library(scran)

## -----------------------------------------------------------------------------
data(PBMC3k_TreeViz)

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  library(Seurat)
#  
#  pbmc.data <- Read10X(data.dir = "../../filtered_gene_bc_matrices/hg19/")
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
#  pbmc <- RunPCA(pbmc)
#  
#  pbmc <- FindNeighbors(pbmc, dims=1:10)
#  
#  pbmc <- FindClusters(pbmc, resolution=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), save.SNN=TRUE)
#  
#  pbmc <- RunTSNE(pbmc)
#  pbmc <- RunUMAP(pbmc, dims=1:5)
#  
#  PBMC3k_TreeViz <- createFromSeurat(pbmc, reduced_dim = c("pca", "tsne", "umap"), check_metadata = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  # to visualize top 100 variable genes
#  app <- startTreeviz(PBMC3k_TreeViz, top_genes = 100, host="http://localhost:8877")
#  
#  # to visualize a genes list
#  # app <- startTreeviz(PBMC3k_TreeViz, genes=c("GBA", "FDPS", "TREM2"))

## ---- eval=FALSE--------------------------------------------------------------
#  sce <- app$extract_SCE()

## ---- eval=FALSE--------------------------------------------------------------
#  library(scater)
#  sce <- logNormCounts(sce)
#  
#  # ignore removed cells
#  sce_rm <- sce[, sce$treeviz_clusters != 'removed']
#  
#  markers <- findMarkers(sce_rm, groups=sce_rm$treeviz_clusters)

## ---- eval=FALSE--------------------------------------------------------------
#  markers_pr <- lapply(markers, function(x) {
#    rownames(head(x[x$p.value < 0.05,]))
#  })
#  markers_pr

## ---- eval=FALSE--------------------------------------------------------------
#  app$plotGene(gene="AIF1")

## ---- eval=FALSE--------------------------------------------------------------
#  app$stop_app()

