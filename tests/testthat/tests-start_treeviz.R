context("create TreeVizApp class")
library(Seurat)
library(scater)

test_that("startTreeViz creates a TreeVizApp Object from Seurat", {
  skip("need data here")
  # pbmc_small
  # treeviz<-createFromSeurat(pbmc_small)
  # app <-
  #   startTreeviz(treeviz,
  #                non_interactive=TRUE)
  # expect_is(app, "TreeVizApp")
  # 
  # expect_is(app$server, "EpivizServer")
  # expect_is(app$chart_mgr, "EpivizChartMgr")
  # expect_is(app$data_mgr, "EpivizDataMgr")
  # 
  # expect_true(app$server$is_closed())
})

test_that("startTreeViz creates a TreeVizApp Object from Single Cell Experiment", {
  skip("need data here")
  # set.seed(1)
  # sce <- mockSCE()
  # treeviz<-createFromSCE(sce) 
  # app <-
  #   startTreeviz(treeviz,
  #                non_interactive=TRUE)
  # expect_is(app, "TreeVizApp")
  # 
  # expect_is(app$server, "EpivizServer")
  # expect_is(app$chart_mgr, "EpivizChartMgr")
  # expect_is(app$data_mgr, "EpivizDataMgr")
  # 
  # expect_true(app$server$is_closed())
})


