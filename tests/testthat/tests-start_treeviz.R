context("create TreeVizApp class")
library(Seurat)

test_that("startTreeViz creates a TreeVizApp Object", {
  #skip("need data here")
  pbmc_small
  treeviz<-createFromSeurat(pbmc_small)
  app <-
    startTreeviz(treeviz,
                 try_ports=TRUE,
                 non_interactive=TRUE)
  expect_is(app, "TreeVizApp")
  
  expect_is(app$server, "EpivizServer")
  expect_is(app$chart_mgr, "EpivizChartMgr")
  expect_is(app$data_mgr, "EpivizDataMgr")
  
  expect_true(app$server$is_closed())
})