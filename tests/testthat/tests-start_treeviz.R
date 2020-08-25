context("create TreVvizApp class")

test_that("startTreeViz creates a TreeVizApp Object", {
  app <- startTreeviz(non_interactive=TRUE)
  expect_is(app, "TreeVizApp")
  
  expect_is(app$server, "EpivizServer")
  expect_is(app$chart_mgr, "EpivizChartMgr")
  expect_is(app$data_mgr, "EpivizDataMgr")
  
  expect_true(app$server$is_closed())
})