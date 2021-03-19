context("create TreeVizApp class")

test_that("startTreeViz creates a TreeVizApp Object", {
  skip("need data here")
  app <- startTreeviz(non_interactive=TRUE)
  expect_is(app, "TreeVizApp")
  
  expect_is(app$server, "EpivizServer")
  expect_is(app$chart_mgr, "EpivizChartMgr")
  expect_is(app$data_mgr, "EpivizDataMgr")
  
  expect_true(app$server$is_closed())
})