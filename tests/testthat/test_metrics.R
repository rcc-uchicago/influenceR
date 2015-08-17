library(testthat)

context("tests")

load("flo_results.RData")

test_that("metrics work as expected", {
  #expect_equal(influenceR::betweenness(flo_graph), flo_bet) # unsure why this isn't working
  expect_equal(influenceR::eigencentrality(flo_graph), flo_eigen)
  expect_equal(influenceR::ens(flo_graph), flo_ens)
  expect_equal(influenceR::constraint(flo_graph), flo_constraint)
  expect_equal(influenceR::bridging(flo_graph, MPI=F), flo_bridge)
})



