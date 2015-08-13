library(testthat)

context("tests")

load("flo_results.RData")

test_that("metrics work as expected", {
  expect_equal(influenceR::betweenness(flo_graph), flo_bet)
  expect_equal(influenceR::eigencentrality(flo_graph), flo_eigen)
  expect_equal(influenceR::ens(flo_graph), flo_ens)
  expect_equal(influenceR::constraint(flo_graph), flo_constraint)
  expect_equal(influenceR::bridging(flo_graph, MPI=F), flo_bridge)
})


ens_test <- function(g) {
  ens <- vector("numeric", length(igraph::V(g)))
  for (i in igraph::V(g)) {
    s <- 0
    ni <- igraph::neighbors(g, i)
    di <- igraph::degree(g, i)
    for (j in ni) {
      t <- 0
      Q <- igraph::intersection(ni, igraph::neighbors(g, j))
      for (q in Q) {
        t <- t + 1/di
      }
      s <- s + 1 - t
    }
    ens[i] <- s
  }
  ens
}


test_that("ens matches simpler function", {
  expect_equal(ens(flo_graph), ens_test(flo_graph))
})
