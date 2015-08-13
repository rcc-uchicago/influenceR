library(testthat)
library(igraph)
library(influenceR)

context("tests")

load("flo_results.RData")

test_that("metrics work as expected", {
  expect_equal(betweenness(flo_graph), flo_bet)
  expect_equal(eigencentrality(flo_graph), flo_eigen)
  expect_equal(ens(flo_graph), flo_ens)
  expect_equal(constraint(flo_graph), flo_constraint)
  expect_equal(bridging(flo_graph, MPI=F), flo_bridge)
})


ens_test <- function(g) {
  ens <- vector("numeric", length(V(g)))
  for (i in V(g)) {
    s <- 0
    ni <- neighbors(g, i)
    di <- degree(g, i)
    for (j in ni) {
      t <- 0
      Q <- intersection(ni, neighbors(g, j))
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
