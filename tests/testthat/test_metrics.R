
library(testthat)
library(ergm)
library(igraph)
library(influenceR)
data(flo)

context("tests")

fe <- which(flo>0, arr.ind=T)
fe_undir <- fe[fe[,1] > fe[,2],]
fg <- igraph::graph(t(fe_undir), directed=F)

load("flo_results.RData")

test_that("metrics work as expected", {
  expect_equal(betweenness(fg), flo_bet*2)
  expect_equal(eigencentrality(fg), flo_eigen)
  expect_equal(ens(fg), flo_ens)
  expect_equal(constraint(fg), flo_constraint)
  expect_equal(bridging(fg, MPI=F), flo_bridge)
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

ens_test2 <- function(g) {
  ens <- vector("numeric", length(V(g)))
  for (i in V(g)) {
    s <- 0
    ni <- neighbors(g, i)
    di <- degree(g, i)
    for (j in ni) {
      Q <- intersection(ni, neighbors(g, j))
      s <- s + 1 - length(Q)/di
    }
    ens[i] <- s
  }
  ens
}

test_that("ens matches simpler function", {
  expect_equal(ens(fg), ens_test(fg))
})
