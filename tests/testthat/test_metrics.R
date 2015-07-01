
library(testthat)
library(ergm)
library(influenceR)
data(flo)

context("tests")

fe <- which(flo>0, arr.ind=T)
fg <- graph(fe, directed=F)

load("flo_results.RData")

test_that("metrics work as expected", {
  expect_equal(betweenness(fg), flo_bet)
  expect_equal(eigencentrality(fg), flo_eigen)
  expect_equal(ens(fg), flo_ens)
  expect_equal(constraint(fg), flo_constraint)
})
