library(intermahpr)
context("Test methods within interface.R file")

test_that("Vector integral fails with different size vector inputs", {
  expect_error(vintegrate(rep(0,2), rep(0,2), rep(0,3)))
  expect_error(vintegrate(rep(0,2), rep(0,3), rep(0,3)))
  expect_error(vintegrate(rep(0,3), rep(0,2), rep(0,3)))
})

test_that("Vector integrate integrates correctly", {
  expect_equal(
    vintegrate(
      list(function(x) 2*x, function(x) 3*(x^2), function(x) 4*(x^3)),
      c(   0,               0,                   0),
      c(   2,               2,                   2)),
    c(     4,               8,                   16)
  )
})

test_that("trigonometric functions match identities", {
  expect_equal(sin(pi / 4), 1 / sqrt(2))
  expect_equal(cos(pi / 4), 1 / sqrt(2))
  expect_equal(tan(pi / 4), 1)
})
