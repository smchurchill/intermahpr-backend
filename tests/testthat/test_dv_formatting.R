library(intermahpr)
context("Data Verification: Input data must be held to a certain standard of formatting")

test_that("errors are thrown when variables are missing", {
  for(i in c(3, 5:12)) {
      expect_error(format_v0_pc(testpc[-i]))
  }
})

test_that("trigonometric functions match identities", {
  expect_equal(sin(pi / 4), 1 / sqrt(2))
  expect_equal(cos(pi / 4), 1 / sqrt(2))
  expect_equal(tan(pi / 4), 1)
})
