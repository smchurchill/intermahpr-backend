library(intermahpr)
context("Data Verification: Input data is converted to tibbles for internal use")

test_that("data frames are converted into tibbles", {
  RR <- data.frame(testrr)
  expect_is(format_v0_rr(RR), 'tbl')
  PC <- data.frame(testpc)
  expect_is(format_v0_pc(PC), 'tbl')
})

test_that("tibbles stay as tibbles", {
  RR <- tibble::as.tibble(testrr)
  expect_is(format_v0_rr(RR), 'tbl')
  PC <- tibble::as.tibble(testpc)
  expect_is(format_v0_pc(PC), 'tbl')
})

test_that("trigonometric functions match identities", {
  expect_equal(sin(pi / 4), 1 / sqrt(2))
  expect_equal(cos(pi / 4), 1 / sqrt(2))
  expect_equal(tan(pi / 4), 1)
})
