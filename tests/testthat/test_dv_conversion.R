library(intermahpr)
context("Data Verification: Input data is converted to tibbles for internal use")

test_that("data frames are converted into tibbles", {
  RR <- data.frame(testrr)
  expect_is(format_rr(RR), 'tbl')
  PC <- data.frame(x=0)
  expect_is(format_pc(PC), 'tbl')
})

test_that("tibbles stay as tibbles", {
  RR <- as.tibble(testrr)
  expect_is(format_rr(RR), 'tbl')
  PC <- data_frame(x=0)
  expect_is(format_pc(PC), 'tbl')
})

test_that("trigonometric functions match identities", {
  expect_equal(sin(pi / 4), 1 / sqrt(2))
  expect_equal(cos(pi / 4), 1 / sqrt(2))
  expect_equal(tan(pi / 4), 1)
})
