# Tests for max_error_summary (tie handling).

test_that("max_error_summary tie_method = 'first' returns the smallest tied index", {
  s <- max_error_summary(c(0.1, 0.5, 0.5, 0.3), tie_method = "first")
  expect_equal(s$max_error_value, 0.5)
  expect_equal(s$argmax_index, 2L)
})

test_that("max_error_summary tie_method = 'last' returns the largest tied index", {
  s <- max_error_summary(c(0.1, 0.5, 0.5, 0.3), tie_method = "last")
  expect_equal(s$max_error_value, 0.5)
  expect_equal(s$argmax_index, 3L)
})

test_that("max_error_summary tie_method = 'random' samples among tied indices only", {
  set.seed(99L)
  idx <- replicate(200L, max_error_summary(c(0.1, 0.5, 0.5, 0.3), tie_method = "random")$argmax_index)
  expect_setequal(unique(idx), c(2L, 3L))
})

test_that("max_error_summary handles a unique maximum", {
  s <- max_error_summary(c(0.1, 0.9, 0.4), tie_method = "first")
  expect_equal(s$argmax_index, 2L)
  expect_equal(s$max_error_value, 0.9)
})
