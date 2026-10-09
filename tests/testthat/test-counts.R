# Tests for count simulators and compute_errors.

test_that("simulate_counts_multinomial returns K nonnegative integers summing to n", {
  p <- generate_proportions_beta(alpha = 2, K = 10)
  set.seed(1L)
  y <- simulate_counts_multinomial(p, n = 200L)
  expect_length(y, 10L)
  expect_equal(sum(y), 200L)
  expect_true(is.integer(y))
  expect_true(all(y >= 0L))
})

test_that("simulate_counts routes the multinomial model to a valid count vector", {
  p <- generate_proportions_beta(alpha = 2, K = 10)
  y <- simulate_counts(p, n = 200L, model = "multinomial")
  expect_length(y, 10L)
  expect_equal(sum(y), 200L)
})

test_that("simulate_counts errors when n is missing for the multinomial model", {
  p <- generate_proportions_beta(alpha = 2, K = 10)
  expect_error(simulate_counts(p, model = "multinomial"), "n must be provided")
})

test_that("compute_errors returns AE and ARE for a known example", {
  p <- c(0.5, 0.3, 0.2)
  phat <- c(0.6, 0.25, 0.15)
  err <- compute_errors(phat, p, metrics = c("AE", "ARE"))
  expect_equal(err$AE, abs(phat - p))
  expect_equal(err$ARE, abs(phat - p) / p)
})
