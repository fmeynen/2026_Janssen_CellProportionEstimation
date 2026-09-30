# Tests for evaluate_thresholds.

toy_max_errors <- function() {
  matrix(
    c(0.05, 0.15, 0.25, 0.10, 0.20, 0.30),
    nrow = 3L, ncol = 2L,
    dimnames = list(NULL, c("AE", "ARE"))
  )
}

rate_at <- function(curves, metric, tau) {
  curves$success_rate[curves$metric == metric & curves$tau == tau]
}

test_that("evaluate_thresholds computes success rates per metric for a shared tau grid", {
  curves <- evaluate_thresholds(toy_max_errors(), taus = c(0.10, 0.20))
  expect_true(all(c("metric", "tau", "success_rate") %in% names(curves)))
  expect_equal(rate_at(curves, "AE", 0.10), 1 / 3)
  expect_equal(rate_at(curves, "AE", 0.20), 2 / 3)
  expect_equal(rate_at(curves, "ARE", 0.10), 1 / 3)
  expect_equal(rate_at(curves, "ARE", 0.20), 2 / 3)
})

test_that("evaluate_thresholds supports a named list of per-metric tau grids", {
  curves <- evaluate_thresholds(
    toy_max_errors(),
    taus = list(AE = c(0.10, 0.20), ARE = c(0.15, 0.25))
  )
  expect_equal(rate_at(curves, "AE", 0.10), 1 / 3)
  expect_equal(rate_at(curves, "AE", 0.20), 2 / 3)
  expect_equal(rate_at(curves, "ARE", 0.15), 1 / 3)
  expect_equal(rate_at(curves, "ARE", 0.25), 2 / 3)
})

test_that("evaluate_thresholds returns one row per metric/tau pair with unequal grid lengths", {
  curves <- evaluate_thresholds(
    toy_max_errors(),
    taus = list(AE = c(0.10, 0.20), ARE = c(0.15, 0.25, 0.35))
  )
  expect_equal(sum(curves$metric == "AE"), 2L)
  expect_equal(sum(curves$metric == "ARE"), 3L)
  expect_equal(nrow(curves), 5L)
})

test_that("evaluate_thresholds errors naming the metric missing from a tau list", {
  expect_error(evaluate_thresholds(toy_max_errors(), taus = list(AE = 0.10)), "ARE")
})
