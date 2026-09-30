# Tests for proportion generators (generate_proportions_beta, generate_props_fixed_*_beta,
# and the generate_proportions dispatcher).

test_that("generate_proportions_beta returns K strictly positive values summing to 1", {
  p <- generate_proportions_beta(alpha = 2, K = 10)
  expect_length(p, 10L)
  expect_equal(sum(p), 1, tolerance = 1e-12)
  expect_true(all(p > 0))
})

test_that("generate_proportions_beta is deterministic given alpha and differs across alpha", {
  p <- generate_proportions_beta(alpha = 2, K = 10)
  expect_identical(p, generate_proportions_beta(alpha = 2, K = 10))
  expect_false(identical(p, generate_proportions_beta(alpha = 5, K = 10)))
})

test_that("the beta dispatcher route matches generate_proportions_beta", {
  expect_identical(
    generate_proportions(alpha = 2, K = 10, method = "beta"),
    generate_proportions_beta(alpha = 2, K = 10)
  )
})

test_that("fixed_max_beta requires p_max", {
  expect_error(
    generate_proportions(alpha = 2, K = 10, method = "fixed_max_beta"),
    "p_max must be provided"
  )
})

test_that("fixed_max_beta puts p_max at the highest index as a strict unique maximum", {
  K <- 10L
  p <- generate_proportions(alpha = 2, K = K, method = "fixed_max_beta", p_max = 0.4)
  expect_length(p, K)
  expect_equal(sum(p), 1, tolerance = 1e-12)
  expect_equal(p[K], 0.4, tolerance = 1e-12)
  expect_equal(which.max(p), K)
  expect_equal(sum(p == max(p)), 1L)
})

test_that("fixed_max_beta returns one valid row per p_max for vector p_max", {
  K <- 10L
  pm <- c(0.3, 0.4)
  m <- generate_props_fixed_max_beta(alpha = 2, K = K, p_max = pm)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, K))
  expect_equal(unname(rowSums(m)), c(1, 1), tolerance = 1e-12)
  expect_equal(unname(m[, K]), pm, tolerance = 1e-12)
  expect_true(all(apply(m, 1, function(r) which.max(r) == K && sum(r == max(r)) == 1L)))
})

test_that("fixed_max_beta warns then errors on impossible alpha/K/p_max combinations", {
  expect_warning(
    expect_error(
      generate_proportions(alpha = 2, K = 2, method = "fixed_max_beta", p_max = 0.4),
      "not strictly unique"
    ),
    "Impossible fixed_max_beta combination"
  )
})

test_that("fixed_min_beta requires p_min", {
  expect_error(
    generate_proportions(alpha = 2, K = 10, method = "fixed_min_beta"),
    "p_min must be provided"
  )
})

test_that("the fixed_min_beta dispatcher route matches generate_props_fixed_min_beta", {
  expect_identical(
    generate_proportions(alpha = 2, K = 10, method = "fixed_min_beta", p_min = 0.01),
    generate_props_fixed_min_beta(alpha = 2, K = 10, p_min = 0.01)
  )
})

test_that("fixed_min_beta puts p_min at the lowest index as the minimum", {
  p <- generate_proportions(alpha = 2, K = 10, method = "fixed_min_beta", p_min = 0.01)
  expect_length(p, 10L)
  expect_equal(sum(p), 1, tolerance = 1e-12)
  expect_equal(p[1L], 0.01, tolerance = 1e-12)
  expect_true(all(p[-1L] >= p[1L]))
})

test_that("fixed_min_beta allows ties with p_min (alpha = 1, p_min = 1/K)", {
  p <- generate_proportions(alpha = 1, K = 10, method = "fixed_min_beta", p_min = 1 / 10)
  expect_equal(p, rep(1 / 10, 10), tolerance = 1e-12)
})

test_that("fixed_min_beta returns one valid row per p_min for vector p_min", {
  pm <- c(0.005, 0.01)
  m <- generate_props_fixed_min_beta(alpha = 2, K = 10, p_min = pm)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, 10L))
  expect_equal(unname(rowSums(m)), c(1, 1), tolerance = 1e-12)
  expect_equal(unname(m[, 1L]), pm, tolerance = 1e-12)
  expect_true(all(apply(m, 1, function(r) all(r[-1L] >= r[1L]))))
})

test_that("fixed_min_beta warns then fails with a classed error on impossible combinations", {
  expect_warning(
    expect_error(
      generate_proportions(alpha = 2, K = 10, method = "fixed_min_beta", p_min = 0.3),
      class = "impossible_fixed_min_error"
    ),
    "Impossible fixed_min_beta combination"
  )
})

test_that("validate_proportions rejects zero proportions", {
  expect_error(validate_proportions(c(0.5, 0.5, 0.0)))
})
