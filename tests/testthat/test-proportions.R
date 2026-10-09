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

test_that("fixed_min_beta requires p_min", {
  expect_error(
    generate_proportions(alpha = 2, K = 10, method = "fixed_min_beta"),
    "p_min must be provided"
  )
})

test_that("fixed_max_beta attains p_max exactly, never exceeds it and sums to 1", {
  for (alpha in c(1, 2, 3, 5)) {
    for (p_max in c(0.12, 0.2, 0.3, 0.5)) {
      p <- generate_props_fixed_max_beta(alpha = alpha, K = 10, p_max = p_max)
      expect_length(p, 10L)
      expect_equal(sum(p), 1, tolerance = 1e-12)
      expect_true(all(p > 0))
      expect_identical(max(p), p_max)
      expect_true(all(p <= p_max))
    }
  }
})

test_that("fixed_min_beta attains p_min exactly, never goes below it and sums to 1", {
  for (alpha in c(1, 2, 3, 5)) {
    for (p_min in c(0.005, 0.01, 0.05, 0.09)) {
      p <- generate_props_fixed_min_beta(alpha = alpha, K = 10, p_min = p_min)
      expect_length(p, 10L)
      expect_equal(sum(p), 1, tolerance = 1e-12)
      expect_true(all(p > 0))
      expect_identical(min(p), p_min)
      expect_true(all(p >= p_min))
    }
  }
})

test_that("fixed_min_beta clips crossing components to p_min and allows ties", {
  p <- generate_props_fixed_min_beta(alpha = 5, K = 10, p_min = 0.01)
  expect_identical(sum(p == 0.01), 4L)
  expect_true(all(p[5:10] > 0.01))
  expect_equal(p[5], 0.020, tolerance = 5e-3)

  expect_identical(sum(generate_props_fixed_min_beta(alpha = 3, K = 10, p_min = 0.01) == 0.01), 2L)
  expect_identical(sum(generate_props_fixed_min_beta(alpha = 4, K = 10, p_min = 0.01) == 0.01), 3L)
})

test_that("fixed_max_beta clips crossing components to p_max and allows ties", {
  p <- generate_props_fixed_max_beta(alpha = 5, K = 10, p_max = 0.2)
  expect_gt(sum(p == 0.2), 1L)
  expect_true(all(p <= 0.2))
  expect_equal(sum(p), 1, tolerance = 1e-12)
})

test_that("fixed_max_beta without crossings moves only the maximum and scales the rest proportionally", {
  beta <- generate_proportions_beta(alpha = 2, K = 10)
  p <- generate_props_fixed_max_beta(alpha = 2, K = 10, p_max = 0.3)
  i <- which.max(beta)
  expect_identical(p[i], 0.3)
  expect_equal(p[-i], beta[-i] * 0.7 / (1 - max(beta)), tolerance = 1e-12)
  expect_identical(sum(p == 0.3), 1L)
})

test_that("fixed_min_beta without crossings moves only the minimum and scales the rest proportionally", {
  beta <- generate_proportions_beta(alpha = 2, K = 10)
  p_min <- min(beta) / 2
  p <- generate_props_fixed_min_beta(alpha = 2, K = 10, p_min = p_min)
  i <- which.min(beta)
  expect_identical(p[i], p_min)
  expect_equal(p[-i], beta[-i] * (1 - p_min) / (1 - min(beta)), tolerance = 1e-12)
})

test_that("fixed_min_beta with p_min equal to the Beta minimum reproduces the Beta vector", {
  beta <- generate_proportions_beta(alpha = 2, K = 10)
  expect_equal(generate_props_fixed_min_beta(alpha = 2, K = 10, p_min = 0.01), beta, tolerance = 1e-12)
})

test_that("p_min = 1/K and p_max = 1/K give the uniform vector", {
  expect_equal(generate_props_fixed_min_beta(alpha = 3, K = 10, p_min = 1 / 10), rep(1 / 10, 10), tolerance = 1e-12)
  expect_equal(generate_props_fixed_max_beta(alpha = 3, K = 10, p_max = 1 / 10), rep(1 / 10, 10), tolerance = 1e-12)
})

test_that("fixed_max_beta returns one valid row per p_max for vector p_max", {
  pm <- c(0.3, 0.4)
  m <- generate_props_fixed_max_beta(alpha = 2, K = 10, p_max = pm)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, 10L))
  expect_equal(colnames(m), paste0("cell_type_", 1:10))
  expect_equal(rownames(m), c("p_max_1_0.3", "p_max_2_0.4"))
  expect_equal(unname(rowSums(m)), c(1, 1), tolerance = 1e-12)
  expect_identical(unname(apply(m, 1, max)), pm)
})

test_that("fixed_min_beta returns one valid row per p_min for vector p_min", {
  pm <- c(0.005, 0.01)
  m <- generate_props_fixed_min_beta(alpha = 2, K = 10, p_min = pm)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, 10L))
  expect_equal(colnames(m), paste0("cell_type_", 1:10))
  expect_equal(rownames(m), c("p_min_1_0.005", "p_min_2_0.010"))
  expect_equal(unname(rowSums(m)), c(1, 1), tolerance = 1e-12)
  expect_identical(unname(apply(m, 1, min)), pm)
})

test_that("fixed_max_beta warns then fails with a classed error when K * p_max < 1", {
  expect_warning(
    expect_error(
      generate_proportions(alpha = 2, K = 10, method = "fixed_max_beta", p_max = 0.05),
      class = "impossible_fixed_max_error"
    ),
    "Impossible fixed_max_beta combination"
  )
  expect_warning(
    expect_error(
      generate_proportions(alpha = 2, K = 2, method = "fixed_max_beta", p_max = 0.4),
      "K * p_max < 1",
      fixed = TRUE
    ),
    "Impossible fixed_max_beta combination"
  )
})

test_that("fixed_min_beta warns then fails with a classed error when K * p_min > 1", {
  expect_warning(
    expect_error(
      generate_proportions(alpha = 2, K = 10, method = "fixed_min_beta", p_min = 0.3),
      class = "impossible_fixed_min_error"
    ),
    "Impossible fixed_min_beta combination"
  )
  expect_warning(
    expect_error(
      generate_proportions(alpha = 2, K = 10, method = "fixed_min_beta", p_min = 0.11),
      "K * p_min > 1",
      fixed = TRUE
    ),
    "Impossible fixed_min_beta combination"
  )
})

test_that("generate_proportions routes to the fixed generators identically to direct calls", {
  expect_identical(
    generate_proportions(alpha = 2, K = 10, method = "fixed_max_beta", p_max = 0.4),
    generate_props_fixed_max_beta(alpha = 2, K = 10, p_max = 0.4)
  )
  expect_identical(
    generate_proportions(alpha = 2, K = 10, method = "fixed_min_beta", p_min = 0.01),
    generate_props_fixed_min_beta(alpha = 2, K = 10, p_min = 0.01)
  )
})

test_that("generate_proportions forwards a custom length-K grid to the fixed methods", {
  grid <- seq(0.2, 0.8, length.out = 10L)
  p_max <- generate_proportions(alpha = 2, K = 10, method = "fixed_max_beta", p_max = 0.4, grid = grid)
  expect_identical(p_max, generate_props_fixed_max_beta(alpha = 2, K = 10, p_max = 0.4, grid = grid))
  expect_false(isTRUE(all.equal(p_max, generate_props_fixed_max_beta(alpha = 2, K = 10, p_max = 0.4))))

  p_min <- generate_proportions(alpha = 2, K = 10, method = "fixed_min_beta", p_min = 0.01, grid = grid)
  expect_identical(p_min, generate_props_fixed_min_beta(alpha = 2, K = 10, p_min = 0.01, grid = grid))
  expect_false(isTRUE(all.equal(p_min, generate_props_fixed_min_beta(alpha = 2, K = 10, p_min = 0.01))))
})

test_that("the fixed generators reject a grid whose length is not K", {
  expect_error(
    generate_proportions(alpha = 2, K = 10, method = "fixed_max_beta", p_max = 0.4, grid = default_beta_grid(9)),
    "grid must have length K"
  )
  expect_error(
    generate_proportions(alpha = 2, K = 10, method = "fixed_min_beta", p_min = 0.01, grid = default_beta_grid(9)),
    "grid must have length K"
  )
})

test_that("validate_proportions rejects zero proportions", {
  expect_error(validate_proportions(c(0.5, 0.5, 0.0)))
})

test_that("the dispatcher errors on a bound the method does not use, for every method and unused bound", {
  expect_error(generate_proportions(alpha = 2, K = 10, method = "beta", p_min = 0.01), "p_min is not used by method = 'beta'")
  expect_error(generate_proportions(alpha = 2, K = 10, method = "beta", p_max = 0.4), "p_max is not used by method = 'beta'")
  expect_error(
    generate_proportions(alpha = 2, K = 10, method = "fixed_max_beta", p_max = 0.4, p_min = 0.01),
    "p_min is not used by method = 'fixed_max_beta'; leave it NULL."
  )
  expect_error(
    generate_proportions(alpha = 2, K = 10, method = "fixed_min_beta", p_min = 0.01, p_max = 0.4),
    "p_max is not used by method = 'fixed_min_beta'; leave it NULL."
  )
})

test_that("the dispatcher still accepts the bound each method uses", {
  expect_length(generate_proportions(alpha = 2, K = 10, method = "beta"), 10L)
  expect_identical(max(generate_proportions(alpha = 2, K = 10, method = "fixed_max_beta", p_max = 0.4)), 0.4)
  expect_identical(min(generate_proportions(alpha = 2, K = 10, method = "fixed_min_beta", p_min = 0.01)), 0.01)
})
