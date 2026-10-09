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

test_that("simulate_counts_dirichlet_multinomial returns matrices of the right shape and totals", {
  p <- c(0.5, 0.3, 0.15, 0.05)
  set.seed(1L)
  out <- simulate_counts_dirichlet_multinomial(p, n_people = 7L, n_per_person = 50L, concentration = 20)
  expect_equal(dim(out$counts), c(7L, 4L))
  expect_equal(dim(out$person_true_proportions), c(7L, 4L))
  expect_true(is.integer(out$counts))
  expect_equal(unname(rowSums(out$counts)), rep(50L, 7L))
  expect_equal(unname(rowSums(out$person_true_proportions)), rep(1, 7L))
  expect_true(all(out$person_true_proportions > 0))
  expect_equal(colnames(out$counts), paste0("cell_type_", 1:4))
  expect_equal(colnames(out$person_true_proportions), colnames(out$counts))
})

test_that("simulate_counts_dirichlet_multinomial is reproducible under the same seed", {
  p <- c(0.5, 0.3, 0.15, 0.05)
  set.seed(42L)
  first <- simulate_counts_dirichlet_multinomial(p, n_people = 20L, n_per_person = 30L, concentration = 5)
  set.seed(42L)
  second <- simulate_counts_dirichlet_multinomial(p, n_people = 20L, n_per_person = 30L, concentration = 5)
  expect_identical(first, second)
})

test_that("simulate_counts_dirichlet_multinomial person proportions average to p", {
  p <- c(0.5, 0.3, 0.15, 0.05)
  set.seed(123L)
  out <- simulate_counts_dirichlet_multinomial(p, n_people = 5000L, n_per_person = 10L, concentration = 10)
  expect_equal(unname(colMeans(out$person_true_proportions)), p, tolerance = 0.01, scale = 1)
})

test_that("simulate_counts_dirichlet_multinomial errors when the gamma draws have an invalid total", {
  # A vanishingly small concentration makes every gamma draw underflow to zero.
  p <- c(0.5, 0.5)
  set.seed(1L)
  expect_error(
    simulate_counts_dirichlet_multinomial(p, n_people = 3L, n_per_person = 10L, concentration = 1e-300),
    "invalid total gamma draw"
  )
})

test_that("compute_errors returns AE and ARE for a known example", {
  p <- c(0.5, 0.3, 0.2)
  phat <- c(0.6, 0.25, 0.15)
  err <- compute_errors(phat, p, metrics = c("AE", "ARE"))
  expect_equal(err$AE, abs(phat - p))
  expect_equal(err$ARE, abs(phat - p) / p)
})
