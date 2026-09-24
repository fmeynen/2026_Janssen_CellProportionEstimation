# Tests for simulate_success_at_n() -----------------------------------------------------------------------------
# Both models must route through the shared replicate_success() success rule (extraction.R); the old
# required_person_fraction / per-person rule has been removed entirely.

multinomial_config <- list(
  K                 = 3L,
  B                 = 5L,
  taus              = list(AE = 0.3, ARE = 2),
  metrics           = c("AE", "ARE"),
  model             = "multinomial",
  tie_method        = "random",
  proportion_method = "beta"
)

dm_config <- list(
  K                 = 3L,
  B                 = 5L,
  taus              = list(AE = 0.3, ARE = 2),
  metrics           = c("AE", "ARE"),
  model             = "dirichlet_multinomial",
  n_people          = 4L,
  concentration     = 5,
  proportion_method = "beta"
)

test_that("multinomial success matches a direct recomputation from rep_out$max_errors", {
  res <- simulate_success_at_n(alpha = 1, n = 10L, config = multinomial_config, seed = 42)

  max_errors <- res$rep_out$max_errors
  metrics    <- names(multinomial_config$taus)
  expected   <- rep(TRUE, nrow(max_errors))
  for (m in metrics) {
    expected <- expected & (max_errors[, m] <= multinomial_config$taus[[m]])
  }

  # The multinomial model has no person structure: simulate_success_at_n() treats each replicate as a single
  # synthetic person, so averaging over persons in replicate_success() is a no-op and this must reduce exactly to
  # "every cell-type error <= tau", i.e. the max-error comparison above.
  expect_identical(res$success, expected)
  expect_identical(res$success_count, sum(expected))
  expect_equal(res$success_rate, mean(expected))
})

test_that("dirichlet_multinomial success matches replicate_success() on rep_out$person_results", {
  res <- simulate_success_at_n(alpha = 1, n = 10L, config = dm_config, seed = 42)

  expected <- replicate_success(res$rep_out$person_results, dm_config$taus)$pass

  expect_identical(res$success, as.logical(expected))
  expect_identical(res$success_count, sum(expected))
  expect_equal(res$success_rate, mean(expected))
})

test_that("results are reproducible for the same seed and change for a different seed", {
  res1 <- simulate_success_at_n(alpha = 1, n = 10L, config = multinomial_config, seed = 42)
  res2 <- simulate_success_at_n(alpha = 1, n = 10L, config = multinomial_config, seed = 42)
  expect_identical(res1$success, res2$success)
  expect_identical(res1$rep_out$phat, res2$rep_out$phat)

  res3 <- simulate_success_at_n(alpha = 1, n = 10L, config = multinomial_config, seed = 999)
  # Compare the underlying draws (phat), not success: with B = 5 the pass/fail pattern could coincidentally tie
  # across seeds, but the simulated proportions themselves are effectively never identical.
  expect_false(isTRUE(all.equal(res1$rep_out$phat, res3$rep_out$phat)))

  res_dm1 <- simulate_success_at_n(alpha = 1, n = 10L, config = dm_config, seed = 42)
  res_dm2 <- simulate_success_at_n(alpha = 1, n = 10L, config = dm_config, seed = 42)
  expect_identical(res_dm1$success, res_dm2$success)
  expect_identical(res_dm1$rep_out$person_results, res_dm2$rep_out$person_results)

  res_dm3 <- simulate_success_at_n(alpha = 1, n = 10L, config = dm_config, seed = 999)
  expect_false(isTRUE(all.equal(res_dm1$rep_out$person_results$error, res_dm3$rep_out$person_results$error)))
})

test_that("common random numbers: the same seed is deterministic at each of two sample sizes", {
  # run_replicates() derives its per-replicate RNG streams from `seed` alone (see replicate_streams()), so the same
  # seed gives common random numbers across n. We don't assert anything about how success compares across n
  # (that would be a stochastic claim); we only check that a fixed seed reproduces exactly at each n in turn.
  res_small_a <- simulate_success_at_n(alpha = 1, n = 5L, config = multinomial_config, seed = 123)
  res_small_b <- simulate_success_at_n(alpha = 1, n = 5L, config = multinomial_config, seed = 123)
  expect_identical(res_small_a$success, res_small_b$success)
  expect_identical(res_small_a$rep_out$phat, res_small_b$rep_out$phat)

  res_large_a <- simulate_success_at_n(alpha = 1, n = 50L, config = multinomial_config, seed = 123)
  res_large_b <- simulate_success_at_n(alpha = 1, n = 50L, config = multinomial_config, seed = 123)
  expect_identical(res_large_a$success, res_large_b$success)
  expect_identical(res_large_a$rep_out$phat, res_large_b$rep_out$phat)

  res_dm_small_a <- simulate_success_at_n(alpha = 1, n = 5L, config = dm_config, seed = 123)
  res_dm_small_b <- simulate_success_at_n(alpha = 1, n = 5L, config = dm_config, seed = 123)
  expect_identical(res_dm_small_a$success, res_dm_small_b$success)

  res_dm_large_a <- simulate_success_at_n(alpha = 1, n = 50L, config = dm_config, seed = 123)
  res_dm_large_b <- simulate_success_at_n(alpha = 1, n = 50L, config = dm_config, seed = 123)
  expect_identical(res_dm_large_a$success, res_dm_large_b$success)
})

test_that("huge taus give all-TRUE and zero taus give all-FALSE for the multinomial model", {
  huge_config <- multinomial_config
  huge_config$taus <- list(AE = 1e6, ARE = 1e6)
  res_huge <- simulate_success_at_n(alpha = 1, n = 10L, config = huge_config, seed = 1)
  expect_true(all(res_huge$success))

  zero_config <- multinomial_config
  zero_config$taus <- list(AE = 0, ARE = 0)
  res_zero <- simulate_success_at_n(alpha = 1, n = 10L, config = zero_config, seed = 1)
  expect_true(all(!res_zero$success))
})

test_that("huge taus give all-TRUE and zero taus give all-FALSE for the dirichlet_multinomial model", {
  huge_config <- dm_config
  huge_config$taus <- list(AE = 1e6, ARE = 1e6)
  res_huge <- simulate_success_at_n(alpha = 1, n = 10L, config = huge_config, seed = 1)
  expect_true(all(res_huge$success))

  zero_config <- dm_config
  zero_config$taus <- list(AE = 0, ARE = 0)
  res_zero <- simulate_success_at_n(alpha = 1, n = 10L, config = zero_config, seed = 1)
  expect_true(all(!res_zero$success))
})

test_that("a metric without a tau is skipped with a warning and does not affect success", {
  partial_config <- multinomial_config
  partial_config$metrics <- c("AE", "ARE")
  partial_config$taus <- list(AE = 0.3)

  expect_warning(
    res <- simulate_success_at_n(alpha = 1, n = 10L, config = partial_config, seed = 42),
    "ARE.*no threshold"
  )

  expected <- res$rep_out$max_errors[, "AE"] <= 0.3
  expect_identical(res$success, expected)
})

test_that("the old required_person_fraction / per-person success rule has been removed", {
  expect_false(exists("validate_required_person_fraction"))
})
