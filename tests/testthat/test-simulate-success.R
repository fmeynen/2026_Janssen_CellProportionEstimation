# Tests for simulate_success_at_n() -----------------------------------------------------------------------------
# Both models must route through the shared pooled-proportion success rule (pooled_error_stat() in calculation.R);
# the old required_person_fraction / per-person rule has been removed entirely. Also tests pass_from_max_errors().

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

test_that("multinomial success matches replicate_success() on the old long-format construction", {
  for (B in c(1L, 40L)) {
    cfg   <- multinomial_config
    cfg$B <- B
    # Tight taus so that both pass and fail occur across replicates.
    cfg$taus <- list(AE = 0.08, ARE = 0.5)
    res <- simulate_success_at_n(alpha = 1, n = 20L, config = cfg, seed = 7)

    phat    <- res$rep_out$phat
    p       <- res$rep_out$inputs$p
    grid <- expand.grid(
      replicate = seq_len(nrow(phat)), cell_type = seq_along(p), metric = cfg$metrics,
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    person_results <- data.frame(
      replicate                  = grid$replicate,
      person_id                  = 1L,
      cell_type                  = grid$cell_type,
      metric                     = grid$metric,
      observed_proportion        = phat[cbind(grid$replicate, grid$cell_type)],
      population_mean_proportion = p[grid$cell_type],
      stringsAsFactors = FALSE
    )
    old <- as.logical(replicate_success(person_results, cfg$taus)$pass)

    expect_identical(res$success, old)
    expect_identical(res$success_count, sum(old))
    expect_equal(res$success_rate, mean(old))
  }
})

test_that("dirichlet_multinomial success matches pooled_error_stat() on rep_out$phat", {
  res <- simulate_success_at_n(alpha = 1, n = 10L, config = dm_config, seed = 42)

  phat     <- res$rep_out$phat
  p        <- res$rep_out$inputs$p
  expected <- rep(TRUE, nrow(phat))
  for (m in names(dm_config$taus)) {
    stat <- pooled_error_stat(phat, p, m)
    expect_identical(unname(res$rep_out$max_errors[, m]), stat)
    expected <- expected & (stat <= dm_config$taus[[m]])
  }

  expect_identical(res$success, expected)
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
  expect_identical(res_dm1$rep_out$phat, res_dm2$rep_out$phat)

  res_dm3 <- simulate_success_at_n(alpha = 1, n = 10L, config = dm_config, seed = 999)
  expect_false(isTRUE(all.equal(res_dm1$rep_out$max_errors, res_dm3$rep_out$max_errors)))
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

test_that("config$p_max is passed to generate_proportions() for proportion_method = 'fixed_max_beta'", {
  fixed_config <- dm_config
  fixed_config$proportion_method <- "fixed_max_beta"
  fixed_config$p_max <- 0.4

  res <- simulate_success_at_n(alpha = 1, n = 10L, config = fixed_config, seed = 42)

  expect_gte(res$success_rate, 0)
  expect_lte(res$success_rate, 1)
  # The true proportions used are exposed as rep_out$inputs$p; the largest must be exactly p_max.
  expect_equal(max(res$rep_out$inputs$p), fixed_config$p_max)
})


# pass_from_max_errors() ------------------------------------------------------------------------------------------

test_that("pass_from_max_errors() ANDs max_errors[, m] <= taus[[m]] across metrics", {
  max_errors <- matrix(
    c(0.01, 0.03, 0.02, 0.01,
      0.5,  0.5,  3,    2),
    ncol = 2, dimnames = list(NULL, c("AE", "ARE"))
  )
  taus <- list(AE = 0.02, ARE = 2)

  # Row 2 fails AE, row 3 fails ARE; the boundary (== tau) passes.
  expect_identical(pass_from_max_errors(max_errors, taus), c(TRUE, FALSE, FALSE, TRUE))
  expect_identical(pass_from_max_errors(max_errors, list(AE = 0.02)), c(TRUE, FALSE, TRUE, TRUE))
})

test_that("pass_from_max_errors() agrees with replicate_success() on the same long-format data", {
  # Two replicates, two persons, three cell types; one AE and one ARE row per (replicate, person, cell type).
  p <- c(0.2, 0.3, 0.5)
  observed <- list(
    rbind(c(0.25, 0.30, 0.45), c(0.15, 0.35, 0.50)),
    rbind(c(0.40, 0.20, 0.40), c(0.30, 0.30, 0.40))
  )
  rows <- list()
  for (b in 1:2) {
    for (person in 1:2) {
      for (m in c("AE", "ARE")) {
        rows[[length(rows) + 1L]] <- data.frame(
          replicate = b, person_id = person, cell_type = 1:3, metric = m,
          observed_proportion = observed[[b]][person, ], population_mean_proportion = p,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  person_results <- do.call(rbind, rows)
  phat <- rbind(colMeans(observed[[1]]), colMeans(observed[[2]]))
  max_errors <- cbind(AE = pooled_error_stat(phat, p, "AE"), ARE = pooled_error_stat(phat, p, "ARE"))

  for (taus in list(list(AE = 0.05, ARE = 0.5), list(AE = 0.2, ARE = 0.1), list(AE = 0.2, ARE = 1))) {
    expect_identical(
      pass_from_max_errors(max_errors, taus),
      replicate_success(person_results, taus)$pass
    )
  }
})

test_that("pass_from_max_errors() skips taus metrics absent from max_errors with a warning, and errors if none remain", {
  max_errors <- matrix(c(0.01, 0.05), ncol = 1, dimnames = list(NULL, "AE"))

  expect_warning(
    pass <- pass_from_max_errors(max_errors, list(AE = 0.02, ARE = 1)),
    "metric 'ARE' is not in max_errors"
  )
  expect_identical(pass, c(TRUE, FALSE))

  expect_error(
    suppressWarnings(pass_from_max_errors(max_errors, list(ARE = 1))),
    "None of the metrics in taus are present in max_errors."
  )
})

test_that("pass_from_max_errors() ignores max_errors columns without a tau", {
  max_errors <- matrix(c(0.01, 0.05, 100, 100), ncol = 2, dimnames = list(NULL, c("AE", "ARE")))
  expect_no_warning(pass <- pass_from_max_errors(max_errors, list(AE = 0.02)))
  expect_identical(pass, c(TRUE, FALSE))
})

test_that("pass_from_max_errors() fails Inf statistics and treats pooled_error_stat()'s NaN -> 0 as passing", {
  # ARE with p_j = 0: pbar_j = 0 gives NaN -> 0 (passes); pbar_j > 0 gives Inf (fails).
  p <- c(0, 0.4, 0.6)
  phat <- rbind(c(0, 0.4, 0.6), c(0.1, 0.3, 0.6))
  max_errors <- cbind(ARE = pooled_error_stat(phat, p, "ARE"))
  expect_identical(max_errors[, "ARE"], c(0, Inf))
  expect_identical(pass_from_max_errors(max_errors, list(ARE = 1e6)), c(TRUE, FALSE))
})

test_that("pass_from_max_errors() validates its inputs", {
  max_errors <- matrix(0.01, dimnames = list(NULL, "AE"))
  expect_error(pass_from_max_errors(unname(max_errors), list(AE = 0.02)), "named column")
  expect_error(pass_from_max_errors(c(AE = 0.01), list(AE = 0.02)), "numeric matrix")
  expect_error(pass_from_max_errors(max_errors, c(AE = 0.02)), "named list")
  expect_error(pass_from_max_errors(max_errors, list(AE = c(0.01, 0.02))), "single numeric threshold")
})
