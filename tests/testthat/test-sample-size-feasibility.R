# Tests for the feasibility check and p_min support of run_sample_size_experiment() with the real Dirichlet-multinomial
# simulator (no fake `simulate`). Under this model n is cells per person: with few people and a low concentration
# the success rate plateaus below 1 however large n gets, because between-person Dirichlet(concentration * p)
# variation does not shrink with n.


#' Build a real-DM experiment config.
#'
#' @param ... Fields overriding the defaults.
make_dm_config <- function(...) {
  cfg <- list(
    alpha                = 2,
    K                    = 4L,
    B                    = 100L,
    taus                 = list(AE = 0.01),
    metrics              = "AE",
    model                = "dirichlet_multinomial",
    tie_method           = "random",
    proportion_method    = "beta",
    n_people             = 1L,
    concentration        = 10,
    seed                 = 7L,
    success_rate_target  = 0.95,
    rel_tol              = 0.05,
    max_iterations       = 10L,
    f0                   = 2,
    f_floor              = 1.1
  )
  utils::modifyList(cfg, list(...))
}


test_that("a low-concentration, single-person DM config is reported infeasible without running the solver", {
  cache_dir <- withr::local_tempdir()
  # concentration 10 with one person: the person's own proportions differ from p by far more than tau_AE = 0.01,
  # so the success rate stays near 0 even at n_max.
  cfg <- make_dm_config()

  expect_warning(
    res <- run_sample_size_experiment(cfg, cache_dir = cache_dir),
    "infeasible"
  )

  expect_identical(res$sample_size$stopping_reason, "infeasible")
  expect_true(is.na(res$sample_size$sample_size))
  expect_identical(res$sample_size$iterations_used, 0L)
  expect_lt(res$sample_size$success_ceiling, 0.5)
  expect_equal(NROW(res$diagnostics), 0L)
})


test_that("a feasible real-DM config converges and the solved n hits the target on a fresh seed", {
  cache_dir <- withr::local_tempdir()
  # Large concentration: the between-person plateau is ~1, so only multinomial noise matters and n* is modest.
  cfg <- make_dm_config(
    n_people = 5L, concentration = 1e4, taus = list(AE = 0.03), B = 200L, rel_tol = 0.05, n_init = 100
  )

  res <- run_sample_size_experiment(cfg, cache_dir = cache_dir)

  expect_identical(res$sample_size$stopping_reason, "tolerance")
  n_star <- res$sample_size$sample_size
  expect_false(is.na(n_star))
  expect_gte(res$sample_size$success_ceiling, 0.95)

  check <- simulate_success_at_n(cfg$alpha, n = n_star, config = cfg, seed = 12345L)
  expect_lte(abs(check$success_rate - cfg$success_rate_target), 0.1)
})


test_that("the real-DM solver converges with fixed_min_beta and the default n_init", {
  cache_dir <- withr::local_tempdir()
  # p_min = 0.05 is feasible for K = 4 (K * p_min < 1). n_init is NULL, so it defaults to the concentration.
  cfg <- make_dm_config(
    n_people = 5L, concentration = 1e4, taus = list(AE = 0.03), B = 200L, rel_tol = 0.05,
    proportion_method = "fixed_min_beta", p_min = 0.05, n_init = NULL
  )

  res <- run_sample_size_experiment(cfg, cache_dir = cache_dir)

  expect_identical(res$sample_size$stopping_reason, "tolerance")
  expect_false(is.na(res$sample_size$sample_size))
})
