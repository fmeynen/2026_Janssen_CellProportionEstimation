# Tests for the iterative sample-size solver in scripts/simulation_layers/calculation.R:
# sample_size_pilots(), fit_success_curve(), solve_success_curve(), estimate_sample_size().


#' Build a solver config with sensible test defaults.
#'
#' @param ... Fields overriding the defaults.
make_solver_config <- function(...) {
  cfg <- list(success_rate_target = 0.95, rel_tol = 0.01, max_iterations = 20L, B = 500L,
              f0 = 2, f_floor = 1.1, seed = 123L)
  utils::modifyList(cfg, list(...))
}

#' Deterministic fake simulator with success rate plogis(a + b * log(n)).
#'
#' @param a,b Logistic intercept and slope on log(n).
make_logistic_sim <- function(a = -20, b = 3) {
  function(alpha, n, config, seed) {
    s <- as.integer(round(config$B * stats::plogis(a + b * log(n))))
    list(success = rep(c(TRUE, FALSE), c(s, config$B - s)), success_count = s, success_rate = s / config$B)
  }
}

#' True n* for the logistic fake at a target success rate.
true_n_star <- function(a = -20, b = 3, target = 0.95) exp((stats::qlogis(target) - a) / b)

#' Fake simulator: 0 successes below `lo`, B at or above `hi`, linear ramp in between.
make_step_sim <- function(lo, hi = lo) {
  function(alpha, n, config, seed) {
    rate <- if (n < lo) 0 else if (n >= hi) 1 else (n - lo) / (hi - lo)
    s <- as.integer(round(config$B * rate))
    list(success_count = s, success_rate = s / config$B)
  }
}


# Helpers ----------------------------------------------------------------------------------------------------------

test_that("sample_size_pilots rounds up, de-duplicates, sorts and floors at 1", {
  expect_identical(sample_size_pilots(100, 2), c(50L, 100L, 200L))
  expect_identical(sample_size_pilots(10, 3), c(4L, 10L, 30L))
  expect_identical(sample_size_pilots(1, 2), c(1L, 2L))
  expect_identical(sample_size_pilots(2, 1.1), c(2L, 3L))
  expect_type(sample_size_pilots(7.2, 1.5), "integer")
  expect_error(sample_size_pilots(10, 1))
})

test_that("fit_success_curve recovers a logistic curve on log(n)", {
  n <- c(100, 200, 400, 800, 1600)
  B <- 1e6
  s <- round(B * stats::plogis(-10 + 2 * log(n)))
  fit <- fit_success_curve(n, s, B)
  expect_s3_class(fit, "glm")
  expect_equal(unname(stats::coef(fit)), c(-10, 2), tolerance = 1e-3)
})

test_that("fit_success_curve suppresses only the separation warnings", {
  expect_no_warning(fit_success_curve(c(10, 20, 40, 80), c(0L, 0L, 100L, 100L), 100L))
  expect_no_warning(fit_success_curve(c(7813, 15625, 31250), c(499L, 500L, 500L), 500L))
  expect_warning(fit_success_curve(c(10, 20, 40), c(10.5, 50, 90), 100L), "non-integer")
})

test_that("solve_success_curve inverts the curve and rounds up", {
  n <- c(100, 200, 400, 800, 1600)
  B <- 1e6
  fit <- fit_success_curve(n, round(B * stats::plogis(-10 + 2 * log(n))), B)
  cf <- stats::coef(fit)
  n_raw <- exp((stats::qlogis(0.95) - cf[[1]]) / cf[[2]])
  out <- solve_success_curve(fit, 0.95)
  expect_type(out, "integer")
  expect_identical(out, as.integer(ceiling(n_raw)))
  expect_gte(out, n_raw)
})

test_that("solve_success_curve refuses a non-positive slope", {
  fit <- fit_success_curve(c(10, 20, 40), c(90L, 50L, 10L), 100L)
  expect_error(solve_success_curve(fit, 0.95), "slope")
})


# Convergence ----------------------------------------------------------------------------------------------------------

test_that("converges to n* from far below and far above", {
  cfg <- make_solver_config()
  n_star <- true_n_star()
  for (n_init in c(10, 1e6)) {
    res <- estimate_sample_size(0.1, n_init, cfg, simulate = make_logistic_sim())
    expect_identical(res$stopping_reason, "tolerance")
    expect_type(res$final_n, "integer")
    expect_lt(abs(res$final_n - n_star) / n_star, 0.03)
  }
})

test_that("converges with a noisy binomial simulator", {
  set.seed(1)
  noisy <- function(alpha, n, config, seed) {
    s <- stats::rbinom(1L, config$B, stats::plogis(-20 + 3 * log(n)))
    list(success_count = s, success_rate = s / config$B)
  }
  res <- estimate_sample_size(0.1, 100, make_solver_config(), simulate = noisy)
  expect_identical(res$stopping_reason, "tolerance")
  expect_lt(abs(res$final_n - true_n_star()) / true_n_star(), 0.1)
})


# Degenerate steps -----------------------------------------------------------------------------------------------------

test_that("expand_up in a flat 0% region jumps by f^2 without clamping", {
  cfg <- make_solver_config()
  res <- estimate_sample_size(0.1, 10, cfg, simulate = make_logistic_sim())
  d <- res$diagnostics
  first <- d[d$iteration == 1L, ]
  expect_true(all(first$step == "expand_up"))
  expect_true(all(first$success_count == 0L))
  # Pilots 5, 10, 20 -> next centre ceiling(2^2 * 20) = 80, not the old 2x clamp (40).
  expect_identical(first$n, c(5L, 10L, 20L))
  expect_identical(unique(first$n_next), 80L)
  expect_identical(unique(d$f[d$iteration == 2L]), cfg$f0)
})

test_that("expand_down in a flat 100% region jumps by 1 / f^2 without clamping", {
  cfg <- make_solver_config()
  res <- estimate_sample_size(0.1, 1e6, cfg, simulate = make_logistic_sim())
  d <- res$diagnostics
  first <- d[d$iteration == 1L, ]
  expect_true(all(first$step == "expand_down"))
  expect_true(all(first$success_count == cfg$B))
  expect_identical(unique(first$n_next), as.integer(ceiling(5e5 / 4)))
  expect_identical(unique(d$f[d$iteration == 2L]), cfg$f0)
})

test_that("expand_down never goes below 1", {
  cfg <- make_solver_config(max_iterations = 3L)
  n_seen <- integer(0L)
  always_pass <- function(alpha, n, config, seed) {
    n_seen <<- c(n_seen, n)
    list(success_count = config$B)
  }
  expect_error(estimate_sample_size(0.1, 2, cfg, simulate = always_pass), "no successful curve fit")
  expect_true(all(n_seen >= 1L))
  expect_true(all(n_seen == as.integer(n_seen)))
})

test_that("bisect step on a step-function simulator uses the geometric mean of the bracket", {
  cfg <- make_solver_config()
  # 0 below 1000, B from 1010 up: pilots almost always land on one side.
  res <- estimate_sample_size(0.1, 10, cfg, simulate = make_step_sim(1000, 1010))
  d <- res$diagnostics
  expect_true("bisect" %in% d$step)
  first_bisect <- min(d$iteration[d$step == "bisect"])
  before <- d[d$iteration <= first_bisect, ]
  n_fail <- max(before$n[before$success_count == 0L])
  n_pass <- min(before$n[before$success_count == cfg$B])
  expect_identical(unique(d$n_next[d$iteration == first_bisect]), as.integer(ceiling(sqrt(n_fail * n_pass))))
  expect_identical(unique(d$f[d$iteration == first_bisect + 1L]), cfg$f0)
  expect_identical(res$stopping_reason, "tolerance")
  expect_gte(res$final_n, 1000L)
  expect_lte(res$final_n, 1010L)
})

test_that("a pure step function that is never resolved errors out", {
  cfg <- make_solver_config(max_iterations = 5L)
  expect_error(estimate_sample_size(0.1, 10, cfg, simulate = make_step_sim(1000)), "no successful curve fit")
})


# Resample ----------------------------------------------------------------------------------------------------------

test_that("a decreasing curve triggers resample with the next seed", {
  cfg <- make_solver_config()
  seeds_seen <- integer(0L)
  flips <- function(alpha, n, config, seed) {
    seeds_seen <<- c(seeds_seen, seed)
    # Seed offset 0: a flat, slightly decreasing curve near the target (as Monte Carlo noise could produce around n*).
    # Afterwards: the true increasing curve.
    rate <- if (seed == config$seed) stats::plogis(3.5 - 0.07 * log(n)) else stats::plogis(-20 + 3 * log(n))
    s <- as.integer(round(config$B * rate))
    list(success_count = s, success_rate = s / config$B)
  }
  res <- estimate_sample_size(0.1, 2000, cfg, simulate = flips)
  d <- res$diagnostics
  expect_identical(d$step[d$iteration == 1L][1], "resample")
  expect_lt(d$glm_slope[d$iteration == 1L][1], 0)
  expect_identical(unique(d$n_next[d$iteration == 1L]), 2000L)
  expect_identical(unique(d$seed_offset[d$iteration == 1L]), 0L)
  expect_identical(unique(d$seed_offset[d$iteration == 2L]), 1L)
  expect_identical(unique(d$f[d$iteration == 2L]), cfg$f0)
  expect_identical(unique(seeds_seen[1:3]), cfg$seed)
  expect_true(all(seeds_seen[-(1:3)] > cfg$seed))
  expect_true("fit" %in% d$step)
  expect_identical(res$stopping_reason, "tolerance")
})

test_that("errors when the slope never becomes positive", {
  cfg <- make_solver_config(max_iterations = 4L)
  decreasing <- make_logistic_sim(a = 20, b = -3)
  expect_error(estimate_sample_size(0.1, 1000, cfg, simulate = decreasing), "no successful curve fit")
})


# Stopping ----------------------------------------------------------------------------------------------------------

test_that("warns and reports max_iterations when too few iterations are allowed", {
  cfg <- make_solver_config(max_iterations = 3L)
  expect_warning(
    res <- estimate_sample_size(0.1, 10, cfg, simulate = make_logistic_sim()),
    "did not converge"
  )
  expect_identical(res$stopping_reason, "max_iterations")
  expect_identical(res$iterations_used, 3L)
  d <- res$diagnostics
  last_fit <- max(d$iteration[d$step == "fit"])
  expect_identical(res$final_n, unique(d$n_next[d$iteration == last_fit]))
})

test_that("f shrinks only after fit steps and never below f_floor", {
  set.seed(2)
  noisy <- function(alpha, n, config, seed) {
    list(success_count = stats::rbinom(1L, config$B, stats::plogis(-20 + 3 * log(n))))
  }
  cfg <- make_solver_config(rel_tol = 0, max_iterations = 12L)
  res <- suppressWarnings(estimate_sample_size(0.1, 10, cfg, simulate = noisy))
  d <- res$diagnostics
  per_iter <- d[!duplicated(d$iteration), c("iteration", "step", "f")]
  expect_identical(per_iter$f[1], cfg$f0)
  expect_true(all(per_iter$f >= cfg$f_floor))
  for (i in seq_len(nrow(per_iter))[-1]) {
    prev <- per_iter[i - 1L, ]
    expected <- if (prev$step == "fit") max(cfg$f_floor, sqrt(prev$f)) else prev$f
    expect_equal(per_iter$f[i], expected)
  }
  expect_true(any(per_iter$f == cfg$f_floor))
  expect_true(any(per_iter$step != "fit"))
})


# Diagnostics ----------------------------------------------------------------------------------------------------------

test_that("diagnostics have the documented shape and final stopping_reason", {
  cfg <- make_solver_config()
  res <- estimate_sample_size(0.25, 10, cfg, simulate = make_logistic_sim())
  d <- res$diagnostics
  expect_named(d, c("alpha", "iteration", "step", "n", "success_count", "success_rate", "f", "seed_offset",
                    "glm_intercept", "glm_slope", "n_next", "stopping_reason"))
  expect_named(res, c("final_n", "stopping_reason", "iterations_used", "diagnostics"))
  expect_true(all(d$alpha == 0.25))
  expect_true(all(d$step %in% c("fit", "expand_up", "expand_down", "bisect", "resample")))
  expect_identical(max(d$iteration), res$iterations_used)
  expect_equal(d$success_rate, d$success_count / cfg$B)
  expect_true(all(is.na(d$glm_slope[d$step %in% c("expand_up", "expand_down", "bisect")])))
  expect_true(all(is.finite(d$glm_slope[d$step == "fit"])))
  expect_true(all(is.na(d$stopping_reason[-nrow(d)])))
  expect_identical(d$stopping_reason[nrow(d)], "tolerance")
  expect_identical(d$step[nrow(d)], "fit")
  expect_identical(d$n_next[nrow(d)], res$final_n)
  # Convergence check holds on the final fit step.
  last <- d[d$iteration == res$iterations_used, ]
  centre <- unique(d$n_next[d$iteration == res$iterations_used - 1L])
  expect_lte(abs(res$final_n - centre) / centre, cfg$rel_tol)
  expect_true(centre %in% last$n)
})


# Validation ----------------------------------------------------------------------------------------------------------

test_that("config is validated up front", {
  sim <- make_logistic_sim()
  expect_error(estimate_sample_size(0.1, 10, make_solver_config(success_rate_target = 1), sim), "success_rate_target")
  expect_error(estimate_sample_size(0.1, 10, make_solver_config(rel_tol = -0.1), sim), "rel_tol")
  expect_error(estimate_sample_size(0.1, 10, make_solver_config(max_iterations = 0), sim), "max_iterations")
  expect_error(estimate_sample_size(0.1, 10, make_solver_config(f0 = 1), sim), "f0")
  expect_error(estimate_sample_size(0.1, 10, make_solver_config(f_floor = 1), sim), "f_floor")
  expect_error(estimate_sample_size(0.1, 10, make_solver_config(f_floor = 3), sim), "f_floor")
  expect_error(estimate_sample_size(0.1, 10, make_solver_config(B = 0), sim), "B")
})
