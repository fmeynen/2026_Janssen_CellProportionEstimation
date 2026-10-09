# Tests for run_sample_size_experiment() in scripts/simulation_layers/orchestration.R:
# grid orchestration over alpha with warm start and per-alpha caching.


#' Build an experiment config with sensible test defaults.
#'
#' @param ... Fields overriding the defaults.
make_experiment_config <- function(...) {
  cfg <- list(
    alpha                = c(2, 3),
    n_init               = 200,
    B                    = 200L,
    success_rate_target  = 0.95,
    rel_tol              = 0.01,
    max_iterations       = 15L,
    f0                   = 2,
    f_floor              = 1.1,
    seed                 = 42L
  )
  utils::modifyList(cfg, list(...))
}

#' Deterministic fake simulator: success rate plogis(-5 + 1.5*log(n) - 0.3*alpha), so the curve (and hence n*)
#' shifts with alpha. Counts calls via a mutable closure so tests can assert whether `simulate` ran at all.
#'
#' @return List with `sim` (the simulator function) and `calls` (a zero-argument function returning the call
#'   count so far).
make_counting_sim <- function() {
  calls <- 0L
  sim <- function(alpha, n, config, seed) {
    calls <<- calls + 1L
    s <- as.integer(round(config$B * stats::plogis(-5 + 1.5 * log(n) - 0.3 * alpha)))
    list(success_count = s, success_rate = s / config$B)
  }
  list(sim = sim, calls = function() calls)
}


# Shape and types ----------------------------------------------------------------------------------------------

test_that("returns the documented shape and types", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config()
  fake <- make_counting_sim()

  res <- run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake$sim)

  expect_named(res, c("sample_size", "diagnostics"))
  expect_s3_class(res$sample_size, "data.frame")
  expect_named(res$sample_size, c("alpha", "sample_size", "stopping_reason", "iterations_used", "success_ceiling"))
  expect_identical(nrow(res$sample_size), length(cfg$alpha))
  expect_identical(res$sample_size$alpha, cfg$alpha)
  expect_type(res$sample_size$sample_size, "integer")
  expect_true(all(res$sample_size$sample_size >= 1L))
  expect_type(res$sample_size$stopping_reason, "character")
  expect_true(all(res$sample_size$stopping_reason %in% c("tolerance", "max_iterations")))

  expect_s3_class(res$diagnostics, "data.frame")
  expect_true(all(c("alpha", "iteration", "step", "n", "success_count", "n_next") %in% names(res$diagnostics)))
  expect_setequal(unique(res$diagnostics$alpha), cfg$alpha)
})


# Warm start -----------------------------------------------------------------------------------------------------

test_that("the second alpha's first-iteration pilots are centred on the first alpha's final_n", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config()
  fake <- make_counting_sim()

  res <- run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake$sim)

  final_n_alpha1 <- res$sample_size$sample_size[res$sample_size$alpha == cfg$alpha[[1]]]
  d2_iter1 <- res$diagnostics[res$diagnostics$alpha == cfg$alpha[[2]] & res$diagnostics$iteration == 1L, ]

  # sample_size_pilots(n, f) always includes ceiling(n) itself as the centre pilot.
  expect_true(final_n_alpha1 %in% d2_iter1$n)
})

test_that("the first alpha's first-iteration pilots are centred on config$n_init", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config(n_init = 777)
  fake <- make_counting_sim()

  res <- run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake$sim)

  d1_iter1 <- res$diagnostics[res$diagnostics$alpha == cfg$alpha[[1]] & res$diagnostics$iteration == 1L, ]
  expect_true(777L %in% d1_iter1$n)
})


# Caching --------------------------------------------------------------------------------------------------------

test_that("per-alpha cache files are created, one per alpha", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config()
  fake <- make_counting_sim()

  run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake$sim)

  cached_files <- list.files(cache_dir, pattern = "^sample_size_.*\\.rds$")
  expect_length(cached_files, length(cfg$alpha))
})

test_that("a second run reads the cache and never calls simulate", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config()
  fake1 <- make_counting_sim()
  res1 <- run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake1$sim)
  expect_gt(fake1$calls(), 0L)

  fake2 <- make_counting_sim()
  res2 <- run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake2$sim)

  expect_identical(fake2$calls(), 0L)
  expect_equal(res1$sample_size, res2$sample_size)
})

test_that("force_recompute ignores the cache and recomputes", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config()
  fake1 <- make_counting_sim()
  run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake1$sim)

  fake2 <- make_counting_sim()
  run_sample_size_experiment(cfg, cache_dir = cache_dir, force_recompute = TRUE, simulate = fake2$sim)

  expect_gt(fake2$calls(), 0L)
})

test_that("cache = FALSE writes nothing", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config()
  fake <- make_counting_sim()

  run_sample_size_experiment(cfg, cache = FALSE, cache_dir = cache_dir, simulate = fake$sim)

  expect_length(list.files(cache_dir), 0L)
})

test_that("changing the first alpha's n_init changes the cache files for later alphas", {
  cache_dir <- withr::local_tempdir()
  fake <- make_counting_sim()

  cfg1 <- make_experiment_config(n_init = 100)
  run_sample_size_experiment(cfg1, cache_dir = cache_dir, simulate = fake$sim)
  files_after_1 <- list.files(cache_dir)
  expect_length(files_after_1, length(cfg1$alpha))

  cfg2 <- make_experiment_config(n_init = 5000)
  run_sample_size_experiment(cfg2, cache_dir = cache_dir, simulate = fake$sim)
  files_after_2 <- list.files(cache_dir)

  new_files <- setdiff(files_after_2, files_after_1)
  # Both alphas' warm-start chains differ from the first run, so both get new cache files.
  expect_length(new_files, length(cfg2$alpha))
})


# Real simulator integration -------------------------------------------------------------------------------------

test_that("runs end-to-end with the real simulator and returns positive integer sample sizes", {
  cache_dir <- withr::local_tempdir()
  cfg <- list(
    alpha                = c(1, 2),
    K                    = 3L,
    n_init               = 200,
    B                    = 50L,
    taus                 = list(AE = 0.2, ARE = 1),
    metrics              = c("AE", "ARE"),
    model                = "multinomial",
    tie_method           = "random",
    proportion_method    = "beta",
    seed                 = 99L,
    success_rate_target  = 0.95,
    rel_tol              = 0.05,
    max_iterations       = 10L,
    f0                   = 2,
    f_floor              = 1.1
  )

  start_time <- Sys.time()
  res <- run_sample_size_experiment(cfg, cache_dir = cache_dir)
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

  expect_identical(nrow(res$sample_size), 2L)
  expect_true(all(res$sample_size$sample_size >= 1L))
  expect_type(res$sample_size$sample_size, "integer")
  message(sprintf("real-simulator integration test elapsed: %.2fs", elapsed))
})


# Bound validation and cache key ---------------------------------------------------------------------------------

make_failing_sim <- function() {
  function(alpha, n, config, seed) stop("should not be called")
}

test_that("run_sample_size_experiment rejects bad bounds before any simulation or cache write", {
  for (case in bound_error_cases) {
    cache_dir <- withr::local_tempdir()
    cfg <- make_experiment_config(K = 10L, proportion_method = case$method)
    cfg["p_min"] <- list(case$p_min)
    cfg["p_max"] <- list(case$p_max)
    expect_error(
      run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = make_failing_sim()),
      case$regex
    )
    expect_length(list.files(cache_dir), 0L)
  }
})

test_that("p_min and p_max each change the sample-size cache key", {
  cache_dir <- withr::local_tempdir()
  fake <- make_counting_sim()
  run <- function(...) {
    cfg <- make_experiment_config(alpha = 2, K = 10L, ...)
    run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake$sim)
  }

  run(proportion_method = "fixed_min_beta", p_min = 0.01)
  expect_length(list.files(cache_dir), 1L)
  calls_before <- fake$calls()
  run(proportion_method = "fixed_min_beta", p_min = 0.02)
  expect_length(list.files(cache_dir), 2L)
  expect_gt(fake$calls(), calls_before)

  run(proportion_method = "fixed_max_beta", p_max = 0.3)
  run(proportion_method = "fixed_max_beta", p_max = 0.4)
  expect_length(list.files(cache_dir), 4L)

  calls_before <- fake$calls()
  run(proportion_method = "fixed_min_beta", p_min = 0.01)
  expect_identical(fake$calls(), calls_before)
  expect_length(list.files(cache_dir), 4L)
})

test_that("run_sample_size_experiment runs with the real simulator for fixed_min_beta", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config(
    alpha = 2, K = 10L, B = 5L, taus = list(AE = 0.3), metrics = "AE", model = "dirichlet_multinomial",
    n_people = 2L, concentration = 50, tie_method = "random", proportion_method = "fixed_min_beta", p_min = 0.01,
    n_init = 50, max_iterations = 2L, rel_tol = 0.5
  )
  res <- suppressWarnings(run_sample_size_experiment(cfg, cache_dir = cache_dir))
  expect_identical(nrow(res$sample_size), 1L)
})


# n_init default -------------------------------------------------------------------------------------------------

#' Fake simulator wrapping the plogis curve that records every `(alpha, n)` call in a mutable log.
#'
#' @param ceiling_count Optional function `(alpha, config)` returning the success count to report at `n_max`.
#' @param flat_for_alpha Optional alpha for which the success count is `round(0.5 * B)` at every n.
#' @return List with `sim` and `log` (zero-argument function returning a data.frame of `alpha`, `n`).
make_logging_sim <- function(ceiling_count = NULL, flat_for_alpha = NULL, n_max = 1e9) {
  log <- data.frame(alpha = numeric(0), n = numeric(0))
  sim <- function(alpha, n, config, seed) {
    log[nrow(log) + 1L, ] <<- list(alpha, n)
    s <- if (!is.null(flat_for_alpha) && alpha == flat_for_alpha) {
      round(0.5 * config$B)
    } else if (!is.null(ceiling_count) && n == n_max) {
      ceiling_count(alpha, config)
    } else {
      round(config$B * stats::plogis(-5 + 1.5 * log(n) - 0.3 * alpha))
    }
    list(success_count = as.integer(s), success_rate = s / config$B)
  }
  list(sim = sim, log = function() log)
}

test_that("a NULL n_init resolves to config$concentration for the Dirichlet-multinomial model", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config(alpha = 2, model = "dirichlet_multinomial", concentration = 321, n_people = 5L)
  cfg$n_init <- NULL
  fake <- make_logging_sim()

  res <- run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake$sim)

  d1 <- res$diagnostics[res$diagnostics$iteration == 1L, ]
  expect_true(321L %in% d1$n)
})

test_that("a NULL n_init errors for the multinomial model", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config(model = "multinomial")
  cfg$n_init <- NULL

  expect_error(
    run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = make_failing_sim()),
    "n_init"
  )
})

test_that("an explicit n_init wins over config$concentration", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config(
    alpha = 2, model = "dirichlet_multinomial", concentration = 321, n_people = 5L, n_init = 444
  )
  fake <- make_logging_sim()

  res <- run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake$sim)

  d1 <- res$diagnostics[res$diagnostics$iteration == 1L, ]
  expect_true(444L %in% d1$n)
  expect_false(321L %in% d1$n)
})


# Feasibility check ----------------------------------------------------------------------------------------------

test_that("the simulate hook is called once per alpha at n_max before the solver", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config(n_max = 100000)
  fake <- make_logging_sim(n_max = 100000)

  res <- run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake$sim)

  log <- fake$log()
  for (a in cfg$alpha) {
    expect_identical(sum(log$alpha == a & log$n == 100000), 1L)
    expect_identical(which(log$alpha == a)[[1L]], which(log$alpha == a & log$n == 100000)[[1L]])
  }
  expect_false(anyNA(res$sample_size$success_ceiling))
})

test_that("a clearly infeasible alpha skips the solver, returns NA and warns", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config(alpha = 2, concentration = 50, n_people = 3L)
  fake <- make_logging_sim(flat_for_alpha = 2)

  expect_warning(
    res <- run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake$sim),
    "infeasible.*cannot be reached"
  )

  n_calls <- nrow(fake$log())
  expect_identical(n_calls, 1L)
  expect_identical(fake$log()$n, 1e9)
  expect_true(is.na(res$sample_size$sample_size))
  expect_type(res$sample_size$sample_size, "integer")
  expect_identical(res$sample_size$stopping_reason, "infeasible")
  expect_identical(res$sample_size$iterations_used, 0L)
  expect_equal(res$sample_size$success_ceiling, 0.5)
  expect_null(res$diagnostics)
})

test_that("a borderline alpha (point estimate below target, upper bound above) still runs the solver", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config(alpha = 2, B = 100L)
  expect_lt(93 / 100, cfg$success_rate_target)
  expect_gte(stats::qbeta(0.95, 94, 7), cfg$success_rate_target)
  fake <- make_logging_sim(ceiling_count = function(alpha, config) 93L)

  res <- run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake$sim)

  expect_false(identical(res$sample_size$stopping_reason, "infeasible"))
  expect_false(is.na(res$sample_size$sample_size))
  expect_equal(res$sample_size$success_ceiling, 0.93)
  expect_gt(nrow(res$diagnostics), 0L)
})

test_that("the warm start skips infeasible alphas and carries the last feasible final_n", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config(alpha = c(2, 3, 4))
  fake <- make_logging_sim(flat_for_alpha = 3)

  res <- expect_warning(
    run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake$sim),
    "infeasible"
  )

  expect_identical(res$sample_size$stopping_reason[[2L]], "infeasible")
  expect_true(is.na(res$sample_size$sample_size[[2L]]))
  expect_false(anyNA(res$sample_size$sample_size[-2L]))
  expect_false(anyNA(res$sample_size$success_ceiling))
  expect_identical(unique(res$diagnostics$alpha), c(2, 4))

  final_n_alpha1 <- res$sample_size$sample_size[[1L]]
  d3_iter1 <- res$diagnostics[res$diagnostics$alpha == 4 & res$diagnostics$iteration == 1L, ]
  expect_true(final_n_alpha1 %in% d3_iter1$n)
})

test_that("a cached infeasible alpha still warns and does not call simulate again", {
  cache_dir <- withr::local_tempdir()
  cfg <- make_experiment_config(alpha = 2)
  fake1 <- make_logging_sim(flat_for_alpha = 2)
  suppressWarnings(run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake1$sim))

  fake2 <- make_logging_sim(flat_for_alpha = 2)
  expect_warning(
    res <- run_sample_size_experiment(cfg, cache_dir = cache_dir, simulate = fake2$sim),
    "infeasible"
  )

  expect_identical(nrow(fake2$log()), 0L)
  expect_identical(res$sample_size$stopping_reason, "infeasible")
})
