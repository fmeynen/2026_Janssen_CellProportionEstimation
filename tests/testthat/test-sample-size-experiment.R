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
  expect_named(res$sample_size, c("alpha", "sample_size", "stopping_reason", "iterations_used"))
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
