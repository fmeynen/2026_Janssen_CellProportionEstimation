source(here::here("scripts", "simulations", "simulation_dm_errorchoice.R"))

small_dm_config <- function(...) {
  utils::modifyList(
    simulation_dm_errorchoice_pmin_defaults(),
    c(
      list(
        B = 5L, alpha = c(2, 4), n_people = 2L, n_per_person_grid = 100L, n_per_person_fixed = 1000L,
        taus = list(AE = NULL, ARE = NULL)
      ),
      list(...)
    )
  )
}

test_that("run_dm_errorchoice_experiment pins the smallest type to p_min with fixed_min_beta", {
  res <- suppressMessages(run_dm_errorchoice_experiment(
    alpha = c(2, 4), K = 10L, B = 5L, metrics = "AE", proportion_method = "fixed_min_beta", p_min = 0.01,
    n_people = 2L, n_per_person = 100L, concentration = 1e4, seed = 1L
  ))
  cell_cols <- paste0("cell_type_", 1:10)
  expect_equal(apply(res$p_table[, cell_cols], 1, min), c(0.01, 0.01))
  expect_equal(unname(rowSums(res$p_table[, cell_cols])), c(1, 1))
  expect_equal(res$inputs$p_min, 0.01)
  expect_null(res$inputs$p_max)
})

test_that("simulation_dm_errorchoice_pmin_defaults differs from the defaults only in method and p_min", {
  base <- simulation_dm_errorchoice_defaults()
  pmin <- simulation_dm_errorchoice_pmin_defaults()
  expect_named(pmin, names(base))
  expect_equal(pmin$proportion_method, "fixed_min_beta")
  expect_equal(pmin$p_min, 0.01)
  same <- setdiff(names(base), c("proportion_method", "p_min"))
  expect_equal(pmin[same], base[same])
})

test_that("run_simulation_dm_errorchoice keys the cache on p_min", {
  dir <- withr::local_tempdir()
  run <- function(p_min) {
    suppressMessages(run_simulation_dm_errorchoice(small_dm_config(p_min = p_min), cache_dir = dir))
  }
  r1 <- run(0.01)
  r2 <- run(0.02)
  expect_length(list.files(dir), 2L)
  expect_false(isTRUE(all.equal(r1$p_table, r2$p_table)))
  expect_equal(apply(r2$p_table[, -1], 1, min), c(0.02, 0.02))
  expect_equal(run(0.01)$p_table, r1$p_table)
  expect_length(list.files(dir), 2L)
})

test_that("run_simulation_dm_errorchoice still works for plain beta with NULL bounds", {
  dir <- withr::local_tempdir()
  cfg <- small_dm_config(proportion_method = "beta", p_min = NULL)
  res <- suppressMessages(run_simulation_dm_errorchoice(cfg, cache_dir = dir))
  expect_true(all(c("p_table", "stats", "curves_tau", "curves_n") %in% names(res)))
  expect_length(list.files(dir), 1L)
})

test_that("run_dm_errorchoice_experiment rejects bad bounds before simulating", {
  for (case in bound_error_cases) {
    expect_error(
      run_dm_errorchoice_experiment(
        alpha = 2, K = 10L, B = 5L, metrics = "AE", proportion_method = case$method, p_min = case$p_min,
        p_max = case$p_max, n_people = 2L, n_per_person = 100L, concentration = 1e4, seed = 1L
      ),
      case$regex
    )
  }
})

# run_dm_errorchoice_samplesize ---------------------------------------------------------------------------------

#' Fake simulator recording what the runner forwards. Success follows plogis in log(n) and reaches 1, except for the
#' metrics in `flat_metrics`, where it stays flat at 0.5 * B (so the target is unreachable).
make_recording_sim <- function(flat_metrics = character()) {
  seen <- list()
  sim <- function(alpha, n, config, seed) {
    seen[[length(seen) + 1L]] <<- list(
      taus = names(config$taus), metrics = config$metrics, n_people = config$n_people, seed = seed
    )
    s <- if (any(names(config$taus) %in% flat_metrics)) {
      as.integer(config$B / 2)
    } else {
      as.integer(round(config$B * stats::plogis(-5 + 1.5 * log(n) - 0.3 * alpha)))
    }
    list(success_count = s, success_rate = s / config$B)
  }
  list(sim = sim, seen = function() seen)
}

samplesize_config <- function() {
  utils::modifyList(
    simulation_dm_errorchoice_defaults(),
    list(alpha = c(2, 3), n_people = c(1L, 2L), B = 20L)
  )
}

test_that("run_dm_errorchoice_samplesize returns one row per (alpha, n_people, metric) and isolates metrics", {
  fake <- make_recording_sim()
  res <- run_dm_errorchoice_samplesize(samplesize_config(), cache_dir = withr::local_tempdir(), simulate = fake$sim)
  expect_named(
    res,
    c("alpha", "n_people", "metric", "sample_size", "stopping_reason", "iterations_used", "success_ceiling")
  )
  expect_equal(nrow(res), 8L)
  expect_equal(nrow(unique(res[, c("alpha", "n_people", "metric")])), 8L)
  expect_setequal(res$metric, c("AE", "ARE"))
  expect_setequal(res$n_people, c(1L, 2L))
  expect_false(anyNA(res$sample_size))

  seen <- fake$seen()
  expect_true(all(vapply(seen, function(s) length(s$metrics) == 1L && identical(s$taus, s$metrics), logical(1))))
  expect_setequal(vapply(seen, function(s) s$metrics, character(1)), c("AE", "ARE"))
  expect_setequal(vapply(seen, function(s) s$n_people, numeric(1)), c(1, 2))
  expect_true(all(vapply(seen, function(s) s$seed == 260926L, logical(1))))
})

test_that("run_dm_errorchoice_samplesize keeps infeasible rows with NA sample_size", {
  fake <- make_recording_sim(flat_metrics = "ARE")
  expect_warning(
    res <- run_dm_errorchoice_samplesize(samplesize_config(), cache_dir = withr::local_tempdir(), simulate = fake$sim),
    "infeasible"
  )
  expect_equal(nrow(res), 8L)
  are <- res[res$metric == "ARE", ]
  ae <- res[res$metric == "AE", ]
  expect_true(all(is.na(are$sample_size)))
  expect_true(all(are$stopping_reason == "infeasible"))
  expect_false(anyNA(ae$sample_size))
  expect_false(any(ae$stopping_reason == "infeasible"))
})
