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
