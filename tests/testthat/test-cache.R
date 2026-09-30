# Tests for cached_result() in scripts/simulation_layers/orchestration.R.

#' Build a zero-argument compute function that counts its calls.
#'
#' @param value Value returned by the compute function.
#' @return List with `compute` and `calls` (a zero-argument function giving the call count).
make_counting_compute <- function(value = list(x = 1)) {
  n_calls <- 0L
  list(
    compute = function() {
      n_calls <<- n_calls + 1L
      value
    },
    calls = function() n_calls
  )
}

test_that("computes and writes the cache file when absent", {
  dir <- withr::local_tempdir()
  cc <- make_counting_compute(list(x = 42))

  res <- cached_result(list(a = 1), "demo", cc$compute, dir = dir)

  expect_identical(res, list(x = 42))
  expect_identical(cc$calls(), 1L)
  expect_length(list.files(dir, pattern = "^demo_.*\\.rds$"), 1L)
})

test_that("reads the cache without calling compute when present", {
  dir <- withr::local_tempdir()
  cached_result(list(a = 1), "demo", make_counting_compute(list(x = 42))$compute, dir = dir)

  cc <- make_counting_compute(list(x = 99))
  res <- cached_result(list(a = 1), "demo", cc$compute, dir = dir)

  expect_identical(cc$calls(), 0L)
  expect_identical(res, list(x = 42))
})

test_that("force_recompute recomputes and overwrites the cache", {
  dir <- withr::local_tempdir()
  cached_result(list(a = 1), "demo", make_counting_compute(list(x = 1))$compute, dir = dir)

  cc <- make_counting_compute(list(x = 2))
  res <- cached_result(list(a = 1), "demo", cc$compute, force_recompute = TRUE, dir = dir)
  expect_identical(cc$calls(), 1L)
  expect_identical(res, list(x = 2))

  cc3 <- make_counting_compute(list(x = 3))
  expect_identical(cached_result(list(a = 1), "demo", cc3$compute, dir = dir), list(x = 2))
  expect_identical(cc3$calls(), 0L)
})

test_that("cache = FALSE never writes and never reads", {
  dir <- withr::local_tempdir()
  cached_result(list(a = 1), "demo", make_counting_compute()$compute, dir = dir)
  files_before <- list.files(dir)

  cc <- make_counting_compute(list(x = 5))
  res <- cached_result(list(a = 1), "demo", cc$compute, cache = FALSE, dir = dir)
  expect_identical(cc$calls(), 1L)
  expect_identical(res, list(x = 5))
  expect_identical(list.files(dir), files_before)

  empty <- withr::local_tempdir()
  cached_result(list(a = 2), "demo", make_counting_compute()$compute, cache = FALSE, dir = empty)
  expect_length(list.files(empty), 0L)
})

test_that("a different key resolves to a different file", {
  dir <- withr::local_tempdir()
  cached_result(list(a = 1), "demo", make_counting_compute()$compute, dir = dir)
  cc <- make_counting_compute()
  cached_result(list(a = 2), "demo", cc$compute, dir = dir)

  expect_identical(cc$calls(), 1L)
  expect_length(list.files(dir), 2L)
})

test_that("CACHE_SCHEMA is part of the hashed key", {
  dir <- withr::local_tempdir()
  key <- list(a = 1)
  cached_result(key, "demo", make_counting_compute()$compute, dir = dir)

  written <- list.files(dir)
  expect_length(written, 1L)
  expect_false(identical(written, basename(simulation_result_path(key, dir, "demo"))))
  expect_identical(
    written,
    basename(simulation_result_path(c(key, list(cache_schema = CACHE_SCHEMA)), dir, "demo"))
  )
})
