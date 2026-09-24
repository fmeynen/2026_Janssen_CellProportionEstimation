# tests/testthat/test-replicates.R
#
# Tests for the replicate-level RNG machinery in scripts/simulation_layers/simulation.R:
# replicate_cores(), replicate_streams(), replicate_apply(), and their use inside run_replicates() /
# run_replicates_dirichlet_multinomial().

testthat::local_edition(3)

# Shared fixtures -----------------------------------------------------------------------------------------------

test_p <- generate_proportions_beta(alpha = 2, K = 5)


# Multinomial reproducibility ------------------------------------------------------------------------------------

test_that("run_replicates (multinomial) is reproducible given the same seed", {
  r1 <- run_replicates(p = test_p, n = 200, B = 15, metrics = c("AE", "ARE"), model = "multinomial", seed = 123)
  r2 <- run_replicates(p = test_p, n = 200, B = 15, metrics = c("AE", "ARE"), model = "multinomial", seed = 123)

  expect_identical(r1$max_errors, r2$max_errors)
  expect_identical(r1$argmax, r2$argmax)
  expect_identical(r1$errors, r2$errors)
  expect_identical(r1$phat, r2$phat)
})

test_that("run_replicates (multinomial) output shape and dimnames are unchanged", {
  metrics <- c("AE", "ARE")
  B <- 10L
  K <- length(test_p)
  out <- run_replicates(p = test_p, n = 200, B = B, metrics = metrics, model = "multinomial", seed = 1)

  expect_equal(dim(out$max_errors), c(B, length(metrics)))
  expect_identical(colnames(out$max_errors), metrics)
  expect_equal(dim(out$argmax), c(B, length(metrics)))
  expect_identical(colnames(out$argmax), metrics)
  expect_equal(dim(out$errors), c(B, K, length(metrics)))
  expect_identical(dimnames(out$errors)[[3]], metrics)
  expect_equal(dim(out$phat), c(B, K))
  expect_true(is.integer(out$argmax))
})


# Dirichlet-multinomial reproducibility --------------------------------------------------------------------------

test_that("run_replicates (dirichlet_multinomial) is reproducible given the same seed", {
  dm1 <- run_replicates(
    p = test_p, B = 6, metrics = c("AE", "ARE"), model = "dirichlet_multinomial",
    n_people = 4, n_per_person = 50, concentration = 10, seed = 55
  )
  dm2 <- run_replicates(
    p = test_p, B = 6, metrics = c("AE", "ARE"), model = "dirichlet_multinomial",
    n_people = 4, n_per_person = 50, concentration = 10, seed = 55
  )

  expect_identical(dm1$person_results, dm2$person_results)
})

test_that("run_replicates (dirichlet_multinomial) output shape and columns are correct", {
  metrics <- c("AE", "ARE")
  B <- 5L
  n_people <- 4L
  K <- length(test_p)

  dm <- run_replicates(
    p = test_p, B = B, metrics = metrics, model = "dirichlet_multinomial",
    n_people = n_people, n_per_person = 50, concentration = 10, seed = 1
  )

  expect_s3_class(dm$person_results, "data.frame")
  expect_equal(nrow(dm$person_results), B * n_people * K * length(metrics))
  expect_identical(
    colnames(dm$person_results),
    c(
      "scenario_id", "n_people", "concentration", "replicate", "person_id", "cell_type", "metric",
      "count", "observed_proportion", "person_true_proportion", "population_mean_proportion", "error"
    )
  )
  expect_true(is.integer(dm$person_results$count) || is.numeric(dm$person_results$count))
  expect_setequal(unique(dm$person_results$metric), metrics)
  expect_setequal(unique(dm$person_results$replicate), seq_len(B))
  expect_setequal(unique(dm$person_results$person_id), seq_len(n_people))
  expect_setequal(unique(dm$person_results$cell_type), seq_len(K))
})


# Stream invariance: replicate b depends only on (seed, b), never on B ------------------------------------------

test_that("multinomial replicate streams do not depend on B (first 10 of B=20 match a B=10 run)", {
  r20 <- run_replicates(p = test_p, n = 200, B = 20, metrics = c("AE", "ARE"), model = "multinomial", seed = 99)
  r10 <- run_replicates(p = test_p, n = 200, B = 10, metrics = c("AE", "ARE"), model = "multinomial", seed = 99)

  expect_identical(r20$max_errors[1:10, , drop = FALSE], r10$max_errors)
  expect_identical(r20$phat[1:10, , drop = FALSE], r10$phat)
  expect_identical(r20$errors[1:10, , , drop = FALSE], r10$errors)
})

test_that("dirichlet_multinomial replicate streams do not depend on B (first 4 of B=8 match a B=4 run)", {
  dm8 <- run_replicates(
    p = test_p, B = 8, metrics = c("AE", "ARE"), model = "dirichlet_multinomial",
    n_people = 4, n_per_person = 50, concentration = 10, seed = 77
  )
  dm4 <- run_replicates(
    p = test_p, B = 4, metrics = c("AE", "ARE"), model = "dirichlet_multinomial",
    n_people = 4, n_per_person = 50, concentration = 10, seed = 77
  )

  sub <- dm8$person_results[dm8$person_results$replicate <= 4, , drop = FALSE]
  rownames(sub) <- NULL
  rownames(dm4$person_results) <- NULL

  expect_identical(sub, dm4$person_results)
})


# Global RNG state is left untouched -----------------------------------------------------------------------------

test_that("run_replicates (multinomial) leaves the global RNG kind and state unchanged", {
  withr::local_preserve_seed()
  RNGkind("Mersenne-Twister")
  set.seed(42)
  kind_before <- RNGkind()
  seed_before <- .Random.seed

  invisible(run_replicates(p = test_p, n = 200, B = 5, metrics = "AE", model = "multinomial", seed = 7))

  expect_identical(RNGkind(), kind_before)
  expect_identical(.Random.seed, seed_before)
})

test_that("run_replicates (dirichlet_multinomial) leaves the global RNG kind and state unchanged", {
  withr::local_preserve_seed()
  RNGkind("Mersenne-Twister")
  set.seed(42)
  kind_before <- RNGkind()
  seed_before <- .Random.seed

  invisible(run_replicates(
    p = test_p, B = 3, metrics = "AE", model = "dirichlet_multinomial",
    n_people = 3, n_per_person = 20, concentration = 5, seed = 7
  ))

  expect_identical(RNGkind(), kind_before)
  expect_identical(.Random.seed, seed_before)
})

test_that("repeated run_replicates() calls with the same seed each leave the RNG state unchanged", {
  # Mirrors how a sample-size solver calls run_replicates() at several n with the same seed
  # (common random numbers): the caller's RNG stream must be identical before and after any number
  # of such calls.
  withr::local_preserve_seed()
  RNGkind("Mersenne-Twister")
  set.seed(1)
  seed_before <- .Random.seed

  invisible(run_replicates(p = test_p, n = 100, B = 5, metrics = "AE", model = "multinomial", seed = 3))
  invisible(run_replicates(p = test_p, n = 150, B = 5, metrics = "AE", model = "multinomial", seed = 3))
  invisible(run_replicates(p = test_p, n = 200, B = 5, metrics = "AE", model = "multinomial", seed = 3))

  expect_identical(.Random.seed, seed_before)
})


# replicate_cores() ------------------------------------------------------------------------------------------------

test_that("replicate_cores() returns 1L on non-unix platforms", {
  old_platform <- .Platform
  withr::defer(assign(".Platform", old_platform, envir = globalenv()))
  assign(".Platform", utils::modifyList(.Platform, list(OS.type = "windows")), envir = globalenv())

  expect_identical(replicate_cores(), 1L)
})

test_that("replicate_cores() uses detectCores() - 1L on unix platforms", {
  old_platform <- .Platform
  withr::defer(assign(".Platform", old_platform, envir = globalenv()))
  assign(".Platform", utils::modifyList(.Platform, list(OS.type = "unix")), envir = globalenv())

  testthat::local_mocked_bindings(detectCores = function(...) 5L, .package = "parallel")
  expect_identical(replicate_cores(), 4L)
})

test_that("replicate_cores() falls back to 1L on unix when detectCores() is NA or < 1", {
  old_platform <- .Platform
  withr::defer(assign(".Platform", old_platform, envir = globalenv()))
  assign(".Platform", utils::modifyList(.Platform, list(OS.type = "unix")), envir = globalenv())

  testthat::local_mocked_bindings(detectCores = function(...) NA_integer_, .package = "parallel")
  expect_identical(replicate_cores(), 1L)

  testthat::local_mocked_bindings(detectCores = function(...) 1L, .package = "parallel")
  expect_identical(replicate_cores(), 1L)
})


# replicate_streams() / replicate_apply() ----------------------------------------------------------------------

test_that("replicate_streams() returns B streams that chain via nextRNGStream() and do not touch caller RNG", {
  withr::local_preserve_seed()
  RNGkind("Mersenne-Twister")
  set.seed(2)
  kind_before <- RNGkind()
  seed_before <- .Random.seed

  streams <- replicate_streams(seed = 10, B = 4)

  expect_length(streams, 4L)
  expect_identical(RNGkind(), kind_before)
  expect_identical(.Random.seed, seed_before)

  # Rebuilding with the same seed reproduces the same streams exactly.
  streams2 <- replicate_streams(seed = 10, B = 4)
  expect_identical(streams, streams2)

  # Chain identity: stream b+1 is nextRNGStream(stream b).
  expect_identical(streams[[2]], parallel::nextRNGStream(streams[[1]]))
  expect_identical(streams[[3]], parallel::nextRNGStream(streams[[2]]))
})

test_that("replicate_streams() draws a working seed when seed is NULL", {
  withr::local_preserve_seed()
  RNGkind("Mersenne-Twister")
  set.seed(5)
  streams <- replicate_streams(seed = NULL, B = 3)
  expect_length(streams, 3L)
  expect_true(all(vapply(streams, is.numeric, logical(1))))
})

test_that("replicate_streams(seed = NULL) advances the caller's RNG instead of rewinding it", {
  # A NULL seed is drawn from the caller's current RNG state via sample.int() before that state is
  # saved/restored, so -- exactly like any other draw -- the caller's RNG must advance permanently.
  # If it didn't (e.g. the seed were drawn after saving/restoring), consecutive unseeded calls would
  # silently draw the *same* seed and return identical streams every time.
  withr::local_preserve_seed()
  RNGkind("Mersenne-Twister")
  set.seed(123)
  seed_before <- .Random.seed

  streams1 <- replicate_streams(seed = NULL, B = 3)
  seed_after_first <- .Random.seed
  streams2 <- replicate_streams(seed = NULL, B = 3)

  expect_false(identical(seed_before, seed_after_first))
  expect_false(identical(streams1, streams2))
})

test_that("two consecutive run_replicates(seed = NULL) calls give different results", {
  withr::local_preserve_seed()
  RNGkind("Mersenne-Twister")
  set.seed(321)

  out1 <- run_replicates(p = test_p, n = 200, B = 10, metrics = "AE", model = "multinomial", seed = NULL)
  out2 <- run_replicates(p = test_p, n = 200, B = 10, metrics = "AE", model = "multinomial", seed = NULL)

  expect_false(identical(out1$max_errors, out2$max_errors))
  expect_false(identical(out1$phat, out2$phat))
})

test_that("replicate_apply() calls FUN once per stream, in order, using each installed stream", {
  # Manually installing a stream below mutates the global RNG kind/seed as a side effect, so save
  # and restore around it (replicate_apply() itself already restores the caller's state internally,
  # but the manual comparison loop below does not).
  withr::local_preserve_seed()

  streams <- replicate_streams(seed = 11, B = 5)
  out <- replicate_apply(streams, function(b) b)
  expect_identical(out, as.list(1:5))

  # Each replicate draws from its own stream: installing the same stream directly and drawing
  # reproduces what replicate_apply() produced for that replicate.
  draws <- replicate_apply(streams, function(b) runif(1))
  manual <- vapply(streams, function(s) {
    assign(".Random.seed", s, envir = globalenv())
    runif(1)
  }, numeric(1))
  expect_equal(unlist(draws), manual)
})


# check_replicate_results() ------------------------------------------------------------------------------------
#
# parallel::mclapply() only runs on unix, but check_replicate_results() is plain list processing with
# no OS dependency, so it is exercised directly here (including on Windows) against hand-built lists
# shaped like what mclapply() would return on worker failure.

test_that("check_replicate_results() passes through a list of successful results unchanged", {
  results <- list(1, 2, 3)
  expect_identical(check_replicate_results(results), results)
})

test_that("check_replicate_results() errors on a NULL element (killed worker)", {
  results <- list(1, NULL, 3)
  expect_error(check_replicate_results(results), "replicate 2", fixed = TRUE)
})

test_that("check_replicate_results() errors on a try-error element (worker error) with its message", {
  failure <- try(stop("boom"), silent = TRUE)
  results <- list(1, failure, 3)
  expect_error(check_replicate_results(results), "replicate 2", fixed = TRUE)
  expect_error(check_replicate_results(results), "boom", fixed = TRUE)
})

test_that("check_replicate_results() reports the first failure when several replicates fail", {
  failure1 <- try(stop("first failure"), silent = TRUE)
  failure2 <- try(stop("second failure"), silent = TRUE)
  results <- list(1, failure1, failure2)
  expect_error(check_replicate_results(results), "replicate 2", fixed = TRUE)
  expect_error(check_replicate_results(results), "first failure", fixed = TRUE)
})
