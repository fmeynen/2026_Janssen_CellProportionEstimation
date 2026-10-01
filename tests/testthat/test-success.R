# Tests for pooled_error_stat(), replicate_pooled_error() and replicate_success() in
# scripts/simulation_layers/calculation.R, and extract_success_rate() in scripts/simulation_layers/extraction.R


#' Build a minimal person_results data.frame.
#'
#' @param scenario_id Character vector (or NULL to omit the column, or NA to include an all-NA column).
#' @param replicate   Integer vector.
#' @param person_id   Integer vector.
#' @param cell_type   Integer vector.
#' @param metric      Character vector.
#' @param observed    Numeric vector; observed (estimated) proportion per row.
#' @param population  Numeric vector; population-level true proportion per row.
make_person_results <- function(scenario_id, replicate, person_id, cell_type, metric, observed, population) {
  df <- data.frame(
    replicate = replicate,
    person_id = person_id,
    cell_type = cell_type,
    metric = metric,
    observed_proportion = observed,
    population_mean_proportion = population,
    stringsAsFactors = FALSE
  )
  if (!is.null(scenario_id)) {
    df <- cbind(scenario_id = scenario_id, df, stringsAsFactors = FALSE)
  }
  df
}


test_that("estimates are pooled over persons before comparing with the population proportion", {
  # Two persons, one cell type, p = 0.5. Person errors are 0.4 each, but the pooled estimate mean(0.9, 0.1) = 0.5
  # hits p exactly, so AE = ARE = 0.
  pr <- make_person_results(
    scenario_id = "scenario_1",
    replicate = c(1L, 1L),
    person_id = c(1L, 2L),
    cell_type = c(1L, 1L),
    metric = c("AE", "AE"),
    observed = c(0.9, 0.1),
    population = c(0.5, 0.5)
  )
  stat <- replicate_pooled_error(pr, "AE")
  expect_equal(stat$stat, 0)
  out <- replicate_success(pr, list(AE = 0))
  expect_equal(nrow(out), 1L)
  expect_true(out$pass_AE)
  expect_true(out$pass)
})


test_that("ARE divides the pooled absolute error by the population proportion", {
  # One cell type, p = 0.2, pooled estimate mean(0.3, 0.2) = 0.25: AE = 0.05, ARE = 0.25.
  pr <- make_person_results(
    scenario_id = "s1",
    replicate = c(1L, 1L),
    person_id = c(1L, 2L),
    cell_type = c(1L, 1L),
    metric = c("ARE", "ARE"),
    observed = c(0.3, 0.2),
    population = c(0.2, 0.2)
  )
  expect_equal(replicate_pooled_error(pr, "ARE")$stat, 0.25)
})


test_that("max over cell types is taken after pooling over persons", {
  # Two cell types (p = 0.5, 0.5), two persons. Pooled estimates: (0.4, 0.6) -> AE = (0.1, 0.1); max = 0.1.
  pr <- make_person_results(
    scenario_id = "scenario_1",
    replicate = c(1L, 1L, 1L, 1L),
    person_id = c(1L, 1L, 2L, 2L),
    cell_type = c(1L, 2L, 1L, 2L),
    metric = rep("AE", 4L),
    observed = c(0.6, 0.4, 0.2, 0.8),
    population = rep(0.5, 4L)
  )
  expect_equal(replicate_pooled_error(pr, "AE")$stat, 0.1)

  out <- replicate_success(pr, list(AE = 0.15))
  expect_equal(out$scenario_id, "scenario_1")
  expect_true(out$pass_AE)

  out2 <- replicate_success(pr, list(AE = 0.05))
  expect_false(out2$pass_AE)
  expect_false(out2$pass)
})


test_that("joint pass requires success on every metric", {
  # One person, one cell type, p = 0.2, observed 0.3: AE = 0.1, ARE = 0.5.
  pr <- make_person_results(
    scenario_id = c("s1", "s1"),
    replicate = c(1L, 1L),
    person_id = c(1L, 1L),
    cell_type = c(1L, 1L),
    metric = c("AE", "ARE"),
    observed = c(0.3, 0.3),
    population = c(0.2, 0.2)
  )
  out <- replicate_success(pr, list(AE = 0.2, ARE = 0.4))
  expect_true(out$pass_AE)
  expect_false(out$pass_ARE)
  expect_false(out$pass)

  out2 <- replicate_success(pr, list(AE = 0.2, ARE = 1))
  expect_true(out2$pass_AE)
  expect_true(out2$pass_ARE)
  expect_true(out2$pass)
})


test_that("ARE with pooled estimate and population proportion both 0 is treated as 0", {
  pr <- make_person_results(
    scenario_id = "s1",
    replicate = c(1L, 1L),
    person_id = c(1L, 2L),
    cell_type = c(1L, 1L),
    metric = c("ARE", "ARE"),
    observed = c(0, 0),
    population = c(0, 0)
  )
  out <- replicate_success(pr, list(ARE = 0))
  expect_true(out$pass_ARE)
  expect_true(out$pass)
})


test_that("metrics other than AE and ARE are rejected", {
  pr <- make_person_results(
    scenario_id = "s1",
    replicate = 1L,
    person_id = 1L,
    cell_type = 1L,
    metric = "TSE",
    observed = 0.3,
    population = 0.2
  )
  expect_error(replicate_pooled_error(pr, "TSE"), "AE and ARE only")
})


test_that("a single-scenario input with scenario_id NA or absent returns B rows", {
  # scenario_id column entirely NA.
  pr_na <- make_person_results(
    scenario_id = rep(NA_character_, 4L),
    replicate = c(1L, 1L, 2L, 2L),
    person_id = c(1L, 2L, 1L, 2L),
    cell_type = c(1L, 1L, 1L, 1L),
    metric = rep("AE", 4L),
    observed = c(0.1, 0.2, 0.3, 0.4),
    population = rep(0.25, 4L)
  )
  out_na <- replicate_success(pr_na, list(AE = 1))
  expect_equal(nrow(out_na), 2L)
  expect_true(all(is.na(out_na$scenario_id)))
  expect_equal(out_na$replicate, c(1L, 2L))

  # scenario_id column entirely absent.
  pr_absent <- pr_na[, setdiff(names(pr_na), "scenario_id")]
  out_absent <- replicate_success(pr_absent, list(AE = 1))
  expect_equal(nrow(out_absent), 2L)
  expect_true(all(is.na(out_absent$scenario_id)))
  expect_equal(out_absent$replicate, c(1L, 2L))
})


test_that("multiple scenarios are handled independently and sorted by scenario_id then replicate", {
  # scenario_2: pooled 0.9 vs p 0.5 -> AE 0.4 (fails 0.2); scenario_1: pooled 0.5 vs p 0.5 -> AE 0 (passes).
  pr <- make_person_results(
    scenario_id = c("scenario_2", "scenario_2", "scenario_1", "scenario_1"),
    replicate = c(1L, 1L, 1L, 1L),
    person_id = c(1L, 2L, 1L, 2L),
    cell_type = c(1L, 1L, 1L, 1L),
    metric = rep("AE", 4L),
    observed = c(0.9, 0.9, 0.4, 0.6),
    population = rep(0.5, 4L)
  )
  out <- replicate_success(pr, list(AE = 0.2))
  expect_equal(out$scenario_id, c("scenario_1", "scenario_2"))
  expect_equal(out$pass, c(TRUE, FALSE))
})


test_that("a missing metric warns and is skipped", {
  pr <- make_person_results(
    scenario_id = "s1",
    replicate = 1L,
    person_id = 1L,
    cell_type = 1L,
    metric = "AE",
    observed = 0.3,
    population = 0.2
  )
  expect_warning(
    out <- replicate_success(pr, list(AE = 0.5, TSE = 0.5)),
    "TSE"
  )
  expect_false("pass_TSE" %in% names(out))
  expect_true("pass_AE" %in% names(out))
  expect_true(out$pass)

  expect_error(
    suppressWarnings(replicate_success(pr, list(TSE = 0.5))),
    "None of the metrics"
  )
})


test_that("extract_success_rate() output format is unchanged", {
  # Replicate 1: pooled 0.5 vs p 0.5 -> AE 0 (passes); replicate 2: pooled 0.9 -> AE 0.4 (fails).
  person_results <- make_person_results(
    scenario_id = c("scenario_1", "scenario_1", "scenario_1", "scenario_1"),
    replicate = c(1L, 1L, 2L, 2L),
    person_id = c(1L, 2L, 1L, 2L),
    cell_type = c(1L, 1L, 1L, 1L),
    metric = rep("AE", 4L),
    observed = c(0.4, 0.6, 0.9, 0.9),
    population = rep(0.5, 4L)
  )
  p_table <- data.frame(
    scenario_id = "scenario_1",
    alpha = 0.05,
    p_max = 0.3,
    n_people = 2L,
    concentration = 10,
    stringsAsFactors = FALSE
  )
  result <- list(person_results = person_results, p_table = p_table)

  out <- extract_success_rate(result, list(AE = 0.2))

  expect_s3_class(out, "data.frame")
  expect_equal(
    names(out),
    c("scenario_id", "alpha", "p_max", "n_people", "concentration", "B", "success_count", "success_rate",
      "success_rate_AE")
  )
  expect_equal(nrow(out), 1L)
  expect_equal(out$scenario_id, "scenario_1")
  expect_equal(out$B, 2L)
  expect_equal(out$success_count, 1L)
  expect_equal(out$success_rate, 0.5)
  expect_equal(out$success_rate_AE, 0.5)
})


test_that("pooled_error_stat() AE is the row-wise max of the absolute errors", {
  # Row 1: errors (0.1, 0.1, 0); row 2: errors (0.2, 0, 0.2) -> 0.1 and 0.2.
  pbar <- rbind(c(0.3, 0.4, 0.3), c(0.4, 0.5, 0.1))
  p <- c(0.2, 0.5, 0.3)
  expect_equal(pooled_error_stat(pbar, p, "AE"), c(0.1, 0.2))
})


test_that("pooled_error_stat() ARE divides each error by its true proportion", {
  # Row 1: relative errors (0.5, 0.2) -> 0.5; row 2: (0, 0.4) -> 0.4.
  pbar <- rbind(c(0.3, 0.4), c(0.2, 0.7))
  p <- c(0.2, 0.5)
  expect_equal(pooled_error_stat(pbar, p, "ARE"), c(0.5, 0.4))
})


test_that("pooled_error_stat() ARE counts 0/0 as 0 and keeps Inf", {
  # Cell type 1 has p = 0. Row 1: pbar_1 = 0 -> NaN -> 0, so the max is cell type 2's 0.25.
  # Row 2: pbar_1 > 0 -> Inf.
  pbar <- rbind(c(0, 0.75), c(0.1, 0.9))
  p <- c(0, 1)
  expect_equal(pooled_error_stat(pbar, p, "ARE"), c(0.25, Inf))
})


test_that("pooled_error_stat() handles a single replicate", {
  pbar <- matrix(c(0.25, 0.75), nrow = 1L)
  expect_equal(pooled_error_stat(pbar, c(0.5, 0.5), "AE"), 0.25)
  expect_equal(pooled_error_stat(pbar, c(0.5, 0.5), "ARE"), 0.5)
})


test_that("pooled_error_stat() rejects metrics other than AE and ARE", {
  pbar <- matrix(c(0.25, 0.75), nrow = 1L)
  expect_error(pooled_error_stat(pbar, c(0.5, 0.5), "TSE"), "AE and ARE only")
})
