# Tests for replicate_success() and extract_success_rate() in scripts/simulation_layers/extraction.R


#' Build a minimal person_results data.frame.
#'
#' @param scenario_id Character vector (or NULL to omit the column, or NA to include an all-NA column).
#' @param replicate   Integer vector.
#' @param person_id   Integer vector.
#' @param cell_type   Integer vector.
#' @param metric      Character vector.
#' @param error       Numeric vector.
make_person_results <- function(scenario_id, replicate, person_id, cell_type, metric, error) {
  df <- data.frame(
    replicate = replicate,
    person_id = person_id,
    cell_type = cell_type,
    metric = metric,
    error = error,
    stringsAsFactors = FALSE
  )
  if (!is.null(scenario_id)) {
    df <- cbind(scenario_id = scenario_id, df, stringsAsFactors = FALSE)
  }
  df
}


test_that("mean over persons then max over cell types disagrees with a per-person rule", {
  # Two persons, one cell type, one replicate, metric AE.
  # Person 1 error = 0.9 (would fail alone at tau = 0.5), person 2 error = 0.1.
  # Mean over persons = 0.5, which passes at tau = 0.5, even though person 1 alone would not.
  pr <- make_person_results(
    scenario_id = "scenario_1",
    replicate = c(1L, 1L),
    person_id = c(1L, 2L),
    cell_type = c(1L, 1L),
    metric = c("AE", "AE"),
    error = c(0.9, 0.1)
  )
  out <- replicate_success(pr, list(AE = 0.5))
  expect_equal(nrow(out), 1L)
  expect_true(out$pass_AE)
  expect_true(out$pass)
})


test_that("max over cell types is taken after the per-cell-type mean over persons", {
  # Two cell types, two persons, one replicate, metric AE.
  # Cell type 1 mean = mean(0.9, 0.1) = 0.5 (passes tau = 0.6)
  # Cell type 2 mean = mean(0.1, 0.9) = 0.5 (passes tau = 0.6)
  # But cell type 2 alone for person 2 is 0.9, so a naive "max raw error <= tau" rule would fail.
  pr <- make_person_results(
    scenario_id = "scenario_1",
    replicate = c(1L, 1L, 1L, 1L),
    person_id = c(1L, 1L, 2L, 2L),
    cell_type = c(1L, 2L, 1L, 2L),
    metric = c("AE", "AE", "AE", "AE"),
    error = c(0.9, 0.1, 0.1, 0.9)
  )
  out <- replicate_success(pr, list(AE = 0.6))
  expect_equal(nrow(out), 1L)
  expect_equal(out$scenario_id, "scenario_1")
  expect_true(out$pass_AE)
  expect_true(out$pass)

  # Tightening tau below the per-cell-type mean (0.5) should fail.
  out2 <- replicate_success(pr, list(AE = 0.4))
  expect_false(out2$pass_AE)
  expect_false(out2$pass)
})


test_that("joint pass requires success on every metric", {
  pr <- make_person_results(
    scenario_id = c("s1", "s1"),
    replicate = c(1L, 1L),
    person_id = c(1L, 1L),
    cell_type = c(1L, 1L),
    metric = c("AE", "ARE"),
    error = c(0.1, 0.9)
  )
  out <- replicate_success(pr, list(AE = 0.5, ARE = 0.5))
  expect_true(out$pass_AE)
  expect_false(out$pass_ARE)
  expect_false(out$pass)

  out2 <- replicate_success(pr, list(AE = 0.5, ARE = 1))
  expect_true(out2$pass_AE)
  expect_true(out2$pass_ARE)
  expect_true(out2$pass)
})


test_that("NaN errors are treated as 0", {
  pr <- make_person_results(
    scenario_id = "s1",
    replicate = c(1L, 1L),
    person_id = c(1L, 2L),
    cell_type = c(1L, 1L),
    metric = c("ARE", "ARE"),
    error = c(NaN, NaN)
  )
  out <- replicate_success(pr, list(ARE = 0))
  expect_true(out$pass_ARE)
  expect_true(out$pass)
})


test_that("a single-scenario input with scenario_id NA or absent returns B rows", {
  # scenario_id column entirely NA.
  pr_na <- make_person_results(
    scenario_id = rep(NA_character_, 4L),
    replicate = c(1L, 1L, 2L, 2L),
    person_id = c(1L, 2L, 1L, 2L),
    cell_type = c(1L, 1L, 1L, 1L),
    metric = rep("AE", 4L),
    error = c(0.1, 0.2, 0.3, 0.4)
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
  pr <- make_person_results(
    scenario_id = c("scenario_2", "scenario_2", "scenario_1", "scenario_1"),
    replicate = c(1L, 1L, 1L, 1L),
    person_id = c(1L, 1L, 1L, 1L),
    cell_type = c(1L, 1L, 1L, 1L),
    metric = rep("AE", 4L),
    error = c(0.9, 0.9, 0.1, 0.1)
  )
  out <- replicate_success(pr, list(AE = 0.5))
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
    error = 0.1
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
  person_results <- make_person_results(
    scenario_id = c("scenario_1", "scenario_1", "scenario_1", "scenario_1"),
    replicate = c(1L, 1L, 2L, 2L),
    person_id = c(1L, 2L, 1L, 2L),
    cell_type = c(1L, 1L, 1L, 1L),
    metric = rep("AE", 4L),
    error = c(0.1, 0.1, 0.9, 0.9)
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

  out <- extract_success_rate(result, list(AE = 0.5))

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
