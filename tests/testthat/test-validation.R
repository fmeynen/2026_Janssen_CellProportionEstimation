# Tests for validate_required_columns() in scripts/simulation_layers/validation_utils.R

test_that("validate_required_columns passes silently when all columns are present", {
  df <- data.frame(a = 1, b = 2, c = 3)
  expect_silent(validate_required_columns(df, c("a", "b")))
  expect_identical(validate_required_columns(df, c("a", "b")), df)
})

test_that("validate_required_columns errors listing the missing columns", {
  df <- data.frame(a = 1)
  expect_error(validate_required_columns(df, c("a", "b", "c")), "missing: b, c", fixed = TRUE)
})

test_that("validate_required_columns errors for non-data.frames and uses the arg name", {
  expect_error(validate_required_columns(list(a = 1), "a", arg = "x"), "x must be a data.frame", fixed = TRUE)
})

test_that("validate_p_max accepts values strictly between 0 and 1 and returns them invisibly", {
  expect_silent(validate_p_max(0.4))
  expect_identical(validate_p_max(c(0.2, 0.5)), c(0.2, 0.5))
})

test_that("validate_p_max errors when p_max is NULL, naming the method argument", {
  expect_error(validate_p_max(NULL), "p_max must be provided when proportion_method = 'fixed_max_beta'", fixed = TRUE)
  expect_error(validate_p_max(NULL, method_arg = "method"), "when method = 'fixed_max_beta'", fixed = TRUE)
})

test_that("validate_p_max rejects non-numeric, non-finite and out-of-range values", {
  for (bad in list("a", NA_real_, Inf, 0, 1, -0.1, 1.5, c(0.3, 1))) {
    expect_error(validate_p_max(bad), "strictly between 0 and 1", fixed = TRUE)
  }
})
