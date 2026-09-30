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
