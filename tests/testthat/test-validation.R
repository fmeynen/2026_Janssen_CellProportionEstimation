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

test_that("is_finite_scalar is TRUE only for a single finite number", {
  expect_true(is_finite_scalar(0.5))
  expect_true(is_finite_scalar(3L))
  for (bad in list(NA_real_, Inf, -Inf, NULL, "a", c(1, 2), numeric(0), TRUE)) {
    expect_false(is_finite_scalar(bad))
  }
})

test_that("validate_positive_integer and validate_positive_numeric coerce and reject", {
  expect_identical(validate_positive_integer(3, "B"), 3L)
  expect_identical(validate_positive_integer(c(1, 2), "K", allow_vector = TRUE), c(1L, 2L))
  expect_error(validate_positive_integer(1.5, "B"), "B must be a positive integer.", fixed = TRUE)
  expect_error(validate_positive_integer(c(1, 2), "B"), "B must be", fixed = TRUE)
  expect_identical(validate_positive_numeric(0.5, "a"), 0.5)
  expect_error(validate_positive_numeric(0, "a"), "a must be a single positive finite number.", fixed = TRUE)
  expect_error(validate_positive_numeric(numeric(0), "a", allow_vector = TRUE), "non-empty vector", fixed = TRUE)
})

test_that("validate_named_matrix requires a matrix with column names", {
  m <- matrix(1:4, 2, dimnames = list(NULL, c("AE", "ARE")))
  expect_identical(validate_named_matrix(m, "max_errors"), m)
  expect_error(validate_named_matrix(matrix(1:4, 2), "max_errors"), "max_errors must be a matrix with column names.", fixed = TRUE)
  expect_error(validate_named_matrix(1:4, "argmax"), "argmax must be a matrix with column names.", fixed = TRUE)
})

test_that("validate_result_fields requires a list containing each field", {
  res <- list(a = 1, b = 2)
  expect_identical(validate_result_fields(res, c("a", "b")), res)
  expect_error(validate_result_fields(1, "a"), "result must be a list.", fixed = TRUE)
  expect_error(validate_result_fields(res, c("a", "z")), "result must contain z.", fixed = TRUE)
})

test_that("validate_open_unit_scalar accepts (0, 1) only", {
  expect_silent(validate_open_unit_scalar(0.95, "target"))
  for (bad in list(0, 1, NA_real_, c(0.1, 0.2), "a", NULL)) {
    expect_error(validate_open_unit_scalar(bad, "target"), "target must be a single number in (0, 1).", fixed = TRUE)
  }
})

test_that("validate_n_max accepts [1, integer.max] only", {
  expect_silent(validate_n_max(1e9))
  for (bad in list(0, 3e9, Inf, NA_real_, c(1, 2))) {
    expect_error(validate_n_max(bad, "config$n_max"), "config$n_max must be a single number in [1, .Machine$integer.max].", fixed = TRUE)
  }
})

test_that("former stopifnot checks now raise explicit call.-free errors", {
  expect_error(normalize_to_simplex(c(-1, 2)), "w must be", fixed = TRUE)
  expect_error(validate_proportions(c(0.5, 0.4)), "p must sum to 1.", fixed = TRUE)
  expect_error(validate_proportions(c(0.5, 0.5, 0)), "strictly positive", fixed = TRUE)
  err <- tryCatch(sample_size_pilots(10, 1), error = identity)
  expect_null(conditionCall(err))
})
