# R/validation_utils.R
#
# Validation helpers and shared constants used across all simulation layers.
# ---------------------------------------------------------------------------


#' Rescale nonnegative weights to the probability simplex.
#'
#' Internal helper used by proportion-generation methods.
normalize_to_simplex <- function(w) {
  if (!is.numeric(w) || any(!is.finite(w)) || any(w < 0) || !(sum(w) > 0)) {
    stop("w must be a numeric vector of finite nonnegative weights with a positive sum.", call. = FALSE)
  }
  w / sum(w)
}

#' Default evaluation grid for Beta-based proportion generators.
#'
#' Returns a shared interior grid of the requested length (K points for every Beta-based generator).
default_beta_grid <- function(K) {
  seq(0.05, 0.95, length.out = K)
}

#' Validate that `p` is a proper proportion vector.
#'
#' Checks: numeric, strictly positive, finite, sums to 1 within `tol`.
validate_proportions <- function(p, tol = 1e-12) {
  if (!is.numeric(p) || any(!is.finite(p)) || any(p <= 0)) {
    stop("p must be a numeric vector of finite, strictly positive values.", call. = FALSE)
  }
  if (!(abs(sum(p) - 1) < tol)) {
    stop("p must sum to 1.", call. = FALSE)
  }
  invisible(p)
}

fixed_max_beta_impossible_error <- paste(
  "method = 'fixed_max_beta' failed because K * p_max < 1, so no proportion vector can have p_max as its maximum."
)

#' Warn and fail when no K-vector can have `p_max` as its maximum.
#'
#' An impossible combination is defined exactly as follows: `K * p_max < 1`, so
#' even K equal proportions of `p_max` would sum to less than 1.
fail_fixed_max_beta_impossible <- function(alpha, K, p_max) {
  warning(
    sprintf(
      paste(
        "Impossible fixed_max_beta combination for alpha=%s, K=%s, p_max=%s:",
        "K * p_max = %s < 1, so p_max cannot be the largest of K proportions summing to 1."
      ),
      format(alpha, trim = TRUE),
      K,
      format(p_max, trim = TRUE),
      format(K * p_max, trim = TRUE)
    ),
    call. = FALSE
  )
  stop(structure(
    list(message = fixed_max_beta_impossible_error, call = NULL),
    class = c("impossible_fixed_max_error", "error", "condition")
  ))
}

fixed_min_beta_impossible_error <- paste(
  "method = 'fixed_min_beta' failed because K * p_min > 1, so no proportion vector can have p_min as its minimum."
)

#' Warn and fail when no K-vector can have `p_min` as its minimum.
#'
#' An impossible combination is defined exactly as follows: `K * p_min > 1`, so
#' even K equal proportions of `p_min` would sum to more than 1.
fail_fixed_min_beta_impossible <- function(alpha, K, p_min) {
  warning(
    sprintf(
      paste(
        "Impossible fixed_min_beta combination for alpha=%s, K=%s, p_min=%s:",
        "K * p_min = %s > 1, so p_min cannot be the smallest of K proportions summing to 1."
      ),
      format(alpha, trim = TRUE),
      K,
      format(p_min, trim = TRUE),
      format(K * p_min, trim = TRUE)
    ),
    call. = FALSE
  )
  stop(structure(
    list(message = fixed_min_beta_impossible_error, call = NULL),
    class = c("impossible_fixed_min_error", "error", "condition")
  ))
}

#' Columns every `person_results` data.frame must contain for extraction.
person_results_required_cols <- c(
  "replicate", "cell_type", "metric", "observed_proportion", "population_mean_proportion"
)

#' Validate that `df` is a data.frame containing all `required` columns.
#'
#' Stops (with `call. = FALSE`) naming the missing columns.
validate_required_columns <- function(df, required, arg = "person_results") {
  if (!is.data.frame(df)) {
    stop(sprintf("%s must be a data.frame.", arg), call. = FALSE)
  }
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "%s must be a data.frame with columns %s; missing: %s.",
        arg, paste(required, collapse = ", "), paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(df)
}

#' Validate `p_max` for `fixed_max_beta` proportions.
#'
#' Stops (with `call. = FALSE`) if `p_max` is NULL or is not a numeric vector of
#' finite values strictly between 0 and 1.
#'
#' @param p_max       Value to check.
#' @param method_arg  Name of the method argument, used in the error message.
validate_p_max <- function(p_max, method_arg = "proportion_method") {
  if (is.null(p_max)) {
    stop(
      sprintf("p_max must be provided when %s = 'fixed_max_beta'.", method_arg),
      call. = FALSE
    )
  }
  if (!is.numeric(p_max) || any(!is.finite(p_max)) || any(p_max <= 0) || any(p_max >= 1)) {
    stop("p_max must contain number(s) strictly between 0 and 1.", call. = FALSE)
  }
  invisible(p_max)
}

#' Validate a positive integer (or vector of them) and return it as integer.
#'
#' Stops (with `call. = FALSE`) unless `x` is numeric, finite, whole and >= 1.
#'
#' @param x            Value to check.
#' @param name         Argument name used in the error message.
#' @param allow_vector If `TRUE`, accept a non-empty vector; otherwise length 1.
validate_positive_integer <- function(x, name, allow_vector = FALSE) {
  valid_length <- if (allow_vector) length(x) >= 1L else length(x) == 1L
  if (!is.numeric(x) || !valid_length || any(!is.finite(x)) || any(x < 1L) || any(x %% 1 != 0)) {
    expected <- if (allow_vector) {
      "a non-empty vector of positive integers"
    } else {
      "a positive integer"
    }
    stop(sprintf("%s must be %s.", name, expected), call. = FALSE)
  }
  as.integer(x)
}

#' Validate a positive finite number (or vector of them) and return it as numeric.
#'
#' Stops (with `call. = FALSE`) unless `x` is numeric, finite and > 0.
#'
#' @param x            Value to check.
#' @param name         Argument name used in the error message.
#' @param allow_vector If `TRUE`, accept a non-empty vector; otherwise length 1.
validate_positive_numeric <- function(x, name, allow_vector = FALSE) {
  valid_length <- if (allow_vector) length(x) >= 1L else length(x) == 1L
  if (!is.numeric(x) || !valid_length || any(!is.finite(x)) || any(x <= 0)) {
    expected <- if (allow_vector) {
      "a non-empty vector of positive finite numbers"
    } else {
      "a single positive finite number"
    }
    stop(sprintf("%s must be %s.", name, expected), call. = FALSE)
  }
  as.numeric(x)
}

#' Is `x` a single finite number?
#'
#' TRUE for a length-1 numeric that is finite; FALSE otherwise (including NA, Inf, NULL, non-numeric).
is_finite_scalar <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

#' Validate that `x` is a matrix with column names (one column per metric).
#'
#' Stops (with `call. = FALSE`) otherwise.
#'
#' @param x    Value to check.
#' @param name Argument name used in the error message.
validate_named_matrix <- function(x, name) {
  if (!is.matrix(x) || is.null(colnames(x))) {
    stop(sprintf("%s must be a matrix with column names.", name), call. = FALSE)
  }
  invisible(x)
}

#' Validate that `result` is a list containing the named `fields`.
#'
#' Stops (with `call. = FALSE`) naming the first missing field.
#'
#' @param result Value to check.
#' @param fields Character vector of required element names.
validate_result_fields <- function(result, fields) {
  if (!is.list(result)) {
    stop("result must be a list.", call. = FALSE)
  }
  for (field in fields) {
    if (!(field %in% names(result))) {
      stop(sprintf("result must contain %s.", field), call. = FALSE)
    }
  }
  invisible(result)
}

#' Validate a single finite number strictly between 0 and 1 (e.g. a target success rate).
#'
#' Stops (with `call. = FALSE`) otherwise.
#'
#' @param x    Value to check.
#' @param name Argument name used in the error message.
validate_open_unit_scalar <- function(x, name) {
  if (!is_finite_scalar(x) || x <= 0 || x >= 1) {
    stop(sprintf("%s must be a single number in (0, 1).", name), call. = FALSE)
  }
  invisible(x)
}

#' Validate an upper cap `n_max`: a single finite number in [1, .Machine$integer.max].
#'
#' Stops (with `call. = FALSE`) otherwise.
#'
#' @param n_max Value to check.
#' @param name  Argument name used in the error message.
validate_n_max <- function(n_max, name = "n_max") {
  if (!is_finite_scalar(n_max) || n_max < 1 || n_max > .Machine$integer.max) {
    stop(sprintf("%s must be a single number in [1, .Machine$integer.max].", name), call. = FALSE)
  }
  invisible(n_max)
}
