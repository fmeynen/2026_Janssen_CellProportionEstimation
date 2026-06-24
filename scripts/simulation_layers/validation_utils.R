# R/validation_utils.R
#
# Validation helpers and shared constants used across all simulation layers.
# ---------------------------------------------------------------------------


#' Rescale nonnegative weights to the probability simplex.
#'
#' Internal helper used by proportion-generation methods.
normalize_to_simplex <- function(w) {
  stopifnot(is.numeric(w), all(is.finite(w)), all(w >= 0), sum(w) > 0)
  w / sum(w)
}

#' Default evaluation grid for Beta-based proportion generators.
#'
#' Returns a shared interior grid of the requested length. The fixed-max Beta
#' generator reuses this helper with `K - 1` points for its non-max remainder.
default_beta_grid <- function(K) {
  seq(0.05, 0.95, length.out = K)
}

#' Validate that `p` is a proper proportion vector.
#'
#' Checks: numeric, strictly positive, finite, sums to 1 within `tol`.
validate_proportions <- function(p, tol = 1e-12) {
  stopifnot(
    is.numeric(p),
    all(is.finite(p)),
    all(p > 0),
    abs(sum(p) - 1) < tol
  )
  invisible(p)
}

fixed_max_beta_impossible_error <- paste(
  "method = 'fixed_max_beta' failed because the fixed largest proportion is not strictly unique."
)

#' Warn and fail when fixed-max Beta proportions cannot keep a unique maximum.
#'
#' An impossible combination is defined exactly as follows: after constructing
#' the Beta-shaped remainder and scaling it to sum to `1 - p_max`, at least one
#' non-max component is `>= p_max`, so the fixed largest component is no longer
#' strictly unique.
fail_fixed_max_beta_impossible <- function(non_max, alpha, K, p_max) {
  warning(
    sprintf(
      paste(
        "Impossible fixed_max_beta combination for alpha=%s, K=%s, p_max=%s:",
        "after scaling the Beta-shaped remainder to sum to 1 - p_max,",
        "at least one non-max component is >= p_max (max non-max = %s)."
      ),
      format(alpha, trim = TRUE),
      K,
      format(p_max, trim = TRUE),
      format(max(non_max), trim = TRUE)
    ),
    call. = FALSE
  )
  stop(structure(
    list(message = fixed_max_beta_impossible_error, call = NULL),
    class = c("impossible_fixed_max_error", "error", "condition")
  ))
}
