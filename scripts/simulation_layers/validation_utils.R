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


#' Validate and normalize isoband triad ranges.
#'
#' @param ranges Named list with numeric length-2 vectors:
#'   `tau_AE`, `tau_ARE`, and `cutoff`.
#'
#' @return Normalized list with elements `tau_AE`, `tau_ARE`, and `cutoff`.
validate_isoband_ranges <- function(ranges) {
  required <- c("tau_AE", "tau_ARE", "cutoff")
  if (!is.list(ranges) || !all(required %in% names(ranges))) {
    stop(
      sprintf("ranges must be a named list containing: %s", paste(required, collapse = ", ")),
      call. = FALSE
    )
  }
  
  normalize_pair <- function(x, name) {
    if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x))) {
      stop(sprintf("ranges$%s must be a numeric length-2 finite vector.", name), call. = FALSE)
    }
    x <- as.numeric(x)
    if (!(x[[1L]] < x[[2L]])) {
      stop(sprintf("ranges$%s must be strictly increasing.", name), call. = FALSE)
    }
    x
  }
  
  tau_AE <- normalize_pair(ranges$tau_AE, "tau_AE")
  tau_ARE <- normalize_pair(ranges$tau_ARE, "tau_ARE")
  cutoff <- normalize_pair(ranges$cutoff, "cutoff")
  
  if (!(tau_AE[[1L]] >= 0 && tau_AE[[2L]] < 1)) {
    stop("tau_AE range must satisfy 0 <= min < max < 1.", call. = FALSE)
  }
  if (!(tau_ARE[[1L]] > 0 && is.finite(tau_ARE[[2L]]) && tau_ARE[[2L]] < Inf)) {
    stop("tau_ARE range must satisfy 0 < min < max < Inf with finite bounds.", call. = FALSE)
  }
  if (!(cutoff[[1L]] >= 0 && cutoff[[2L]] <= 1)) {
    stop("cutoff range must satisfy 0 <= min < max <= 1.", call. = FALSE)
  }
  
  list(tau_AE = tau_AE, tau_ARE = tau_ARE, cutoff = cutoff)
}