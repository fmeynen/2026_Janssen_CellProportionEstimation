# Calculation Layer -----------------------------------------------------------------------------------------------
# Calculation layer: error metrics and threshold evaluation.
# NOTE: the old iterative-GLM sample-size helpers (fit_success_glm,
# solve_sample_size_from_glm, iterate_sample_size_for_alpha) have been moved to
# scripts/deprecated/deprecated.R and are being replaced by a new solver below.


# Compute Success Rate --------------------------------------------------------------------------------------------

#' Compute success rates for each metric across a grid of thresholds.
#'
#' Runs post-hoc on stored max errors; no re-simulation needed.
#'
#' @param max_errors B x M matrix of max error values (from `run_replicates()`).
#' @param taus       Either a numeric vector of threshold values applied to every metric, or a named list with one
#'    numeric vector per metric (e.g. `list(AE = c(...), ARE = c(...))`).  When a list is supplied every metric present
#'    in `max_errors` must have an entry.
#' @param errors     Optional B x K x M array of per-cell-type errors
#'   (from `run_replicates()`), used to compute `mean_n_above`.
#'
#' @return Tidy data.frame with columns: metric, tau, success_rate, mean_n_above.
evaluate_thresholds <- function(max_errors, taus, errors = NULL) {
  stopifnot(is.matrix(max_errors), !is.null(colnames(max_errors)))
  metrics <- colnames(max_errors)
  has_error_array <- !is.null(errors)
  if (has_error_array) {
    stopifnot(length(dim(errors)) == 3L)
    stopifnot(dim(errors)[1] == nrow(max_errors))
    stopifnot(dim(errors)[3] == length(metrics))
  }

  # Normalise taus: plain vector -> same grid for all metrics;
  # named list     -> per-metric grids.
  if (is.numeric(taus)) {
    taus_list <- setNames(rep(list(taus), length(metrics)), metrics)
  } else if (is.list(taus)) {
    missing_metrics <- setdiff(metrics, names(taus))
    if (length(missing_metrics) > 0L) {
      stop(sprintf(
        "`taus` list is missing entries for metric(s): %s",
        paste(missing_metrics, collapse = ", ")
      ))
    }
    taus_list <- taus[metrics]
  } else {
    stop("`taus` must be a numeric vector or a named list of numeric vectors.")
  }

  rows <- vector("list", length(metrics))
  for (i in seq_along(metrics)) {
    m     <- metrics[[i]]
    tau_m <- taus_list[[m]]
    rates <- vapply(tau_m, function(tau) mean(max_errors[, m] <= tau, na.rm = TRUE), numeric(1L))
    if (has_error_array) {
      errors_m <- errors[, , m, drop = TRUE]
      if (is.null(dim(errors_m))) errors_m <- matrix(errors_m, nrow = nrow(max_errors))
      mean_n_above <- vapply(
        tau_m,
        function(tau) mean(rowSums(errors_m > tau), na.rm = TRUE),
        numeric(1L)
      )
    } else {
      mean_n_above <- rep(NA_real_, length(tau_m))
    }
    rows[[i]] <- data.frame(metric = m, tau = tau_m, success_rate = rates,
                            mean_n_above = mean_n_above,
                            stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}


# Sample-Size Estimation --------------------------------------------------------------------------------------

#' Pilot sample sizes around a centre on a multiplicative grid.
#'
#' @param n Numeric scalar, the current centre (> 0).
#' @param f Numeric scalar spread factor (> 1). Pilots are placed at `n / f`, `n` and `n * f`.
#'
#' @return Sorted, unique integer vector of pilot sizes, each rounded up and at least 1.
sample_size_pilots <- function(n, f) {
  stopifnot(is.numeric(n), length(n) == 1L, is.finite(n), n > 0,
            is.numeric(f), length(f) == 1L, is.finite(f), f > 1)
  sort(unique(pmax(1L, as.integer(ceiling(c(n / f, n, n * f))))))
}


#' Fit a logistic success curve on log(n).
#'
#' Fits `glm(cbind(s, B - s) ~ log(n), family = binomial)`. Only the two warnings caused by (quasi-)separation are
#' suppressed: "fitted probabilities numerically 0 or 1" and "algorithm did not converge". The latter is routine when
#' a single pilot is interior and all others are saturated (e.g. just after an `expand_down`); the resulting steep fit
#' still gives a usable next centre. Any other warning is propagated.
#'
#' @param n_values      Numeric vector of pilot sample sizes (>= 1).
#' @param success_count Integer vector of successful replicates per pilot (0..B).
#' @param B             Number of replicates per pilot.
#'
#' @return The fitted `glm` object.
fit_success_curve <- function(n_values, success_count, B) {
  stopifnot(is.numeric(n_values), length(n_values) >= 1L, all(n_values >= 1),
            is.numeric(success_count), length(success_count) == length(n_values),
            is.numeric(B), length(B) == 1L, B >= 1,
            all(success_count >= 0), all(success_count <= B))
  dat <- data.frame(n = n_values, s = success_count, fails = B - success_count)
  withCallingHandlers(
    stats::glm(cbind(s, fails) ~ log(n), family = stats::binomial(), data = dat),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("fitted probabilities numerically 0 or 1", msg, fixed = TRUE) ||
          grepl("algorithm did not converge", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}


#' Invert a fitted success curve at a target success rate.
#'
#' Solves `qlogis(target) = intercept + slope * log(n)` for n.
#'
#' @param fit    A `glm` from `fit_success_curve()`; its slope must be finite and positive.
#' @param target Target success rate in (0, 1).
#'
#' @return Integer n, rounded up, at least 1.
solve_success_curve <- function(fit, target) {
  stopifnot(is.numeric(target), length(target) == 1L, target > 0, target < 1)
  coefs     <- stats::coef(fit)
  intercept <- unname(coefs[[1L]])
  slope     <- unname(coefs[[2L]])
  if (!is.finite(intercept) || !is.finite(slope) || slope <= 0) {
    stop("Cannot invert the success curve: the slope must be finite and positive.", call. = FALSE)
  }
  n_raw <- exp((stats::qlogis(target) - intercept) / slope)
  if (!is.finite(n_raw) || n_raw > .Machine$integer.max) {
    stop(sprintf("Solved sample size is not representable as an integer (n = %g).", n_raw), call. = FALSE)
  }
  max(1L, as.integer(ceiling(n_raw)))
}


#' Estimate the smallest sample size reaching a target success rate for one alpha.
#'
#' Iterative solver. Each iteration simulates pilots at `n / f`, `n` and `n * f` (see `sample_size_pilots()`) with
#' seed `config$seed + seed_offset` and accumulates them. The step type is chosen in this order:
#' \enumerate{
#'   \item No accumulated pilot has 0 < success_count < B (degenerate):
#'     `expand_up` (all of this iteration's pilots failed; centre <- ceiling(f^2 * max(pilots))),
#'     `expand_down` (all succeeded; centre <- max(1, ceiling(min(pilots) / f^2))), or
#'     `bisect` (the accumulated data brackets the answer; centre <- ceiling(sqrt(largest all-fail n *
#'     smallest all-success n))). `f` is unchanged.
#'   \item Otherwise fit `fit_success_curve()` on all accumulated pilots. A non-finite or non-positive slope gives
#'     `resample`: same centre, `seed_offset + 1`, `f` unchanged.
#'   \item Otherwise `fit`: centre <- `solve_success_curve()`; converged when `|new_n - n| / n <= rel_tol`; then
#'     `f <- max(f_floor, sqrt(f))`. Convergence is only checked on `fit` steps.
#' }
#' No clamping is applied to any step.
#'
#' @param alpha    Numeric scalar passed through to `simulate`.
#' @param n_init   Starting centre (>= 1); rounded up.
#' @param config   List with `success_rate_target`, `rel_tol`, `max_iterations`, `B`, `f0`, `f_floor`, `seed`.
#' @param simulate Function `(alpha, n, config, seed)` returning a list with at least `success_count`.
#'   Defaults to `simulate_success_at_n()`; injectable for testing.
#'
#' @return List with `final_n` (integer), `stopping_reason` ("tolerance" or "max_iterations"), `iterations_used`,
#'   and `diagnostics`: a data.frame with one row per pilot per iteration and columns alpha, iteration, step, n,
#'   success_count, success_rate, f, seed_offset, glm_intercept, glm_slope, n_next, stopping_reason (NA except in
#'   the final row). `f` and `seed_offset` are the values used for that iteration's pilots. Warns when stopping at
#'   max_iterations after at least one `fit` step; errors when no `fit` step ever happened.
estimate_sample_size <- function(alpha, n_init, config, simulate = simulate_success_at_n) {
  stopifnot(is.numeric(n_init), length(n_init) == 1L, is.finite(n_init), n_init >= 1, is.function(simulate))
  target    <- config$success_rate_target
  rel_tol   <- config$rel_tol
  max_iter  <- config$max_iterations
  B         <- config$B
  f0        <- config$f0
  f_floor   <- config$f_floor
  base_seed <- config$seed
  is_scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!is_scalar(target) || target <= 0 || target >= 1) {
    stop("`config$success_rate_target` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is_scalar(rel_tol) || rel_tol < 0) stop("`config$rel_tol` must be a single number >= 0.", call. = FALSE)
  if (!is_scalar(max_iter) || max_iter < 1) {
    stop("`config$max_iterations` must be a single number >= 1.", call. = FALSE)
  }
  if (!is_scalar(B) || B < 1) stop("`config$B` must be a single number >= 1.", call. = FALSE)
  if (!is_scalar(f0) || f0 <= 1) stop("`config$f0` must be a single number > 1.", call. = FALSE)
  if (!is_scalar(f_floor) || f_floor <= 1 || f_floor > f0) {
    stop("`config$f_floor` must be a single number with 1 < f_floor <= f0.", call. = FALSE)
  }
  if (!is_scalar(base_seed)) stop("`config$seed` must be a single finite number.", call. = FALSE)
  max_iter <- as.integer(max_iter)

  n           <- as.integer(ceiling(n_init))
  f           <- f0
  seed_offset <- 0L
  acc_n       <- integer(0L)
  acc_s       <- integer(0L)
  diag_rows   <- vector("list", max_iter)
  last_fit_n  <- NA_integer_
  converged   <- FALSE
  iter        <- 0L

  while (iter < max_iter && !converged) {
    iter   <- iter + 1L
    pilots <- sample_size_pilots(n, f)
    seed   <- base_seed + seed_offset
    s <- vapply(pilots, function(p) {
      as.integer(simulate(alpha = alpha, n = p, config = config, seed = seed)$success_count)
    }, integer(1L))
    if (anyNA(s) || any(s < 0L) || any(s > B)) {
      stop(sprintf("`simulate` returned a success_count outside [0, B] at iteration %d.", iter), call. = FALSE)
    }
    acc_n <- c(acc_n, pilots)
    acc_s <- c(acc_s, s)

    intercept   <- NA_real_
    slope       <- NA_real_
    f_used      <- f
    offset_used <- seed_offset

    if (!any(acc_s > 0L & acc_s < B)) {
      if (all(s == 0L)) {
        step   <- "expand_up"
        n_next <- as.integer(ceiling(f^2 * max(pilots)))
      } else if (all(s == B)) {
        step   <- "expand_down"
        n_next <- max(1L, as.integer(ceiling(min(pilots) / f^2)))
      } else {
        step   <- "bisect"
        n_fail <- max(acc_n[acc_s == 0L])
        n_pass <- min(acc_n[acc_s == B])
        n_next <- as.integer(ceiling(sqrt(n_fail * n_pass)))
      }
    } else {
      fit       <- fit_success_curve(acc_n, acc_s, B)
      coefs     <- stats::coef(fit)
      intercept <- unname(coefs[[1L]])
      slope     <- unname(coefs[[2L]])
      if (!is.finite(intercept) || !is.finite(slope) || slope <= 0) {
        step        <- "resample"
        n_next      <- n
        seed_offset <- seed_offset + 1L
      } else {
        step       <- "fit"
        n_next     <- solve_success_curve(fit, target)
        converged  <- abs(n_next - n) / n <= rel_tol
        f          <- max(f_floor, sqrt(f))
        last_fit_n <- n_next
      }
    }

    diag_rows[[iter]] <- data.frame(
      alpha            = alpha,
      iteration        = iter,
      step             = step,
      n                = pilots,
      success_count    = s,
      success_rate     = s / B,
      f                = f_used,
      seed_offset      = offset_used,
      glm_intercept    = intercept,
      glm_slope        = slope,
      n_next           = n_next,
      stopping_reason  = NA_character_,
      stringsAsFactors = FALSE
    )
    n <- n_next
  }

  diagnostics <- do.call(rbind, diag_rows[seq_len(iter)])
  rownames(diagnostics) <- NULL
  if (converged) {
    stopping_reason <- "tolerance"
  } else if (!is.na(last_fit_n)) {
    stopping_reason <- "max_iterations"
    warning(sprintf(
      "Sample-size solver for alpha = %s did not converge within %d iterations; returning last fitted n = %d.",
      format(alpha), max_iter, last_fit_n
    ), call. = FALSE)
  } else {
    stop(sprintf(paste0(
      "Sample-size solver for alpha = %s made no successful curve fit in %d iterations (last step: %s). ",
      "The success curve never showed a positive slope in log(n); increase max_iterations or check the simulator."
    ), format(alpha), max_iter, diagnostics$step[nrow(diagnostics)]), call. = FALSE)
  }
  diagnostics$stopping_reason[nrow(diagnostics)] <- stopping_reason

  list(
    final_n         = last_fit_n,
    stopping_reason = stopping_reason,
    iterations_used = iter,
    diagnostics     = diagnostics
  )
}
