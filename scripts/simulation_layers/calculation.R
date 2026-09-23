# Calculation Layer -----------------------------------------------------------------------------------------------
# Calculation layer: error metrics, threshold evaluation, and sample size estimation helpers


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


# Sample-Size Estimation Helpers--------------------------------------------------------------------------------------

#' Fit a binomial GLM to relate sample size to observed success counts.
#'
#' @param n_values      Integer or numeric vector of pilot sample sizes.
#' @param success_count Integer vector of success counts (same length as `n_values`).
#' @param B             Integer; number of replicates per pilot size (the
#'   denominator for the proportion).
#'
#' @return A fitted `glm` object (`family = binomial`).
fit_success_glm <- function(n_values, success_count, B) {
  stopifnot(length(n_values) == length(success_count))
  stopifnot(is.numeric(B), length(B) == 1L, B >= 1L)
  failure_count <- B - success_count
  stats::glm(
    cbind(success_count, failure_count) ~ n_values,
    family = stats::binomial(link = "logit")
  )
}


#' Solve for the sample size that achieves a target success rate from a GLM.
#'
#' Inverts the logistic link: n = (logit(target) - intercept) / slope.
#' Always rounds the result **up** to the nearest integer.
#'
#' @param glm_fit            A fitted `glm` with a single predictor `n_values`.
#' @param target_success_rate Numeric scalar in (0, 1).
#'
#' @return List with:
#'   \describe{
#'     \item{n_raw}{Numeric raw (unrounded) estimate.}
#'     \item{n_rounded}{Integer; `ceiling(n_raw)`.}
#'   }
solve_sample_size_from_glm <- function(glm_fit, target_success_rate) {
  stopifnot(
    is.numeric(target_success_rate),
    length(target_success_rate) == 1L,
    target_success_rate > 0,
    target_success_rate < 1
  )
  coefs     <- stats::coef(glm_fit)
  intercept <- coefs[[1L]]
  slope     <- coefs[[2L]]
  if (!is.finite(slope) || slope == 0) {
    stop("GLM slope is zero or non-finite; cannot solve for sample size.",
         call. = FALSE)
  }
  logit_target <- log(target_success_rate / (1 - target_success_rate))
  n_raw     <- (logit_target - intercept) / slope
  list(
    n_raw     = n_raw,
    n_rounded = as.integer(ceiling(n_raw))
  )
}


#' Iteratively estimate the required sample size for one alpha.
#'
#' At each iteration three pilot sample sizes are evaluated: 95%, 100%, and 105% of the current estimate.
#' A binomial GLM is fitted to the resulting success counts and inverted to obtain the next estimate.
#' Iteration stops when `abs(new_n - old_n) <= config$sample_size_tolerance` or `config$max_iterations` is reached.
#' Sample sizes are always rounded up.
#'
#' @param alpha   Positive scalar Beta shape parameter.
#' @param n_init  Initial sample-size estimate (positive integer).
#' @param config  Named list; must contain all fields required by `simulate_success_at_n()` plus:
#'   \describe{
#'     \item{success_rate_target}{Target success probability in (0, 1).}
#'     \item{sample_size_tolerance}{Stopping tolerance (non-negative integer or numeric).}
#'     \item{max_iterations}{Maximum number of iterations (positive integer).}
#'   }
#'
#' @return List with:
#'   \describe{
#'     \item{final_n}{Final integer sample-size estimate.}
#'     \item{stopping_reason}{Character; `"tolerance"` or `"max_iterations"`.}
#'     \item{iterations_used}{Integer; number of iterations performed.}
#'     \item{diagnostics}{Long-format `data.frame` with one row per pilot
#'       point per iteration.}
#'   }
iterate_sample_size_for_alpha <- function(alpha, n_init, config) {
  stopifnot(
    is.numeric(alpha), length(alpha) == 1L, alpha > 0,
    is.numeric(n_init), length(n_init) == 1L, n_init >= 1L
  )
  target    <- config$success_rate_target
  tolerance <- config$sample_size_tolerance
  max_iter  <- config$max_iterations
  B         <- config$B
  stopifnot(
    is.numeric(target),    length(target)    == 1L, target > 0,    target < 1,
    is.numeric(tolerance), length(tolerance) == 1L, tolerance >= 0,
    is.numeric(max_iter) || is.integer(max_iter),
    length(max_iter) == 1L, max_iter >= 1L
  )
  max_iter <- as.integer(max_iter)

  current_n      <- as.integer(ceiling(n_init))
  diag_rows      <- vector("list", max_iter * 3L)
  diag_idx       <- 0L
  stopping_reason <- "max_iterations"

  # Pre-allocate accumulated pilot data across all iterations (for GLM fitting)
  all_pilot_ns       <- integer(max_iter * 3L)
  all_success_counts <- integer(max_iter * 3L)
  n_accumulated      <- 0L

  for (iter in seq_len(max_iter)) {
    pilot_ns <- as.integer(ceiling(c(0.95, 1.00, 1.05) * current_n))
    pilot_ns <- pmax(pilot_ns, 1L)   # guard against n < 1

    success_counts <- integer(3L)
    success_rates  <- numeric(3L)

    for (j in seq_along(pilot_ns)) {
      sim_j            <- simulate_success_at_n(alpha, pilot_ns[j], config)
      success_counts[j] <- sim_j$success_count
      success_rates[j]  <- sim_j$success_rate
    }

    # Accumulate evidence: add this iteration's pilot points to the history
    idx <- n_accumulated + seq_len(3L)
    all_pilot_ns[idx]       <- pilot_ns
    all_success_counts[idx] <- success_counts
    n_accumulated           <- n_accumulated + 3L

    glm_fit  <- fit_success_glm(
      all_pilot_ns[seq_len(n_accumulated)],
      all_success_counts[seq_len(n_accumulated)],
      B
    )
    solved   <- tryCatch(
      solve_sample_size_from_glm(glm_fit, target),
      error = function(e) list(n_raw = current_n, n_rounded = current_n)
    )
    new_n <- as.integer(ceiling(solved$n_rounded))

    # Clamp the new estimate to prevent extreme jumps when the pilot success
    # rates are all near 0% or all near 100%.
    # - If mean success rate < target (need more n): cap at 2x the largest pilot.
    # - If mean success rate >= target (need less n): floor at 0.5x the smallest pilot.
    mean_success_rate <- mean(success_rates)
    if (mean_success_rate < target) {
      upper_bound <- as.integer(ceiling(2.0 * max(pilot_ns)))
      new_n <- min(new_n, upper_bound)
      if(is.na(new_n)){new_n <- upper_bound}
    } else {
      lower_bound <- as.integer(ceiling(0.5 * min(pilot_ns)))
      new_n <- max(new_n, lower_bound)
      if(is.na(new_n)){new_n <- lower_bound}
    }
    new_n <- pmax(new_n, 1L)

    coefs <- stats::coef(glm_fit)
    n_unclamped <- as.integer(ceiling(solved$n_rounded))
    clamped     <- (new_n != n_unclamped)

    for (j in seq_along(pilot_ns)) {
      diag_idx <- diag_idx + 1L
      diag_rows[[diag_idx]] <- data.frame(
        alpha               = alpha,
        iteration           = iter,
        pilot_index         = j,
        n                   = pilot_ns[j],
        success_count       = success_counts[j],
        success_rate        = success_rates[j],
        target_success_rate = target,
        glm_intercept       = coefs[[1L]],
        glm_slope           = coefs[[2L]],
        n_raw               = solved$n_raw,
        n_rounded           = new_n,
        n_clamped           = clamped,
        stopping_reason     = NA_character_,
        stringsAsFactors    = FALSE
      )
    }

    if (abs(new_n - current_n) <= tolerance) {
      stopping_reason <- "tolerance"
      current_n       <- new_n
      break
    }
    current_n <- new_n
  }

  diagnostics <- do.call(rbind, diag_rows[seq_len(diag_idx)])
  diagnostics$stopping_reason[nrow(diagnostics)] <- stopping_reason

  list(
    final_n         = current_n,
    stopping_reason = stopping_reason,
    iterations_used = as.integer(diag_idx / 3L),
    diagnostics     = diagnostics
  )
}
