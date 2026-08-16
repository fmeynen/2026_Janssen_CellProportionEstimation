# R/calculation.R
#
# Calculation layer: error metrics, threshold evaluation, and hybrid cutoff logic.
#
# Depends on: (none from this package; uses base R only)
# ---------------------------------------------------------------------------


#' Compute per-cell-type error vectors for each requested metric.
#'
#' @param phat    Observed proportion vector (length K).
#' @param p       True proportion vector (length K).
#' @param metrics Character vector; any subset of `c("AE", "ARE", "TSE", "LAE")`.
#' @param n       Total sample size; required only when `"TSE"` is in `metrics`.
#'
#' @return Named list with one numeric vector per metric (length K).
#'
#' @details
#' AE  = abs(phat - p)
#' ARE = abs(phat - p) / p  (no epsilon stabilisation; NaN/Inf for p == 0 is expected)
#' TSE = asinh(sqrt(2 * n^2 * (phat - p)^2))  (requires n)
#' LAE = log(abs(phat - p))
compute_errors <- function(phat, p, metrics = c("AE", "ARE"), n = NULL) {
  stopifnot(length(phat) == length(p))
  if ("TSE" %in% metrics && is.null(n)) {
    stop("n must be provided when metric 'TSE' is requested.", call. = FALSE)
  }
  result <- list()
  if ("AE" %in% metrics) {
    result[["AE"]] <- abs(phat - p)
  }
  if ("ARE" %in% metrics) {
    result[["ARE"]] <- abs(phat - p) / p   # ARE: NaN when p == 0 and phat == 0 (0/0); Inf when p == 0 and phat != 0; no stabilisation by design
  }
  if ("TSE" %in% metrics) {
    result[["TSE"]] <- asinh(sqrt(2 * n^2 * (phat - p)^2))
  }
  if ("LAE" %in% metrics) {
    result[["LAE"]] <- log(abs(phat - p))
  }
  result
}


#' Compute success rates for each metric across a grid of thresholds.
#'
#' Runs post-hoc on stored max errors; no re-simulation needed.
#'
#' @param max_errors B x M matrix of max error values (from `run_replicates()`).
#' @param taus       Either a numeric vector of threshold values applied to
#'   every metric, or a named list with one numeric vector per metric
#'   (e.g. `list(AE = c(...), ARE = c(...))`).  When a list is supplied every
#'   metric present in `max_errors` must have an entry.
#'
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


#' Evaluate hybrid AE/ARE success at the cell-type level for one cutoff.
#'
#' @param phat_mat B x K numeric matrix of observed proportions.
#' @param p        True proportion vector of length K.
#' @param cutoff   Numeric scalar cutoff on observed proportions.
#' @param tau_AE   Numeric scalar AE success threshold.
#' @param tau_ARE  Numeric scalar ARE success threshold.
#'
#' @return List with elements:
#'   \describe{
#'     \item{success_matrix}{B x K logical matrix of hybrid success indicators.}
#'     \item{use_AE_matrix}{B x K logical matrix; TRUE when AE is selected.}
#'     \item{success_rate_cell}{Mean success across all replicate-index pairs.}
#'     \item{success_rate_replicate}{Fraction of replicates where all K indices succeed.}
#'   }
evaluate_hybrid_success_cell_level <- function(phat_mat, p, cutoff, tau_AE, tau_ARE) {
  stopifnot(
    is.matrix(phat_mat),
    is.numeric(phat_mat),
    is.numeric(p),
    length(p) == ncol(phat_mat),
    all(is.finite(p)),
    all(p > 0),
    is.numeric(cutoff),
    length(cutoff) == 1L,
    is.finite(cutoff),
    is.numeric(tau_AE),
    length(tau_AE) == 1L,
    is.finite(tau_AE),
    is.numeric(tau_ARE),
    length(tau_ARE) == 1L,
    is.finite(tau_ARE)
  )

  p_mat <- matrix(p, nrow = nrow(phat_mat), ncol = ncol(phat_mat), byrow = TRUE)
  ae_mat <- abs(phat_mat - p_mat)
  are_mat <- ae_mat / p_mat
  use_AE_matrix <- phat_mat <= cutoff
  success_matrix <- ifelse(use_AE_matrix, ae_mat <= tau_AE, are_mat <= tau_ARE)

  list(
    success_matrix = success_matrix,
    use_AE_matrix = use_AE_matrix,
    success_rate_cell = mean(success_matrix, na.rm = TRUE),
    success_rate_replicate = mean(rowSums(success_matrix, na.rm = TRUE) == ncol(success_matrix))
  )
}

#' Sweep hybrid AE/ARE cutoffs and summarize success rates.
#'
#' @param phat_mat B x K numeric matrix of observed proportions.
#' @param p        True proportion vector of length K.
#' @param cutoffs  Numeric vector of candidate cutoffs.
#' @param tau_AE   Numeric scalar AE success threshold.
#' @param tau_ARE  Numeric scalar ARE success threshold.
#'
#' @return Tidy data.frame with one row per cutoff and columns:
#'   cutoff, success_rate_cell, success_rate_replicate,
#'   prop_using_AE, prop_using_ARE.
sweep_hybrid_cutoffs_cell_level <- function(phat_mat, p, cutoffs, tau_AE, tau_ARE) {
  stopifnot(is.numeric(cutoffs), length(cutoffs) >= 1L, all(is.finite(cutoffs)))

  rows <- lapply(cutoffs, function(cutoff_i) {
    eval_i <- evaluate_hybrid_success_cell_level(
      phat_mat = phat_mat,
      p = p,
      cutoff = cutoff_i,
      tau_AE = tau_AE,
      tau_ARE = tau_ARE
    )
    prop_using_AE_i <- mean(eval_i$use_AE_matrix, na.rm = TRUE)
    data.frame(
      cutoff = cutoff_i,
      success_rate_cell = eval_i$success_rate_cell,
      success_rate_replicate = eval_i$success_rate_replicate,
      prop_using_AE = prop_using_AE_i,
      prop_using_ARE = 1 - prop_using_AE_i,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' Find the best hybrid cutoff from a sweep curve.
#'
#' @param curve_df  Data.frame returned by `sweep_hybrid_cutoffs_cell_level()`.
#' @param maximize  Which success rate to maximize: `"cell"` or `"replicate"`.
#' @param tie_break Tie-breaker among maximizing cutoffs:
#'   `"smallest"`, `"largest"`, or `"median"`.
#'
#' @return List with elements:
#'   best_cutoff, best_success_rate_cell, best_success_rate_replicate,
#'   best_prop_using_AE, n_tied_maximizers.
find_best_hybrid_cutoff <- function(curve_df,
                                    maximize = c("cell", "replicate"),
                                    tie_break = c("smallest", "largest", "median")) {
  maximize <- match.arg(maximize)
  tie_break <- match.arg(tie_break)

  required_cols <- c(
    "cutoff", "success_rate_cell", "success_rate_replicate",
    "prop_using_AE", "prop_using_ARE"
  )
  if (!all(required_cols %in% names(curve_df))) {
    stop(
      sprintf(
        "curve_df must contain columns: %s",
        paste(required_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (nrow(curve_df) == 0L) {
    stop("curve_df must contain at least one row.", call. = FALSE)
  }

  target_col <- if (identical(maximize, "cell")) "success_rate_cell" else "success_rate_replicate"
  target_values <- curve_df[[target_col]]
  best_value <- max(target_values, na.rm = TRUE)
  tied_idx <- which(target_values == best_value)
  tied_cutoffs <- curve_df$cutoff[tied_idx]

  best_cutoff <- switch(
    tie_break,
    smallest = min(tied_cutoffs, na.rm = TRUE),
    largest = max(tied_cutoffs, na.rm = TRUE),
    median = stats::median(tied_cutoffs, na.rm = TRUE)
  )

  chosen_row <- curve_df[curve_df$cutoff == best_cutoff, , drop = FALSE][1L, , drop = FALSE]

  list(
    best_cutoff = best_cutoff,
    best_success_rate_cell = chosen_row$success_rate_cell[[1L]],
    best_success_rate_replicate = chosen_row$success_rate_replicate[[1L]],
    best_prop_using_AE = chosen_row$prop_using_AE[[1L]],
    n_tied_maximizers = length(tied_idx)
  )
}

#' Run hybrid cutoff analysis for one simulation scenario.
#'
#' @param phat_mat  B x K numeric matrix of observed proportions.
#' @param p         True proportion vector of length K.
#' @param cutoffs   Numeric vector of candidate cutoffs.
#' @param tau_AE    Numeric scalar AE success threshold.
#' @param tau_ARE   Numeric scalar ARE success threshold.
#' @param maximize  Which success rate to maximize: `"cell"` or `"replicate"`.
#' @param alpha     Optional scalar metadata carried into outputs.
#' @param p_max     Optional scalar metadata carried into outputs.
#'
#' @return List with two elements:
#'   \describe{
#'     \item{cutoff_curve}{Data.frame from `sweep_hybrid_cutoffs_cell_level()`.}
#'     \item{best_summary}{One-row data.frame with the selected best cutoff summary.}
#'   }
run_hybrid_cutoff_analysis <- function(phat_mat, p, cutoffs, tau_AE, tau_ARE,
                                       maximize = c("cell", "replicate"),
                                       alpha = NULL, p_max = NULL) {
  maximize <- match.arg(maximize)
  cutoff_curve <- sweep_hybrid_cutoffs_cell_level(
    phat_mat = phat_mat,
    p = p,
    cutoffs = cutoffs,
    tau_AE = tau_AE,
    tau_ARE = tau_ARE
  )
  best <- find_best_hybrid_cutoff(
    curve_df = cutoff_curve,
    maximize = maximize,
    tie_break = "smallest"
  )

  best_summary <- data.frame(
    alpha = if (is.null(alpha)) NA_real_ else alpha,
    p_max = if (is.null(p_max)) NA_real_ else p_max,
    maximize = maximize,
    tau_AE = tau_AE,
    tau_ARE = tau_ARE,
    best_cutoff = best$best_cutoff,
    best_success_rate_cell = best$best_success_rate_cell,
    best_success_rate_replicate = best$best_success_rate_replicate,
    best_prop_using_AE = best$best_prop_using_AE,
    n_tied_maximizers = best$n_tied_maximizers,
    stringsAsFactors = FALSE
  )

  list(
    cutoff_curve = cutoff_curve,
    best_summary = best_summary
  )
}


# ---------------------------------------------------------------------------
# Sample-size estimation helpers
# ---------------------------------------------------------------------------

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
#' At each iteration three pilot sample sizes are evaluated: 95%, 100%, and
#' 105% of the current estimate.  A binomial GLM is fitted to the resulting
#' success counts and inverted to obtain the next estimate.  Iteration stops
#' when `abs(new_n - old_n) <= config$sample_size_tolerance` or
#' `config$max_iterations` is reached.  Sample sizes are always rounded up.
#'
#' @param alpha   Positive scalar Beta shape parameter.
#' @param n_init  Initial sample-size estimate (positive integer).
#' @param config  Named list; must contain all fields required by
#'   `simulate_success_at_n()` plus:
#'   \describe{
#'     \item{success_rate_target}{Target success probability in (0, 1).}
#'     \item{sample_size_tolerance}{Stopping tolerance (non-negative integer
#'       or numeric).}
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

    glm_fit  <- fit_success_glm(pilot_ns, success_counts, B)
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
    } else {
      lower_bound <- as.integer(ceiling(0.5 * min(pilot_ns)))
      new_n <- max(new_n, lower_bound)
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
