# Calculation Layer -----------------------------------------------------------------------------------------------
# Calculation layer: error metrics (compute_errors), threshold evaluation and success rates,
# the success rule (pooled-proportion statistic), and the sample-size solver.

# Compute Errors --------------------------------------------------------------------------------------------------

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
  if (length(phat) != length(p)) {
    stop("phat and p must have the same length.", call. = FALSE)
  }
  if ("TSE" %in% metrics && is.null(n)) {
    stop("n must be provided when metric 'TSE' is requested.", call. = FALSE)
  }
  result <- list()
  if ("AE" %in% metrics) {
    result[["AE"]] <- abs(phat - p)
  }
  if ("ARE" %in% metrics) {
    result[["ARE"]] <- abs(phat - p) / p   # ARE: NaN when p == 0 and phat == 0 (0/0);
                                           # Inf when p == 0 and phat != 0; no stabilisation by design
  }
  if ("TSE" %in% metrics) {
    result[["TSE"]] <- asinh(sqrt(2 * n^2 * (phat - p)^2))
  }
  if ("LAE" %in% metrics) {
    result[["LAE"]] <- log(abs(phat - p))
  }
  result
}


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
  validate_named_matrix(max_errors, "max_errors")
  metrics <- colnames(max_errors)
  has_error_array <- !is.null(errors)
  if (has_error_array) {
    if (length(dim(errors)) != 3L) {
      stop("errors must be a 3-dimensional array.", call. = FALSE)
    }
    if (dim(errors)[1] != nrow(max_errors)) {
      stop("errors must have the same number of replicates (first dimension) as max_errors has rows.", call. = FALSE)
    }
    if (dim(errors)[3] != length(metrics)) {
      stop("errors must have one slice (third dimension) per column of max_errors.", call. = FALSE)
    }
  }

  # Normalise taus: plain vector -> same grid for all metrics;
  # named list     -> per-metric grids.
  if (is.numeric(taus)) {
    taus_list <- setNames(rep(list(taus), length(metrics)), metrics)
  } else if (is.list(taus)) {
    missing_metrics <- setdiff(metrics, names(taus))
    if (length(missing_metrics) > 0L) {
      stop(
        sprintf(
          "`taus` list is missing entries for metric(s): %s",
          paste(missing_metrics, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    taus_list <- taus[metrics]
  } else {
    stop("`taus` must be a numeric vector or a named list of numeric vectors.", call. = FALSE)
  }

  rows <- vector("list", length(metrics))
  for (i in seq_along(metrics)) {
    m <- metrics[[i]]
    tau_m <- taus_list[[m]]
    rates <- vapply(
      tau_m,
      function(tau) mean(max_errors[, m] <= tau, na.rm = TRUE),
      numeric(1L)
    )
    if (has_error_array) {
      errors_m <- errors[,, m, drop = TRUE]
      if (is.null(dim(errors_m))) {
        errors_m <- matrix(errors_m, nrow = nrow(max_errors))
      }
      mean_n_above <- vapply(
        tau_m,
        function(tau) mean(rowSums(errors_m > tau), na.rm = TRUE),
        numeric(1L)
      )
    } else {
      mean_n_above <- rep(NA_real_, length(tau_m))
    }
    rows[[i]] <- data.frame(
      metric = m,
      tau = tau_m,
      success_rate = rates,
      mean_n_above = mean_n_above,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


# Success Rate From Stat Vector --------------------------------------------------------------------------------------

#' Compute success rate at each threshold from a per-replicate statistic.
#'
#' A replicate succeeds at threshold `tau` iff `stat <= tau`. Vectorised over `taus`. `stat` may contain `Inf` (e.g.
#' relative error when a true proportion is 0), which never counts as a success at any finite `tau`.
#'
#' @param stat Numeric vector (length >= 1, no `NA`) of per-replicate values, e.g. the per-replicate max-over-cell-
#'   types error.
#' @param taus Numeric vector (no `NA`) of threshold values.
#'
#' @return Numeric vector, same length as `taus`, giving `mean(stat <= tau)` for each `tau`. Unnamed.
success_rate_from_stat <- function(stat, taus) {
  if (!is.numeric(stat) || length(stat) < 1L || anyNA(stat)) {
    stop("`stat` must be a numeric vector of length >= 1 with no NA values.", call. = FALSE)
  }
  if (!is.numeric(taus) || anyNA(taus)) {
    stop("`taus` must be a numeric vector with no NA values.", call. = FALSE)
  }
  unname(vapply(taus, function(t) mean(stat <= t), numeric(1L)))
}


#' Build a data-driven grid of thresholds for plotting success rate vs tau.
#'
#' The grid spans `[0, upper]` where `upper` is the `prob`-quantile of `stat`. If that quantile is not finite (e.g.
#' `stat` contains `Inf` at high `prob`), it is recomputed using only the finite values of `stat`.
#'
#' @param stat     Numeric vector (length >= 1, no `NA`) of per-replicate values; may contain `Inf`.
#' @param n_points Number of grid points (positive integer, >= 2). Default 200.
#' @param prob     Quantile probability in `(0, 1]` used to set the grid's upper end. Default 0.999.
#'
#' @return Numeric vector of length `n_points`, `seq(0, upper, length.out = n_points)`.
default_tau_grid <- function(stat, n_points = 200, prob = 0.999) {
  if (!is.numeric(stat) || length(stat) < 1L || anyNA(stat)) {
    stop("`stat` must be a numeric vector of length >= 1 with no NA values.", call. = FALSE)
  }
  n_points <- validate_positive_integer(n_points, "n_points")
  if (n_points < 2L) {
    stop("`n_points` must be >= 2.", call. = FALSE)
  }
  if (!is.numeric(prob) || length(prob) != 1L || !is.finite(prob) || prob <= 0 || prob > 1) {
    stop("`prob` must be a single number in (0, 1].", call. = FALSE)
  }

  upper <- stats::quantile(stat, probs = prob, names = FALSE, na.rm = TRUE)
  if (!is.finite(upper)) {
    finite_stat <- stat[is.finite(stat)]
    if (length(finite_stat) == 0L) {
      stop("`stat` has no finite values; cannot build a tau grid.", call. = FALSE)
    }
    upper <- stats::quantile(finite_stat, probs = prob, names = FALSE, na.rm = TRUE)
  }
  if (!is.finite(upper) || upper <= 0) {
    stop("Computed an invalid (non-finite or non-positive) upper bound for the tau grid.", call. = FALSE)
  }
  seq(0, upper, length.out = n_points)
}


# Success Rule ----------------------------------------------------------------------------------------------------

#' Identifier of the current replicate success rule.
#'
#' Included in the cache keys of experiments whose cached results depend on the rule (`run_sample_size_experiment()`,
#' `run_dm_errorchoice_experiment()`), so that changing the rule never silently reuses results computed under an
#' earlier one. Change this string whenever `replicate_pooled_error()` changes what it computes.
#'
#' @return Character scalar.
success_rule_id <- function() {
  "pooled_proportion_v1"
}


#' Compute the per-replicate pooled-proportion error statistic for one metric.
#'
#' Single source of truth for the success statistic used by `replicate_success()` (and by any caller that needs the
#' raw per-replicate statistic, e.g. to evaluate many tau thresholds via `mean(stat <= tau)` without
#' recomputation). Per (scenario, replicate, cell type j), the estimated proportions are first averaged over persons,
#' giving the pooled estimate `pbar_j`, which is compared with the population-level true proportion `p_j`:
#'   * AE:  `|pbar_j - p_j|`
#'   * ARE: `|pbar_j - p_j| / p_j`
#' The largest of these cell-type errors is the per-replicate "stat"; a replicate passes the metric iff
#' `stat <= tau`. `NaN` (ARE with `pbar_j` and `p_j` both 0) counts as 0; `Inf` (ARE with `p_j = 0` and
#' `pbar_j > 0`) is kept. With a single person per replicate this reduces to the per-draw error against `p`.
#'
#' The computation is fully vectorised: integer group keys, `rowsum()` for the per-cell-type means, and a single
#' `order()` + `duplicated()` pass for the max-by-group step. There is no `stats::aggregate()` call, no
#' `tapply()`, and no per-row loop.
#'
#' @param person_results Data.frame with at least the columns `replicate`, `cell_type`, `metric`,
#'   `observed_proportion` and `population_mean_proportion`, one row per (replicate, person, cell type, metric). A
#'   `scenario_id` column is optional: when it is absent, or present but entirely `NA`, all rows are treated as a
#'   single scenario and the output's `scenario_id` is `NA_character_`.
#' @param metric Single metric name, `"AE"` or `"ARE"`. Rows of `person_results` are filtered to this metric (so
#'   each person's proportions are counted once); an error is raised if none match.
#'
#' @return Data.frame with one row per (scenario_id, replicate), sorted by scenario_id then replicate, with
#'   columns scenario_id, replicate, stat (max over cell types of the pooled-proportion error).
replicate_pooled_error <- function(person_results, metric) {
  validate_required_columns(person_results, person_results_required_cols)
  if (!is.character(metric) || length(metric) != 1L || is.na(metric)) {
    stop("metric must be a single character string.", call. = FALSE)
  }
  if (!metric %in% c("AE", "ARE")) {
    stop(sprintf("The pooled success rule is defined for AE and ARE only, not '%s'.", metric), call. = FALSE)
  }

  idx <- person_results$metric == metric
  if (!any(idx)) {
    stop(sprintf("metric '%s' is not present in person_results.", metric), call. = FALSE)
  }

  # A missing/all-NA scenario_id is replaced by a non-NA sentinel for grouping purposes (rowsum() drops NA
  # groups), and mapped back to NA_character_ in the output.
  has_scenario <- "scenario_id" %in% names(person_results)
  scenario_raw <- if (has_scenario) {
    as.character(person_results$scenario_id)
  } else {
    rep(NA_character_, nrow(person_results))
  }
  no_scenario <- all(is.na(scenario_raw))
  scenario_key <- if (no_scenario) rep("__single_scenario__", length(scenario_raw)) else scenario_raw

  scenario_i   <- scenario_key[idx]
  replicate_i  <- person_results$replicate[idx]
  cell_type_i  <- person_results$cell_type[idx]
  observed_i   <- person_results$observed_proportion[idx]
  population_i <- person_results$population_mean_proportion[idx]

  scenario_lv  <- sort(unique(scenario_i))
  replicate_lv <- sort(unique(replicate_i))
  cell_type_lv <- sort(unique(cell_type_i))
  nS <- length(scenario_lv)
  nR <- length(replicate_lv)
  nC <- length(cell_type_lv)

  s_idx <- match(scenario_i, scenario_lv)
  r_idx <- match(replicate_i, replicate_lv)
  c_idx <- match(cell_type_i, cell_type_lv)

  # Stage 1: pooled estimate (mean over persons) and population proportion per (scenario, replicate, cell_type) via
  # one integer group key and rowsum(); then the per-cell-type error of the pooled estimate.
  key1 <- (s_idx - 1L) * nR * nC + (r_idx - 1L) * nC + c_idx
  counts1     <- rowsum(rep(1L, length(observed_i)), key1)[, 1L]
  pooled1     <- rowsum(observed_i, key1)[, 1L] / counts1
  population1 <- rowsum(population_i, key1)[, 1L] / counts1
  error1 <- abs(pooled1 - population1)
  if (identical(metric, "ARE")) {
    error1 <- error1 / population1
  }
  error1[is.nan(error1)] <- 0
  key1_sorted <- as.numeric(names(counts1))

  tmp1 <- (key1_sorted - 1L) %/% nC
  r1 <- (tmp1 %% nR) + 1L
  s1 <- (tmp1 %/% nR) + 1L

  # Stage 2: max over cell types per (scenario, replicate). Sort by group then by value descending and keep the
  # first row of each group -- fully vectorised, no tapply()/aggregate() and no per-row loop.
  key2 <- (s1 - 1L) * nR + r1
  ord  <- order(key2, -error1)
  keep <- !duplicated(key2[ord])
  max_key2 <- key2[ord][keep]
  stat     <- error1[ord][keep]

  r2 <- ((max_key2 - 1L) %% nR) + 1L
  s2 <- ((max_key2 - 1L) %/% nR) + 1L

  out <- data.frame(
    scenario_id = scenario_lv[s2],
    replicate   = replicate_lv[r2],
    stat        = stat,
    stringsAsFactors = FALSE
  )
  if (no_scenario) out$scenario_id <- NA_character_

  out <- out[order(out$scenario_id, out$replicate), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Determine per-replicate success against a set of metric thresholds.
#'
#' Replicate success rule shared by `extract_success_rate()` and `simulate_success_at_n()` (both models). For each
#' metric, the per-replicate statistic comes from `replicate_pooled_error()`, the single source of truth for the rule:
#' the estimated proportions are averaged over persons per cell type, compared with the population-level true
#' proportion, and the largest cell-type error is taken. A replicate passes that metric if this value is `<= tau`,
#' and passes overall if it passes every metric in `taus`. This function only turns each metric's `stat` into a
#' `pass_<metric>` flag and ANDs them together.
#'
#' @param person_results Data.frame with at least the columns `replicate`, `cell_type`, `metric`,
#'   `observed_proportion`, `population_mean_proportion`. A `scenario_id` column is optional: when it is absent, or present but entirely `NA`, all rows are treated as a
#'   single scenario and the output's `scenario_id` is `NA_character_`.
#' @param taus   Named list with one scalar threshold per metric (e.g. `list(AE = 0.02, ARE = 0.5)`). Metrics missing
#'   from `person_results$metric` are skipped with a warning; an error is raised if none of the metrics remain.
#'
#' @return Data.frame with one row per (scenario_id, replicate), sorted by scenario_id then replicate, with columns
#'   scenario_id, replicate, `pass_<metric>` for each metric used, and `pass` (logical AND across all `pass_<metric>`
#'   columns).
replicate_success <- function(person_results, taus) {
  validate_required_columns(person_results, person_results_required_cols)
  if (!is.list(taus) || is.null(names(taus)) || any(names(taus) == "")) {
    stop("taus must be a named list with one scalar threshold per metric.", call. = FALSE)
  }
  for (m in names(taus)) {
    if (!is.numeric(taus[[m]]) || length(taus[[m]]) != 1L || is.na(taus[[m]])) {
      stop(sprintf("taus$%s must be a single numeric threshold.", m), call. = FALSE)
    }
  }

  metrics <- names(taus)
  missing_metrics <- setdiff(metrics, unique(person_results$metric))
  for (m in missing_metrics) {
    warning(sprintf("replicate_success: metric '%s' is not in person_results; it is skipped.", m), call. = FALSE)
  }
  metrics <- setdiff(metrics, missing_metrics)
  if (length(metrics) == 0L) {
    stop("None of the metrics in taus are present in person_results.", call. = FALSE)
  }

  replicate_pass <- NULL
  for (m in metrics) {
    stat_df <- replicate_pooled_error(person_results, m)
    stat_df[[paste0("pass_", m)]] <- stat_df$stat <= taus[[m]]
    stat_df$stat <- NULL
    replicate_pass <- if (is.null(replicate_pass)) {
      stat_df
    } else {
      merge(replicate_pass, stat_df, by = c("scenario_id", "replicate"))
    }
  }
  pass_cols <- paste0("pass_", metrics)
  replicate_pass$pass <- Reduce(`&`, replicate_pass[pass_cols])

  replicate_pass <- replicate_pass[order(replicate_pass$scenario_id, replicate_pass$replicate), , drop = FALSE]
  rownames(replicate_pass) <- NULL
  replicate_pass[, c("scenario_id", "replicate", pass_cols, "pass")]
}


# Sample-Size Estimation --------------------------------------------------------------------------------------

#' Pilot sample sizes around a centre on a multiplicative grid.
#'
#' @param n Numeric scalar, the current centre (> 0).
#' @param f Numeric scalar spread factor (> 1). Pilots are placed at `n / f`, `n` and `n * f`.
#' @param n_max Upper cap on any pilot size (<= `.Machine$integer.max`).
#'
#' @return Sorted, unique integer vector of pilot sizes, each rounded up and within `[1, n_max]`.
sample_size_pilots <- function(n, f, n_max = 1e9) {
  if (!is_finite_scalar(n) || n <= 0) {
    stop("n must be a single finite number > 0.", call. = FALSE)
  }
  if (!is_finite_scalar(f) || f <= 1) {
    stop("f must be a single finite number > 1.", call. = FALSE)
  }
  validate_n_max(n_max)
  sort(unique(as.integer(pmin(n_max, pmax(1, ceiling(c(n / f, n, n * f)))))))
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
  if (!is.numeric(n_values) || length(n_values) < 1L || !isTRUE(all(n_values >= 1))) {
    stop("n_values must be a non-empty numeric vector of values >= 1.", call. = FALSE)
  }
  if (!is.numeric(success_count) || length(success_count) != length(n_values)) {
    stop("success_count must be a numeric vector with the same length as n_values.", call. = FALSE)
  }
  if (!is_finite_scalar(B) || B < 1) {
    stop("B must be a single number >= 1.", call. = FALSE)
  }
  if (!isTRUE(all(success_count >= 0)) || !isTRUE(all(success_count <= B))) {
    stop("success_count must lie between 0 and B.", call. = FALSE)
  }
  dat <- data.frame(n = n_values, s = success_count, fails = B - success_count)
  withCallingHandlers(
    stats::glm(
      cbind(s, fails) ~ log(n),
      family = stats::binomial(),
      data = dat
    ),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (
        grepl("fitted probabilities numerically 0 or 1", msg, fixed = TRUE) ||
          grepl("algorithm did not converge", msg, fixed = TRUE)
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )
}


#' Invert a fitted success curve at a target success rate.
#'
#' Solves `qlogis(target) = intercept + slope * log(n)` for n. A non-finite or negative slope is set to 0: the curve
#' is then flat, so the target is `n_max` when the fitted success rate is below the target and 1 otherwise. A solved n
#' that is not finite or exceeds `n_max` falls back to `n_max`.
#'
#' @param fit    A `glm` from `fit_success_curve()`; its intercept must be finite.
#' @param target Target success rate in (0, 1).
#' @param n_max  Fallback and upper cap for the solved n (<= `.Machine$integer.max`).
#'
#' @return Integer n, rounded up, within `[1, n_max]`.
solve_success_curve <- function(fit, target, n_max = 1e9) {
  validate_open_unit_scalar(target, "target")
  validate_n_max(n_max)
  coefs <- stats::coef(fit)
  intercept <- unname(coefs[[1L]])
  slope <- unname(coefs[[2L]])
  if (!is.finite(intercept)) {
    stop("Cannot invert the success curve: the intercept is not finite.", call. = FALSE)
  }
  gap <- stats::qlogis(target) - intercept
  n_raw <- if (is.finite(slope) && slope > 0) {
    exp(gap / slope)
  } else if (gap > 0) {
    Inf
  } else {
    0
  }
  if (!is.finite(n_raw) || n_raw > n_max) {
    return(as.integer(n_max))
  }
  max(1L, as.integer(ceiling(n_raw)))
}


#' Estimate the smallest sample size reaching a target success rate for one alpha.
#'
#' Iterative solver. Each iteration simulates pilots at `n / f`, `n` and `n * f` (see `sample_size_pilots()`) with
#' seed `config$seed` and accumulates them. The step type is chosen in this order:
#' \enumerate{
#'   \item No accumulated pilot has 0 < success_count < B (degenerate):
#'     `expand_up` (all of this iteration's pilots failed; centre <- min(n_max, ceiling(f^2 * max(pilots)))),
#'     `expand_down` (all succeeded; centre <- max(1, ceiling(min(pilots) / f^2))), or
#'     `bisect` (the accumulated data brackets the answer; centre <- ceiling(sqrt(largest all-fail n *
#'     smallest all-success n))). `f` is unchanged.
#'   \item Otherwise fit `fit_success_curve()` on all accumulated pilots. A non-finite or non-positive slope gives
#'     `flat`: the slope is taken as 0, so centre <- `n_max` (fitted success below target) or 1 (above), via
#'     `solve_success_curve()`; `f` unchanged, no convergence check.
#'   \item Otherwise `fit`: centre <- `solve_success_curve()`; converged when `|new_n - n| / n <= rel_tol`; then
#'     `f <- max(f_floor, sqrt(f))`. Convergence is only checked on `fit` steps.
#' }
#' No clamping is applied to any step other than the `n_max` cap.
#'
#' @param alpha    Numeric scalar passed through to `simulate`.
#' @param n_init   Starting centre (>= 1); rounded up and capped at `n_max`.
#' @param config   List with `success_rate_target`, `rel_tol`, `max_iterations`, `B`, `f0`, `f_floor`, `seed`, and
#'   optionally `n_max` (default 1e9): the cap on every centre and pilot, and the fallback target when the solved n
#'   is not finite or too large.
#' @param simulate Function `(alpha, n, config, seed)` returning a list with at least `success_count`.
#'   Defaults to `simulate_success_at_n()`; injectable for testing.
#'
#' @return List with `final_n` (integer), `stopping_reason` ("tolerance", "max_iterations" or "n_max"),
#'   `iterations_used`, and `diagnostics`: a data.frame with one row per pilot per iteration and columns alpha,
#'   iteration, step, n, success_count, success_rate, f, glm_intercept, glm_slope, n_next, stopping_reason (NA except
#'   in the final row). `f` is the value used for that iteration's pilots. When the final centre is at `n_max`, returns
#'   `n_max` with stopping_reason "n_max" and a warning. Otherwise warns when stopping at max_iterations after at
#'   least one `fit` step, and errors when no `fit` step ever happened.
estimate_sample_size <- function(
  alpha,
  n_init,
  config,
  simulate = simulate_success_at_n
) {
  if (!is_finite_scalar(n_init) || n_init < 1) {
    stop("n_init must be a single finite number >= 1.", call. = FALSE)
  }
  if (!is.function(simulate)) {
    stop("simulate must be a function.", call. = FALSE)
  }
  target <- config$success_rate_target
  rel_tol <- config$rel_tol
  max_iter <- config$max_iterations
  B <- config$B
  f0 <- config$f0
  f_floor <- config$f_floor
  base_seed <- config$seed
  validate_open_unit_scalar(target, "`config$success_rate_target`")
  if (!is_finite_scalar(rel_tol) || rel_tol < 0) {
    stop("`config$rel_tol` must be a single number >= 0.", call. = FALSE)
  }
  if (!is_finite_scalar(max_iter) || max_iter < 1) {
    stop("`config$max_iterations` must be a single number >= 1.", call. = FALSE)
  }
  if (!is_finite_scalar(B) || B < 1) {
    stop("`config$B` must be a single number >= 1.", call. = FALSE)
  }
  if (!is_finite_scalar(f0) || f0 <= 1) {
    stop("`config$f0` must be a single number > 1.", call. = FALSE)
  }
  if (!is_finite_scalar(f_floor) || f_floor <= 1 || f_floor > f0) {
    stop(
      "`config$f_floor` must be a single number with 1 < f_floor <= f0.",
      call. = FALSE
    )
  }
  if (!is_finite_scalar(base_seed)) {
    stop("`config$seed` must be a single finite number.", call. = FALSE)
  }
  n_max <- if (is.null(config$n_max)) 1e9 else config$n_max
  validate_n_max(n_max, "`config$n_max`")
  n_max <- as.integer(floor(n_max))
  max_iter <- as.integer(max_iter)

  n <- as.integer(min(n_max, ceiling(n_init)))
  f <- f0
  acc_n <- integer(0L)
  acc_s <- integer(0L)
  diag_rows <- vector("list", max_iter)
  last_fit_n <- NA_integer_
  converged <- FALSE
  iter <- 0L

  while (iter < max_iter && !converged) {
    iter <- iter + 1L
    pilots <- sample_size_pilots(n, f, n_max)
    s <- vapply(
      pilots,
      function(p) {
        as.integer(
          simulate(
            alpha = alpha,
            n = p,
            config = config,
            seed = base_seed
          )$success_count
        )
      },
      integer(1L)
    )
    if (anyNA(s) || any(s < 0L) || any(s > B)) {
      stop(
        sprintf(
          "`simulate` returned a success_count outside [0, B] at iteration %d.",
          iter
        ),
        call. = FALSE
      )
    }
    acc_n <- c(acc_n, pilots)
    acc_s <- c(acc_s, s)

    intercept <- NA_real_
    slope <- NA_real_
    f_used <- f

    if (!any(acc_s > 0L & acc_s < B)) {
      if (all(s == 0L)) {
        step <- "expand_up"
        n_next <- as.integer(min(n_max, ceiling(f^2 * max(pilots))))
      } else if (all(s == B)) {
        step <- "expand_down"
        n_next <- max(1L, as.integer(ceiling(min(pilots) / f^2)))
      } else {
        step <- "bisect"
        n_fail <- max(acc_n[acc_s == 0L])
        n_pass <- min(acc_n[acc_s == B])
        n_next <- as.integer(ceiling(sqrt(n_fail * n_pass)))
      }
    } else {
      fit <- fit_success_curve(acc_n, acc_s, B)
      coefs <- stats::coef(fit)
      intercept <- unname(coefs[[1L]])
      slope <- unname(coefs[[2L]])
      if (!is.finite(slope) || slope <= 0) {
        step <- "flat"
        n_next <- solve_success_curve(fit, target, n_max)
      } else {
        step <- "fit"
        n_next <- solve_success_curve(fit, target, n_max)
        converged <- abs(n_next - n) / n <= rel_tol
        f <- max(f_floor, sqrt(f))
        last_fit_n <- n_next
      }
    }

    diag_rows[[iter]] <- data.frame(
      alpha = alpha,
      iteration = iter,
      step = step,
      n = pilots,
      success_count = s,
      success_rate = s / B,
      f = f_used,
      glm_intercept = intercept,
      glm_slope = slope,
      n_next = n_next,
      stopping_reason = NA_character_,
      stringsAsFactors = FALSE
    )
    n <- n_next
  }

  diagnostics <- do.call(rbind, diag_rows[seq_len(iter)])
  rownames(diagnostics) <- NULL
  if (n == n_max) {
    stopping_reason <- "n_max"
    last_fit_n <- n_max
    warning(
      sprintf(
        "Sample-size solver for alpha = %s ended at the cap n_max = %d; the required n is likely larger.",
        format(alpha),
        n_max
      ),
      call. = FALSE
    )
  } else if (converged) {
    stopping_reason <- "tolerance"
  } else if (!is.na(last_fit_n)) {
    stopping_reason <- "max_iterations"
    warning(
      sprintf(
        "Sample-size solver for alpha = %s did not converge within %d iterations; returning last fitted n = %d.",
        format(alpha),
        max_iter,
        last_fit_n
      ),
      call. = FALSE
    )
  } else {
    stop(
      sprintf(
        paste0(
          "Sample-size solver for alpha = %s made no successful curve fit in %d iterations (last step: %s). ",
          "The success curve never showed a positive slope in log(n); increase max_iterations or check the simulator."
        ),
        format(alpha),
        max_iter,
        diagnostics$step[nrow(diagnostics)]
      ),
      call. = FALSE
    )
  }
  diagnostics$stopping_reason[nrow(diagnostics)] <- stopping_reason

  list(
    final_n = last_fit_n,
    stopping_reason = stopping_reason,
    iterations_used = iter,
    diagnostics = diagnostics
  )
}
