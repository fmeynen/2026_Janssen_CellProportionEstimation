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

#' Evaluate hybrid AE/ARE success at the cell-type level with strict cutoff split.
#'
#' Uses the strict rule requested for isoband automation:
#' values with `obs < cutoff` are checked against `tau_AE`,
#' and values with `obs >= cutoff` are checked against `tau_ARE`.
#'
#' @param phat_mat B x K numeric matrix of observed proportions.
#' @param p        True proportion vector of length K.
#' @param cutoff   Numeric scalar cutoff on observed proportions.
#' @param tau_AE   Numeric scalar AE success threshold.
#' @param tau_ARE  Numeric scalar ARE success threshold.
#'
#' @return List with the same structure as `evaluate_hybrid_success_cell_level()`.
evaluate_hybrid_success_cell_level_strict <- function(phat_mat, p, cutoff, tau_AE, tau_ARE) {
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
  use_AE_matrix <- phat_mat < cutoff
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

#' Fit hard-coded isoband GAM surrogate for replicate-level success.
#'
#' @param design_results_df Data frame containing:
#'   `tau_AE`, `tau_ARE`, `cutoff`, `n_success`, and `n_total`.
#'
#' @return List with elements `model`, `fit_data`, and `converged`.
fit_isoband_gam <- function(design_results_df) {
  required_cols <- c("tau_AE", "tau_ARE", "cutoff", "n_success", "n_total")
  if (!is.data.frame(design_results_df) || !all(required_cols %in% names(design_results_df))) {
    stop(
      sprintf("design_results_df must contain columns: %s", paste(required_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  if (nrow(design_results_df) == 0L) {
    stop("design_results_df must contain at least one row.", call. = FALSE)
  }

  fit_data <- design_results_df[, required_cols, drop = FALSE]
  if (any(!is.finite(as.matrix(fit_data)))) {
    stop("design_results_df contains non-finite values in required columns.", call. = FALSE)
  }
  if (any(fit_data$n_total <= 0) || any(fit_data$n_success < 0) || any(fit_data$n_success > fit_data$n_total)) {
    stop("n_success and n_total must satisfy 0 <= n_success <= n_total and n_total > 0.", call. = FALSE)
  }

  gam_model <- mgcv::gam(
    cbind(n_success, n_total - n_success) ~ mgcv::te(tau_AE, tau_ARE, cutoff, k = c(10, 10, 10)),
    family = stats::binomial(link = "logit"),
    method = "REML",
    data = fit_data
  )
  converged <- isTRUE(gam_model$converged)
  if (!converged && !is.null(gam_model$outer.info$converged)) {
    converged <- isTRUE(gam_model$outer.info$converged)
  }

  list(model = gam_model, fit_data = fit_data, converged = converged)
}

#' Predict isoband GAM surrogate values for candidate triads.
#'
#' @param gam_fit Output list from `fit_isoband_gam()` or a `gam` model object.
#' @param newdata Data frame containing `tau_AE`, `tau_ARE`, and `cutoff`.
#' @param se_fit  Logical; whether to include standard error on link scale.
#'
#' @return Data frame containing `newdata`, `p_hat`, and optional `se_link`.
predict_isoband_gam <- function(gam_fit, newdata, se_fit = TRUE) {
  required_cols <- c("tau_AE", "tau_ARE", "cutoff")
  if (!is.data.frame(newdata) || !all(required_cols %in% names(newdata))) {
    stop(
      sprintf("newdata must contain columns: %s", paste(required_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  if (nrow(newdata) == 0L) {
    out <- newdata
    out$p_hat <- numeric(0L)
    if (isTRUE(se_fit)) out$se_link <- numeric(0L)
    return(out)
  }

  model <- if (inherits(gam_fit, "gam")) gam_fit else gam_fit$model
  if (is.null(model) || !inherits(model, "gam")) {
    stop("gam_fit must be a gam model or the output list from fit_isoband_gam().", call. = FALSE)
  }

  pred <- stats::predict(model, newdata = newdata, type = "link", se.fit = isTRUE(se_fit))
  if (isTRUE(se_fit)) {
    out <- newdata
    out$p_hat <- stats::plogis(as.numeric(pred$fit))
    out$se_link <- as.numeric(pred$se.fit)
    return(out)
  }

  out <- newdata
  out$p_hat <- stats::plogis(as.numeric(pred))
  out
}

#' Extract point-cloud isoband from surrogate predictions.
#'
#' @param pred_df Prediction data frame containing at least `p_hat`.
#' @param p0      Target success probability in (0, 1).
#' @param eps     Tolerance around `p0`.
#'
#' @return Filtered data frame with added `delta = abs(p_hat - p0)`.
extract_isoband_points <- function(pred_df, p0, eps) {
  if (!is.data.frame(pred_df) || !("p_hat" %in% names(pred_df))) {
    stop("pred_df must be a data.frame containing a p_hat column.", call. = FALSE)
  }
  if (!is.numeric(p0) || length(p0) != 1L || !is.finite(p0) || p0 <= 0 || p0 >= 1) {
    stop("p0 must be a finite scalar strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("eps must be a finite positive scalar.", call. = FALSE)
  }

  out <- pred_df
  out$delta <- abs(out$p_hat - p0)
  out[out$delta <= eps, , drop = FALSE]
}

#' Score candidate triads for adaptive isoband sampling.
#'
#' @param pred_df Data frame containing `p_hat` and optional `se_link`.
#' @param p0      Target success probability in (0, 1).
#'
#' @return Input data frame with `score_primary` and `score_secondary`.
score_candidates_isoband <- function(pred_df, p0) {
  if (!is.data.frame(pred_df) || !("p_hat" %in% names(pred_df))) {
    stop("pred_df must be a data.frame containing a p_hat column.", call. = FALSE)
  }
  if (!is.numeric(p0) || length(p0) != 1L || !is.finite(p0) || p0 <= 0 || p0 >= 1) {
    stop("p0 must be a finite scalar strictly between 0 and 1.", call. = FALSE)
  }

  out <- pred_df
  out$score_primary <- abs(out$p_hat - p0)
  out$score_secondary <- if ("se_link" %in% names(out)) -out$se_link else 0
  out
}
