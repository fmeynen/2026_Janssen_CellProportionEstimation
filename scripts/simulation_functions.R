# simulation_functions.R
#
# Functions for multinomial sampling error simulation.
# Implements the MVP workflow:
#   1. Generate true proportions from a Beta curve
#   2. Simulate multinomial counts (with an extension placeholder dispatcher)
#   3. Compute AE / ARE error metrics
#   4. Extract max error with configurable tie handling
#   5. Run B replicates efficiently (simulate once, evaluate thresholds post-hoc)
#   6. Evaluate success-rate curves over a threshold grid
#   7. Summarize which cell type most often drives the max error
#   8. Orchestrate the whole experiment in one call
#
# Future extensions (dirichlet-multinomial, logistic-normal) only need:
#   - A concrete `simulate_counts_<model>()` function
#   - A new branch in the `simulate_counts()` dispatcher
# ---------------------------------------------------------------------------


# ---- A) Proportion-generation layer ----------------------------------------

#' Rescale nonnegative weights to the probability simplex.
#'
#' Internal helper used by `generate_proportions_beta()`.
normalize_to_simplex <- function(w) {
  stopifnot(is.numeric(w), all(is.finite(w)), all(w >= 0), sum(w) > 0)
  w / sum(w)
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

#' Generate K true proportions from dbeta(grid, alpha, 1), always rescaled.
#'
#' @param alpha  shape1 parameter of Beta(alpha, 1); controls curve steepness.
#' @param K      number of cell types (default 10).
#' @param grid   evaluation points in (0,1) (default K equidistant points from
#'               0.1 to 0.9, avoiding boundary values where dbeta returns 0).
#'
#' @return Numeric vector of length K, all strictly positive, summing to 1.
#'
#' @details
#' Unscaled weights: w = dbeta(grid, shape1 = alpha, shape2 = 1).
#' The default grid avoids 0 and 1 so that all weights — and therefore all
#' proportions — are strictly positive.  This prevents ARE from producing
#' NaN or Inf values for any cell type under the default parameters.
generate_proportions_beta <- function(alpha, K = 10, grid = seq(0.05, 0.95, length.out = K)) {
  stopifnot(is.numeric(alpha), length(alpha) == 1L, alpha > 0)
  stopifnot(length(grid) == K)
  w <- dbeta(grid, shape1 = alpha, shape2 = 1)
  p <- normalize_to_simplex(w)
  validate_proportions(p)
  p
}

generate_proportions_curve <- function(alpha, K = 10, grid = seq(0.05, 0.95, length.out = K)) {
  w <- dbeta(grid, shape1 = alpha, shape2 = 1)
  p <- normalize_to_simplex(w)
  s <- seq(0, 1, length.out = 1000)
  f <- dbeta(s, shape1 = alpha, shape2 = 1) / sum(w)
  ggplot2::ggplot(mapping = ggplot2::aes(x = s, y = f)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(ggplot2::aes(x = grid, y = p, color = "red")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x     = "x",
      y     = "p",
      title = paste0("Proportions taken with alpha = ", alpha)
    )
}


# ---- B) Count-simulation layer ---------------------------------------------

#' Simulate one vector of counts from Multinomial(n, p).
#'
#' @param p  True proportion vector (length K, sums to 1).
#' @param n  Total sample size (positive integer).
#'
#' @return Integer vector of length K summing to n.
simulate_counts_multinomial <- function(p, n) {
  stopifnot(is.numeric(p), all(p >= 0), abs(sum(p) - 1) < 1e-10)
  stopifnot(is.numeric(n), length(n) == 1L, n >= 1L)
  as.integer(rmultinom(1L, size = n, prob = p))
}

# Placeholder: to be implemented when overdispersion support is added.
# simulate_counts_dirichlet_multinomial <- function(p, n, concentration, ...) {
#   stop("Dirichlet-multinomial not yet implemented.")
# }

# Placeholder: to be implemented when correlation support is added.
# simulate_counts_logistic_normal_multinomial <- function(p, n, Sigma, ...) {
#   stop("Logistic-normal multinomial not yet implemented.")
# }

#' Dispatcher: simulate counts from the requested model.
#'
#' @param p      True proportion vector.
#' @param n      Total sample size.
#' @param model  Sampling model; currently only "multinomial" is implemented.
#' @param ...    Additional arguments forwarded to the concrete simulator
#'               (reserved for future overdispersed / correlated models).
#'
#' @return Integer vector of length K summing to n.
simulate_counts <- function(p, n,
                            model = c("multinomial",
                                      "dirichlet_multinomial",
                                      "logistic_normal_multinomial"),
                            ...) {
  model <- match.arg(model)
  switch(model,
    multinomial = simulate_counts_multinomial(p, n),
    dirichlet_multinomial = stop(
      "model = 'dirichlet_multinomial' is not yet implemented."
    ),
    logistic_normal_multinomial = stop(
      "model = 'logistic_normal_multinomial' is not yet implemented."
    )
  )
}

#' Convert count vector to observed proportions.
#'
#' @param y  Integer count vector (length K, sums to n).
#' @param n  Total sample size; defaults to sum(y).
#'
#' @return Numeric vector of observed proportions summing to 1.
counts_to_proportions <- function(y, n = sum(y)) {
  stopifnot(is.numeric(y) || is.integer(y), n > 0)
  y / n
}


# ---- C) Error-metric layer -------------------------------------------------

#' Compute per-cell-type error vectors for each requested metric.
#'
#' @param phat    Observed proportion vector (length K).
#' @param p       True proportion vector (length K).
#' @param metrics Character vector; any subset of c("AE", "ARE").
#'
#' @return Named list with one numeric vector per metric (length K).
#'
#' @details
#' AE  = abs(phat - p)
#' ARE = abs(phat - p) / p  (no epsilon stabilisation; NaN/Inf for p == 0 is expected)
compute_errors <- function(phat, p, metrics = c("AE", "ARE")) {
  stopifnot(length(phat) == length(p))
  result <- list()
  if ("AE" %in% metrics) {
    result[["AE"]] <- abs(phat - p)
  }
  if ("ARE" %in% metrics) {
    result[["ARE"]] <- abs(phat - p) / p   # NaN when p == 0 and phat == 0; Inf when p == 0 and phat != 0; no stabilisation by design
  }
  result
}


# ---- D) Max-error extraction with configurable tie handling ----------------

#' Extract the maximum error value and its index from one error vector.
#'
#' @param error_vec  Numeric error vector (length K).
#' @param tie_method How to break ties among indices sharing the maximum:
#'   \describe{
#'     \item{"random"}{sample uniformly among tied indices (default)}
#'     \item{"first"}{smallest tied index}
#'     \item{"last"}{largest tied index}
#'   }
#'
#' @return Named list: `max_error_value` (numeric scalar) and `argmax_index` (integer scalar).
max_error_summary <- function(error_vec, tie_method = c("random", "first", "last")) {
  tie_method <- match.arg(tie_method)
  max_val <- max(error_vec, na.rm = TRUE)
  tied <- which(error_vec == max_val)          # always a proper integer vector
  argmax_idx <- switch(tie_method,
    random = tied[sample.int(length(tied), 1L)],   # safe even when length(tied)==1
    first  = tied[[1L]],
    last   = tied[[length(tied)]]
  )
  list(max_error_value = max_val, argmax_index = argmax_idx)
}


# ---- E) Replicate runner ---------------------------------------------------

#' Run B simulation replicates and store max errors + argmax indices.
#'
#' Efficiency strategy: simulate B times once, store only the per-replicate
#' max error values and argmax indices.  Threshold evaluation is done
#' post-hoc by `evaluate_thresholds()` without re-simulating.
#'
#' @param p          True proportion vector (length K).
#' @param n          Total sample size.
#' @param B          Number of replicates.
#' @param metrics    Error metrics to compute (subset of c("AE", "ARE")).
#' @param model      Sampling model passed to `simulate_counts()`.
#' @param tie_method Tie-breaking rule passed to `max_error_summary()`.
#' @param seed       Optional integer seed for reproducibility.
#' @param ...        Additional arguments forwarded to `simulate_counts()`.
#'
#' @return List with elements:
#'   \describe{
#'     \item{max_errors}{B x M numeric matrix of max error values.}
#'     \item{argmax}{B x M integer matrix of argmax indices.}
#'     \item{errors}{B x K x M numeric array of per-cell-type errors.}
#'     \item{inputs}{Copy of all input arguments (including seed used).}
#'   }
run_replicates <- function(p, n, B,
                           metrics = c("AE", "ARE"),
                           model = "multinomial",
                           tie_method = "random",
                           seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  K <- length(p)
  stopifnot(K >= 1L, B >= 1L)

  max_errors <- matrix(NA_real_,    nrow = B, ncol = length(metrics),
                       dimnames = list(NULL, metrics))
  argmax     <- matrix(NA_integer_, nrow = B, ncol = length(metrics),
                       dimnames = list(NULL, metrics))
  errors     <- array(NA_real_, dim = c(B, K, length(metrics)),
                      dimnames = list(NULL, seq_len(K), metrics))

  for (b in seq_len(B)) {
    y_b      <- simulate_counts(p, n, model = model, ...)
    phat_b   <- counts_to_proportions(y_b, n)
    errors_b <- compute_errors(phat_b, p, metrics = metrics)

    for (m in metrics) {
      s <- max_error_summary(errors_b[[m]], tie_method = tie_method)
      max_errors[b, m] <- s$max_error_value
      argmax[b, m]     <- s$argmax_index
      errors[b, , m]   <- errors_b[[m]]
    }
  }

  list(
    max_errors = max_errors,
    argmax     = argmax,
    errors     = errors,
    inputs     = list(p = p, n = n, B = B, metrics = metrics,
                      model = model, tie_method = tie_method, seed = seed)
  )
}


# ---- F) Threshold evaluation -----------------------------------------------

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
      if (is.null(dim(errors_m))) errors_m <- matrix(errors_m, ncol = 1L)
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


# ---- G) Argmax summarisation -----------------------------------------------

#' Tabulate which cell-type index most often achieved the max error.
#'
#' @param argmax  B x M integer matrix of argmax indices (from `run_replicates()`).
#' @param p       True proportion vector (length K).
#'
#' @return Tidy data.frame with columns: metric, index, count, fraction, p_value.
summarize_argmax <- function(argmax, p) {
  stopifnot(is.matrix(argmax), !is.null(colnames(argmax)))
  K <- length(p)
  metrics <- colnames(argmax)
  rows <- vector("list", length(metrics))
  for (i in seq_along(metrics)) {
    m <- metrics[[i]]
    counts <- tabulate(argmax[, m], nbins = K)
    rows[[i]] <- data.frame(
      metric   = m,
      index    = seq_len(K),
      count    = counts,
      fraction = counts / sum(counts),
      p_value  = p,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


# ---- H) Visualisation helpers ----------------------------------------------

#' Plot success-rate curves from a simulation result object.
#'
#' @param result  Output list from `run_simulation_experiment()`.
#' @param metric  Character scalar. If NULL, plot all metrics (faceted).
#' @param alphas  Optional numeric vector; subset of alpha values to plot.
#' @param target  Success-rate reference line drawn as a horizontal dotted line
#'                (default 0.95).
#'
#' @return A ggplot object.
plot_success_rate_curve <- function(result, metric = NULL, alphas = NULL, target = 0.95) {
  stopifnot(is.list(result), "curves" %in% names(result))
  df <- result$curves
  if (!is.null(metric)) {
    df <- df[df$metric %in% metric, , drop = FALSE]
  }
  if (!is.null(alphas)) {
    df <- df[df$alpha %in% alphas, , drop = FALSE]
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = tau, y = success_rate, color = factor(alpha))) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = target, linetype = "dotted") +
    ggplot2::labs(
      x     = "Threshold (tau)",
      y     = "Success rate",
      color = "alpha",
      title = "Success-rate curve(s)"
    ) +
    ggplot2::theme_bw()
  if (is.null(metric) || length(unique(df$metric)) > 1L) {
    p <- p + ggplot2::facet_wrap(~metric, scales = "free_x")
  }
  p
}

#' Plot a histogram of which cell-type proportion drives the maximum error.
#'
#' Maps argmax indices to the true proportions so that the x-axis shows
#' actual proportion values rather than integer indices.
#'
#' @param argmax  B x M integer matrix of argmax indices (from `run_replicates()`).
#' @param p       True proportion vector (length K).
#' @param metric  Character scalar; which metric to plot ("AE" or "ARE").
#'
#' @return A ggplot object.
plot_argmax_histogram <- function(argmax, p, metric) {
  idx_values <- argmax[, metric]
  df <- data.frame(index = factor(idx_values, levels = seq_along(p)))
  ggplot2::ggplot(df, ggplot2::aes(x = index)) +
    ggplot2::geom_bar() +
    ggplot2::scale_x_discrete(labels = round(p, 4)) +
    ggplot2::labs(
      x     = "True proportion",
      y     = "Count",
      title = paste0("Argmax distribution: ", metric)
    )
}


# ---- I) High-level orchestration -------------------------------------------

#' Run the full simulation experiment end-to-end.
#'
#' @param alpha      Numeric vector; one or more shape1 values for
#'   Beta(alpha, 1) proportion generation.
#' @param K          Number of cell types (default 10).
#' @param n          Total sample size per replicate.
#' @param B          Number of replicates.
#' @param taus       Numeric vector of thresholds (same for all metrics) or a
#'   named list with one numeric vector per metric
#'   (e.g. `list(AE = c(...), ARE = c(...))`).
#' @param metrics    Error metrics; subset of c("AE", "ARE").
#' @param model      Sampling model (currently only "multinomial").
#' @param tie_method Tie-breaking rule for max-error argmax.
#' @param seed       Optional integer seed for reproducibility.
#' @param ...        Additional arguments forwarded to `simulate_counts()`.
#'
#' @return List with elements:
#'   \describe{
#'     \item{inputs}{All input arguments.}
#'     \item{p_table}{Tidy data.frame: alpha, index, p_value.}
#'     \item{replicate_summaries}{Tidy data.frame:
#'       alpha, replicate, metric, max_error, argmax_index.}
#'     \item{errors_long}{Tidy data.frame:
#'       alpha, replicate, metric, index, error.}
#'     \item{curves}{Tidy data.frame:
#'       alpha, metric, tau, success_rate, mean_n_above.}
#'     \item{argmax_summary}{Tidy data.frame:
#'       alpha, metric, index, count, fraction, p_value.}
#'   }
run_simulation_experiment <- function(alpha, K = 10, n, B, taus,
                                      metrics = c("AE", "ARE"),
                                      model = "multinomial",
                                      tie_method = "random",
                                      seed = NULL, ...) {
  stopifnot(is.numeric(alpha), length(alpha) >= 1L, all(alpha > 0))

  p_table_list <- vector("list", length(alpha))
  replicate_summaries_list <- vector("list", length(alpha))
  errors_long_list <- vector("list", length(alpha))
  curves_list <- vector("list", length(alpha))
  argmax_summary_list <- vector("list", length(alpha))

  for (i in seq_along(alpha)) {
    alpha_i <- alpha[[i]]
    seed_i <- if (is.null(seed)) NULL else seed + i - 1L
    p <- generate_proportions_beta(alpha = alpha_i, K = K)
    rep_out <- run_replicates(
      p, n, B, metrics = metrics, model = model,
      tie_method = tie_method, seed = seed_i, ...
    )

    p_table_list[[i]] <- data.frame(
      alpha = alpha_i,
      index = seq_len(K),
      p_value = p,
      stringsAsFactors = FALSE
    )

    replicate_summaries_list[[i]] <- data.frame(
      alpha = alpha_i,
      replicate = rep(seq_len(B), times = length(metrics)),
      metric = rep(metrics, each = B),
      max_error = as.vector(rep_out$max_errors[, metrics, drop = FALSE]),
      argmax_index = as.vector(rep_out$argmax[, metrics, drop = FALSE]),
      stringsAsFactors = FALSE
    )

    errors_m_list <- vector("list", length(metrics))
    for (j in seq_along(metrics)) {
      m <- metrics[[j]]
      errors_m <- rep_out$errors[, , m, drop = TRUE]
      if (is.null(dim(errors_m))) errors_m <- matrix(errors_m, ncol = 1L)
      errors_m_list[[j]] <- data.frame(
        alpha = alpha_i,
        replicate = rep(seq_len(B), times = ncol(errors_m)),
        metric = m,
        index = rep(seq_len(ncol(errors_m)), each = B),
        error = as.vector(errors_m),
        stringsAsFactors = FALSE
      )
    }
    errors_long_list[[i]] <- do.call(rbind, errors_m_list)

    curves_i <- evaluate_thresholds(rep_out$max_errors, taus, errors = rep_out$errors)
    curves_i$alpha <- alpha_i
    curves_list[[i]] <- curves_i[, c("alpha", "metric", "tau", "success_rate", "mean_n_above")]

    argmax_i <- summarize_argmax(rep_out$argmax, p)
    argmax_i$alpha <- alpha_i
    argmax_summary_list[[i]] <- argmax_i[, c("alpha", "metric", "index", "count", "fraction", "p_value")]
  }

  list(
    inputs = list(alpha = alpha, K = K, n = n, B = B, taus = taus,
                  metrics = metrics, model = model,
                  tie_method = tie_method, seed = seed),
    p_table = do.call(rbind, p_table_list),
    replicate_summaries = do.call(rbind, replicate_summaries_list),
    errors_long = do.call(rbind, errors_long_list),
    curves = do.call(rbind, curves_list),
    argmax_summary = do.call(rbind, argmax_summary_list)
  )
}
