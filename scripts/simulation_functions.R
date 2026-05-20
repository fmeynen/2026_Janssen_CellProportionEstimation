# simulation_functions.R
#
# Functions for multinomial sampling error simulation.
# Implements the MVP workflow:
#   1. Generate true proportions from a configurable method-based generator
#      (currently Beta and fixed-max Beta, with placeholders for future extensions)
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

#' Generate K true proportions from dbeta(grid, alpha, 1), always rescaled.
#'
#' @param alpha  shape1 parameter of Beta(alpha, 1); controls curve steepness.
#' @param K      number of cell types (default 10).
#' @param grid   evaluation points in (0,1) (default K equidistant points from
#'               0.05 to 0.95, avoiding boundary values where dbeta returns 0).
#'
#' @return Numeric vector of length K, all strictly positive, summing to 1.
#'
#' @details
#' Unscaled weights: w = dbeta(grid, shape1 = alpha, shape2 = 1).
#' The default grid avoids 0 and 1 so that all weights — and therefore all
#' proportions — are strictly positive.  This prevents ARE from producing
#' NaN or Inf values for any cell type under the default parameters.
generate_proportions_beta <- function(alpha, K = 10, grid = default_beta_grid(K)) {
  stopifnot(is.numeric(alpha), length(alpha) == 1L, alpha > 0)
  stopifnot(length(grid) == K)
  w <- dbeta(grid, shape1 = alpha, shape2 = 1)
  p <- normalize_to_simplex(w)
  validate_proportions(p)
  p
}

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
  stop(
    fixed_max_beta_impossible_error,
    call. = FALSE
  )
}

#' Generate K true proportions with a fixed maximum at the highest index.
#'
#' @param alpha  shape1 parameter of the Beta(alpha, 1) remainder curve.
#' @param K      number of cell types (default 10; must be at least 2).
#' @param p_max  fixed largest true proportion, placed at the highest index.
#' @param grid   evaluation points in (0,1), length K - 1, used to construct
#'               the Beta-shaped remainder over the first K - 1 indices.
#'
#' @return Numeric vector of length K, all strictly positive, summing to 1,
#'   with a strictly unique largest value at index K.
#'
#' @details
#' The first `K - 1` proportions are built from Beta(alpha, 1) weights,
#' normalized and then rescaled to sum to `1 - p_max`. The final proportion is
#' set to `p_max`, so the largest true proportion is fixed at the highest
#' index. If the rescaled remainder contains any value `>= p_max`, the
#' combination of `alpha`, `K`, and `p_max` is impossible for a strictly unique
#' fixed maximum; the function warns and then fails.
generate_proportions_fixed_max_beta <- function(alpha, K = 10, p_max,
                                                grid = default_beta_grid(K - 1L)) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
    stop("alpha must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K %% 1 != 0 || K < 2L) {
    stop("K must be a single integer >= 2 for method = 'fixed_max_beta'.", call. = FALSE)
  }
  if (is.null(p_max)) {
    stop("p_max must be provided when method = 'fixed_max_beta'.", call. = FALSE)
  }
  if (!is.numeric(p_max) || length(p_max) != 1L || !is.finite(p_max) || p_max <= 0 || p_max >= 1) {
    stop("p_max must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  if (length(grid) != K - 1L) {
    stop("grid must have length K - 1 for method = 'fixed_max_beta'.", call. = FALSE)
  }

  remainder_weights <- dbeta(grid, shape1 = alpha, shape2 = 1)
  remainder <- (1 - p_max) * normalize_to_simplex(remainder_weights)

  if (any(remainder >= p_max)) {
    fail_fixed_max_beta_impossible(
      non_max = remainder,
      alpha = alpha,
      K = K,
      p_max = p_max
    )
  }

  p <- c(remainder, p_max)
  validate_proportions(p)
  p
}

#' Dispatcher: generate true proportions from the requested method.
#'
#' @param alpha   Shape parameter used by the requested generation method.
#' @param K       Number of cell types (default 10).
#' @param method  Proportion-generation method: `"beta"` or `"fixed_max_beta"`.
#' @param p_max   Fixed largest true proportion for `"fixed_max_beta"`. The
#'   largest value is always placed at the highest index and must remain
#'   strictly unique; impossible combinations warn and fail.
#' @param grid    Evaluation points in (0,1), length K (used by `"beta"`).
#'
#' @return Numeric vector of length K, all strictly positive, summing to 1.
generate_proportions <- function(alpha, K = 10,
                                 method = c("beta", "fixed_max_beta"),
                                 p_max = NULL,
                                 grid = default_beta_grid(K)) {
  method <- match.arg(method)
  switch(method,
    beta = generate_proportions_beta(alpha = alpha, K = K, grid = grid),
    fixed_max_beta = generate_proportions_fixed_max_beta(
      alpha = alpha,
      K = K,
      p_max = p_max
    )
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

#' Plot normalized Beta(alpha, 1) proportion curves.
#'
#' @param alpha Numeric vector of one or more positive alpha values.
#' @param K Number of cell types (default 10).
#' @param grid Evaluation points in (0,1), length K.
#'
#' @return A ggplot object with one facet per alpha value, using `facet_grid`.
plot_proportions_curve <- function(alpha, K = 10, grid = seq(0.05, 0.95, length.out = K)) {
  stopifnot(is.numeric(alpha), length(alpha) >= 1L, all(alpha > 0))
  stopifnot(length(grid) == K)
  
  s <- seq(0, 1, length.out = 1000)
  
  curve_rows <- vector("list", length(alpha))
  point_rows <- vector("list", length(alpha))
  
  for (i in seq_along(alpha)) {
    alpha_i <- alpha[[i]]
    w <- dbeta(grid, shape1 = alpha_i, shape2 = 1)
    p <- normalize_to_simplex(w)
    
    curve_rows[[i]] <- data.frame(
      alpha = alpha_i,
      s = s,
      f = dbeta(s, shape1 = alpha_i, shape2 = 1) / sum(w),
      stringsAsFactors = FALSE
    )
    
    point_rows[[i]] <- data.frame(
      alpha = alpha_i,
      grid = grid,
      p = p,
      stringsAsFactors = FALSE
    )
  }
  
  curve_df <- do.call(rbind, curve_rows)
  point_df <- do.call(rbind, point_rows)
  alpha_levels <- alpha[!duplicated(alpha)]
  curve_df$alpha <- factor(curve_df$alpha, levels = alpha_levels)
  point_df$alpha <- factor(point_df$alpha, levels = alpha_levels)
  
  ggplot2::ggplot(curve_df, ggplot2::aes(x = s, y = f)) +
    ggplot2::geom_line(linewidth = 0.7, color = "#2C3E50") +
    ggplot2::geom_point(
      data = point_df,
      ggplot2::aes(x = grid, y = p),
      inherit.aes = FALSE,
      size = 1.8,
      color = "#D62728"
    ) +
    ggplot2::facet_grid(cols = ggplot2::vars(alpha), scales = "fixed") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(
      x = "x",
      y = "p",
      title = "Proportions taken with Beta(alpha, 1)"
    )
}


#' Plot success-rate curves from a simulation result object.
#'
#' @param result  Output list from `run_simulation_experiment()`.
#' @param metric  Character scalar. If NULL, plot all metrics (faceted).
#' @param alphas  Optional numeric vector; subset of alpha values to plot.
#' @param p_maxs  Optional numeric vector; subset of p_max values to plot.
#' @param target  Success-rate reference line drawn as a horizontal dotted line
#'                (default 0.95).
#'
#' @return A ggplot object.
plot_success_rate_curve <- function(result, metric = NULL, alphas = NULL,
                                    p_maxs = NULL, target = 0.95) {
  stopifnot(is.list(result), "curves" %in% names(result))
  df <- result$curves
  if (!is.null(metric)) {
    df <- df[df$metric %in% metric, , drop = FALSE]
  }
  if (!is.null(alphas)) {
    df <- df[df$alpha %in% alphas, , drop = FALSE]
  }
  if (!is.null(p_maxs)) {
    if (!("p_max" %in% names(df))) {
      stop("No p_max metadata found in result$curves.", call. = FALSE)
    }
    df <- df[df$p_max %in% p_maxs, , drop = FALSE]
  }
  if (nrow(df) == 0L) {
    stop("No rows match the requested metric/alpha/p_max filter(s).", call. = FALSE)
  }

  has_p_max <- "p_max" %in% names(df) && any(!is.na(df$p_max))
  aes_mapping <- ggplot2::aes(x = tau, y = success_rate, color = factor(alpha))
  if (has_p_max) {
    aes_mapping <- ggplot2::aes(
      x = tau, y = success_rate, color = factor(alpha),
      linetype = factor(p_max), group = interaction(alpha, p_max)
    )
  }

  p <- ggplot2::ggplot(df, aes_mapping) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = target, linetype = "dotted") +
    ggplot2::labs(
      x     = "Threshold (tau)",
      y     = "Success rate",
      color = "alpha",
      title = "Success-rate curve(s)"
    ) +
    ggplot2::theme_bw()
  if (has_p_max) {
    p <- p + ggplot2::labs(linetype = "p_max")
  }
  if (is.null(metric) || length(unique(df$metric)) > 1L) {
    p <- p + ggplot2::facet_wrap(~metric, scales = "free_x")
  }
  p
}

#' Plot a histogram of which cell-type proportion drives the maximum error.
#'
#' @param result  Output list from `run_simulation_experiment()`.
#' @param metric  Character scalar. If NULL, plot all metrics (faceted by row).
#' @param alphas  Optional numeric vector; subset of alpha values to plot.
#' @param p_maxs  Optional numeric vector; subset of p_max values to plot.
#'
#' @return A ggplot object.
plot_argmax_histogram <- function(result, metric = NULL, alphas = NULL, p_maxs = NULL) {
  stopifnot(
    "result must be a list" = is.list(result),
    "result must contain replicate_summaries" = "replicate_summaries" %in% names(result)
  )
  df <- result$replicate_summaries
  if (!is.null(metric)) {
    if (!any(df$metric %in% metric)) {
      stop("No rows match the requested metric value(s).", call. = FALSE)
    }
    df <- df[df$metric %in% metric, , drop = FALSE]
  }
  if (!is.null(alphas)) {
    if (!any(df$alpha %in% alphas)) {
      stop("No rows match the requested alpha value(s).", call. = FALSE)
    }
    df <- df[df$alpha %in% alphas, , drop = FALSE]
  }
  if (!is.null(p_maxs)) {
    if (!("p_max" %in% names(df))) {
      stop("No p_max metadata found in result$replicate_summaries.", call. = FALSE)
    }
    if (!any(df$p_max %in% p_maxs)) {
      stop("No rows match the requested p_max value(s).", call. = FALSE)
    }
    df <- df[df$p_max %in% p_maxs, , drop = FALSE]
  }
  has_p_max <- "p_max" %in% names(df) && any(!is.na(df$p_max))

  argmax_plot <- ggplot2::ggplot(df, ggplot2::aes(x = argmax_index)) +
    ggplot2::geom_bar() +
    ggplot2::labs(
      x     = "Cell-type index (argmax)",
      y     = "Count",
      title = "Distribution of maximum-error indices"
    ) +
    ggplot2::theme_bw()

  facet_cols <- if (has_p_max) ggplot2::vars(alpha, p_max) else ggplot2::vars(alpha)
  if (length(unique(df$metric)) > 1L) {
    argmax_plot <- argmax_plot + ggplot2::facet_grid(rows = ggplot2::vars(metric), cols = facet_cols)
  } else {
    argmax_plot <- argmax_plot + ggplot2::facet_grid(cols = facet_cols)
  }
  argmax_plot
}


# ---- I) High-level orchestration -------------------------------------------

#' Run the full simulation experiment end-to-end.
#'
#' @param alpha      Numeric vector; one or more positive shape values used by
#'   the selected proportion-generation method (default method is Beta-based).
#' @param K          Number of cell types (default 10).
#' @param n          Total sample size per replicate.
#' @param B          Number of replicates.
#' @param taus       Numeric vector of thresholds (same for all metrics) or a
#'   named list with one numeric vector per metric
#'   (e.g. `list(AE = c(...), ARE = c(...))`).
#' @param metrics    Error metrics; subset of c("AE", "ARE").
#' @param proportion_method Proportion-generation method (`"beta"` or
#'   `"fixed_max_beta"`). The fixed-max Beta method places `p_max` at the
#'   highest index and warns then fails for impossible combinations.
#' @param p_max      Fixed largest true proportion(s) used by
#'   `proportion_method = "fixed_max_beta"`. Can be a numeric vector.
#'   When multiple values are provided, all alpha × p_max combinations are
#'   attempted; impossible fixed-max combinations are warned and skipped.
#' @param model      Sampling model (currently only "multinomial").
#' @param tie_method Tie-breaking rule for max-error argmax.
#' @param seed       Optional integer seed for reproducibility.
#' @param ...        Additional arguments forwarded to `simulate_counts()`.
#'
#' @return List with elements:
#'   \describe{
#'     \item{inputs}{All input arguments.}
#'     \item{p_table}{Data.frame with one row per simulated alpha/p_max combination and one column per index
#'       (`index_1`, ..., `index_K`) containing the corresponding p values,
#'       plus a `p_max` column.}
#'     \item{replicate_summaries}{Tidy data.frame:
#'       alpha, p_max, replicate, metric, max_error, argmax_index.}
#'     \item{errors_long}{Tidy data.frame:
#'       alpha, p_max, replicate, metric, index, error.}
#'     \item{curves}{Tidy data.frame:
#'       alpha, p_max, metric, tau, success_rate, mean_n_above.}
#'     \item{argmax_summary}{Tidy data.frame:
#'       alpha, p_max, metric, index, count, fraction, p_value.}
#'   }
run_simulation_experiment <- function(alpha, K = 10, n, B, taus,
                                      metrics = c("AE", "ARE"),
                                      proportion_method = "beta",
                                      p_max = NULL,
                                      model = "multinomial",
                                      tie_method = "random",
                                      seed = NULL, ...) {
  stopifnot(is.numeric(alpha), length(alpha) >= 1L, all(alpha > 0))

  if (identical(proportion_method, "fixed_max_beta")) {
    if (is.null(p_max)) {
      stop("p_max must be provided when proportion_method = 'fixed_max_beta'.", call. = FALSE)
    }
    if (!is.numeric(p_max) || any(!is.finite(p_max)) || any(p_max <= 0) || any(p_max >= 1)) {
      stop("p_max must contain numbers strictly between 0 and 1.", call. = FALSE)
    }
    p_max_values <- as.numeric(p_max)
  } else {
    p_max_values <- NA_real_
  }

  combinations <- expand.grid(
    alpha = alpha,
    p_max = p_max_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  n_combinations <- nrow(combinations)
  p_table_list <- vector("list", n_combinations)
  replicate_summaries_list <- vector("list", n_combinations)
  errors_long_list <- vector("list", n_combinations)
  curves_list <- vector("list", n_combinations)
  argmax_summary_list <- vector("list", n_combinations)
  keep <- logical(n_combinations)
  skip_impossible <- identical(proportion_method, "fixed_max_beta") &&
    length(p_max_values) > 1L

  for (i in seq_len(n_combinations)) {
    alpha_i <- combinations$alpha[[i]]
    p_max_i <- combinations$p_max[[i]]
    seed_i <- if (is.null(seed)) NULL else seed + i - 1L
    p <- tryCatch(
      generate_proportions(
        alpha = alpha_i,
        K = K,
        method = proportion_method,
        p_max = if (is.na(p_max_i)) NULL else p_max_i
      ),
      error = function(e) {
        is_impossible_fixed_max <- grepl(
          fixed_max_beta_impossible_error,
          conditionMessage(e),
          fixed = TRUE
        )
        if (skip_impossible && is_impossible_fixed_max) {
          return(NULL)
        }
        stop(e)
      }
    )
    if (is.null(p)) next # Skip impossible alpha/p_max combinations in vectorized p_max runs.

    rep_out <- run_replicates(
      p, n, B, metrics = metrics, model = model,
      tie_method = tie_method, seed = seed_i, ...
    )
    keep[[i]] <- TRUE

    p_table_list[[i]] <- data.frame(
      alpha = alpha_i,
      p_max = p_max_i,
      as.list(stats::setNames(as.numeric(p), paste0("index_", seq_len(K)))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    replicate_summaries_list[[i]] <- data.frame(
      alpha = alpha_i,
      p_max = p_max_i,
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
      if (is.null(dim(errors_m))) errors_m <- matrix(errors_m, nrow = B)
      errors_m_list[[j]] <- data.frame(
        alpha = alpha_i,
        p_max = p_max_i,
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
    curves_i$p_max <- p_max_i
    curves_list[[i]] <- curves_i[, c("alpha", "p_max", "metric", "tau", "success_rate", "mean_n_above")]

    argmax_i <- summarize_argmax(rep_out$argmax, p)
    argmax_i$alpha <- alpha_i
    argmax_i$p_max <- p_max_i
    argmax_summary_list[[i]] <- argmax_i[, c("alpha", "p_max", "metric", "index", "count", "fraction", "p_value")]
  }

  if (!any(keep)) {
    stop("No feasible alpha/p_max combinations produced simulation output.", call. = FALSE)
  }

  p_table_list <- p_table_list[keep]
  replicate_summaries_list <- replicate_summaries_list[keep]
  errors_long_list <- errors_long_list[keep]
  curves_list <- curves_list[keep]
  argmax_summary_list <- argmax_summary_list[keep]

  list(
    inputs = list(alpha = alpha, K = K, n = n, B = B, taus = taus,
                  metrics = metrics, proportion_method = proportion_method,
                  p_max = p_max,
                  model = model,
                  tie_method = tie_method, seed = seed),
    p_table = do.call(rbind, p_table_list),
    replicate_summaries = do.call(rbind, replicate_summaries_list),
    errors_long = do.call(rbind, errors_long_list),
    curves = do.call(rbind, curves_list),
    argmax_summary = do.call(rbind, argmax_summary_list)
  )
}
