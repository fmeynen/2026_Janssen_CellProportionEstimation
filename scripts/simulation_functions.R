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

#' Generate K true proportions from dbeta(grid, 1, alpha), always rescaled.
#'
#' @param alpha  shape2 parameter of Beta(1, alpha); controls curve steepness.
#' @param K      number of cell types (default 10).
#' @param grid   evaluation points in (0,1) (default K equidistant points from
#'               0.1 to 0.9, avoiding boundary values where dbeta returns 0).
#'
#' @return Numeric vector of length K, all strictly positive, summing to 1.
#'
#' @details
#' Unscaled weights: w = dbeta(grid, shape1 = 1, shape2 = alpha).
#' The default grid avoids 0 and 1 so that all weights — and therefore all
#' proportions — are strictly positive.  This prevents ARE from producing
#' NaN or Inf values for any cell type under the default parameters.
generate_proportions_beta <- function(alpha, K = 10, grid = seq(0.1, 0.9, length.out = K)) {
  stopifnot(is.numeric(alpha), length(alpha) == 1L, alpha > 0)
  stopifnot(length(grid) == K)
  w <- dbeta(grid, shape1 = 1, shape2 = alpha)
  p <- normalize_to_simplex(w)
  validate_proportions(p)
  p
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

  for (b in seq_len(B)) {
    y_b      <- simulate_counts(p, n, model = model, ...)
    phat_b   <- counts_to_proportions(y_b, n)
    errors_b <- compute_errors(phat_b, p, metrics = metrics)

    for (m in metrics) {
      s <- max_error_summary(errors_b[[m]], tie_method = tie_method)
      max_errors[b, m] <- s$max_error_value
      argmax[b, m]     <- s$argmax_index
    }
  }

  list(
    max_errors = max_errors,
    argmax     = argmax,
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
#' @param taus       Numeric vector of threshold values.
#'
#' @return Tidy data.frame with columns: metric, tau, success_rate.
evaluate_thresholds <- function(max_errors, taus) {
  stopifnot(is.matrix(max_errors), !is.null(colnames(max_errors)))
  metrics <- colnames(max_errors)
  rows <- vector("list", length(metrics))
  for (i in seq_along(metrics)) {
    m <- metrics[[i]]
    rates <- vapply(taus, function(tau) mean(max_errors[, m] <= tau, na.rm = TRUE), numeric(1L))
    rows[[i]] <- data.frame(metric = m, tau = taus, success_rate = rates,
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


# ---- H) High-level orchestration -------------------------------------------

#' Run the full simulation experiment end-to-end.
#'
#' @param alpha      shape2 parameter of Beta(1, alpha) for proportion generation.
#' @param K          Number of cell types (default 10).
#' @param n          Total sample size per replicate.
#' @param B          Number of replicates.
#' @param taus       Numeric vector of thresholds for success-rate curves.
#' @param metrics    Error metrics; subset of c("AE", "ARE").
#' @param model      Sampling model (currently only "multinomial").
#' @param tie_method Tie-breaking rule for max-error argmax.
#' @param seed       Optional integer seed for reproducibility.
#' @param ...        Additional arguments forwarded to `simulate_counts()`.
#'
#' @return List with elements:
#'   \describe{
#'     \item{inputs}{All input arguments.}
#'     \item{p}{True proportion vector (length K).}
#'     \item{max_errors}{B x M matrix of max error values.}
#'     \item{argmax}{B x M matrix of argmax indices.}
#'     \item{curves}{Tidy data.frame: metric, tau, success_rate.}
#'     \item{argmax_summary}{Tidy data.frame: metric, index, count, fraction, p_value.}
#'   }
run_simulation_experiment <- function(alpha, K = 10, n, B, taus,
                                      metrics = c("AE", "ARE"),
                                      model = "multinomial",
                                      tie_method = "random",
                                      seed = NULL, ...) {
  p               <- generate_proportions_beta(alpha = alpha, K = K)
  rep_out         <- run_replicates(p, n, B, metrics = metrics, model = model,
                                    tie_method = tie_method, seed = seed, ...)
  curves          <- evaluate_thresholds(rep_out$max_errors, taus)
  argmax_summary  <- summarize_argmax(rep_out$argmax, p)

  list(
    inputs         = list(alpha = alpha, K = K, n = n, B = B, taus = taus,
                          metrics = metrics, model = model,
                          tie_method = tie_method, seed = seed),
    p              = p,
    max_errors     = rep_out$max_errors,
    argmax         = rep_out$argmax,
    curves         = curves,
    argmax_summary = argmax_summary
  )
}
