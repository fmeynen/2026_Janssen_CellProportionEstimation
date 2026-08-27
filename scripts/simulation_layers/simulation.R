# Simulation layer ---------------------------------------------------------------------------
# Simulation layer: generate true proportions, simulate counts, run replicates.

# Generate proportions --------------------------------------------------------------------------------------------

#' Generate K true proportions from dbeta(grid, alpha, 1), always rescaled.
#'
#' @param alpha  shape1 parameter of Beta(alpha, 1); controls curve steepness.
#' @param K      number of cell types (default 10).
#' @param grid   evaluation points in (0,1) (default K equidistant points from 0.05 to 0.95, avoiding boundary values
#'               where dbeta returns 0).
#'
#' @return Numeric vector of length K, all strictly positive, summing to 1.
#'
#' @details
#' Unscaled weights: w = dbeta(grid, shape1 = alpha, shape2 = 1). The default grid avoids 0 and 1 so that all weights —
#' and therefore all proportions — are strictly positive.
generate_proportions_beta <- function(alpha, K = 10, grid = default_beta_grid(K)) {
  stopifnot(is.numeric(alpha), length(alpha) == 1L, alpha > 0)
  stopifnot(length(grid) == K)
  w <- dbeta(grid, shape1 = alpha, shape2 = 1)
  p <- normalize_to_simplex(w)
  validate_proportions(p)
  p
}

#' Generate K true proportions with a fixed maximum at the highest index.
#'
#' @param alpha  shape1 parameter of the Beta(alpha, 1) remainder curve.
#' @param K      number of cell types (default 10; must be at least 2).
#' @param p_max  fixed largest true proportion(s), placed at the highest index.
#' @param grid   evaluation points in (0,1), length K - 1, used to construct the Beta-shaped remainder over the first
#'               K - 1 indices.
#'
#' @return If `length(p_max) == 1`, a numeric vector of length K, all strictly positive, summing to 1, with a strictly
#'         unique largest value at index K.
#'         If `length(p_max) > 1`, a numeric matrix with one row per `p_max` value and K columns
#'         (`index_1`, ..., `index_K`).
#'
#' @details
#' The first `K - 1` proportions are built from Beta(alpha, 1) weights, normalized and then rescaled to sum to
#' `1 - p_max`. The final proportion is set to `p_max`, so the largest true proportion is fixed at the highest index.
#' If the rescaled remainder contains any value `>= p_max`, the combination of `alpha`, `K`, and `p_max` is impossible
#' for a strictly unique fixed maximum; the function warns and then fails. When `p_max` contains multiple values, this
#' construction is applied independently per value.
generate_props_fixed_max_beta <- function(alpha, K = 10, p_max,
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
  if (!is.numeric(p_max) || any(!is.finite(p_max)) || any(p_max <= 0) || any(p_max >= 1)) {
    stop("p_max must contain number(s) strictly between 0 and 1.", call. = FALSE)
  }
  if (length(grid) != K - 1L) {
    stop("grid must have length K - 1 for method = 'fixed_max_beta'.", call. = FALSE)
  }
  if (length(p_max) > 1L) {
    p_mat <- t(vapply(
      p_max,
      function(p_max_i) {
        generate_props_fixed_max_beta(alpha = alpha, K = K, p_max = p_max_i, grid = grid)
      },
      FUN.VALUE = numeric(K)
    ))
    colnames(p_mat) <- paste0("index_", seq_len(K))
    rownames(p_mat) <- paste0("p_max_", seq_along(p_max), "_", format(p_max, trim = TRUE))
    return(p_mat)
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
#' @param p_max   Fixed largest true proportion for `"fixed_max_beta"`. The largest value is always placed at the
#'                highest index and must remain strictly unique; impossible combinations warn and fail. May be a numeric
#'                 vector when calling the fixed-max generator directly.
#' @param grid    Evaluation points in (0,1), length K (used by `"beta"`).
#'
#' @return Numeric vector of length K, all strictly positive, summing to 1.
#'   For method `"fixed_max_beta"` with vector `p_max`, returns a numeric matrix with one row per `p_max` value and
#'   K columns.
generate_proportions <- function(alpha, K = 10,
                                 method = c("beta", "fixed_max_beta"),
                                 p_max = NULL,
                                 grid = default_beta_grid(K)) {
  method <- match.arg(method)
  switch(method,
    beta = generate_proportions_beta(alpha = alpha, K = K, grid = grid),
    fixed_max_beta = generate_props_fixed_max_beta(
      alpha = alpha,
      K = K,
      p_max = p_max
    )
  )
}


# Simulate counts -------------------------------------------------------------------------------------------------

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
#' @param ...    Additional arguments forwarded to the concrete simulator (reserved for future overdispersed /
#'               correlated models).
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


# Convert counts to proportions -----------------------------------------------------------------------------------

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
  stopifnot(length(phat) == length(p))
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

# Coordinate Simulation --------------------------------------------------------------------------------------------


#' Run B simulation replicates and store max errors + argmax indices.
#'
#' Efficiency strategy: simulate B times once, store only the per-replicate max error values and argmax indices.
#' Threshold evaluation is done post-hoc by `evaluate_thresholds()` without re-simulating.
#'
#' @param p          True proportion vector (length K).
#' @param n          Total sample size.
#' @param B          Number of replicates.
#' @param metrics    Error metrics to compute; any subset of
#'   `c("AE", "ARE", "TSE", "LAE")`.
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
#'     \item{phat}{B x K numeric matrix of observed proportions.}
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
  phat       <- matrix(NA_real_, nrow = B, ncol = K)

  for (b in seq_len(B)) {
    y_b      <- simulate_counts(p, n, model = model, ...)
    phat_b   <- counts_to_proportions(y_b, n)
    phat[b, ] <- phat_b
    errors_b <- compute_errors(phat_b, p, metrics = metrics, n = n)

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
    phat       = phat,
    inputs     = list(p = p, n = n, B = B, metrics = metrics,
                      model = model, tie_method = tie_method, seed = seed)
  )
}


# Simulate success at sample size ----------------------------------------------------------------------------------



#' Simulate replicates at one sample size and derive per-replicate success.
#'
#' Success is defined as: for every metric that has a corresponding threshold in `taus`, the max error across all K cell
#' types must be at or below that threshold.  When both AE and ARE are requested and both have thresholds, a replicate
#' is successful only if *both* conditions hold simultaneously.
#'
#' @param alpha          Positive scalar; Beta shape parameter used to generate the true proportions.
#' @param n              Total sample size (positive integer).
#' @param config         Named list; must contain at minimum:
#'   \describe{
#'     \item{K}{Number of cell types.}
#'     \item{B}{Number of replicates.}
#'     \item{taus}{Named list with one scalar threshold per metric (e.g.
#'       `list(AE = 0.02, ARE = 0.5)`); scalar thresholds only.}
#'     \item{metrics}{Character vector of metric names to simulate.}
#'     \item{model}{Sampling model (currently `"multinomial"`).}
#'     \item{tie_method}{Tie-breaking rule for max-error argmax.}
#'     \item{proportion_method}{Proportion-generation method.}
#'     \item{seed}{Optional integer RNG seed.}
#'   }
#'
#' @return List with elements:
#'   \describe{
#'     \item{success}{Logical vector of length B.}
#'     \item{success_count}{Integer; number of successful replicates.}
#'     \item{success_rate}{Numeric; fraction of successful replicates.}
#'     \item{rep_out}{Raw output of `run_replicates()`.}
#'   }
simulate_success_at_n <- function(alpha, n, config) {
  p <- generate_proportions(
    alpha  = alpha,
    K      = config$K,
    method = config$proportion_method
  )
  rep_out <- run_replicates(
    p          = p,
    n          = n,
    B          = config$B,
    metrics    = config$metrics,
    model      = config$model,
    tie_method = config$tie_method,
    seed       = config$seed
  )
  max_errors <- rep_out$max_errors
  metrics    <- config$metrics
  taus       <- config$taus

  # Build a B-length logical success vector: replicate passes iff every metric
  # with a threshold has its max error <= that threshold.
  success <- rep(TRUE, config$B)
  for (m in metrics) {
    if (is.null(taus[[m]])) {
      warning(sprintf(
        "simulate_success_at_n: metric '%s' has no threshold in config$taus; it will not contribute to the success
        criterion.",
        m
      ))
    } else if (length(taus[[m]]) == 1L) {
      success <- success & (max_errors[, m] <= taus[[m]])
    }
  }

  list(
    success       = success,
    success_count = sum(success),
    success_rate  = mean(success),
    rep_out       = rep_out
  )
}


simulate_success_curve_for_alpha <- function(alpha, n_values, config, seed_offset = 0L) {
  rows <- vector("list", length(n_values))
  for (i in seq_along(n_values)) {
    config_i <- config
    if (!is.null(config$seed)) {
      config_i$seed <- as.integer(config$seed + seed_offset + i - 1L)
    }
    
    sim_i <- simulate_success_at_n(alpha = alpha, n = n_values[[i]], config = config_i)
    rows[[i]] <- data.frame(
      alpha = alpha,
      n = as.integer(n_values[[i]]),
      success_rate = sim_i$success_rate,
      success_count = sim_i$success_count,
      B = as.integer(config$B),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


# Multi-patient success simulation -----------------------------------------------------------------------------------

#' Draw one patient-specific alpha constrained to be above a minimum.
#'
#' @param alpha      Baseline alpha value.
#' @param sigma      Standard deviation of the patient-specific Normal draw.
#' @param alpha_min  Strict lower bound for accepted draws.
#' @param max_tries  Maximum number of rejection-sampling attempts.
#'
#' @return Numeric scalar greater than `alpha_min`.
draw_alpha_j <- function(alpha, sigma, alpha_min = 1 + 1e-8, max_tries = 100000L) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha)) {
    stop("alpha must be a single finite number.", call. = FALSE)
  }
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma < 0) {
    stop("sigma must be a single finite number >= 0.", call. = FALSE)
  }
  if (!is.numeric(alpha_min) || length(alpha_min) != 1L || !is.finite(alpha_min) || alpha_min <= 1) {
    stop("alpha_min must be a single finite number > 1.", call. = FALSE)
  }
  if (!is.numeric(max_tries) || length(max_tries) != 1L || !is.finite(max_tries) || max_tries %% 1 != 0 || max_tries < 1L) {
    stop("max_tries must be a single integer >= 1.", call. = FALSE)
  }

  if (sigma == 0) {
    return(max(alpha, alpha_min))
  }

  for (attempt in seq_len(max_tries)) {
    alpha_j <- stats::rnorm(1L, mean = alpha, sd = sigma)
    if (is.finite(alpha_j) && alpha_j > alpha_min) {
      return(alpha_j)
    }
  }

  stop(
    sprintf(
      "Failed to draw alpha_j > %s from N(%s, %s) after %d attempts.",
      format(alpha_min, trim = TRUE),
      format(alpha, trim = TRUE),
      format(sigma, trim = TRUE),
      as.integer(max_tries)
    ),
    call. = FALSE
  )
}


#' Simulate multi-patient replicate success at one fixed sample size.
#'
#' @param alpha       Positive scalar baseline alpha.
#' @param n           Total sample size per patient.
#' @param n_patients  Number of patients per replicate.
#' @param config      Named list with the standard simulation fields plus
#'   `alpha_sigma` and `alpha_min`.
#'
#' @return List with elements `success`, `success_count`, and `success_rate`.
simulate_success_at_n_multi_patient <- function(alpha, n, n_patients, config) {
  if (!is.null(config$seed)) set.seed(config$seed)

  B <- as.integer(config$B)
  K <- as.integer(config$K)
  n_patients <- as.integer(n_patients)
  metrics <- config$metrics
  taus <- config$taus
  active_metrics <- intersect(metrics, names(taus))

  if (length(active_metrics) < 1L) {
    stop("config$taus must provide at least one threshold for a metric in config$metrics.", call. = FALSE)
  }

  success <- rep(TRUE, B)

  for (b in seq_len(B)) {
    metric_sums <- lapply(active_metrics, function(m) numeric(K))
    names(metric_sums) <- active_metrics

    for (j in seq_len(n_patients)) {
      alpha_j <- draw_alpha_j(
        alpha = alpha,
        sigma = config$alpha_sigma,
        alpha_min = config$alpha_min
      )
      p_j <- generate_proportions(
        alpha = alpha_j,
        K = K,
        method = config$proportion_method
      )
      rep_out_j <- run_replicates(
        p = p_j,
        n = n,
        B = 1L,
        metrics = active_metrics,
        model = config$model,
        tie_method = config$tie_method,
        seed = NULL
      )

      for (m in active_metrics) {
        metric_sums[[m]] <- metric_sums[[m]] + rep_out_j$errors[1, , m]
      }
    }

    success_b <- TRUE
    for (m in active_metrics) {
      mean_errors_m <- metric_sums[[m]] / n_patients
      success_b <- success_b & all(mean_errors_m <= taus[[m]])
    }
    success[[b]] <- success_b
  }

  list(
    success = success,
    success_count = sum(success),
    success_rate = mean(success)
  )
}


#' Simulate success-rate curves over sample size and patient counts for one alpha.
#'
#' @param alpha        Positive scalar baseline alpha.
#' @param n_values     Numeric vector of sample sizes per patient; values are
#'   rounded to integers before simulation.
#' @param n_patients   Integer vector of patient counts.
#' @param config       Multi-patient simulation configuration.
#' @param seed_offset  Optional integer seed offset used when looping over
#'   multiple scenarios externally.
#'
#' @return Data.frame with one row per `n` x `n_patients` combination.
simulate_success_curve_for_alpha_over_n_and_patients <- function(alpha, n_values, n_patients, config, seed_offset = 0L) {
  n_rows <- length(n_values) * length(n_patients)
  rows <- vector("list", n_rows)
  row_i <- 1L

  for (i in seq_along(n_values)) {
    for (j in seq_along(n_patients)) {
      n_i <- as.integer(round(n_values[[i]]))
      config_i <- config
      if (!is.null(config$seed)) {
        config_i$seed <- as.integer(config$seed + seed_offset + row_i - 1L)
      }

      sim_i <- simulate_success_at_n_multi_patient(
        alpha = alpha,
        n = n_i,
        n_patients = n_patients[[j]],
        config = config_i
      )
      rows[[row_i]] <- data.frame(
        alpha = alpha,
        n = n_i,
        n_patients = as.integer(n_patients[[j]]),
        success_rate = sim_i$success_rate,
        success_count = sim_i$success_count,
        B = as.integer(config$B),
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }
  do.call(rbind, rows)
}

simulate_success_curve_for_alpha_over_patients <- function(alpha, n_values, n_patients, config, seed_offset = 0L) {
  simulate_success_curve_for_alpha_over_n_and_patients(
    alpha = alpha,
    n_values = n_values,
    n_patients = n_patients,
    config = config,
    seed_offset = seed_offset
  )
}