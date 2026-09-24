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

validate_positive_integer <- function(x, name, allow_vector = FALSE) {
  valid_length <- if (allow_vector) length(x) >= 1L else length(x) == 1L
  if (!is.numeric(x) || !valid_length || any(!is.finite(x)) || any(x < 1L) || any(x %% 1 != 0)) {
    expected <- if (allow_vector) {
      "a non-empty vector of positive integers"
    } else {
      "a positive integer"
    }
    stop(sprintf("%s must be %s.", name, expected), call. = FALSE)
  }
  as.integer(x)
}

validate_positive_numeric <- function(x, name, allow_vector = FALSE) {
  valid_length <- if (allow_vector) length(x) >= 1L else length(x) == 1L
  if (!is.numeric(x) || !valid_length || any(!is.finite(x)) || any(x <= 0)) {
    expected <- if (allow_vector) {
      "a non-empty vector of positive finite numbers"
    } else {
      "a single positive finite number"
    }
    stop(sprintf("%s must be %s.", name, expected), call. = FALSE)
  }
  as.numeric(x)
}

validate_required_person_fraction <- function(required_person_fraction) {
  if (!is.numeric(required_person_fraction) ||
      length(required_person_fraction) != 1L ||
      !is.finite(required_person_fraction) ||
      required_person_fraction <= 0 ||
      required_person_fraction > 1) {
    stop("required_person_fraction must be a single number in (0, 1].", call. = FALSE)
  }
  as.numeric(required_person_fraction)
}

#' Draw one composition from a Dirichlet distribution.
#'
#' @param concentration_parameters Positive Dirichlet concentration parameters.
#'
#' @return Numeric vector on the simplex with the same length as
#'   `concentration_parameters`.
sample_dirichlet <- function(concentration_parameters) {
  concentration_parameters <- validate_positive_numeric(
    concentration_parameters,
    "concentration_parameters",
    allow_vector = TRUE
  )
  draws <- stats::rgamma(length(concentration_parameters), shape = concentration_parameters, rate = 1)
  total <- sum(draws)
  if (!is.finite(total) || total <= 0) {
    stop("Dirichlet sampling produced an invalid total gamma draw.", call. = FALSE)
  }
  draws / total
}

#' Simulate person-level counts from a Dirichlet-multinomial hierarchy.
#'
#' @param p Population mean proportion vector.
#' @param n_people Number of people to sample.
#' @param n_per_person Number of cells sampled for each person.
#' @param concentration Positive Dirichlet concentration parameter.
#'
#' @return List containing an `n_people` by `K` count matrix and a matching
#'   matrix of person-specific latent true proportions.
simulate_counts_dirichlet_multinomial <- function(p, n_people, n_per_person, concentration) {
  validate_proportions(p)
  n_people      <- validate_positive_integer(n_people, "n_people")
  n_per_person  <- validate_positive_integer(n_per_person, "n_per_person")
  concentration <- validate_positive_numeric(concentration, "concentration")

  K <- length(p)
  person_true_proportions <- t(vapply(
    seq_len(n_people),
    function(person_id) sample_dirichlet(concentration * p),
    FUN.VALUE = numeric(K)
  ))
  counts <- t(vapply(
    seq_len(n_people),
    function(person_id) simulate_counts_multinomial(person_true_proportions[person_id, ], n_per_person),
    FUN.VALUE = integer(K)
  ))
  colnames(counts) <- paste0("cell_type_", seq_len(K))
  colnames(person_true_proportions) <- colnames(counts)

  list(
    counts = counts,
    person_true_proportions = person_true_proportions
  )
}

# Placeholder: to be implemented when correlation support is added.
# simulate_counts_logistic_normal_multinomial <- function(p, n, Sigma, ...) {
#   stop("Logistic-normal multinomial not yet implemented.")
# }

#' Dispatcher: simulate counts from the requested model.
#'
#' @param p      True proportion vector.
#' @param n      Total sample size for the multinomial model.
#' @param model  Sampling model.
#' @param n_people Number of people for the Dirichlet-multinomial model.
#' @param n_per_person Number of cells sampled for each person in the
#'   Dirichlet-multinomial model.
#' @param concentration Dirichlet concentration parameter for the
#'   Dirichlet-multinomial model.
#' @param ...    Additional arguments forwarded to the concrete simulator (reserved for future overdispersed /
#'               correlated models).
#'
#' @return For `"multinomial"`, an integer vector of length K summing to n.
#'   For `"dirichlet_multinomial"`, a list with person-level count and latent
#'   proportion matrices.
simulate_counts <- function(p, n = NULL,
                            model = c("multinomial",
                                      "dirichlet_multinomial",
                                      "logistic_normal_multinomial"),
                            n_people = NULL,
                            n_per_person = NULL,
                            concentration = NULL,
                            ...) {
  model <- match.arg(model)
  switch(model,
    multinomial = {
      if (is.null(n)) {
        stop("n must be provided for model = 'multinomial'.", call. = FALSE)
      }
      simulate_counts_multinomial(p, n)
    },
    dirichlet_multinomial = simulate_counts_dirichlet_multinomial(
      p = p,
      n_people = n_people,
      n_per_person = n_per_person,
      concentration = concentration
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


#' Run B simulation replicates and store error results.
#'
#' Efficiency strategy: simulate B times once, store only the per-replicate max error values and argmax indices.
#' Threshold evaluation is done post-hoc by `evaluate_thresholds()` without re-simulating.
#'
#' @param p          True proportion vector (length K).
#' @param n          Total sample size for the multinomial model.
#' @param B          Number of replicates.
#' @param metrics    Error metrics to compute; any subset of
#'   `c("AE", "ARE", "TSE", "LAE")`.
#' @param model      Sampling model passed to `simulate_counts()`.
#' @param n_people   Number of people for the Dirichlet-multinomial model.
#' @param n_per_person Number of cells sampled for each person in the
#'   Dirichlet-multinomial model.
#' @param concentration Positive Dirichlet concentration parameter.
#' @param scenario_id Optional scenario identifier included in person-level
#'   Dirichlet-multinomial output.
#' @param tie_method Tie-breaking rule passed to `max_error_summary()`.
#' @param seed       Optional integer seed for reproducibility.
#' @param ...        Additional arguments forwarded to `simulate_counts()`.
#'
#' @return For `"multinomial"`, the existing list with max-error matrices and
#'   arrays. For `"dirichlet_multinomial"`, a list with `person_results`, a
#'   tidy data.frame containing one row per replicate, person, cell type, and
#'   metric, plus `inputs`.
#'   \describe{
#'     \item{max_errors}{B x M numeric matrix of max error values.}
#'     \item{argmax}{B x M integer matrix of argmax indices.}
#'     \item{errors}{B x K x M numeric array of per-cell-type errors.}
#'     \item{phat}{B x K numeric matrix of observed proportions.}
#'     \item{inputs}{Copy of all input arguments (including seed used).}
#'   }
run_replicates_dirichlet_multinomial <- function(p, B, metrics,
                                                 n_people, n_per_person,
                                                 concentration, scenario_id = NA_character_,
                                                 seed = NULL) {
  B             <- validate_positive_integer(B, "B")
  n_people      <- validate_positive_integer(n_people, "n_people")
  n_per_person  <- validate_positive_integer(n_per_person, "n_per_person")
  concentration <- validate_positive_numeric(concentration, "concentration")
  validate_proportions(p)

  K <- length(p)
  rows <- vector("list", B * n_people * K * length(metrics))
  row_index <- 0L
  for (replicate_id in seq_len(B)) {
    draw <- simulate_counts(
      p = p,
      model = "dirichlet_multinomial",
      n_people = n_people,
      n_per_person = n_per_person,
      concentration = concentration
    )
    for (person_id in seq_len(n_people)) {
      counts <- draw$counts[person_id, ]
      person_p <- draw$person_true_proportions[person_id, ]
      observed_p <- counts_to_proportions(counts, n_per_person)
      errors <- compute_errors(observed_p, person_p, metrics = metrics, n = n_per_person)
      for (metric in metrics) {
        for (cell_type in seq_len(K)) {
          row_index <- row_index + 1L
          rows[[row_index]] <- data.frame(
            scenario_id = scenario_id,
            n_people = n_people,
            concentration = concentration,
            replicate = replicate_id,
            person_id = person_id,
            cell_type = cell_type,
            metric = metric,
            count = as.integer(counts[[cell_type]]),
            observed_proportion = observed_p[[cell_type]],
            person_true_proportion = person_p[[cell_type]],
            population_mean_proportion = p[[cell_type]],
            error = errors[[metric]][[cell_type]],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  list(
    person_results = do.call(rbind, rows),
    inputs = list(
      p = p,
      B = B,
      metrics = metrics,
      model = "dirichlet_multinomial",
      n_people = n_people,
      n_per_person = n_per_person,
      concentration = concentration,
      scenario_id = scenario_id,
      seed = seed
    )
  )
}

run_replicates <- function(p, n = NULL, B,
                           metrics = c("AE", "ARE"),
                           model = "multinomial",
                           tie_method = "random",
                           seed = NULL,
                           n_people = NULL,
                           n_per_person = NULL,
                           concentration = NULL,
                           scenario_id = NA_character_,
                           ...) {
  if (!is.null(seed)) set.seed(seed)
  if (identical(model, "dirichlet_multinomial")) {
    return(run_replicates_dirichlet_multinomial(
      p = p,
      B = B,
      metrics = metrics,
      n_people = n_people,
      n_per_person = n_per_person,
      concentration = concentration,
      scenario_id = scenario_id,
      seed = seed
    ))
  }
  if (!identical(model, "multinomial")) {
    stop(sprintf("model = '%s' is not supported by run_replicates().", model), call. = FALSE)
  }
  if (is.null(n)) {
    stop("n must be provided for model = 'multinomial'.", call. = FALSE)
  }
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
#' @param n              Total sample size (positive integer) for the
#'   multinomial model. For the Dirichlet-multinomial model, use
#'   `n_per_person`.
#' @param n_per_person   Number of sampled cells for each person in the
#'   Dirichlet-multinomial model.
#' @param config         Named list; must contain at minimum:
#'   \describe{
#'     \item{K}{Number of cell types.}
#'     \item{B}{Number of replicates.}
#'     \item{taus}{Named list with one scalar threshold per metric (e.g.
#'       `list(AE = 0.02, ARE = 0.5)`); scalar thresholds only.}
#'     \item{metrics}{Character vector of metric names to simulate.}
#'     \item{model}{Sampling model (`"multinomial"` or
#'       `"dirichlet_multinomial"`).}
#'     \item{n_people}{Required for `"dirichlet_multinomial"`; number of
#'       people per replicate.}
#'     \item{concentration}{Required for `"dirichlet_multinomial"`; positive
#'       Dirichlet concentration parameter.}
#'     \item{required_person_fraction}{Optional fraction in `(0, 1]` of
#'       people who must pass all active thresholds; defaults to `1`.}
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
# TODO: the Dirichlet-multinomial branch still uses the old per-person rule (a person passes if all cell-type errors are
#   within tau; a replicate succeeds if `required_person_fraction` of people pass). Update it to the current definition
#   used by `extract_success_rate()`: mean error over persons per cell type, max over cell types, <= tau for every metric.
simulate_success_at_n <- function(alpha, n = NULL, config, n_per_person = NULL) {
  p <- generate_proportions(
    alpha  = alpha,
    K      = config$K,
    method = config$proportion_method
  )
  if (identical(config$model, "dirichlet_multinomial")) {
    if (is.null(n_per_person)) {
      n_per_person <- n
    }
    if (is.null(n_per_person)) {
      stop("n_per_person must be provided for model = 'dirichlet_multinomial'.", call. = FALSE)
    }
    required_person_fraction <- if (is.null(config$required_person_fraction)) {
      1
    } else {
      config$required_person_fraction
    }
    required_person_fraction <- validate_required_person_fraction(required_person_fraction)
    rep_out <- run_replicates(
      p = p,
      B = config$B,
      metrics = config$metrics,
      model = config$model,
      seed = config$seed,
      n_people = config$n_people,
      n_per_person = n_per_person,
      concentration = config$concentration
    )
    person_results <- rep_out$person_results
    person_pass <- rep(TRUE, config$B * config$n_people)
    person_keys <- expand.grid(
      replicate = seq_len(config$B),
      person_id = seq_len(config$n_people),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    for (metric in config$metrics) {
      threshold <- config$taus[[metric]]
      if (is.null(threshold)) {
        warning(sprintf(
          "simulate_success_at_n: metric '%s' has no threshold in config$taus; it will not contribute to the success criterion.",
          metric
        ))
      } else if (length(threshold) == 1L) {
        metric_results <- person_results[person_results$metric == metric, , drop = FALSE]
        metric_pass <- vapply(seq_len(nrow(person_keys)), function(i) {
          rows <- metric_results[
            metric_results$replicate == person_keys$replicate[[i]] &
              metric_results$person_id == person_keys$person_id[[i]],
            ,
            drop = FALSE
          ]
          nrow(rows) > 0L && all(rows$error <= threshold)
        }, logical(1L))
        person_pass <- person_pass & metric_pass
      }
    }

    person_keys$pass <- person_pass
    people_passing <- tapply(
      person_keys$pass,
      person_keys$replicate,
      sum
    )
    people_required <- ceiling(required_person_fraction * config$n_people)
    success <- people_passing >= people_required

    return(list(
      success = as.logical(success),
      success_count = sum(success),
      success_rate = mean(success),
      person_success = person_keys,
      required_people_passing = people_required,
      rep_out = rep_out
    ))
  }

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
