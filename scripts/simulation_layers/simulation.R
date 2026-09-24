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

# Replicate RNG helpers ---------------------------------------------------------------------------------------

#' Number of cores to use for replicate-level parallelism.
#'
#' Returns `1L` on any non-unix platform (in particular Windows), so replicates always run serially
#' there. On unix platforms, returns `parallel::detectCores() - 1L`, falling back to `1L` when that
#' value is not available (`NA`) or less than `1`. There is no configuration override: this is a
#' deliberate, fixed policy.
#'
#' @return A single positive integer.
replicate_cores <- function() {
  if (!identical(.Platform$OS.type, "unix")) {
    return(1L)
  }
  cores <- parallel::detectCores() - 1L
  if (is.na(cores) || cores < 1L) {
    return(1L)
  }
  as.integer(cores)
}

#' Build B independent, reproducible L'Ecuyer-CMRG RNG streams.
#'
#' Each returned element is a `.Random.seed` vector that, once installed as the active RNG state
#' (under `RNGkind("L'Ecuyer-CMRG")`), starts an independent substream. Stream `b` is derived by
#' chaining `parallel::nextRNGStream()` `b - 1` times from the stream seeded directly from `seed`, so
#' it depends only on `(seed, b)` — never on `B` or on how many cores later consume the streams. As a
#' result, the first `k` streams of a `B`-replicate call are identical to the first `k` streams of any
#' larger `B' > k` call made with the same `seed`.
#'
#' The caller's RNG kind is always restored on exit, and `.Random.seed` is restored to whatever it
#' was immediately *after* resolving `seed` (see below), so calling this function with an explicit
#' `seed` has no visible effect on the global RNG state.
#'
#' @param seed Optional single integer seed. If `NULL`, a seed is drawn from the caller's current RNG
#'   state via `sample.int()` *before* that state is saved, so — exactly like any other call that
#'   consumes randomness — the caller's RNG advances and is not rewound afterwards; consecutive
#'   unseeded calls therefore draw different seeds and produce different streams.
#' @param B    Number of streams to build (positive integer).
#'
#' @return A list of length `B` of `.Random.seed` vectors, one per replicate.
replicate_streams <- function(seed, B) {
  B <- validate_positive_integer(B, "B")

  # Draw a NULL seed *before* saving/restoring the caller's RNG state below, so an unseeded call
  # advances the caller's RNG exactly like any other call to sample.int() would -- otherwise the
  # restore would silently undo the draw and every unseeded call would return identical streams.
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  }

  old_kind <- RNGkind()
  has_old_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (has_old_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    # RNGkind() must be restored *before* re-assigning .Random.seed: assigning .Random.seed while
    # the active kind still differs re-derives/mutates the seed instead of reinstating it exactly.
    suppressWarnings(RNGkind(old_kind[[1L]], old_kind[[2L]], old_kind[[3L]]))
    if (has_old_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  stream  <- .Random.seed
  streams <- vector("list", B)
  for (b in seq_len(B)) {
    streams[[b]] <- stream
    stream       <- parallel::nextRNGStream(stream)
  }
  streams
}

#' Check `parallel::mclapply()` results for failed or missing replicates.
#'
#' `parallel::mclapply()` does not stop the parent process on a worker error: a replicate whose call
#' raised an error is returned as a `"try-error"` object instead of its real result, and a replicate
#' whose worker was killed (e.g. by the OS, out of memory) is returned as `NULL`. Left unchecked,
#' either failure mode would silently propagate a bogus value into downstream aggregation. This
#' helper scans `results` and stops with the first failure's message if any replicate failed; it does
#' not itself depend on `parallel::mclapply()` or on unix, so it can be unit-tested on any platform.
#'
#' @param results List of per-replicate results, as returned by `parallel::mclapply()` (or any list
#'   that may contain `NULL` or `"try-error"` elements).
#'
#' @return `results`, invisibly and unchanged, if every replicate succeeded.
check_replicate_results <- function(results) {
  for (i in seq_along(results)) {
    res <- results[[i]]
    if (is.null(res)) {
      stop(
        sprintf("replicate_apply: worker for replicate %d was killed or returned no result.", i),
        call. = FALSE
      )
    }
    if (inherits(res, "try-error")) {
      stop(
        sprintf(
          "replicate_apply: replicate %d failed: %s",
          i,
          conditionMessage(attr(res, "condition"))
        ),
        call. = FALSE
      )
    }
  }
  invisible(results)
}

#' Apply FUN once per replicate, each with its own independent RNG stream.
#'
#' Before calling `FUN(b)`, installs `streams[[b]]` as the active `.Random.seed` (under
#' `RNGkind("L'Ecuyer-CMRG")`), so replicate `b` always draws from the same stream no matter which
#' worker processes it or in what order. On unix platforms with more than one core available
#' (`replicate_cores()`), replicates are distributed across cores with `parallel::mclapply()`
#' (`mc.set.seed = FALSE`, since seeding is handled explicitly per replicate); on Windows, or when
#' only one core is available, replicates run serially via `lapply()`. Either way, results are
#' returned in replicate order and do not depend on the number of cores used.
#'
#' `parallel::mclapply()` does not raise an error in the parent process when a worker fails or is
#' killed; `check_replicate_results()` is used to detect and re-raise those failures so they are not
#' silently swallowed.
#'
#' The caller's RNG kind and `.Random.seed` are saved and restored on exit.
#'
#' @param streams List of `.Random.seed` vectors (as produced by `replicate_streams()`), one per
#'   replicate.
#' @param FUN     Function of one argument, the replicate index `b` (the position of the
#'   corresponding stream in `streams`), returning that replicate's result.
#'
#' @return A list of length `length(streams)`, in replicate order, of `FUN`'s return values.
replicate_apply <- function(streams, FUN) {
  old_kind <- RNGkind()
  has_old_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (has_old_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    # RNGkind() must be restored *before* re-assigning .Random.seed: assigning .Random.seed while
    # the active kind still differs re-derives/mutates the seed instead of reinstating it exactly.
    suppressWarnings(RNGkind(old_kind[[1L]], old_kind[[2L]], old_kind[[3L]]))
    if (has_old_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")

  run_one <- function(b) {
    assign(".Random.seed", streams[[b]], envir = globalenv())
    FUN(b)
  }

  cores <- replicate_cores()
  if (cores > 1L) {
    results <- parallel::mclapply(seq_along(streams), run_one, mc.cores = cores, mc.set.seed = FALSE)
    check_replicate_results(results)
  } else {
    lapply(seq_along(streams), run_one)
  }
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
  M <- length(metrics)

  # Row order within one replicate's data.frame: person_id (outer), metric, cell_type (inner) --
  # matching expand.grid()'s fastest-first variation of its first argument.
  grid <- expand.grid(
    cell_type = seq_len(K),
    metric    = metrics,
    person_id = seq_len(n_people),
    KEEP.OUT.ATTRS   = FALSE,
    stringsAsFactors = FALSE
  )
  idx_person_cell <- cbind(grid$person_id, grid$cell_type)
  metric_index    <- match(grid$metric, metrics)

  streams <- replicate_streams(seed, B)
  replicate_frames <- replicate_apply(streams, function(b) {
    draw <- simulate_counts(
      p = p,
      model = "dirichlet_multinomial",
      n_people = n_people,
      n_per_person = n_per_person,
      concentration = concentration
    )
    counts      <- draw$counts
    person_p    <- draw$person_true_proportions
    observed_p  <- counts_to_proportions(counts, n_per_person)
    errors_list <- compute_errors(observed_p, person_p, metrics = metrics, n = n_per_person)

    errors_3d <- array(NA_real_, dim = c(n_people, K, M))
    for (mi in seq_len(M)) {
      errors_3d[, , mi] <- errors_list[[metrics[[mi]]]]
    }

    data.frame(
      scenario_id                = scenario_id,
      n_people                   = n_people,
      concentration              = concentration,
      replicate                  = b,
      person_id                  = grid$person_id,
      cell_type                  = grid$cell_type,
      metric                     = grid$metric,
      count                      = as.integer(counts[idx_person_cell]),
      observed_proportion        = observed_p[idx_person_cell],
      person_true_proportion     = person_p[idx_person_cell],
      population_mean_proportion = p[grid$cell_type],
      error                      = errors_3d[cbind(grid$person_id, grid$cell_type, metric_index)],
      stringsAsFactors = FALSE
    )
  })

  list(
    person_results = do.call(rbind, replicate_frames),
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

  streams <- replicate_streams(seed, B)
  replicate_results <- replicate_apply(streams, function(b) {
    y_b      <- simulate_counts(p, n, model = model, ...)
    phat_b   <- counts_to_proportions(y_b, n)
    errors_b <- compute_errors(phat_b, p, metrics = metrics, n = n)

    max_error_b <- stats::setNames(rep(NA_real_,    length(metrics)), metrics)
    argmax_b    <- stats::setNames(rep(NA_integer_, length(metrics)), metrics)
    for (m in metrics) {
      s <- max_error_summary(errors_b[[m]], tie_method = tie_method)
      max_error_b[[m]] <- s$max_error_value
      argmax_b[[m]]    <- s$argmax_index
    }

    list(phat = phat_b, errors = errors_b, max_error = max_error_b, argmax = argmax_b)
  })

  for (b in seq_len(B)) {
    res <- replicate_results[[b]]
    phat[b, ] <- res$phat
    for (m in metrics) {
      max_errors[b, m] <- res$max_error[[m]]
      argmax[b, m]     <- res$argmax[[m]]
      errors[b, , m]   <- res$errors[[m]]
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
#' Both sampling models share a single success rule, `replicate_success()`: for each metric in `config$taus`, the
#' error is averaged over persons per (replicate, cell type), and the largest of these per-cell-type means must be
#' `<= tau`; a replicate succeeds jointly only if it succeeds for every metric in `config$taus`.
#'
#' The multinomial model has no person structure, so its per-cell-type errors (one draw per replicate) are treated
#' as a single synthetic person (`person_id = 1`) before being handed to `replicate_success()` — averaging over one
#' person is a no-op, so this reduces to "every cell-type error <= tau" for that model, matching its previous
#' behaviour. The Dirichlet-multinomial model's `person_results` (one row per replicate, person, cell type, metric)
#' is passed to `replicate_success()` directly.
#'
#' Metrics in `config$metrics` that have no entry in `config$taus` are simulated but do not contribute to the
#' success criterion: `simulate_success_at_n()` warns about them once per call, and — because such metrics are
#' simply absent from `config$taus` — `replicate_success()` never sees them and so never emits its own "metric not
#' in person_results" warning for the same metric. That second warning path only fires for a metric that has a tau
#' but was not simulated at all, a genuinely different (misconfiguration) case.
#'
#' @param alpha  Positive scalar; Beta shape parameter used to generate the true proportions.
#' @param n      Positive integer sample size. For `config$model == "multinomial"`, the total sample size (cells)
#'   per replicate. For `config$model == "dirichlet_multinomial"`, the number of cells sampled per person
#'   (`n_per_person`); the number of people per replicate is fixed by `config$n_people`, not by `n`.
#' @param config Named list; must contain at minimum:
#'   \describe{
#'     \item{K}{Number of cell types.}
#'     \item{B}{Number of replicates.}
#'     \item{taus}{Named list with one scalar threshold per metric (e.g. `list(AE = 0.02, ARE = 0.5)`); scalar
#'       thresholds only. Metrics simulated but absent here are skipped (see Details).}
#'     \item{metrics}{Character vector of metric names to simulate.}
#'     \item{model}{Sampling model (`"multinomial"` or `"dirichlet_multinomial"`).}
#'     \item{n_people}{Required for `"dirichlet_multinomial"`; number of people per replicate.}
#'     \item{concentration}{Required for `"dirichlet_multinomial"`; positive Dirichlet concentration parameter.}
#'     \item{tie_method}{Tie-breaking rule for max-error argmax (multinomial only; unused by
#'       `replicate_success()`, but still forwarded to `run_replicates()`).}
#'     \item{proportion_method}{Proportion-generation method.}
#'   }
#' @param seed   Optional integer RNG seed forwarded to `run_replicates()`; defaults to `config$seed`. Because
#'   `run_replicates()` derives per-replicate RNG streams from `seed` alone (see `replicate_streams()`), calling
#'   this function with the same `seed` at different `n` draws from the same per-replicate streams (common random
#'   numbers across `n`), which is what the sample-size solver relies on when comparing pilots.
#'
#' @return List with elements:
#'   \describe{
#'     \item{success}{Logical vector of length B, in replicate order.}
#'     \item{success_count}{Integer; number of successful replicates.}
#'     \item{success_rate}{Numeric; fraction of successful replicates.}
#'     \item{rep_out}{Raw output of `run_replicates()`.}
#'   }
simulate_success_at_n <- function(alpha, n = NULL, config, seed = config$seed) {
  p <- generate_proportions(
    alpha  = alpha,
    K      = config$K,
    method = config$proportion_method
  )

  missing_tau_msg <- paste0(
    "simulate_success_at_n: metric '%s' has no threshold in config$taus; ",
    "it will not contribute to the success criterion."
  )
  for (m in setdiff(config$metrics, names(config$taus))) {
    warning(sprintf(missing_tau_msg, m), call. = FALSE)
  }

  if (identical(config$model, "dirichlet_multinomial")) {
    rep_out <- run_replicates(
      p             = p,
      B             = config$B,
      metrics       = config$metrics,
      model         = config$model,
      seed          = seed,
      n_people      = config$n_people,
      n_per_person  = n,
      concentration = config$concentration
    )
    pass <- replicate_success(rep_out$person_results, config$taus)

    return(list(
      success       = as.logical(pass$pass),
      success_count = sum(pass$pass),
      success_rate  = mean(pass$pass),
      rep_out       = rep_out
    ))
  }

  rep_out <- run_replicates(
    p          = p,
    n          = n,
    B          = config$B,
    metrics    = config$metrics,
    model      = config$model,
    tie_method = config$tie_method,
    seed       = seed
  )

  # rep_out$errors is a B x K x M array (dimnames list(NULL, 1..K, metrics)); flatten it into a long data.frame
  # with a single synthetic person_id = 1 per replicate so replicate_success() (which expects one row per
  # replicate/person/cell_type/metric) can be reused for the multinomial model too. expand.grid()'s default
  # variation order (first argument fastest) matches the array's column-major storage order (dim 1 fastest, then
  # dim 2, then dim 3), so `error = as.vector(rep_out$errors)` lines up exactly with `grid`.
  errors_dim <- dim(rep_out$errors)
  metrics    <- dimnames(rep_out$errors)[[3L]]
  grid <- expand.grid(
    replicate = seq_len(errors_dim[[1L]]),
    cell_type = seq_len(errors_dim[[2L]]),
    metric    = metrics,
    KEEP.OUT.ATTRS   = FALSE,
    stringsAsFactors = FALSE
  )
  person_results <- data.frame(
    replicate = grid$replicate,
    person_id = 1L,
    cell_type = grid$cell_type,
    metric    = grid$metric,
    error     = as.vector(rep_out$errors),
    stringsAsFactors = FALSE
  )

  pass <- replicate_success(person_results, config$taus)

  list(
    success       = as.logical(pass$pass),
    success_count = sum(pass$pass),
    success_rate  = mean(pass$pass),
    rep_out       = rep_out
  )
}
