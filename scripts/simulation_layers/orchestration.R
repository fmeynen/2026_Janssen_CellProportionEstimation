# Orchestration ---------------------------------------------------------------------------------------------------
# Orchestration layer: high-level experiment runners that coordinate all other layers.
#
# Depends on: simulation.R, calculation.R, extraction.R

# Simulation - I/O -------------------------------------------------------------------------------------------------

# I/O helpers for persisting simulation results.
#
# Each unique combination of simulation parameters is saved as a separate .rds file whose name contains a short MD5 hash
# of those parameters.

#' Compute an MD5 hash of an R object.
#'
#' Serialises `x` to a temporary file, computes its MD5 checksum, and returns the hex string.
#' Only `tools` (part of base R) is required.
#'
#' @param x Any R object.
#' @return A 32-character hexadecimal string.
hash_config <- function(x) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  saveRDS(x, tmp)
  unname(tools::md5sum(tmp))
}


#' Build the file path for a simulation result.
#'
#' The file name is `<name>_<md5>.rds`, where the MD5 is derived from the serialised `config`.
#' Identical parameters always resolve to the same path; any change produces a different path.
#'
#' @param config  Named list of simulation parameters.
#' @param dir     Directory that will hold the result files.  Created automatically if it does not yet exist.
#' @param name    Short label prefixed to the file name
#'   (e.g. `"simulation_errorchoice"`).
#'
#' @return Absolute path string.
simulation_result_path <- function(config, dir, name) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(name, "_", hash_config(config), ".rds"))
}


# Schema version of cached results. Bump this whenever the shape of a cached result changes (columns, list elements,
# types), so that stale files on disk are never read back as if they had the new shape.
CACHE_SCHEMA <- 3L


#' Return a cached result, or compute it and cache it.
#'
#' The cache file is `simulation_result_path(c(key, list(cache_schema = CACHE_SCHEMA)), dir, name)`. Callers should
#' pass only the fields that influence the computed result in `key`, so that e.g. plotting-only settings never
#' trigger a recomputation.
#'
#' @param key             Named list of the simulation-relevant inputs that identify the result.
#' @param name            Short label prefixed to the file name (e.g. `"errorchoice"`).
#' @param compute         Zero-argument function returning the result to cache.
#' @param cache           Logical; read/write the `.rds` cache file under `dir`. If `FALSE`, always computes and
#'   never touches the disk.
#' @param force_recompute Logical; ignore an existing cache file and recompute (still writes the new result when
#'   `cache` is `TRUE`).
#' @param dir             Directory holding the cache files. Created automatically if it does not yet exist.
#'
#' @return The cached or freshly computed result.
cached_result <- function(key, name, compute, cache = TRUE, force_recompute = FALSE, dir) {
  result_file <- simulation_result_path(c(key, list(cache_schema = CACHE_SCHEMA)), dir, name)

  if (cache && !force_recompute && file.exists(result_file)) {
    return(readRDS(result_file))
  }

  result <- compute()
  if (cache) {
    saveRDS(result, result_file)
  }
  result
}

#' Run the sample-size solver over a grid of alphas with warm start and per-alpha caching.
#'
#' For each alpha in `config$alpha`, in the given order, runs `estimate_sample_size()` to find the smallest sample
#' size reaching `config$success_rate_target`. Alphas are chained: the first alpha's solver starts from
#' `config$n_init`, and every later alpha's solver starts from the previous alpha's `final_n` (warm start). Each
#' alpha's result is cached in its own file under `cache_dir`, keyed by a config that includes that alpha's
#' `n_init` — so changing an earlier alpha (and hence a later alpha's warm-start value) invalidates that later
#' alpha's cache entry, while leaving unaffected entries untouched.
#'
#' @param config          List with `alpha` (numeric vector, the grid) and `n_init` (the first alpha's starting
#'   centre), plus every field required by `estimate_sample_size()` (`success_rate_target`, `rel_tol`,
#'   `max_iterations`, `B`, `f0`, `f_floor`, `seed`) and by `simulate` (by default `simulate_success_at_n()`, which
#'   also needs `K`, `taus`, `metrics`, `model`, `tie_method`, `proportion_method`, and — for
#'   `model == "dirichlet_multinomial"` — `n_people` and `concentration`).
#' @param cache           Logical; read/write per-alpha `.rds` cache files under `cache_dir`.
#' @param force_recompute Logical; ignore any existing cache file and recompute (still writes the new result when
#'   `cache` is `TRUE`).
#' @param cache_dir       Directory holding the per-alpha cache files.
#' @param simulate        Function `(alpha, n, config, seed)` forwarded to `estimate_sample_size()`. Defaults to
#'   `simulate_success_at_n()`.
#'
#' @return List with elements:
#'   \describe{
#'     \item{sample_size}{Data.frame with one row per alpha and columns `alpha`, `sample_size` (integer,
#'       `final_n`), `stopping_reason`, `iterations_used`.}
#'     \item{diagnostics}{Data.frame; `rbind()` of every alpha's `estimate_sample_size()` diagnostics, in grid
#'       order.}
#'   }
run_sample_size_experiment <- function(
  config,
  cache = TRUE,
  force_recompute = FALSE,
  cache_dir = here::here("results", "simresults"),
  simulate = simulate_success_at_n
) {
  alphas <- validate_positive_numeric(
    config$alpha,
    "config$alpha",
    allow_vector = TRUE
  )
  n_init <- validate_positive_numeric(config$n_init, "config$n_init")

  n_alpha <- length(alphas)
  sample_size_rows <- vector("list", n_alpha)
  diag_list <- vector("list", n_alpha)

  for (i in seq_len(n_alpha)) {
    alpha_i <- alphas[[i]]
    alpha_config <- config
    alpha_config$alpha <- alpha_i
    alpha_config$n_init <- n_init
    # Only the fields read by estimate_sample_size() / simulate_success_at_n() go into the cache key.
    key <- list(
      alpha = alpha_i,
      n_init = n_init,
      success_rate_target = config$success_rate_target,
      rel_tol = config$rel_tol,
      max_iterations = config$max_iterations,
      f0 = config$f0,
      f_floor = config$f_floor,
      n_max = config$n_max,
      B = config$B,
      seed = config$seed,
      K = config$K,
      taus = config$taus,
      metrics = config$metrics,
      model = config$model,
      tie_method = config$tie_method,
      proportion_method = config$proportion_method,
      p_max = config$p_max,
      n_people = config$n_people,
      concentration = config$concentration,
      success_rule = success_rule_id()
    )
    alpha_result <- cached_result(
      key = key,
      name = "sample_size",
      compute = function() {
        estimate_sample_size(alpha_i, n_init, alpha_config, simulate = simulate)
      },
      cache = cache,
      force_recompute = force_recompute,
      dir = cache_dir
    )

    sample_size_rows[[i]] <- data.frame(
      alpha = alpha_i,
      sample_size = as.integer(alpha_result$final_n),
      stopping_reason = alpha_result$stopping_reason,
      iterations_used = alpha_result$iterations_used,
      stringsAsFactors = FALSE
    )
    diag_list[[i]] <- alpha_result$diagnostics
    n_init <- alpha_result$final_n
  }

  list(
    sample_size = do.call(rbind, sample_size_rows),
    diagnostics = do.call(rbind, diag_list)
  )
}


# Run individual experiments --------------------------------------------------------------------------------------
#' Generate true proportions for every alpha x p_max scenario, skipping impossible fixed-max combinations.
#'
#' Validates `p_max` (for `proportion_method = "fixed_max_beta"`), builds the alpha x p_max grid (alpha varies
#' fastest) and calls `generate_proportions()` once per row. `generate_proportions()` is RNG-free, so calling it
#' up front does not change any random draws made later by the simulation. When several `p_max` values are given,
#' fixed-max combinations with `K * p_max < 1` are skipped (`generate_props_fixed_max_beta()` warns once per
#' combination);
#' otherwise the error is propagated. Stops if no combination is feasible.
#'
#' @return A list with `grid` (data.frame with `alpha`, `p_max` (NA unless fixed-max), `population_id`),
#'   `p` (list of proportion vectors, one per grid row; NULL for skipped rows) and `feasible` (logical vector).
feasible_scenarios <- function(alpha, K, proportion_method, p_max) {
  if (identical(proportion_method, "fixed_max_beta")) {
    validate_p_max(p_max)
    p_max_values <- as.numeric(p_max)
  } else {
    p_max_values <- NA_real_
  }

  grid <- expand.grid(
    alpha = alpha,
    p_max = p_max_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$population_id <- seq_len(nrow(grid))

  skip_impossible <- identical(proportion_method, "fixed_max_beta") &&
    length(p_max_values) > 1L

  p_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    p_max_i <- grid$p_max[[i]]
    p <- tryCatch(
      generate_proportions(
        alpha = grid$alpha[[i]],
        K = K,
        method = proportion_method,
        p_max = if (is.na(p_max_i)) NULL else p_max_i
      ),
      error = function(e) {
        if (skip_impossible && inherits(e, "impossible_fixed_max_error")) {
          return(NULL)
        }
        stop(e)
      }
    )
    if (!is.null(p)) {
      p_list[[i]] <- p
    }
  }

  feasible <- !vapply(p_list, is.null, logical(1L))
  if (!any(feasible)) {
    stop(
      "No feasible alpha/p_max combinations produced simulation output.",
      call. = FALSE
    )
  }
  list(grid = grid, p = p_list, feasible = feasible)
}


#' Run the Dirichlet-multinomial branch of `run_simulation_experiment()`.
#'
#' One scenario per feasible (alpha, p_max) x (n_people, concentration) combination; each scenario calls
#' `run_replicates(..., keep_person_results = TRUE)` and keeps its long person-level output.
#'
#' @param metrics Error metrics; any subset of `c("AE", "ARE")`. Other arguments as in `run_simulation_experiment()`.
#'
#' @return List with elements:
#'   \describe{
#'     \item{inputs}{All input arguments, with `model = "dirichlet_multinomial"`.}
#'     \item{p_table}{Data.frame with one row per scenario: `scenario_id`, `alpha`, `p_max`, `n_people`,
#'       `concentration`, `cell_type_1`, ..., `cell_type_K`.}
#'     \item{person_results}{Tidy data.frame with one row per scenario x replicate x person x cell type x metric and
#'       columns `scenario_id`, `alpha`, `p_max`, `n_people`, `concentration`, `replicate`, `person_id`, `cell_type`,
#'       `metric`, `count`, `observed_proportion`, `person_true_proportion`, `population_mean_proportion`.}
#'   }
run_dirichlet_multinomial_experiment <- function(
  alpha,
  K,
  B,
  metrics,
  proportion_method = "beta",
  p_max = NULL,
  n_people,
  n_per_person,
  concentration,
  seed
) {
  n_people <- validate_positive_integer(
    n_people,
    "n_people",
    allow_vector = TRUE
  )
  n_per_person <- validate_positive_integer(n_per_person, "n_per_person")
  concentration <- validate_positive_numeric(
    concentration,
    "concentration",
    allow_vector = TRUE
  )

  feasible <- feasible_scenarios(alpha, K, proportion_method, p_max)
  population_scenarios <- feasible$grid
  person_scenarios <- expand.grid(
    n_people = n_people,
    concentration = concentration,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  scenario_grid <- merge(
    population_scenarios,
    person_scenarios,
    by = NULL,
    sort = FALSE
  )
  scenario_grid$scenario_id <- paste0("scenario_", seq_len(nrow(scenario_grid)))

  p_table_list <- vector("list", nrow(scenario_grid))
  person_results_list <- vector("list", nrow(scenario_grid))
  keep <- logical(nrow(scenario_grid))

  for (i in seq_len(nrow(scenario_grid))) {
    scenario <- scenario_grid[i, , drop = FALSE]
    p <- feasible$p[[scenario$population_id[[1L]]]]
    if (is.null(p)) {
      next
    }

    rep_out <- run_replicates(
      p = p,
      B = B,
      metrics = metrics,
      model = "dirichlet_multinomial",
      seed = if (is.null(seed)) NULL else seed + i - 1L,
      n_people = scenario$n_people[[1L]],
      n_per_person = n_per_person,
      concentration = scenario$concentration[[1L]],
      scenario_id = scenario$scenario_id[[1L]],
      keep_person_results = TRUE
    )
    person_results <- rep_out$person_results
    person_results$alpha <- scenario$alpha[[1L]]
    person_results$p_max <- scenario$p_max[[1L]]
    person_results_list[[i]] <- person_results[, c(
      "scenario_id",
      "alpha",
      "p_max",
      "n_people",
      "concentration",
      "replicate",
      "person_id",
      "cell_type",
      "metric",
      "count",
      "observed_proportion",
      "person_true_proportion",
      "population_mean_proportion"
    )]
    p_table_list[[i]] <- extract_p_table_row(
      scenario$alpha[[1L]],
      scenario$p_max[[1L]],
      p,
      K,
      before = list(scenario_id = scenario$scenario_id[[1L]]),
      after = list(
        n_people = scenario$n_people[[1L]],
        concentration = scenario$concentration[[1L]]
      )
    )
    keep[[i]] <- TRUE
  }

  list(
    inputs = list(
      alpha = alpha,
      K = K,
      B = B,
      metrics = metrics,
      proportion_method = proportion_method,
      p_max = p_max,
      model = "dirichlet_multinomial",
      n_people = n_people,
      n_per_person = n_per_person,
      concentration = concentration,
      seed = seed
    ),
    p_table = do.call(rbind, p_table_list[keep]),
    person_results = do.call(rbind, person_results_list[keep])
  )
}


#' Run the Dirichlet-multinomial "errorchoice" experiment over an (alpha, n_people, n_per_person) grid.
#'
#' AE and ARE are studied separately: for each scenario (one combination of `alpha`, `n_people`, `n_per_person`)
#' and each metric, every replicate is reduced immediately to one scalar "stat" -- the max, over cell types, of the
#' error of the person-pooled proportion estimate against the population proportion -- which `run_replicates()`
#' returns directly as `max_errors` (computed by `pooled_error_stat()`). No per-person data are kept, so memory use
#' stays small even for large `B` and `n_people`.
#'
#' Common random numbers: the seed depends only on the (alpha, n_people) pair. Pairs are enumerated in the order
#' of `expand.grid(n_people = n_people, alpha = alpha)`; the j-th pair uses `seed + j - 1L` (or `NULL` when `seed`
#' is `NULL`), identical for every `n_per_person` value, so scenarios that share an (alpha, n_people) pair but
#' differ only in `n_per_person` draw from the same underlying Dirichlet-multinomial streams.
#'
#' @param alpha  Numeric vector; positive Beta shape parameter(s) used by `proportion_method` to generate each
#'   population mean proportion vector (one vector per alpha, computed once and reused across the grid).
#' @param K      Positive integer; number of cell types.
#' @param B      Positive integer; number of replicates per scenario.
#' @param metrics Character vector; error metrics, any of `"AE"`, `"ARE"` (forwarded to `run_replicates()`),
#'   studied independently of one another.
#' @param proportion_method Proportion-generation method forwarded to `generate_proportions()` (default
#'   `"beta"`).
#' @param n_people Positive integer vector; number(s) of people per replicate.
#' @param n_per_person Positive integer vector; number(s) of cells sampled per person.
#' @param concentration Positive numeric scalar; Dirichlet concentration parameter (shared by every scenario).
#' @param seed   Optional single integer seed; see Details for how it is combined with the (alpha, n_people) pair
#'   index. `NULL` means every scenario draws an unseeded (non-reproducible) stream.
#'
#' @return List with elements:
#'   \describe{
#'     \item{inputs}{All input arguments (`alpha`, `K`, `B`, `metrics`, `proportion_method`, `n_people`,
#'       `n_per_person`, `concentration`, `seed`).}
#'     \item{p_table}{Data.frame with one row per alpha and columns `alpha`, `cell_type_1`, ..., `cell_type_K`.}
#'     \item{stats}{Data.frame with one row per (alpha, n_people, n_per_person, metric, replicate) and columns
#'       `alpha`, `n_people`, `concentration`, `n_per_person`, `metric`, `replicate`, `stat`.}
#'   }
run_dm_errorchoice_experiment <- function(
  alpha,
  K,
  B,
  metrics,
  proportion_method = "beta",
  n_people,
  n_per_person,
  concentration,
  seed
) {
  alpha <- validate_positive_numeric(alpha, "alpha", allow_vector = TRUE)
  K <- validate_positive_integer(K, "K")
  B <- validate_positive_integer(B, "B")
  n_people <- validate_positive_integer(n_people, "n_people", allow_vector = TRUE)
  n_per_person <- validate_positive_integer(n_per_person, "n_per_person", allow_vector = TRUE)
  concentration <- validate_positive_numeric(concentration, "concentration")

  # p is a deterministic function of alpha (and K, proportion_method), so it is computed once per alpha and
  # reused across every n_people / n_per_person scenario for that alpha.
  p_list <- lapply(alpha, function(a) {
    generate_proportions(alpha = a, K = K, method = proportion_method)
  })

  pairs <- expand.grid(
    n_people = n_people,
    alpha = alpha,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  n_pairs <- nrow(pairs)
  n_values <- length(n_per_person)
  total_scenarios <- n_pairs * n_values

  stats_list <- vector("list", total_scenarios)
  scenario_counter <- 0L

  for (j in seq_len(n_pairs)) {
    alpha_j <- pairs$alpha[[j]]
    n_people_j <- pairs$n_people[[j]]
    p_j <- p_list[[match(alpha_j, alpha)]]
    seed_j <- if (is.null(seed)) NULL else seed + j - 1L

    for (n in n_per_person) {
      scenario_counter <- scenario_counter + 1L

      rep_out <- run_replicates(
        p = p_j,
        B = B,
        metrics = metrics,
        model = "dirichlet_multinomial",
        seed = seed_j,
        n_people = n_people_j,
        n_per_person = n,
        concentration = concentration
      )

      metric_rows <- vector("list", length(metrics))
      for (mi in seq_along(metrics)) {
        m <- metrics[[mi]]
        metric_rows[[mi]] <- data.frame(
          alpha = alpha_j,
          n_people = n_people_j,
          concentration = concentration,
          n_per_person = n,
          metric = m,
          replicate = seq_len(B),
          stat = unname(rep_out$max_errors[, m]),
          stringsAsFactors = FALSE
        )
      }
      stats_list[[scenario_counter]] <- do.call(rbind, metric_rows)

      message(sprintf(
        "[%d/%d] alpha=%s n_people=%d n_per_person=%d",
        scenario_counter,
        total_scenarios,
        format(alpha_j, trim = TRUE),
        n_people_j,
        n
      ))
    }
  }

  p_table_list <- lapply(seq_along(alpha), function(i) {
    data.frame(
      alpha = alpha[[i]],
      as.list(stats::setNames(as.numeric(p_list[[i]]), paste0("cell_type_", seq_len(K)))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  list(
    inputs = list(
      alpha = alpha,
      K = K,
      B = B,
      metrics = metrics,
      proportion_method = proportion_method,
      n_people = n_people,
      n_per_person = n_per_person,
      concentration = concentration,
      seed = seed
    ),
    p_table = do.call(rbind, p_table_list),
    stats = do.call(rbind, stats_list)
  )
}


# Original Simulation ---------------------------------------------------------------------------------------------

#' Run the original full simulation experiment end-to-end.
#'
#' @param alpha      Numeric vector; one or more positive shape values used by the selected proportion-generation method
#'      (default method is Beta-based).
#' @param K          Number of cell types (default 10).
#' @param n          Total sample size per replicate for the multinomial
#'   model.
#' @param B          Number of replicates.
#' @param taus       Numeric vector of thresholds (same for all metrics) or a named list with one numeric vector per
#'      metric (e.g. `list(AE = c(...), ARE = c(...))`).
#' @param metrics    Error metrics; any subset of `c("AE", "ARE", "TSE", "LAE")` (only `"AE"` and `"ARE"` for
#'      `model = "dirichlet_multinomial"`).
#' @param proportion_method Proportion-generation method (`"beta"` or `"fixed_max_beta"`).
#'    The fixed-max Beta method places `p_max` at the highest index and warns then fails for impossible combinations.
#' @param p_max      Fixed largest true proportion(s) used by `proportion_method = "fixed_max_beta"`.
#'      Can be a numeric vector.
#'      When multiple values are provided, all alpha × p_max combinations are attempted; impossible fixed-max
#'      combinations are warned and skipped.
#' @param model      Sampling model (`"multinomial"` or
#'   `"dirichlet_multinomial"`).
#' @param n_people   Positive integer vector of people per replicate for the
#'   Dirichlet-multinomial model.
#' @param n_per_person Positive integer sampled-cell count for every person
#'   in the Dirichlet-multinomial model.
#' @param concentration Positive numeric vector of Dirichlet concentration
#'   values for the Dirichlet-multinomial model.
#' @param tie_method Tie-breaking rule for max-error argmax.
#' @param seed       Optional integer seed for reproducibility.
#' @param ...        Additional arguments forwarded to `simulate_counts()`.
#'
#' @return For `model = "dirichlet_multinomial"`, the return value of `run_dirichlet_multinomial_experiment()`
#'   (`inputs`, `p_table`, `person_results`). For `model = "multinomial"`, a list with elements:
#'   \describe{
#'     \item{inputs}{All input arguments.}
#'     \item{p_table}{Data.frame with one row per simulated alpha/p_max combination,
#'       an `alpha` column, a `p_max` column, and one column per cell type
#'       (`cell_type_1`, ..., `cell_type_K`) containing the corresponding p values.}
#'     \item{replicate_summaries}{Tidy data.frame:
#'       alpha, p_max, replicate, metric, max_error, argmax_index.}
#'     \item{errors_long}{Tidy data.frame:
#'       alpha, p_max, replicate, metric, cell_type, error.}
#'     \item{phat_long}{Tidy data.frame:
#'       alpha, p_max, replicate, cell_type, phat.}
#'     \item{curves}{Tidy data.frame:
#'       alpha, p_max, metric, tau, success_rate, mean_n_above.}
#'     \item{argmax_summary}{Tidy data.frame:
#'       alpha, p_max, metric, cell_type, count, fraction, true_proportion.}
#'   }
run_simulation_experiment <- function(
  alpha,
  K = 10,
  n = NULL,
  B,
  taus,
  metrics = c("AE", "ARE"),
  proportion_method = "beta",
  p_max = NULL,
  model = "multinomial",
  tie_method = "random",
  seed = NULL,
  n_people = NULL,
  n_per_person = NULL,
  concentration = NULL,
  ...
) {
  validate_positive_numeric(alpha, "alpha", allow_vector = TRUE)
  model <- match.arg(model, c("multinomial", "dirichlet_multinomial"))

  if (identical(model, "dirichlet_multinomial")) {
    return(run_dirichlet_multinomial_experiment(
      alpha = alpha,
      K = K,
      B = B,
      metrics = metrics,
      proportion_method = proportion_method,
      p_max = p_max,
      n_people = n_people,
      n_per_person = n_per_person,
      concentration = concentration,
      seed = seed
    ))
  }
  if (is.null(n)) {
    stop("n must be provided for model = 'multinomial'.", call. = FALSE)
  }

  feasible <- feasible_scenarios(alpha, K, proportion_method, p_max)
  combinations <- feasible$grid
  n_combinations <- nrow(combinations)
  p_table_list <- vector("list", n_combinations)
  replicate_summaries_list <- vector("list", n_combinations)
  errors_long_list <- vector("list", n_combinations)
  phat_long_list <- vector("list", n_combinations)
  curves_list <- vector("list", n_combinations)
  argmax_summary_list <- vector("list", n_combinations)
  keep <- logical(n_combinations)

  for (i in seq_len(n_combinations)) {
    alpha_i <- combinations$alpha[[i]]
    p_max_i <- combinations$p_max[[i]]
    seed_i <- if (is.null(seed)) NULL else seed + i - 1L
    p <- feasible$p[[i]]
    if (is.null(p)) {
      next
    } # Skip impossible alpha/p_max combinations (already warned in feasible_scenarios()).

    rep_out <- run_replicates(
      p,
      n,
      B,
      metrics = metrics,
      model = model,
      tie_method = tie_method,
      seed = seed_i,
      ...
    )
    keep[[i]] <- TRUE

    p_table_list[[i]] <- extract_p_table_row(alpha_i, p_max_i, p, K)

    replicate_summaries_list[[i]] <- extract_replicate_summaries(
      rep_out,
      alpha_i,
      p_max_i,
      B,
      metrics
    )

    errors_long_list[[i]] <- extract_errors_long(
      rep_out,
      alpha_i,
      p_max_i,
      B,
      metrics
    )

    phat_long_list[[i]] <- extract_phat_long(rep_out, alpha_i, p_max_i, B)

    curves_i <- evaluate_thresholds(
      rep_out$max_errors,
      taus,
      errors = rep_out$errors
    )
    curves_i$alpha <- alpha_i
    curves_i$p_max <- p_max_i
    curves_list[[i]] <- curves_i[, c(
      "alpha",
      "p_max",
      "metric",
      "tau",
      "success_rate",
      "mean_n_above"
    )]

    argmax_i <- summarize_argmax(rep_out$argmax, p)
    argmax_i$alpha <- alpha_i
    argmax_i$p_max <- p_max_i
    argmax_summary_list[[i]] <- argmax_i[, c(
      "alpha",
      "p_max",
      "metric",
      "cell_type",
      "count",
      "fraction",
      "true_proportion"
    )]
  }

  p_table_list <- p_table_list[keep]
  replicate_summaries_list <- replicate_summaries_list[keep]
  errors_long_list <- errors_long_list[keep]
  phat_long_list <- phat_long_list[keep]
  curves_list <- curves_list[keep]
  argmax_summary_list <- argmax_summary_list[keep]

  list(
    inputs = list(
      alpha = alpha,
      K = K,
      n = n,
      B = B,
      taus = taus,
      metrics = metrics,
      proportion_method = proportion_method,
      p_max = p_max,
      model = model,
      tie_method = tie_method,
      seed = seed
    ),
    p_table = do.call(rbind, p_table_list),
    replicate_summaries = do.call(rbind, replicate_summaries_list),
    errors_long = do.call(rbind, errors_long_list),
    phat_long = do.call(rbind, phat_long_list),
    curves = do.call(rbind, curves_list),
    argmax_summary = do.call(rbind, argmax_summary_list)
  )
}
