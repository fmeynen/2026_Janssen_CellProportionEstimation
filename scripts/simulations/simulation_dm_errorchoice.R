# simulation_dm_errorchoice.R
#
# Dirichlet-multinomial "errorchoice" simulation.
#
# Research question: for a population of persons whose latent cell-type composition varies around a population
# mean (Dirichlet-multinomial hierarchy), how often does the *maximum* (over cell types) absolute (AE) or relative
# (ARE) error of the person-pooled proportion estimate against the population proportion exceed a threshold tau --
# as a function of alpha (population-composition skew), n_people (persons per replicate) and n_per_person (cells
# sampled per person)? For cell type j: AE_j = |mean_i phat_ij - p_j|, ARE_j = AE_j / p_j.
#
# AE and ARE are studied separately throughout (never combined into a joint success criterion).
#
# Workflow
#   1. Generate population mean proportions from a monotone Beta curve (deterministic), once per alpha. Optionally
#      pin the smallest (largest) types to p_min (p_max) via proportion_method "fixed_min_beta" ("fixed_max_beta").
#   2. For each (alpha, n_people, n_per_person) scenario, simulate B Dirichlet-multinomial replicates and reduce
#      each replicate immediately to one scalar "stat" per metric (max over cell types of the pooled-estimate
#      error; see replicate_pooled_error()).
#      Common random numbers: the seed depends only on (alpha, n_people), identical across n_per_person.
#   3. curves_tau: with n_per_person held fixed at `n_per_person_fixed`, sweep a tau grid to get success rate vs
#      tau, per (alpha, n_people, metric).
#   4. curves_n: with tau held fixed at `taus_fixed[[metric]]`, sweep n_per_person (`n_per_person_grid`) to get
#      success rate vs cells-per-person, per (alpha, n_people, metric).
#   5. run_dm_errorchoice_samplesize(): with tau held fixed at `taus_fixed[[metric]]`, solve for the smallest
#      cells-per-person (n_per_person) reaching the success-rate `target`, per (alpha, n_people, metric), using the
#      sample-size solver (one solver run per n_people x metric, since the solver's success rule is joint over metrics).
#
# Possible future changes:
#   * Additional error metrics
# ---------------------------------------------------------------------------
source(here::here("scripts", "load_layers.R"))

# ---- Parameters ------------------------------------------------------------

#' Default configuration for the Dirichlet-multinomial errorchoice simulation.
#'
#' @return Named list of simulation parameters (see individual fields below).
#'   \describe{
#'     \item{alpha}{Beta shape parameter grid.}
#'     \item{K}{Number of cell types.}
#'     \item{B}{Replicates per scenario.}
#'     \item{metrics}{Error metrics studied separately (`"AE"`, `"ARE"`).}
#'     \item{proportion_method}{Proportion-generation method (`"beta"`, `"fixed_min_beta"`, `"fixed_max_beta"`).}
#'     \item{p_min}{Lower bound pinned by `"fixed_min_beta"`; `NULL` for plain beta.}
#'     \item{p_max}{Upper bound pinned by `"fixed_max_beta"`; `NULL` for plain beta.}
#'     \item{n_people}{Grid of persons per replicate.}
#'     \item{concentration}{Dirichlet concentration parameter.}
#'     \item{n_per_person_fixed}{Cells per person held fixed for the tau-sweep (`curves_tau`).}
#'     \item{n_per_person_grid}{Grid of cells per person for the n_per_person-sweep (`curves_n`).}
#'     \item{taus_fixed}{Named list; per-metric tau held fixed for `curves_n`.}
#'     \item{taus}{Named list; per-metric tau grid override for `curves_tau`. A non-`NULL` `taus[[metric]]`
#'       overrides the data-driven tau grid (`default_tau_grid()`) for that metric.}
#'     \item{tau_grid_points}{Number of points in the data-driven tau grid.}
#'     \item{tau_grid_prob}{Quantile probability used to set the data-driven tau grid's upper end.}
#'     \item{target}{Success-rate reference line used by the `curves_n` plot.}
#'     \item{seed}{Base seed; see `run_dm_errorchoice_experiment()` for how it is combined with (alpha, n_people).}
#'     \item{tie_method}{Tie handling in the sample-size solver's success rule (used only by
#'       `run_dm_errorchoice_samplesize()`).}
#'     \item{rel_tol}{Sample-size solver relative tolerance on the success rate.}
#'     \item{max_iterations}{Sample-size solver iteration cap.}
#'     \item{f0}{Sample-size solver initial bracket factor.}
#'     \item{f_floor}{Sample-size solver minimum bracket factor.}
#'     \item{n_max}{Largest cells-per-person value the sample-size solver considers; also used for the feasibility
#'       check.}
#'     \item{n_init}{Sample-size solver starting value; `NULL` starts from `concentration`.}
#'   }
simulation_dm_errorchoice_defaults <- function() {
  list(
    alpha = c(2, 3, 4, 5),
    K = 10L,
    B = 1000L,
    metrics = c("AE", "ARE"),
    proportion_method = "beta",
    p_min = NULL,
    p_max = NULL,
    n_people = c(1L, 2L, 3L, 5L, 10L),
    concentration = 1e4,
    n_per_person_fixed = 200000L,
    n_per_person_grid = as.integer(round(10^seq(2, 8, by = 0.25))),
    taus_fixed = list(AE = 0.01, ARE = 0.5),
    taus = list(AE = NULL, ARE = NULL),
    tau_grid_points = 200L,
    tau_grid_prob = 0.95,
    target = 0.95,
    seed = 260926L,
    tie_method = "random",
    rel_tol = 0.01,
    max_iterations = 30L,
    f0 = 2,
    f_floor = 1.1,
    n_max = 1e9,
    n_init = NULL
  )
}

#' Errorchoice configuration with the smallest cell types pinned to `p_min`.
#'
#' Same as `simulation_dm_errorchoice_defaults()` except `proportion_method = "fixed_min_beta"` and `p_min = 0.01`.
#'
#' @return Named list of simulation parameters; see `simulation_dm_errorchoice_defaults()`.
simulation_dm_errorchoice_pmin_defaults <- function() {
  utils::modifyList(
    simulation_dm_errorchoice_defaults(),
    list(proportion_method = "fixed_min_beta", p_min = 0.01)
  )
}

# ---- Run experiment --------------------------------------------------------

#' Run (or load from cache) the Dirichlet-multinomial errorchoice simulation.
#'
#' Runs `run_dm_errorchoice_experiment()` over `n_values <- sort(unique(c(config$n_per_person_grid,
#' config$n_per_person_fixed)))`, then derives two tidy success-rate views from the resulting `stats` table:
#' `curves_tau` (tau-sweep at `n_per_person_fixed`) and `curves_n` (n_per_person-sweep at `taus_fixed`). The cache
#' key only includes simulation-relevant fields (`alpha`, `K`, `B`, `metrics`, `proportion_method`, `p_min`, `p_max`,
#' `n_people`, `concentration`, `n_per_person = n_values`, `seed`), so changing `taus`, `taus_fixed`, `target`,
#' `tau_grid_points` or `tau_grid_prob` reuses the same cached simulation result.
#'
#' @param config    List as returned by `simulation_dm_errorchoice_defaults()`.
#' @param cache     Logical; read/write the cached simulation result under `cache_dir`.
#' @param force_recompute Logical; ignore any existing cache file and recompute (still writes the new result when
#'   `cache` is `TRUE`).
#' @param cache_dir Directory holding the cached `.rds` file.
#'
#' @return List with elements:
#'   \describe{
#'     \item{inputs}{`config`, as passed in.}
#'     \item{p_table}{From `run_dm_errorchoice_experiment()`: one row per alpha, columns `alpha`, `cell_type_1`, ...,
#'       `cell_type_K`.}
#'     \item{stats}{From `run_dm_errorchoice_experiment()`: one row per (alpha, n_people, n_per_person, metric,
#'       replicate).}
#'     \item{curves_tau}{Data.frame: `alpha`, `n_people`, `metric`, `tau`, `success_rate`, at
#'       `n_per_person = config$n_per_person_fixed`.}
#'     \item{curves_n}{Data.frame: `alpha`, `n_people`, `n_per_person`, `metric`, `tau`, `success_rate`, over
#'       `config$n_per_person_grid` at each metric's `config$taus_fixed` threshold.}
#'   }
run_simulation_dm_errorchoice <- function(config = simulation_dm_errorchoice_defaults(),
                                          cache = TRUE,
                                          force_recompute = FALSE,
                                          cache_dir = here::here("results", "simresults")) {
  n_values <- sort(unique(c(config$n_per_person_grid, config$n_per_person_fixed)))

  # Only simulation-relevant fields go into the cache key: changing taus/taus_fixed/target/tau_grid_* below re-uses
  # the same cached run_dm_errorchoice_experiment() result.
  sim_config <- list(
    alpha = config$alpha,
    K = config$K,
    B = config$B,
    metrics = config$metrics,
    proportion_method = config$proportion_method,
    p_min = config$p_min,
    p_max = config$p_max,
    n_people = config$n_people,
    concentration = config$concentration,
    n_per_person = n_values,
    seed = config$seed,
    success_rule = success_rule_id()
  )
  sim_result <- cached_result(
    key = sim_config,
    name = "dm_errorchoice",
    compute = function() {
      run_dm_errorchoice_experiment(
        alpha = config$alpha,
        K = config$K,
        B = config$B,
        metrics = config$metrics,
        proportion_method = config$proportion_method,
        p_min = config$p_min,
        p_max = config$p_max,
        n_people = config$n_people,
        n_per_person = n_values,
        concentration = config$concentration,
        seed = config$seed
      )
    },
    cache = cache,
    force_recompute = force_recompute,
    dir = cache_dir
  )

  stats <- sim_result$stats

  # ---- curves_tau: tau-sweep at n_per_person_fixed --------------------------
  fixed_stats <- stats[stats$n_per_person == config$n_per_person_fixed, , drop = FALSE]
  curves_tau_list <- lapply(config$metrics, function(m) {
    fixed_stats_m <- fixed_stats[fixed_stats$metric == m, , drop = FALSE]
    tau_grid <- config$taus[[m]]
    if (is.null(tau_grid)) {
      tau_grid <- default_tau_grid(
        fixed_stats_m$stat,
        n_points = config$tau_grid_points,
        prob = config$tau_grid_prob
      )
    }
    groups <- split(fixed_stats_m, list(fixed_stats_m$alpha, fixed_stats_m$n_people), drop = TRUE)
    group_rows <- lapply(groups, function(g) {
      data.frame(
        alpha = g$alpha[[1L]],
        n_people = g$n_people[[1L]],
        metric = m,
        tau = tau_grid,
        success_rate = success_rate_from_stat(g$stat, tau_grid),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, group_rows)
  })
  curves_tau <- do.call(rbind, curves_tau_list)
  rownames(curves_tau) <- NULL

  # ---- curves_n: n_per_person-sweep at taus_fixed ----------------------------
  grid_stats <- stats[stats$n_per_person %in% config$n_per_person_grid, , drop = FALSE]
  curves_n_list <- lapply(config$metrics, function(m) {
    grid_stats_m <- grid_stats[grid_stats$metric == m, , drop = FALSE]
    tau_m <- config$taus_fixed[[m]]
    groups <- split(
      grid_stats_m,
      list(grid_stats_m$alpha, grid_stats_m$n_people, grid_stats_m$n_per_person),
      drop = TRUE
    )
    group_rows <- lapply(groups, function(g) {
      data.frame(
        alpha = g$alpha[[1L]],
        n_people = g$n_people[[1L]],
        n_per_person = g$n_per_person[[1L]],
        metric = m,
        tau = tau_m,
        success_rate = success_rate_from_stat(g$stat, tau_m),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, group_rows)
  })
  curves_n <- do.call(rbind, curves_n_list)
  rownames(curves_n) <- NULL

  list(
    inputs = config,
    p_table = sim_result$p_table,
    stats = stats,
    curves_tau = curves_tau,
    curves_n = curves_n
  )
}

# ---- Required cells per person (sample-size solver) ------------------------

#' Solve for the cells per person needed to reach the target success rate, per (alpha, n_people, metric).
#'
#' Runs `run_sample_size_experiment()` once per (`n_people`, metric) pair, because the solver's success rule is joint
#' over all metrics in its `taus`; each run therefore gets exactly one metric, with `taus = config$taus_fixed[metric]`.
#' Every run uses `config$seed` (common random numbers across n_people and metrics). Alphas whose target cannot be
#' reached emit the solver's warning (not suppressed here) and are kept as rows with `sample_size` `NA` and
#' `stopping_reason` `"infeasible"`. Caching is per alpha inside `run_sample_size_experiment()`.
#'
#' @param config    List as returned by `simulation_dm_errorchoice_defaults()`.
#' @param cache     Logical; read/write the per-alpha cached results under `cache_dir`.
#' @param force_recompute Logical; ignore existing cache files and recompute.
#' @param cache_dir Directory holding the cached `.rds` files.
#' @param simulate  Function `(alpha, n, config, seed)` forwarded to `run_sample_size_experiment()`. Defaults to
#'   `simulate_success_at_n()`.
#'
#' @return Data.frame with columns `alpha`, `n_people`, `metric`, `sample_size` (cells per person; `NA` if
#'   infeasible), `stopping_reason`, `iterations_used`, `success_ceiling`.
run_dm_errorchoice_samplesize <- function(config = simulation_dm_errorchoice_defaults(),
                                          cache = TRUE,
                                          force_recompute = FALSE,
                                          cache_dir = here::here("results", "simresults"),
                                          simulate = simulate_success_at_n) {
  rows <- list()
  for (n_people in config$n_people) {
    for (metric in config$metrics) {
      solver_config <- list(
        alpha = config$alpha,
        K = config$K,
        B = config$B,
        taus = config$taus_fixed[metric],
        metrics = metric,
        model = "dirichlet_multinomial",
        tie_method = config$tie_method,
        proportion_method = config$proportion_method,
        p_min = config$p_min,
        p_max = config$p_max,
        n_people = n_people,
        concentration = config$concentration,
        seed = config$seed,
        success_rate_target = config$target,
        rel_tol = config$rel_tol,
        max_iterations = config$max_iterations,
        f0 = config$f0,
        f_floor = config$f_floor,
        n_max = config$n_max,
        n_init = config$n_init
      )
      res <- run_sample_size_experiment(
        solver_config,
        cache = cache,
        force_recompute = force_recompute,
        cache_dir = cache_dir,
        simulate = simulate
      )
      ss <- res$sample_size
      rows[[length(rows) + 1L]] <- data.frame(
        alpha = ss$alpha,
        n_people = n_people,
        metric = metric,
        sample_size = ss$sample_size,
        stopping_reason = ss$stopping_reason,
        iterations_used = ss$iterations_used,
        success_ceiling = ss$success_ceiling,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
