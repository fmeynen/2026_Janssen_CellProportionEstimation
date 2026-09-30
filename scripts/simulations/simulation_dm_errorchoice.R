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
#   1. Generate population mean proportions from a monotone Beta curve (deterministic), once per alpha.
#   2. For each (alpha, n_people, n_per_person) scenario, simulate B Dirichlet-multinomial replicates and reduce
#      each replicate immediately to one scalar "stat" per metric (max over cell types of the pooled-estimate
#      error; see replicate_pooled_error()).
#      Common random numbers: the seed depends only on (alpha, n_people), identical across n_per_person.
#   3. curves_tau: with n_per_person held fixed at `n_per_person_fixed`, sweep a tau grid to get success rate vs
#      tau, per (alpha, n_people, metric).
#   4. curves_n: with tau held fixed at `taus_fixed[[metric]]`, sweep n_per_person (`n_per_person_grid`) to get
#      success rate vs cells-per-person, per (alpha, n_people, metric).
#
# Possible future changes:
#   * Additional error metrics
#   * Non-beta proportion-generation methods
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
#'     \item{proportion_method}{Proportion-generation method.}
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
#'   }
simulation_dm_errorchoice_defaults <- function() {
  list(
    alpha = c(2, 3, 4, 5),
    K = 10L,
    B = 1000L,
    metrics = c("AE", "ARE"),
    proportion_method = "beta",
    n_people = c(1L, 2L, 3L, 5L, 10L),
    concentration = 1e4,
    n_per_person_fixed = 200000L,
    n_per_person_grid = as.integer(round(10^seq(2, 8, by = 0.25))),
    taus_fixed = list(AE = 0.005, ARE = 0.5),
    taus = list(AE = NULL, ARE = NULL),
    tau_grid_points = 200L,
    tau_grid_prob = 0.95,
    target = 0.95,
    seed = 260926L
  )
}

# ---- Run experiment --------------------------------------------------------

#' Run (or load from cache) the Dirichlet-multinomial errorchoice simulation.
#'
#' Runs `run_dm_errorchoice_experiment()` over `n_values <- sort(unique(c(config$n_per_person_grid,
#' config$n_per_person_fixed)))`, then derives two tidy success-rate views from the resulting `stats` table:
#' `curves_tau` (tau-sweep at `n_per_person_fixed`) and `curves_n` (n_per_person-sweep at `taus_fixed`). The cache
#' key only includes simulation-relevant fields (`alpha`, `K`, `B`, `metrics`, `proportion_method`, `n_people`,
#' `concentration`, `n_per_person = n_values`, `seed`), so changing `taus`, `taus_fixed`, `target`,
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
