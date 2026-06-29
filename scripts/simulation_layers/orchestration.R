# R/orchestration.R
#
# Orchestration layer: high-level experiment runners that coordinate all other layers.
#
# Depends on: R/validation_utils.R, R/simulation.R, R/calculation.R, R/extraction.R
# ---------------------------------------------------------------------------


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
#' @param metrics    Error metrics; any subset of `c("AE", "ARE", "TSE", "LAE")`.
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
#'     \item{p_table}{Data.frame with one row per simulated alpha/p_max combination,
#'       an `alpha` column, a `p_max` column, and one column per index
#'       (`index_1`, ..., `index_K`) containing the corresponding p values.}
#'     \item{replicate_summaries}{Tidy data.frame:
#'       alpha, p_max, replicate, metric, max_error, argmax_index.}
#'     \item{errors_long}{Tidy data.frame:
#'       alpha, p_max, replicate, metric, index, error.}
#'     \item{phat_long}{Tidy data.frame:
#'       alpha, p_max, replicate, index, phat.}
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
  phat_long_list <- vector("list", n_combinations)
  curves_list <- vector("list", n_combinations)
  argmax_summary_list <- vector("list", n_combinations)
  keep <- logical(n_combinations)
  should_skip_impossible_combinations <- identical(proportion_method, "fixed_max_beta") &&
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
        if (should_skip_impossible_combinations && inherits(e, "impossible_fixed_max_error")) {
          return(NULL)
        }
        stop(e)
      }
    )
    if (is.null(p)) next # Skip impossible alpha/p_max combinations caught by the error handler.

    rep_out <- run_replicates(
      p, n, B, metrics = metrics, model = model,
      tie_method = tie_method, seed = seed_i, ...
    )
    keep[[i]] <- TRUE

    p_table_list[[i]] <- extract_p_table_row(alpha_i, p_max_i, p, K)

    replicate_summaries_list[[i]] <- extract_replicate_summaries(
      rep_out, alpha_i, p_max_i, B, metrics
    )

    errors_long_list[[i]] <- extract_errors_long(rep_out, alpha_i, p_max_i, B, metrics)

    phat_long_list[[i]] <- extract_phat_long(rep_out, alpha_i, p_max_i, B)

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
  phat_long_list <- phat_long_list[keep]
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
    phat_long = do.call(rbind, phat_long_list),
    curves = do.call(rbind, curves_list),
    argmax_summary = do.call(rbind, argmax_summary_list)
  )
}


#' Run simulation and hybrid cutoff analysis across all alpha/p_max scenarios.
#'
#' @param alpha             Numeric vector of positive Beta shape values.
#' @param K                 Number of cell types.
#' @param n                 Total sample size per replicate.
#' @param B                 Number of replicates.
#' @param cutoffs           Numeric vector of candidate cutoffs.
#' @param tau_AE            Numeric scalar AE success threshold.
#' @param tau_ARE           Numeric scalar ARE success threshold.
#' @param proportion_method Proportion-generation method (`"beta"` or
#'   `"fixed_max_beta"`).
#' @param p_max             Optional fixed largest true proportion(s).
#' @param model             Sampling model passed to `run_simulation_experiment()`.
#' @param tie_method        Tie-breaking rule for max-error argmax.
#' @param maximize          Which success rate to maximize (`"cell"` or
#'   `"replicate"`).
#' @param seed              Optional integer seed for reproducibility.
#' @param ...               Additional arguments forwarded to
#'   `run_simulation_experiment()`.
#'
#' @return List with elements:
#'   \describe{
#'     \item{inputs}{Copy of all function inputs.}
#'     \item{p_table}{Scenario-level true proportions table.}
#'     \item{phat_long}{Tidy replicate-level observed proportions.}
#'     \item{cutoff_curves}{Combined hybrid cutoff curves across scenarios.}
#'     \item{best_cutoff_summary}{One-row-per-scenario best-cutoff summary.}
#'   }
run_simulation_hybrid_cutoff_experiment <- function(alpha, K, n, B, cutoffs,
                                                     tau_AE, tau_ARE,
                                                     proportion_method = "beta",
                                                     p_max = NULL,
                                                     model = "multinomial",
                                                     tie_method = "random",
                                                     maximize = c("cell", "replicate"),
                                                     seed = NULL, ...) {
  maximize <- match.arg(maximize)
  sim_out <- run_simulation_experiment(
    alpha = alpha,
    K = K,
    n = n,
    B = B,
    taus = list(AE = tau_AE, ARE = tau_ARE),
    metrics = c("AE", "ARE"),
    proportion_method = proportion_method,
    p_max = p_max,
    model = model,
    tie_method = tie_method,
    seed = seed,
    ...
  )

  p_table <- sim_out$p_table
  index_cols <- grep("^index_", names(p_table), value = TRUE)
  cutoff_curves_list <- vector("list", nrow(p_table))
  best_summary_list <- vector("list", nrow(p_table))

  for (i in seq_len(nrow(p_table))) {
    alpha_i <- p_table$alpha[[i]]
    p_max_i <- if ("p_max" %in% names(p_table)) p_table$p_max[[i]] else NA_real_
    p_i <- as.numeric(p_table[i, index_cols, drop = FALSE])

    same_p_max <- if (is.na(p_max_i)) {
      is.na(sim_out$phat_long$p_max)
    } else {
      sim_out$phat_long$p_max == p_max_i
    }
    phat_subset <- sim_out$phat_long[
      sim_out$phat_long$alpha == alpha_i & same_p_max,
      c("replicate", "index", "phat"),
      drop = FALSE
    ]
    phat_subset <- phat_subset[order(phat_subset$index, phat_subset$replicate), , drop = FALSE]
    phat_mat_i <- matrix(phat_subset$phat, nrow = B, ncol = length(p_i))

    hybrid_i <- run_hybrid_cutoff_analysis(
      phat_mat = phat_mat_i,
      p = p_i,
      cutoffs = cutoffs,
      tau_AE = tau_AE,
      tau_ARE = tau_ARE,
      maximize = maximize,
      alpha = alpha_i,
      p_max = p_max_i
    )

    curve_i <- hybrid_i$cutoff_curve
    curve_i$alpha <- alpha_i
    curve_i$p_max <- p_max_i
    curve_i$tau_AE <- tau_AE
    curve_i$tau_ARE <- tau_ARE
    cutoff_curves_list[[i]] <- curve_i[, c(
      "alpha", "p_max", "cutoff", "tau_AE", "tau_ARE",
      "success_rate_cell", "success_rate_replicate",
      "prop_using_AE", "prop_using_ARE"
    )]

    best_summary_list[[i]] <- hybrid_i$best_summary
  }

  list(
    inputs = list(
      alpha = alpha, K = K, n = n, B = B, cutoffs = cutoffs,
      tau_AE = tau_AE, tau_ARE = tau_ARE,
      proportion_method = proportion_method, p_max = p_max,
      model = model, tie_method = tie_method, maximize = maximize, seed = seed
    ),
    p_table = sim_out$p_table,
    phat_long = sim_out$phat_long,
    cutoff_curves = do.call(rbind, cutoff_curves_list),
    best_cutoff_summary = do.call(rbind, best_summary_list)
  )
}


# ---------------------------------------------------------------------------
# Iso-probability pipeline
# ---------------------------------------------------------------------------


#' Build a Cartesian candidate design over (tau_AE, tau_ARE, cutoff) grids.
#'
#' @param tau_AE_grid  Numeric vector of candidate AE thresholds.
#' @param tau_ARE_grid Numeric vector of candidate ARE thresholds.
#' @param cutoff_grid  Numeric vector of candidate cutoffs.
#'
#' @return data.frame with columns:
#'   `candidate_id` (integer), `tau_AE` (numeric), `tau_ARE` (numeric),
#'   `cutoff` (numeric).
make_iso_parameter_design <- function(tau_AE_grid, tau_ARE_grid, cutoff_grid) {
  stopifnot(
    is.numeric(tau_AE_grid),  length(tau_AE_grid)  >= 1L, all(is.finite(tau_AE_grid)),
    is.numeric(tau_ARE_grid), length(tau_ARE_grid) >= 1L, all(is.finite(tau_ARE_grid)),
    is.numeric(cutoff_grid),  length(cutoff_grid)  >= 1L, all(is.finite(cutoff_grid))
  )
  design <- expand.grid(
    tau_AE  = tau_AE_grid,
    tau_ARE = tau_ARE_grid,
    cutoff  = cutoff_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  design <- design[, c("tau_AE", "tau_ARE", "cutoff"), drop = FALSE]
  design$candidate_id <- seq_len(nrow(design))
  design[, c("candidate_id", "tau_AE", "tau_ARE", "cutoff")]
}


#' Run the iso-probability success search over a candidate parameter design.
#'
#' Estimates the baseline joint success probability for a fixed baseline triple
#' (tau_AE, tau_ARE, cutoff), then evaluates all candidate triples in the
#' design and flags those that are equivalent to the baseline (point rule:
#' |p_hat - p0_hat| <= delta).
#'
#' **Single-scenario constraint**: this function accepts exactly one alpha value
#' and at most one p_max value so that a single true-proportion vector `p` is
#' unambiguously defined.  To search over multiple proportion scenarios, call
#' this function once per scenario and combine results externally.
#'
#' @param alpha              Numeric scalar (> 0); Beta shape for proportion generation.
#' @param K                  Integer; number of cell types (default 10).
#' @param n                  Integer; total sample size per replicate.
#' @param baseline_tau_AE    Numeric scalar; baseline AE threshold.
#' @param baseline_tau_ARE   Numeric scalar; baseline ARE threshold.
#' @param baseline_cutoff    Numeric scalar; baseline cutoff value.
#' @param B_baseline         Integer; replicates used to estimate baseline p0.
#' @param design             data.frame from `make_iso_parameter_design()`.
#' @param B_screen           Integer; replicates used to evaluate each candidate.
#' @param delta              Numeric scalar >= 0; absolute equivalence tolerance
#'   for the point rule `|p_hat - p0_hat| <= delta`.
#' @param proportion_method  Proportion-generation method (`"beta"` or
#'   `"fixed_max_beta"`).
#' @param p_max              Numeric scalar (or NULL); fixed largest proportion
#'   for `"fixed_max_beta"`.
#' @param model              Sampling model passed to `run_replicates()`.
#' @param tie_method         Tie-breaking rule passed to `run_replicates()`.
#' @param seed               Optional integer seed; if provided, baseline uses
#'   `seed` and candidate i uses `seed + i`.
#' @param conf_level         Numeric scalar in (0, 1); CI confidence level.
#' @param ci_method          Character scalar: `"wilson"` or `"jeffreys"`.
#' @param ...                Additional arguments forwarded to `run_replicates()`.
#'
#' @details
#' Runtime scales as O(nrow(design) * B_screen) replicates plus O(B_baseline)
#' for the baseline.  For a grid of G candidates with B_screen replicates each,
#' total replicates simulated is B_baseline + G * B_screen.  Keep grid sizes
#' manageable (e.g., a few hundred candidates) to avoid long runtimes.
#'
#' @return Named list with elements:
#'   \describe{
#'     \item{inputs}{List of all input arguments.}
#'     \item{baseline}{One-row data.frame:
#'       tau_AE, tau_ARE, cutoff, n_success, n_total, p_hat (= p0_hat),
#'       se_mc, ci_low, ci_high, conf_level, ci_method.}
#'     \item{screen_results}{data.frame with one row per candidate:
#'       candidate_id, tau_AE, tau_ARE, cutoff, n_success, n_total,
#'       p_hat, se_mc, ci_low, ci_high, abs_diff_p0, pass.}
#'     \item{iso_candidates}{Subset of screen_results where pass == TRUE.}
#'     \item{diagnostics}{List: n_candidates_total, n_pass_screen, runtime_sec.}
#'   }
run_iso_success_search <- function(alpha, K = 10L, n,
                                   baseline_tau_AE, baseline_tau_ARE, baseline_cutoff,
                                   B_baseline, design, B_screen,
                                   delta = 0.01,
                                   proportion_method = "beta",
                                   p_max = NULL,
                                   model = "multinomial",
                                   tie_method = "random",
                                   seed = NULL,
                                   conf_level = 0.95,
                                   ci_method = c("wilson", "jeffreys"),
                                   ...) {
  ci_method <- match.arg(ci_method)

  stopifnot(
    is.numeric(alpha),  length(alpha) == 1L,  is.finite(alpha),  alpha > 0,
    is.numeric(K),      length(K) == 1L,      K %% 1 == 0,       K >= 1L,
    is.numeric(n),      length(n) == 1L,      is.finite(n),      n >= 1L,
    is.numeric(baseline_tau_AE),   length(baseline_tau_AE)   == 1L, is.finite(baseline_tau_AE),
    is.numeric(baseline_tau_ARE),  length(baseline_tau_ARE)  == 1L, is.finite(baseline_tau_ARE),
    is.numeric(baseline_cutoff),   length(baseline_cutoff)   == 1L, is.finite(baseline_cutoff),
    is.numeric(B_baseline), length(B_baseline) == 1L, B_baseline >= 1L,
    is.data.frame(design),
    all(c("candidate_id", "tau_AE", "tau_ARE", "cutoff") %in% names(design)),
    is.numeric(B_screen), length(B_screen) == 1L, B_screen >= 1L,
    is.numeric(delta), length(delta) == 1L, is.finite(delta), delta >= 0
  )

  t_start <- proc.time()[["elapsed"]]

  # Generate a single true-proportion vector for this scenario.
  p <- generate_proportions(
    alpha  = alpha,
    K      = as.integer(K),
    method = proportion_method,
    p_max  = p_max
  )

  # --- Baseline ---
  seed_baseline <- if (is.null(seed)) NULL else seed
  rep_baseline  <- run_replicates(
    p, n, B_baseline,
    metrics    = c("AE", "ARE"),
    model      = model,
    tie_method = tie_method,
    seed       = seed_baseline,
    ...
  )
  success_baseline <- compute_joint_success(
    phat_mat = rep_baseline$phat,
    p        = p,
    tau_AE   = baseline_tau_AE,
    tau_ARE  = baseline_tau_ARE,
    cutoff   = baseline_cutoff
  )
  summ_baseline <- summarize_joint_success(success_baseline, conf_level, ci_method)

  baseline_df <- data.frame(
    tau_AE     = baseline_tau_AE,
    tau_ARE    = baseline_tau_ARE,
    cutoff     = baseline_cutoff,
    n_success  = summ_baseline$n_success,
    n_total    = summ_baseline$n_total,
    p_hat      = summ_baseline$p_hat,
    se_mc      = summ_baseline$se_mc,
    ci_low     = summ_baseline$ci_low,
    ci_high    = summ_baseline$ci_high,
    conf_level = summ_baseline$conf_level,
    ci_method  = summ_baseline$ci_method,
    stringsAsFactors = FALSE
  )

  p0_hat <- summ_baseline$p_hat

  # --- Screen candidates ---
  n_candidates <- nrow(design)
  screen_rows  <- vector("list", n_candidates)

  for (i in seq_len(n_candidates)) {
    seed_i <- if (is.null(seed)) NULL else seed + i
    rep_i  <- run_replicates(
      p, n, B_screen,
      metrics    = c("AE", "ARE"),
      model      = model,
      tie_method = tie_method,
      seed       = seed_i,
      ...
    )
    success_i <- compute_joint_success(
      phat_mat = rep_i$phat,
      p        = p,
      tau_AE   = design$tau_AE[[i]],
      tau_ARE  = design$tau_ARE[[i]],
      cutoff   = design$cutoff[[i]]
    )
    summ_i    <- summarize_joint_success(success_i, conf_level, ci_method)
    equiv_i   <- classify_iso_equivalence(summ_i$p_hat, p0_hat, delta)

    screen_rows[[i]] <- data.frame(
      candidate_id = design$candidate_id[[i]],
      tau_AE       = design$tau_AE[[i]],
      tau_ARE      = design$tau_ARE[[i]],
      cutoff       = design$cutoff[[i]],
      n_success    = summ_i$n_success,
      n_total      = summ_i$n_total,
      p_hat        = summ_i$p_hat,
      se_mc        = summ_i$se_mc,
      ci_low       = summ_i$ci_low,
      ci_high      = summ_i$ci_high,
      abs_diff_p0  = equiv_i$abs_diff_p0,
      pass         = equiv_i$pass,
      stringsAsFactors = FALSE
    )
  }

  screen_results <- do.call(rbind, screen_rows)
  iso_candidates <- screen_results[screen_results$pass, , drop = FALSE]

  runtime_sec <- proc.time()[["elapsed"]] - t_start

  list(
    inputs = list(
      alpha             = alpha,
      K                 = K,
      n                 = n,
      baseline_tau_AE   = baseline_tau_AE,
      baseline_tau_ARE  = baseline_tau_ARE,
      baseline_cutoff   = baseline_cutoff,
      B_baseline        = B_baseline,
      B_screen          = B_screen,
      delta             = delta,
      proportion_method = proportion_method,
      p_max             = p_max,
      model             = model,
      tie_method        = tie_method,
      seed              = seed,
      conf_level        = conf_level,
      ci_method         = ci_method
    ),
    baseline       = baseline_df,
    screen_results = screen_results,
    iso_candidates = iso_candidates,
    diagnostics    = list(
      n_candidates_total = n_candidates,
      n_pass_screen      = sum(screen_results$pass),
      runtime_sec        = runtime_sec
    )
  )
}
