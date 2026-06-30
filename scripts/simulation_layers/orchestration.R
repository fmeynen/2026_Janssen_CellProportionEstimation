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

#' Generate initial LHS design for isoband automation.
#'
#' @param ranges Named list of validated ranges for `tau_AE`, `tau_ARE`, `cutoff`.
#' @param n_init Number of initial design points.
#' @param seed   Optional integer seed for reproducibility.
#'
#' @return Data frame with columns `tau_AE`, `tau_ARE`, `cutoff`, `round`, `stage`.
generate_isoband_initial_design <- function(ranges, n_init = 200, seed = NULL) {
  ranges <- validate_isoband_ranges(ranges)
  if (!is.numeric(n_init) || length(n_init) != 1L || !is.finite(n_init) || n_init < 1 || n_init %% 1 != 0) {
    stop("n_init must be a positive integer.", call. = FALSE)
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or a finite numeric scalar.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  lhs_unit <- lhs::randomLHS(n = as.integer(n_init), k = 3L)
  design_df <- data.frame(
    tau_AE = ranges$tau_AE[[1L]] + lhs_unit[, 1L] * diff(ranges$tau_AE),
    tau_ARE = ranges$tau_ARE[[1L]] + lhs_unit[, 2L] * diff(ranges$tau_ARE),
    cutoff = ranges$cutoff[[1L]] + lhs_unit[, 3L] * diff(ranges$cutoff),
    round = 0L,
    stage = "init",
    stringsAsFactors = FALSE
  )
  design_df
}

#' Generate candidate pool for isoband surrogate prediction.
#'
#' @param ranges       Named list of validated ranges.
#' @param n_candidates Number of candidate points.
#' @param seed         Optional integer seed for reproducibility.
#'
#' @return Data frame with columns `tau_AE`, `tau_ARE`, and `cutoff`.
generate_isoband_candidate_pool <- function(ranges, n_candidates = 5000, seed = NULL) {
  ranges <- validate_isoband_ranges(ranges)
  if (!is.numeric(n_candidates) || length(n_candidates) != 1L ||
    !is.finite(n_candidates) || n_candidates < 1 || n_candidates %% 1 != 0) {
    stop("n_candidates must be a positive integer.", call. = FALSE)
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or a finite numeric scalar.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  lhs_unit <- lhs::randomLHS(n = as.integer(n_candidates), k = 3L)
  data.frame(
    tau_AE = ranges$tau_AE[[1L]] + lhs_unit[, 1L] * diff(ranges$tau_AE),
    tau_ARE = ranges$tau_ARE[[1L]] + lhs_unit[, 2L] * diff(ranges$tau_ARE),
    cutoff = ranges$cutoff[[1L]] + lhs_unit[, 3L] * diff(ranges$cutoff),
    stringsAsFactors = FALSE
  )
}

#' Simulate replicate-level hybrid success for a set of isoband design points.
#'
#' @param design_df Data frame containing `tau_AE`, `tau_ARE`, and `cutoff`.
#' @param alpha Scalar alpha used to generate true proportions.
#' @param K Number of cell types.
#' @param n Total sample size per replicate.
#' @param B_point Number of replicates per design point.
#' @param proportion_method Proportion-generation method (`"beta"` or `"fixed_max_beta"`).
#' @param p_max Optional fixed largest true proportion.
#' @param model Sampling model passed to `run_replicates()`.
#' @param seed Optional integer seed.
#' @param ... Additional arguments passed to `run_replicates()`.
#'
#' @return Data frame with one row per design point and replicate-level success summaries.
simulate_isoband_design_points <- function(design_df, alpha, K, n, B_point,
                                           proportion_method = "beta",
                                           p_max = NULL,
                                           model = "multinomial",
                                           seed = NULL, ...) {
  required_cols <- c("tau_AE", "tau_ARE", "cutoff")
  if (!is.data.frame(design_df) || !all(required_cols %in% names(design_df))) {
    stop(
      sprintf("design_df must contain columns: %s", paste(required_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
    stop("alpha must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 1 || K %% 1 != 0) {
    stop("K must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n %% 1 != 0) {
    stop("n must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(B_point) || length(B_point) != 1L || !is.finite(B_point) ||
    B_point < 1 || B_point %% 1 != 0) {
    stop("B_point must be a positive integer.", call. = FALSE)
  }

  p <- generate_proportions(
    alpha = alpha,
    K = as.integer(K),
    method = proportion_method,
    p_max = p_max
  )
  n_design <- nrow(design_df)
  if (n_design == 0L) {
    out <- design_df
    out$n_success <- integer(0L)
    out$n_total <- integer(0L)
    out$success_rate <- numeric(0L)
    return(out)
  }

  rows <- lapply(seq_len(n_design), function(i) {
    point_i <- design_df[i, , drop = FALSE]
    seed_i <- if (is.null(seed)) NULL else as.integer(seed) + i - 1L
    rep_out <- run_replicates(
      p = p,
      n = as.integer(n),
      B = as.integer(B_point),
      metrics = c("AE", "ARE"),
      model = model,
      seed = seed_i,
      ...
    )
    eval_i <- evaluate_hybrid_success_cell_level_strict(
      phat_mat = rep_out$phat,
      p = p,
      cutoff = point_i$cutoff[[1L]],
      tau_AE = point_i$tau_AE[[1L]],
      tau_ARE = point_i$tau_ARE[[1L]]
    )
    n_success_i <- sum(rowSums(eval_i$success_matrix, na.rm = TRUE) == ncol(eval_i$success_matrix))
    n_total_i <- nrow(eval_i$success_matrix)

    out_i <- point_i
    out_i$n_success <- n_success_i
    out_i$n_total <- n_total_i
    out_i$success_rate <- n_success_i / n_total_i
    out_i
  })

  do.call(rbind, rows)
}

#' Propose adaptive isoband points from surrogate predictions.
#'
#' @param design_results_df Existing simulated design results.
#' @param ranges Triad ranges.
#' @param p0 Target replicate-level success probability.
#' @param n_add Number of points to propose.
#' @param seed Optional integer seed.
#' @param n_candidates Candidate pool size.
#'
#' @return List containing `next_points_df`, `pred_candidates_df`, and `gam_fit`.
propose_isoband_adaptive_points <- function(design_results_df, ranges, p0, n_add = 200,
                                            seed = NULL, n_candidates = 5000) {
  required_cols <- c("tau_AE", "tau_ARE", "cutoff", "n_success", "n_total")
  if (!is.data.frame(design_results_df) || !all(required_cols %in% names(design_results_df))) {
    stop(
      sprintf("design_results_df must contain columns: %s", paste(required_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!is.numeric(n_add) || length(n_add) != 1L || !is.finite(n_add) || n_add < 1 || n_add %% 1 != 0) {
    stop("n_add must be a positive integer.", call. = FALSE)
  }

  gam_fit <- fit_isoband_gam(design_results_df = design_results_df)
  if (!isTRUE(gam_fit$converged)) {
    stop("fit_isoband_gam did not converge for adaptive proposal step.", call. = FALSE)
  }
  candidates_df <- generate_isoband_candidate_pool(
    ranges = ranges,
    n_candidates = n_candidates,
    seed = seed
  )
  pred_candidates_df <- predict_isoband_gam(gam_fit = gam_fit, newdata = candidates_df, se_fit = TRUE)
  scored_df <- score_candidates_isoband(pred_df = pred_candidates_df, p0 = p0)
  scored_df <- scored_df[order(scored_df$score_primary, scored_df$score_secondary), , drop = FALSE]
  scored_df <- scored_df[!duplicated(scored_df[, c("tau_AE", "tau_ARE", "cutoff")]), , drop = FALSE]
  n_take <- min(as.integer(n_add), nrow(scored_df))
  next_points_df <- scored_df[seq_len(n_take), c("tau_AE", "tau_ARE", "cutoff"), drop = FALSE]
  next_round <- if ("round" %in% names(design_results_df)) max(design_results_df$round, na.rm = TRUE) + 1L else 1L
  next_points_df$round <- as.integer(next_round)
  next_points_df$stage <- "adaptive"

  list(
    next_points_df = next_points_df,
    pred_candidates_df = scored_df,
    gam_fit = gam_fit
  )
}

#' Re-simulate isoband seed points with higher replication.
#'
#' @param band_seed_points_df Data frame of unrefined band points.
#' @param alpha Scalar alpha value.
#' @param K Number of cell types.
#' @param n Total sample size per replicate.
#' @param R_final Number of final refinement replicates.
#' @param proportion_method Proportion-generation method.
#' @param p_max Optional fixed largest true proportion.
#' @param model Sampling model.
#' @param seed Optional integer seed.
#' @param ... Additional arguments passed to `run_replicates()`.
#'
#' @return Data frame of band points with refined simulation metrics attached.
refine_isoband_band_points <- function(band_seed_points_df, alpha, K, n, R_final = 1000,
                                       proportion_method = "beta",
                                       p_max = NULL,
                                       model = "multinomial",
                                       seed = NULL, ...) {
  if (!is.data.frame(band_seed_points_df)) {
    stop("band_seed_points_df must be a data.frame.", call. = FALSE)
  }
  required_cols <- c("tau_AE", "tau_ARE", "cutoff")
  if (!all(required_cols %in% names(band_seed_points_df))) {
    stop(
      sprintf("band_seed_points_df must contain columns: %s", paste(required_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  if (nrow(band_seed_points_df) == 0L) {
    return(band_seed_points_df)
  }

  unique_points <- unique(band_seed_points_df[, required_cols, drop = FALSE])
  unique_points$stage <- "refine"
  refined_df <- simulate_isoband_design_points(
    design_df = unique_points,
    alpha = alpha,
    K = K,
    n = n,
    B_point = R_final,
    proportion_method = proportion_method,
    p_max = p_max,
    model = model,
    seed = seed,
    ...
  )
  refined_df <- refined_df[, c(required_cols, "n_success", "n_total", "success_rate"), drop = FALSE]
  names(refined_df)[4:6] <- c("n_success_refined", "n_total_refined", "success_rate_refined")
  merge(
    band_seed_points_df,
    refined_df,
    by = required_cols,
    all.x = TRUE,
    sort = FALSE
  )
}

#' Run the automated isoband pipeline for one fixed alpha scenario.
#'
#' @param alpha Scalar alpha value for true-proportion generation.
#' @param p0 Target replicate-level success probability in (0, 1).
#' @param ranges Named list of triad ranges (`tau_AE`, `tau_ARE`, `cutoff`).
#' @param eps Iso-band tolerance (`|p_hat - p0| <= eps`).
#' @param seed Optional RNG seed.
#' @param n Total sample size per replicate.
#' @param K Number of cell types (default 10).
#' @param proportion_method Proportion-generation method (default `"beta"`).
#' @param p_max Optional fixed largest true proportion.
#' @param model Sampling model passed to `run_replicates()`.
#' @param n_init Initial number of design points.
#' @param R_init Initial replicates per design point.
#' @param n_rounds_max Maximum adaptive rounds.
#' @param n_add Points added per adaptive round.
#' @param R_final Replicates for final band refinement.
#' @param n_candidates Candidate pool size for proposal/prediction.
#' @param ... Additional arguments passed to simulation calls.
#'
#' @return List with isoband design history, model, and final band point cloud.
run_isoband_pipeline <- function(alpha, p0, ranges, eps, seed = NULL, n,
                                 K = 10,
                                 proportion_method = "beta",
                                 p_max = NULL,
                                 model = "multinomial",
                                 n_init = 200,
                                 R_init = 200,
                                 n_rounds_max = 5,
                                 n_add = 200,
                                 R_final = 1000,
                                 n_candidates = 5000,
                                 ...) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
    stop("alpha must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.numeric(p0) || length(p0) != 1L || !is.finite(p0) || p0 <= 0 || p0 >= 1) {
    stop("p0 must be a finite scalar strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("eps must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n %% 1 != 0) {
    stop("n must be a positive integer.", call. = FALSE)
  }
  is_positive_integer <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x %% 1 == 0
  if (!is_positive_integer(n_init)) stop("n_init must be a positive integer.", call. = FALSE)
  if (!is_positive_integer(R_init)) stop("R_init must be a positive integer.", call. = FALSE)
  if (!is.numeric(n_rounds_max) || length(n_rounds_max) != 1L || !is.finite(n_rounds_max) ||
    n_rounds_max < 0 || n_rounds_max %% 1 != 0) {
    stop("n_rounds_max must be a non-negative integer.", call. = FALSE)
  }
  if (!is_positive_integer(n_add)) stop("n_add must be a positive integer.", call. = FALSE)
  if (!is_positive_integer(R_final)) stop("R_final must be a positive integer.", call. = FALSE)
  if (!is_positive_integer(n_candidates)) stop("n_candidates must be a positive integer.", call. = FALSE)

  ranges <- validate_isoband_ranges(ranges)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))) {
    stop("seed must be NULL or a finite numeric scalar.", call. = FALSE)
  }

  inputs <- list(
    alpha = alpha,
    p0 = p0,
    ranges = ranges,
    eps = eps,
    seed = seed,
    n = as.integer(n),
    K = as.integer(K),
    proportion_method = proportion_method,
    p_max = p_max,
    model = model,
    n_init = as.integer(n_init),
    R_init = as.integer(R_init),
    n_rounds_max = as.integer(n_rounds_max),
    n_add = as.integer(n_add),
    R_final = as.integer(R_final),
    n_candidates = as.integer(n_candidates)
  )
  result <- new_isoband_result_container(inputs = inputs)

  init_design <- generate_isoband_initial_design(
    ranges = ranges,
    n_init = n_init,
    seed = seed
  )
  design_history <- simulate_isoband_design_points(
    design_df = init_design,
    alpha = alpha,
    K = K,
    n = n,
    B_point = R_init,
    proportion_method = proportion_method,
    p_max = p_max,
    model = model,
    seed = if (is.null(seed)) NULL else as.integer(seed) + 1000L,
    ...
  )
  round_summaries <- data.frame(
    round = 0L,
    stage = "init",
    n_points = nrow(init_design),
    n_total_points = nrow(design_history),
    stringsAsFactors = FALSE
  )

  if (n_rounds_max > 0L) {
    for (round_i in seq_len(as.integer(n_rounds_max))) {
      proposal_i <- propose_isoband_adaptive_points(
        design_results_df = design_history,
        ranges = ranges,
        p0 = p0,
        n_add = n_add,
        seed = if (is.null(seed)) NULL else as.integer(seed) + round_i,
        n_candidates = n_candidates
      )
      next_points <- proposal_i$next_points_df
      if (nrow(next_points) == 0L) {
        round_summaries <- rbind(
          round_summaries,
          data.frame(
            round = round_i,
            stage = "adaptive",
            n_points = 0L,
            n_total_points = nrow(design_history),
            stringsAsFactors = FALSE
          )
        )
        next
      }

      sim_i <- simulate_isoband_design_points(
        design_df = next_points,
        alpha = alpha,
        K = K,
        n = n,
        B_point = R_init,
        proportion_method = proportion_method,
        p_max = p_max,
        model = model,
        seed = if (is.null(seed)) NULL else as.integer(seed) + 2000L + round_i * 10000L,
        ...
      )
      design_history <- rbind(design_history, sim_i)
      round_summaries <- rbind(
        round_summaries,
        data.frame(
          round = round_i,
          stage = "adaptive",
          n_points = nrow(sim_i),
          n_total_points = nrow(design_history),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  final_gam_fit <- fit_isoband_gam(design_results_df = design_history)
  if (!isTRUE(final_gam_fit$converged)) {
    stop("fit_isoband_gam did not converge on final design history.", call. = FALSE)
  }
  dense_candidates <- generate_isoband_candidate_pool(
    ranges = ranges,
    n_candidates = n_candidates,
    seed = if (is.null(seed)) NULL else as.integer(seed) + 500000L
  )
  dense_pred <- predict_isoband_gam(gam_fit = final_gam_fit, newdata = dense_candidates, se_fit = TRUE)
  band_points_unrefined <- extract_isoband_points(pred_df = dense_pred, p0 = p0, eps = eps)
  final_band_points <- refine_isoband_band_points(
    band_seed_points_df = band_points_unrefined,
    alpha = alpha,
    K = K,
    n = n,
    R_final = R_final,
    proportion_method = proportion_method,
    p_max = p_max,
    model = model,
    seed = if (is.null(seed)) NULL else as.integer(seed) + 600000L,
    ...
  )

  result$design_history <- design_history
  result$round_summaries <- round_summaries
  result$final_band_points <- final_band_points
  result$final_band_points_unrefined <- band_points_unrefined
  result$final_model <- final_gam_fit
  result
}
