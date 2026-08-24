
# Orchestration ---------------------------------------------------------------------------------------------------
# Orchestration layer: high-level experiment runners that coordinate all other layers.
#
# Depends on: simulation.R, calculation.R, extraction.R



# Original Simulation ---------------------------------------------------------------------------------------------

#' Run the original full simulation experiment end-to-end.
#'
#' @param alpha      Numeric vector; one or more positive shape values used by the selected proportion-generation method
#'      (default method is Beta-based).
#' @param K          Number of cell types (default 10).
#' @param n          Total sample size per replicate.
#' @param B          Number of replicates.
#' @param taus       Numeric vector of thresholds (same for all metrics) or a named list with one numeric vector per
#'      metric (e.g. `list(AE = c(...), ARE = c(...))`).
#' @param metrics    Error metrics; any subset of `c("AE", "ARE", "TSE", "LAE")`.
#' @param proportion_method Proportion-generation method (`"beta"` or `"fixed_max_beta"`).
#'    The fixed-max Beta method places `p_max` at the highest index and warns then fails for impossible combinations.
#' @param p_max      Fixed largest true proportion(s) used by `proportion_method = "fixed_max_beta"`.
#'      Can be a numeric vector.
#'      When multiple values are provided, all alpha × p_max combinations are attempted; impossible fixed-max
#'      combinations are warned and skipped.
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


# Hybrid Cutoff Simulation ----------------------------------------------------------------------------------------

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
