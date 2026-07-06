# simulation_cutoff_threshold_contours.R
#
# Simulation orchestrator: hybrid AE/ARE cutoff threshold-grid evaluation
# with success-probability contour visualisation.
#
# For a grid of (tau_AE, tau_ARE) pairs and a set of candidate cutoff values,
# estimates the probability that every cell type's selected error metric is
# below its corresponding threshold.
#
# Routing rule (applied per cell type per replicate):
#   use AE  when  phat  < cutoff   (strictly below)
#   use ARE when  phat >= cutoff   (at or above)
#
# Success (per replicate): every cell type satisfies its selected threshold
#   selected AE  < tau_AE
#   selected ARE < tau_ARE
# ---------------------------------------------------------------------------

simulation_helper_files <- list.files(here::here("scripts", "simulation_layers"))
lapply(simulation_helper_files, function(f) {
  source(here::here("scripts", "simulation_layers", f))
})


# ---- Default parameters ----------------------------------------------------

#' Default configuration for the cutoff threshold-grid contour simulation.
#'
#' @return Named list of simulation parameters.
simulation_contour_defaults <- function() {
  list(
    alpha          = c(2, 3, 5),
    K              = 10L,
    n              = 1000L,
    B              = 500L,
    cutoffs        = c(0.05, 0.10, 0.15),
    tau_AE_values  = seq(0.02, 0.20, by = 0.02),
    tau_ARE_values = seq(0.10, 1.20, by = 0.10),
    contour_levels = c(0.5, 0.6, 0.7, 0.8, 0.9),
    proportion_method = "beta",
    p_max          = NULL,
    model          = "multinomial",
    seed           = 260926L
  )
}


# ---- Orchestrator ----------------------------------------------------------

#' Run the cutoff threshold-grid contour simulation.
#'
#' Simulates observed cell-type proportions for each alpha/p_max scenario,
#' then evaluates the hybrid AE/ARE success probability at every point of the
#' (tau_AE, tau_ARE) grid for each cutoff value.  One contour plot per
#' alpha/p_max scenario is produced, faceted by cutoff.
#'
#' @param config          Named list of simulation parameters; see
#'   [simulation_contour_defaults()] for the full set of required fields.
#' @param cache           Logical; when `TRUE` (default) the result is saved to
#'   `cache_dir` and returned from cache on subsequent calls with the same
#'   config.
#' @param force_recompute Logical; when `TRUE` re-runs even if a cached result
#'   exists.
#' @param cache_dir       Directory for cached `.rds` files.
#'
#' @return Named list with elements:
#'   \describe{
#'     \item{config}{The configuration list passed in.}
#'     \item{p_table}{Scenario-level true proportions table.}
#'     \item{phat_long}{Tidy replicate-level observed proportions.}
#'     \item{grid_results}{Tidy data.frame with columns alpha, p_max, cutoff,
#'       tau_AE, tau_ARE, p_success, n_sim — one row per (scenario, cutoff,
#'       threshold pair).}
#'     \item{contour_plots}{Named list of ggplot objects, one per alpha/p_max
#'       scenario.}
#'   }
run_simulation_cutoff_threshold_contours <- function(
    config            = simulation_contour_defaults(),
    cache             = TRUE,
    force_recompute   = FALSE,
    cache_dir         = here::here("results", "simresults")) {

  # ---- Input validation ----------------------------------------------------
  required_fields <- c(
    "alpha", "K", "n", "B", "cutoffs",
    "tau_AE_values", "tau_ARE_values", "contour_levels",
    "proportion_method", "seed"
  )
  missing_fields <- setdiff(required_fields, names(config))
  if (length(missing_fields) > 0L) {
    stop(
      sprintf("config is missing required field(s): %s",
              paste(missing_fields, collapse = ", ")),
      call. = FALSE
    )
  }
  stopifnot(
    is.numeric(config$alpha),    length(config$alpha) >= 1L,
    is.numeric(config$cutoffs),  length(config$cutoffs) >= 1L,
    is.numeric(config$tau_AE_values),  length(config$tau_AE_values) >= 1L,
    is.numeric(config$tau_ARE_values), length(config$tau_ARE_values) >= 1L,
    is.numeric(config$contour_levels),
    all(config$contour_levels > 0 & config$contour_levels < 1),
    config$B >= 1L,
    config$n >= 1L
  )

  # ---- Cache check ---------------------------------------------------------
  result_file <- simulation_result_path(
    config = config,
    dir    = cache_dir,
    name   = "cutoff_threshold_contours"
  )

  if (cache && !force_recompute && file.exists(result_file)) {
    return(readRDS(result_file))
  }

  # ---- Simulate observed proportions ---------------------------------------
  sim_out <- run_simulation_experiment(
    alpha             = config$alpha,
    K                 = config$K,
    n                 = config$n,
    B                 = config$B,
    taus              = list(AE = config$tau_AE_values, ARE = config$tau_ARE_values),
    metrics           = c("AE", "ARE"),
    proportion_method = config$proportion_method,
    p_max             = config$p_max,
    model             = config$model,
    seed              = config$seed
  )

  p_table   <- sim_out$p_table
  index_cols <- grep("^index_", names(p_table), value = TRUE)

  # ---- Evaluate threshold grid per scenario --------------------------------
  grid_results_list  <- vector("list", nrow(p_table))
  contour_plots_list <- vector("list", nrow(p_table))
  plot_names         <- character(nrow(p_table))

  for (i in seq_len(nrow(p_table))) {
    alpha_i <- p_table$alpha[[i]]
    p_max_i <- if ("p_max" %in% names(p_table)) p_table$p_max[[i]] else NA_real_
    p_i     <- as.numeric(p_table[i, index_cols, drop = FALSE])

    phat_long_i <- sim_out$phat_long[sim_out$phat_long$alpha == alpha_i, , drop = FALSE]
    same_p_max <- if (is.na(p_max_i)) {
      is.na(phat_long_i$p_max)
    } else {
      phat_long_i$p_max == p_max_i
    }
    phat_subset <- phat_long_i[same_p_max, c("replicate", "index", "phat"), drop = FALSE]
    phat_subset <- phat_subset[
      order(phat_subset$index, phat_subset$replicate), , drop = FALSE
    ]
    phat_mat_i <- matrix(phat_subset$phat, nrow = config$B, ncol = length(p_i))

    grid_i <- evaluate_hybrid_success_grid_multi_cutoff(
      phat_mat       = phat_mat_i,
      p              = p_i,
      cutoffs        = config$cutoffs,
      tau_AE_values  = config$tau_AE_values,
      tau_ARE_values = config$tau_ARE_values
    )
    grid_i$alpha <- alpha_i
    grid_i$p_max <- p_max_i
    grid_results_list[[i]] <- grid_i[, c(
      "alpha", "p_max", "cutoff", "tau_AE", "tau_ARE", "p_success", "n_sim"
    )]

    scenario_label <- if (is.na(p_max_i)) {
      paste0("alpha = ", alpha_i)
    } else {
      paste0("alpha = ", alpha_i, ", p_max = ", p_max_i)
    }
    plot_names[[i]] <- gsub("[^A-Za-z0-9_]", "_", scenario_label)

    contour_plots_list[[i]] <- plot_success_contours(
      grid_df         = grid_i,
      levels          = config$contour_levels,
      facet_by_cutoff = TRUE,
      fill_contours   = FALSE
    ) + ggplot2::labs(title = paste0("Hybrid cutoff contours (", scenario_label, ")"))
  }

  names(contour_plots_list) <- plot_names

  result <- list(
    config         = config,
    p_table        = p_table,
    phat_long      = sim_out$phat_long,
    grid_results   = do.call(rbind, grid_results_list),
    contour_plots  = contour_plots_list
  )

  if (cache) {
    saveRDS(result, result_file)
  }

  result
}
