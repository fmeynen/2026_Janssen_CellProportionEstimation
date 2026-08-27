# Success is defined jointly across all active scalar thresholds in `config$taus`,
# after averaging patient-level errors within each replicate for every cell type.

rm(list = ls())
# Source Functions ------------------------------------------------------------------------------------------------
simulation_helper_files <- list.files(here::here("scripts", "simulation_layers"))
lapply(simulation_helper_files, function(f) {
  source(here::here("scripts", "simulation_layers", f))
})


# Helper functions ------------------------------------------------------------------------------------------------

sim_success_curve_multiple_patients_defaults <- function() {
  list(
    alpha = c(2, 2.5, 3, 4, 5),
    n_values = 10^seq(from = 3, to = 5, by = 0.5), # sample size (# of cells) per patient
    n_patients = c(1L, 2L, 3L, 5L),
    alpha_sigma = 0.1,
    alpha_min = 1 + 1e-8,
    K = 10L, # number of cell types
    B = 500L, # number of simulations to estimate success rate
    taus = list(AE = 0.02, ARE = 0.5), # AE = absolute error, ARE = absolute relative error
    metrics = c("AE", "ARE"),
    model = "multinomial",
    tie_method = "random",
    proportion_method = "beta",
    seed = 260925L,
    success_rate_target = 0.95
  )
}

validate_sim_success_curve_multiple_patients_config <- function(config) {
  required_fields <- c(
    "alpha", "n_values", "n_patients", "alpha_sigma", "alpha_min", "K", "B", "taus", "metrics",
    "model", "tie_method", "proportion_method", "seed"
  )
  missing_fields <- setdiff(required_fields, names(config))
  if (length(missing_fields) > 0L) {
    stop(
      sprintf("config is missing required fields: %s", paste(missing_fields, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!is.numeric(config$alpha) || length(config$alpha) < 1L || any(!is.finite(config$alpha)) || any(config$alpha <= 0)) {
    stop("config$alpha must be a non-empty numeric vector of positive values.", call. = FALSE)
  }
  if (!is.numeric(config$n_values) || length(config$n_values) < 1L ||
      any(!is.finite(config$n_values)) || any(config$n_values < 1L)) {
    stop("config$n_values must be a non-empty numeric vector of values >= 1.", call. = FALSE)
  }
  if (!is.numeric(config$n_patients) || length(config$n_patients) < 1L || any(!is.finite(config$n_patients)) ||
      any(config$n_patients %% 1 != 0) || any(config$n_patients < 1L)) {
    stop("config$n_patients must be a non-empty numeric vector of integers >= 1.", call. = FALSE)
  }
  if (!is.numeric(config$alpha_sigma) || length(config$alpha_sigma) != 1L ||
      !is.finite(config$alpha_sigma) || config$alpha_sigma < 0) {
    stop("config$alpha_sigma must be a single finite numeric value >= 0.", call. = FALSE)
  }
  if (!is.numeric(config$alpha_min) || length(config$alpha_min) != 1L ||
      !is.finite(config$alpha_min) || config$alpha_min <= 1) {
    stop("config$alpha_min must be a single finite numeric value > 1.", call. = FALSE)
  }
  if (!is.numeric(config$K) || length(config$K) != 1L || !is.finite(config$K) || config$K %% 1 != 0 || config$K < 1L) {
    stop("config$K must be a single integer >= 1.", call. = FALSE)
  }
  if (!is.numeric(config$B) || length(config$B) != 1L || !is.finite(config$B) || config$B %% 1 != 0 || config$B < 1L) {
    stop("config$B must be a single integer >= 1.", call. = FALSE)
  }
  if (!is.character(config$metrics) || length(config$metrics) < 1L || any(is.na(config$metrics))) {
    stop("config$metrics must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.list(config$taus) || is.null(names(config$taus))) {
    stop("config$taus must be a named list.", call. = FALSE)
  }

  active_thresholds <- intersect(config$metrics, names(config$taus))
  if (length(active_thresholds) != length(config$metrics)) {
    stop("config$taus must provide one scalar threshold for every metric in config$metrics.", call. = FALSE)
  }
  for (metric in active_thresholds) {
    tau <- config$taus[[metric]]
    if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau)) {
      stop(sprintf("config$taus[['%s']] must be a finite numeric scalar.", metric), call. = FALSE)
    }
  }

  if (!is.null(config$seed)) {
    if (!is.numeric(config$seed) || length(config$seed) != 1L || !is.finite(config$seed) || config$seed %% 1 != 0) {
      stop("config$seed must be NULL or a single finite integer value.", call. = FALSE)
    }
  }

  config
}

run_simulation_success_curve_multiple_patients <- function(config = sim_success_curve_multiple_patients_defaults(),
                                                           overwrite = FALSE,
                                                           dir = "results/simresults",
                                                           name = "alpha_curve_multiple_patients") {
  config <- validate_sim_success_curve_multiple_patients_config(config)

  path <- simulation_result_path(config, dir, name)
  if (!overwrite && file.exists(path)) {
    return(readRDS(path))
  }

  alpha_values <- config$alpha
  n_values <- config$n_values
  n_patients <- as.integer(config$n_patients)
  curves_list <- vector("list", length(alpha_values))

  for (a in seq_along(alpha_values)) {
    seed_offset <- (a - 1L) * length(n_values) * length(n_patients)
    curves_list[[a]] <- simulate_success_curve_for_alpha_over_n_and_patients(
      alpha = alpha_values[[a]],
      n_values = n_values,
      n_patients = n_patients,
      config = config,
      seed_offset = seed_offset
    )
  }
  curves <- do.call(rbind, curves_list)
  rownames(curves) <- NULL

  if ("success_rate_target" %in% names(config)) {
    curves$success_rate_target <- config$success_rate_target
  }

  res <- list(
    inputs = config,
    curves = curves
  )
  saveRDS(res, path)
  res
}


# Perform simulation ----------------------------------------------------------------------------------------------
## Configuration
cfg <- sim_success_curve_multiple_patients_defaults()
cfg$n_values <- 10^seq(from = 3, to = 4, by = 0.05)
cfg$B <- 1000L
## Simulation
result <- run_simulation_success_curve_multiple_patients(cfg)

plot_success_rate_vs_n(result, target = 0.95, smooth = FALSE, color_by = "n_patients", facet_by_alpha = TRUE)
