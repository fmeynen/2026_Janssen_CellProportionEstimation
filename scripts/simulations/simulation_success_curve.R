# Success is defined *jointly* across all active scalar thresholds in `config$taus`.

rm(list = ls())
# Source Functions ------------------------------------------------------------------------------------------------
simulation_helper_files <- list.files(here::here("scripts", "simulation_layers"))
lapply(simulation_helper_files, function(f) {
  source(here::here("scripts", "simulation_layers", f))
})


# Helper functions ------------------------------------------------------------------------------------------------

sim_success_curve_defaults <- function() {
  list(
    alpha = c(2, 2.5, 3, 4, 5),
    n_values = 10^seq(from = 3, to = 7, by = 0.5), #sample size (# of cells)
    K = 10L, #number of cell types
    B = 500L, #number of simulations to estimate success rate
    taus = list(AE = 0.02, ARE = 0.5),
    metrics = c("AE", "ARE"), #absolute error, absolute relative error
    model = "multinomial",
    tie_method = "random",
    proportion_method = "beta",
    seed = 260925L,
    success_rate_target = 0.95
  )
}

validate_sim_success_curve_config <- function(config) {
  required_fields <- c(
    "alpha", "n_values", "K", "B", "taus", "metrics",
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
  if (!is.numeric(config$n_values) || length(config$n_values) < 1L || any(!is.finite(config$n_values)) || any(config$n_values < 1L)) {
    stop("config$n_values must be a non-empty numeric vector of values >= 1.", call. = FALSE)
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
  if (length(active_thresholds) < 1L) {
    stop("config$taus must provide at least one threshold for a metric in config$metrics.", call. = FALSE)
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

run_simulation_success_curve <- function(config = sim_success_curve_defaults(),
                                         overwrite = FALSE, dir = "results/simresults", name = "alpha_curve") {
  config <- validate_sim_success_curve_config(config)

  path <- simulation_result_path(config, dir, name)
  if(file.exists(path)){
    return(readRDS(path))
  }
  
  alpha_values <- config$alpha
  n_values     <- as.integer(config$n_values)
  curves_list  <- vector("list", length(alpha_values))

  for (a in seq_along(alpha_values)) {
    seed_offset <- (a - 1L) * length(n_values)
    curves_list[[a]] <- simulate_success_curve_for_alpha(
      alpha = alpha_values[[a]],
      n_values = n_values,
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
cfg          <- sim_success_curve_defaults()
cfg$alpha    <- seq(from = 2, to = 3.5, by = 0.5)
cfg$n_values <- 10^seq(from = 3, to = 5, by = 0.05)
cfg$B        <- 10000
## Simulation
result <- run_simulation_success_curve(cfg)

plot_success_rate_vs_n(result, target = 0.95, smooth = FALSE)

