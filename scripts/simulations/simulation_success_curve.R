# simulation_success_curve.R
#
# MVP driver for success-rate curves vs sample size.
#
# Success is defined jointly across all active scalar thresholds in `config$taus`.

simulation_helper_files <- list.files(here::here("scripts", "simulation_layers"))
lapply(simulation_helper_files, function(f) {
  source(here::here("scripts", "simulation_layers", f))
})


simulation_success_curve_defaults <- function() {
  list(
    alpha = c(2, 2.5, 3, 4, 5),
    n_values = seq(500L, 5000L, by = 500L),
    K = 10L,
    B = 500L,
    taus = list(AE = 0.02, ARE = 0.5),
    metrics = c("AE", "ARE"),
    model = "multinomial",
    tie_method = "random",
    proportion_method = "beta",
    seed = 260925L,
    success_rate_target = 0.95
  )
}


validate_simulation_success_curve_config <- function(config) {
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


simulate_success_curve_for_alpha <- function(alpha, n_values, config, seed_offset = 0L) {
  rows <- vector("list", length(n_values))
  for (i in seq_along(n_values)) {
    config_i <- config
    if (!is.null(config$seed)) {
      config_i$seed <- as.integer(config$seed + seed_offset + i - 1L)
    }

    sim_i <- simulate_success_at_n(alpha = alpha, n = n_values[[i]], config = config_i)
    rows[[i]] <- data.frame(
      alpha = alpha,
      n = as.integer(n_values[[i]]),
      success_rate = sim_i$success_rate,
      success_count = sim_i$success_count,
      B = as.integer(config$B),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


run_simulation_success_curve <- function(config = simulation_success_curve_defaults()) {
  config <- validate_simulation_success_curve_config(config)

  alpha_values <- config$alpha
  n_values <- as.integer(config$n_values)
  curves_list <- vector("list", length(alpha_values))

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

  list(
    inputs = config,
    curves = curves
  )
}


plot_success_rate_vs_n <- function(result, target = NULL) {
  if (is.list(result) && "curves" %in% names(result)) {
    df <- result$curves
    if (is.null(target) && "inputs" %in% names(result) && "success_rate_target" %in% names(result$inputs)) {
      target <- result$inputs$success_rate_target
    }
  } else {
    df <- result
  }

  required_cols <- c("alpha", "n", "success_rate")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("result is missing required columns: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = n,
      y = success_rate,
      color = factor(alpha),
      group = factor(alpha)
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Sample size (n)",
      y = "Success rate",
      color = "alpha",
      title = "Success-rate curves vs sample size"
    ) +
    ggplot2::theme_bw()

  if (!is.null(target)) {
    p <- p + ggplot2::geom_hline(yintercept = target, linetype = "dotted")
  }

  p
}


run_simulation_success_curve_example <- function() {
  config <- simulation_success_curve_defaults()
  config$alpha <- c(2, 3)
  config$n_values <- c(500L, 1000L, 1500L)
  config$B <- 100L

  result <- run_simulation_success_curve(config)
  print(head(result$curves, 10))
  plot_success_rate_vs_n(result)
}
