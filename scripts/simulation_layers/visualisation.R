
# Visualisation Layer ---------------------------------------------------------------------------------------------

# Visualisation layer: plotting functions for simulation results.
#
# Depends on ggplot2


#' Plot true-proportion points with their underlying Beta-shaped curves.
#'
#' @param result Output list from `run_simulation_experiment()`.
#'
#' @return A ggplot object. For Beta proportions, facets are by alpha.
#'    For fixed-max Beta proportions,facets are by p_max (rows) and alpha (columns).
plot_proportions_curve <- function(result) {
  validate_result_fields(result, c("inputs", "p_table"))
  if (!("proportion_method" %in% names(result$inputs))) {
    stop("result$inputs must contain proportion_method.", call. = FALSE)
  }

  p_table <- result$p_table
  if (!is.data.frame(p_table) || nrow(p_table) == 0L) {
    stop("result$p_table must be a non-empty data.frame.", call. = FALSE)
  }
  cell_type_cols <- grep("^cell_type_", names(p_table), value = TRUE)
  if (length(cell_type_cols) == 0L) {
    stop("result$p_table must contain cell_type_1 ... cell_type_K columns.", call. = FALSE)
  }
  if (!("alpha" %in% names(p_table))) {
    stop("result$p_table must contain an alpha column.", call. = FALSE)
  }

  method <- result$inputs$proportion_method
  if (!method %in% c("beta", "fixed_max_beta")) {
    stop("Unsupported proportion_method in result$inputs.", call. = FALSE)
  }

  K <- length(cell_type_cols)
  n_rows <- nrow(p_table)
  curve_rows <- vector("list", n_rows)
  point_rows <- vector("list", n_rows)
  curve_upper_bound <- NA_real_
  if (identical(method, "fixed_max_beta")) {
    curve_upper_bound <- default_beta_grid(K - 1L)[K - 1L]
  }

  for (i in seq_len(n_rows)) {
    alpha_i <- as.numeric(p_table$alpha[[i]])
    p_i <- as.numeric(p_table[i, cell_type_cols, drop = FALSE])
    p_max_i <- if ("p_max" %in% names(p_table)) as.numeric(p_table$p_max[[i]]) else NA_real_

    if (identical(method, "beta")) {
      grid_i <- default_beta_grid(K)
      w <- dbeta(grid_i, shape1 = alpha_i, shape2 = 1)
      s <- seq(0, 1, length.out = 1000)
      f <- dbeta(s, shape1 = alpha_i, shape2 = 1) / sum(w)
      x_points <- grid_i
    } else {
      fixed_max_point_x <- 1
      if (!is.finite(p_max_i)) {
        stop("result$p_table must contain finite p_max values for fixed_max_beta.", call. = FALSE)
      }
      grid_i <- default_beta_grid(K - 1L)
      w <- dbeta(grid_i, shape1 = alpha_i, shape2 = 1)
      s <- seq(0, curve_upper_bound, length.out = 1000)
      f <- (1 - p_max_i) * dbeta(s, shape1 = alpha_i, shape2 = 1) / sum(w)
      x_points <- c(grid_i, fixed_max_point_x)
    }

    curve_rows[[i]] <- data.frame(
      alpha = alpha_i,
      p_max = p_max_i,
      s = s,
      f = f,
      stringsAsFactors = FALSE
    )
    point_rows[[i]] <- data.frame(
      alpha = alpha_i,
      p_max = p_max_i,
      x = x_points,
      p = p_i,
      stringsAsFactors = FALSE
    )
  }

  curve_df <- do.call(rbind, curve_rows)
  point_df <- do.call(rbind, point_rows)
  alpha_levels <- unique(p_table$alpha)
  curve_df$alpha <- factor(curve_df$alpha, levels = alpha_levels)
  point_df$alpha <- factor(point_df$alpha, levels = alpha_levels)
  if (identical(method, "fixed_max_beta") && "p_max" %in% names(p_table)) {
    p_max_levels <- unique(p_table$p_max)
    curve_df$p_max <- factor(curve_df$p_max, levels = p_max_levels)
    point_df$p_max <- factor(point_df$p_max, levels = p_max_levels)
  }

  proportions_plot <- ggplot2::ggplot(curve_df, ggplot2::aes(x = s, y = f)) +
    ggplot2::geom_line(linewidth = 0.7, color = "#2C3E50") +
    ggplot2::geom_point(
      data = point_df,
      ggplot2::aes(x = x, y = p),
      inherit.aes = FALSE,
      size = 1.8,
      color = "#D62728"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(
      x = "x",
      y = "p",
      title = "True proportions and Beta-shaped curves"
    )

  if (identical(method, "fixed_max_beta")) {
    proportions_plot <- proportions_plot + ggplot2::facet_grid(rows = ggplot2::vars(p_max), cols = ggplot2::vars(alpha))
  } else {
    proportions_plot <- proportions_plot + ggplot2::facet_grid(cols = ggplot2::vars(alpha))
  }
  proportions_plot
}


#' Plot success-rate curves from a simulation result object.
#'
#' @param result  Output list from `run_simulation_experiment()`.
#' @param metric  Character scalar. If NULL, plot all metrics (faceted).
#' @param alphas  Optional numeric vector; subset of alpha values to plot.
#' @param p_maxs  Optional numeric vector; subset of p_max values to plot.
#' @param target  Success-rate reference line drawn as a horizontal dotted line (default 0.95);
#'   `NULL` hides the line.
#'
#' @return A ggplot object.
plot_success_rate_curve <- function(result, metric = NULL, alphas = NULL,
                                    p_maxs = NULL, target = 0.95) {
  validate_result_fields(result, "curves")
  df <- result$curves
  if (!is.null(metric)) {
    df <- df[df$metric %in% metric, , drop = FALSE]
  }
  if (!is.null(alphas)) {
    df <- df[df$alpha %in% alphas, , drop = FALSE]
  }
  if (!is.null(p_maxs)) {
    if (!("p_max" %in% names(df))) {
      stop("No p_max metadata found in result$curves.", call. = FALSE)
    }
    df <- df[df$p_max %in% p_maxs, , drop = FALSE]
  }
  if (nrow(df) == 0L) {
    stop("No rows match the requested metric/alpha/p_max filter(s).", call. = FALSE)
  }

  has_p_max <- "p_max" %in% names(df) && any(!is.na(df$p_max))
  aes_mapping <- ggplot2::aes(x = tau, y = success_rate, color = factor(alpha))
  if (has_p_max) {
    aes_mapping <- ggplot2::aes(
      x = tau, y = success_rate, color = factor(alpha),
      linetype = factor(p_max), group = interaction(alpha, p_max)
    )
  }

  p <- ggplot2::ggplot(df, aes_mapping) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x     = "Threshold (tau)",
      y     = "Success rate",
      color = "alpha",
      title = "Success-rate curve(s)"
    ) +
    ggplot2::theme_bw()
  if (!is.null(target)) {
    p <- p + ggplot2::geom_hline(yintercept = target, linetype = "dotted")
  }
  if (has_p_max) {
    p <- p + ggplot2::labs(linetype = "p_max")
  }
  if (is.null(metric) || length(unique(df$metric)) > 1L) {
    p <- p + ggplot2::facet_wrap(~metric, scales = "free_x")
  }
  p
}

#' Plot a histogram of which cell-type proportion drives the maximum error.
#'
#' @param result  Output list from `run_simulation_experiment()`.
#' @param metric  Character scalar; one of `"AE"`, `"ARE"`, `"TSE"`, or `"LAE"`.
#' @param alphas  Optional numeric vector; subset of alpha values to plot.
#' @param p_maxs  Optional numeric vector; subset of p_max values to plot.
#'
#' @return A ggplot object.
plot_argmax_histogram <- function(result, metric, alphas = NULL, p_maxs = NULL) {
  validate_result_fields(result, "replicate_summaries")
  df <- result$replicate_summaries
  metric_is_scalar_character <- is.character(metric) && length(metric) == 1L && !is.na(metric)
  metric_is_supported <- metric_is_scalar_character && metric %in% c("AE", "ARE", "TSE", "LAE")
  if (!metric_is_supported) {
    stop("metric must be a single value: 'AE', 'ARE', 'TSE', or 'LAE'.", call. = FALSE)
  }
  if (!any(df$metric %in% metric)) {
    stop("No rows match the requested metric value.", call. = FALSE)
  }
  df <- df[df$metric %in% metric, , drop = FALSE]
  if (!is.null(alphas)) {
    if (!any(df$alpha %in% alphas)) {
      stop("No rows match the requested alpha value(s).", call. = FALSE)
    }
    df <- df[df$alpha %in% alphas, , drop = FALSE]
  }
  if (!is.null(p_maxs)) {
    if (!("p_max" %in% names(df))) {
      stop("No p_max metadata found in result$replicate_summaries.", call. = FALSE)
    }
    if (!any(df$p_max %in% p_maxs)) {
      stop("No rows match the requested p_max value(s).", call. = FALSE)
    }
    df <- df[df$p_max %in% p_maxs, , drop = FALSE]
  }
  has_p_max <- "p_max" %in% names(df) && any(!is.na(df$p_max))

  argmax_plot <- ggplot2::ggplot(df, ggplot2::aes(x = argmax_index)) +
    ggplot2::geom_bar() +
    ggplot2::labs(
      x     = "Cell-type index (argmax)",
      y     = "Count",
      title = "Distribution of maximum-error indices"
    ) +
    ggplot2::theme_bw()

  if (has_p_max) {
    argmax_plot <- argmax_plot + ggplot2::facet_grid(rows = ggplot2::vars(p_max), cols = ggplot2::vars(alpha))
  } else {
    argmax_plot <- argmax_plot + ggplot2::facet_grid(cols = ggplot2::vars(alpha))
  }
  argmax_plot
}

#' Plot success rate against sample size.
#'
#' @param result  Output list with a `curves` data.frame, or such a data.frame directly (columns
#'   `alpha`, `n`, `success_rate`).
#' @param target  Success-rate reference line drawn as a horizontal dotted line (default 0.95). If
#'   not supplied and `result$inputs$success_rate_target` exists, that value is used. Pass `NULL`
#'   to hide the line.
#' @param smooth  Logical; draw a smoothed curve instead of a line.
#'
#' @return A ggplot object.
plot_success_rate_vs_n <- function(result, target = 0.95, smooth = FALSE) {
  if (is.list(result) && "curves" %in% names(result)) {
    df <- result$curves
    if (missing(target) && "inputs" %in% names(result) && "success_rate_target" %in% names(result$inputs)) {
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
    (if (smooth) ggplot2::geom_smooth() else ggplot2::geom_line()) +
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


# Dirichlet-multinomial errorchoice ---------------------------------------------------------------------------------


#' Plot Dirichlet-multinomial success rate vs threshold (tau).
#'
#' One panel per alpha; within a panel, one curve per n_people (colour).
#'
#' @param curves_tau data.frame with columns `alpha`, `n_people`, `metric`,
#'   `tau`, `success_rate` (as produced by the errorchoice tau-sweep, with
#'   `n_per_person` held fixed).
#' @param metric     Character scalar; the metric to plot (e.g. `"AE"` or
#'   `"ARE"`).
#' @param subtitle   Optional character scalar; passed through to `labs()`
#'   (e.g. to state the fixed `n_per_person`).
#' @param target     Numeric scalar (default `0.95`); if non-`NULL`, draws a
#'   dotted horizontal reference line at this success rate.
#'
#' @return A ggplot object.
plot_success_vs_tau <- function(curves_tau, metric, subtitle = NULL, target = 0.95) {
  if (!is.data.frame(curves_tau)) {
    stop("curves_tau must be a data.frame.", call. = FALSE)
  }
  required_cols <- c("alpha", "n_people", "metric", "tau", "success_rate")
  missing_cols <- setdiff(required_cols, names(curves_tau))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("curves_tau is missing required columns: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!is.character(metric) || length(metric) != 1L || is.na(metric)) {
    stop("metric must be a single character string.", call. = FALSE)
  }

  df <- curves_tau[curves_tau$metric == metric, , drop = FALSE]
  if (nrow(df) == 0L) {
    stop(sprintf("No rows in curves_tau match metric '%s'.", metric), call. = FALSE)
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = tau,
      y = success_rate,
      color = factor(n_people),
      group = factor(n_people)
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 0.6) +
    ggplot2::facet_wrap(~alpha, labeller = ggplot2::label_both) +
    ggplot2::scale_color_viridis_d(end = 0.85) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x        = "Threshold (tau)",
      y        = "Success rate",
      color    = "n_people",
      title    = sprintf("%s: success rate vs threshold", metric),
      subtitle = subtitle
    ) +
    ggplot2::theme_bw()

  if (!is.null(target)) {
    p <- p + ggplot2::geom_hline(yintercept = target, linetype = "dotted")
  }

  p
}


#' Plot Dirichlet-multinomial success rate vs cells per person.
#'
#' One panel per alpha; within a panel, one curve per n_people (colour).
#'
#' @param curves_n  data.frame with columns `alpha`, `n_people`,
#'   `n_per_person`, `metric`, `tau`, `success_rate` (as produced by the
#'   errorchoice n_per_person-sweep, with `tau` held fixed per metric).
#' @param metric    Character scalar; the metric to plot (e.g. `"AE"` or
#'   `"ARE"`).
#' @param subtitle  Optional character scalar; passed through to `labs()`
#'   (e.g. to state the fixed `tau`).
#' @param target    Numeric scalar (default `0.95`); if non-`NULL`, draws a
#'   dotted horizontal reference line at this success rate.
#'
#' @return A ggplot object.
plot_success_vs_n <- function(curves_n, metric, subtitle = NULL, target = 0.95) {
  if (!is.data.frame(curves_n)) {
    stop("curves_n must be a data.frame.", call. = FALSE)
  }
  required_cols <- c("alpha", "n_people", "n_per_person", "metric", "tau", "success_rate")
  missing_cols <- setdiff(required_cols, names(curves_n))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("curves_n is missing required columns: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!is.character(metric) || length(metric) != 1L || is.na(metric)) {
    stop("metric must be a single character string.", call. = FALSE)
  }

  df <- curves_n[curves_n$metric == metric, , drop = FALSE]
  if (nrow(df) == 0L) {
    stop(sprintf("No rows in curves_n match metric '%s'.", metric), call. = FALSE)
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = n_per_person,
      y = success_rate,
      color = factor(n_people),
      group = factor(n_people)
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~alpha, labeller = ggplot2::label_both) +
    ggplot2::scale_color_viridis_d(end = 0.85) +
    ggplot2::scale_x_log10() +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x        = "Cells per person (n_per_person, log scale)",
      y        = "Success rate",
      color    = "n_people",
      title    = sprintf("%s: success rate vs cells per person", metric),
      subtitle = subtitle
    ) +
    ggplot2::theme_bw()

  if (!is.null(target)) {
    p <- p + ggplot2::geom_hline(yintercept = target, linetype = "dotted")
  }

  p
}
