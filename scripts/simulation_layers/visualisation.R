
# Visualization Layer ---------------------------------------------------------------------------------------------

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
  stopifnot(
    "result must be a list" = is.list(result),
    "result must contain inputs" = "inputs" %in% names(result),
    "result must contain p_table" = "p_table" %in% names(result)
  )
  if (!("proportion_method" %in% names(result$inputs))) {
    stop("result$inputs must contain proportion_method.", call. = FALSE)
  }

  p_table <- result$p_table
  if (!is.data.frame(p_table) || nrow(p_table) == 0L) {
    stop("result$p_table must be a non-empty data.frame.", call. = FALSE)
  }
  index_cols <- grep("^index_", names(p_table), value = TRUE)
  if (length(index_cols) == 0L) {
    stop("result$p_table must contain index_1 ... index_K columns.", call. = FALSE)
  }
  if (!("alpha" %in% names(p_table))) {
    stop("result$p_table must contain an alpha column.", call. = FALSE)
  }

  method <- result$inputs$proportion_method
  if (!method %in% c("beta", "fixed_max_beta")) {
    stop("Unsupported proportion_method in result$inputs.", call. = FALSE)
  }

  K <- length(index_cols)
  n_rows <- nrow(p_table)
  curve_rows <- vector("list", n_rows)
  point_rows <- vector("list", n_rows)
  curve_upper_bound <- NA_real_
  if (identical(method, "fixed_max_beta")) {
    curve_upper_bound <- default_beta_grid(K - 1L)[K - 1L]
  }

  for (i in seq_len(n_rows)) {
    alpha_i <- as.numeric(p_table$alpha[[i]])
    p_i <- as.numeric(p_table[i, index_cols, drop = FALSE])
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
#' @param target  Success-rate reference line drawn as a horizontal dotted line (default 0.95).
#'
#' @return A ggplot object.
plot_success_rate_curve <- function(result, metric = NULL, alphas = NULL,
                                    p_maxs = NULL, target = 0.95) {
  stopifnot(is.list(result), "curves" %in% names(result))
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
    ggplot2::geom_hline(yintercept = target, linetype = "dotted") +
    ggplot2::labs(
      x     = "Threshold (tau)",
      y     = "Success rate",
      color = "alpha",
      title = "Success-rate curve(s)"
    ) +
    ggplot2::theme_bw()
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
  stopifnot(
    "result must be a list" = is.list(result),
    "result must contain replicate_summaries" = "replicate_summaries" %in% names(result)
  )
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

plot_success_rate_vs_n <- function(result, target = NULL, smooth = FALSE) {
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


# Hybrid Cutoff ---------------------------------------------------------------------------------------------------


#' Plot a heatmap of best hybrid cutoffs across AE/ARE threshold pairs.
#'
#' @param phat_mat       B x K numeric matrix of observed proportions (B
#'   replicates by K cell types).
#' @param p              True proportion vector of length K.
#' @param cutoffs        Numeric vector of candidate cutoffs.
#' @param tau_AE_values  Numeric vector of AE thresholds for heatmap rows.
#' @param tau_ARE_values Numeric vector of ARE thresholds for heatmap columns.
#' @param maximize       Which success rate to maximize (`"cell"` or
#'   `"replicate"`).
#' @param tie_break      Tie-breaker among maximizing cutoffs:
#'   `"smallest"`, `"largest"`, or `"median"`.
#' @param label_digits   Number of digits displayed in cell labels.
#'
#' @return A ggplot heatmap object with fill and text labels equal to the selected best cutoff for each
#'  (tau_AE, tau_ARE) pair.
plot_hybrid_best_cutoff_heatmap <- function(phat_mat, p, cutoffs,
                                            tau_AE_values, tau_ARE_values,
                                            maximize = c("cell", "replicate"),
                                            tie_break = c("smallest", "largest", "median"),
                                            label_digits = 3L) {
  maximize <- match.arg(maximize)
  tie_break <- match.arg(tie_break)

  stopifnot(
    is.matrix(phat_mat),
    is.numeric(p),
    length(p) == ncol(phat_mat),
    is.numeric(cutoffs),
    length(cutoffs) >= 1L,
    all(is.finite(cutoffs)),
    is.numeric(tau_AE_values),
    length(tau_AE_values) >= 1L,
    all(is.finite(tau_AE_values)),
    is.numeric(tau_ARE_values),
    length(tau_ARE_values) >= 1L,
    all(is.finite(tau_ARE_values)),
    is.numeric(label_digits),
    length(label_digits) == 1L,
    is.finite(label_digits),
    label_digits %% 1 == 0,
    label_digits >= 0
  )

  threshold_grid <- expand.grid(
    tau_AE = tau_AE_values,
    tau_ARE = tau_ARE_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  rows <- lapply(seq_len(nrow(threshold_grid)), function(i) {
    tau_AE_i <- threshold_grid$tau_AE[[i]]
    tau_ARE_i <- threshold_grid$tau_ARE[[i]]
    curve_i <- sweep_hybrid_cutoffs_cell_level(
      phat_mat = phat_mat,
      p = p,
      cutoffs = cutoffs,
      tau_AE = tau_AE_i,
      tau_ARE = tau_ARE_i
    )
    best_i <- find_best_hybrid_cutoff(
      curve_df = curve_i,
      maximize = maximize,
      tie_break = tie_break
    )
    data.frame(
      tau_AE = tau_AE_i,
      tau_ARE = tau_ARE_i,
      best_cutoff = best_i$best_cutoff,
      stringsAsFactors = FALSE
    )
  })

  heatmap_df <- do.call(rbind, rows)
  heatmap_df$tau_AE <- factor(heatmap_df$tau_AE, levels = sort(tau_AE_values, decreasing = TRUE))
  heatmap_df$tau_ARE <- factor(heatmap_df$tau_ARE, levels = sort(tau_ARE_values))
  heatmap_df$label <- formatC(heatmap_df$best_cutoff, digits = label_digits, format = "f")

  ggplot2::ggplot(
    heatmap_df,
    ggplot2::aes(x = tau_ARE, y = tau_AE, fill = best_cutoff)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08306B") +
    ggplot2::labs(
      x = "ARE threshold (tau_ARE)",
      y = "AE threshold (tau_AE)",
      fill = "Best cutoff",
      title = "Best hybrid cutoff heatmap"
    ) +
    ggplot2::theme_bw()
}


