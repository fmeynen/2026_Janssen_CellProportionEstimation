# Tests for plotting helpers (p_max filters, facet layout). Small B to stay fast.

skip_if_not_installed("ggplot2")

plot_res_beta <- function() {
  run_simulation_experiment(
    alpha = 2, K = 10L, n = 100L, B = 10L,
    taus = c(0.05, 0.10, 0.20), metrics = c("AE", "ARE"), seed = 7L
  )
}

plot_res_fixed <- function() {
  run_simulation_experiment(
    alpha = 2, K = 10L, n = 100L, B = 5L, taus = c(0.05, 0.10),
    proportion_method = "fixed_max_beta", p_max = 0.4, seed = 9L
  )
}

plot_res_multi_pmax <- function() {
  suppressWarnings(run_simulation_experiment(
    alpha = c(2, 3), K = 2L, n = 100L, B = 5L, taus = c(0.05, 0.10),
    proportion_method = "fixed_max_beta", p_max = c(0.4, 0.8), seed = 11L
  ))
}

test_that("plot_proportions_curve returns a ggplot for beta results", {
  expect_s3_class(plot_proportions_curve(plot_res_beta()), "ggplot")
})

test_that("plot_proportions_curve returns a ggplot for fixed_max_beta results", {
  expect_s3_class(plot_proportions_curve(plot_res_fixed()), "ggplot")
})

test_that("fixed_max_beta proportions plot puts points on the beta grid with unpinned points on the curve", {
  res <- plot_res_fixed()
  b <- ggplot2::ggplot_build(plot_proportions_curve(res))
  curve <- b$data[[1]]
  pts <- b$data[[2]]
  expect_equal(pts$x, default_beta_grid(10L))
  p_max <- res$p_table$p_max[[1]]
  unpinned <- abs(pts$y - p_max) > 1e-9
  expect_true(any(!unpinned))
  expect_true(any(unpinned))
  on_curve <- stats::approx(curve$x, curve$y, xout = pts$x[unpinned])$y
  expect_equal(on_curve, pts$y[unpinned], tolerance = 1e-3)
})

test_that("fixed_max_beta proportions plot draws a dashed hline at the bound", {
  res <- plot_res_fixed()
  b <- ggplot2::ggplot_build(plot_proportions_curve(res))
  hline <- Filter(function(d) "yintercept" %in% names(d), b$data)
  expect_length(hline, 1L)
  expect_equal(unique(hline[[1]]$yintercept), res$p_table$p_max[[1]])
})

test_that("fixed_min_beta proportions plot works from a minimal dm_errorchoice-style result", {
  props <- generate_proportions(alpha = 3, K = 10L, method = "fixed_min_beta", p_min = 0.01)
  props5 <- generate_proportions(alpha = 5, K = 10L, method = "fixed_min_beta", p_min = 0.01)
  p_table <- data.frame(alpha = c(3, 5), rbind(props, props5))
  names(p_table) <- c("alpha", paste0("cell_type_", 1:10))
  res <- list(inputs = list(proportion_method = "fixed_min_beta", p_min = 0.01), p_table = p_table)
  p <- plot_proportions_curve(res)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  pts <- b$data[[2]]
  expect_equal(unique(pts$x), default_beta_grid(10L))
  hline <- Filter(function(d) "yintercept" %in% names(d), b$data)
  expect_equal(unique(hline[[1]]$yintercept), 0.01)
  expect_error(
    plot_proportions_curve(list(inputs = list(proportion_method = "fixed_min_beta"), p_table = p_table)),
    "p_min"
  )
})

test_that("plot_success_rate_curve supports p_max filtering", {
  p <- plot_success_rate_curve(plot_res_multi_pmax(), metric = "AE", alphas = c(2, 3), p_maxs = 0.8)
  expect_s3_class(p, "ggplot")
})

test_that("plot_argmax_histogram supports p_max filtering", {
  p <- plot_argmax_histogram(plot_res_multi_pmax(), metric = "AE", alphas = c(2, 3), p_maxs = 0.8)
  expect_s3_class(p, "ggplot")
})

test_that("plot_argmax_histogram rejects an unsupported metric", {
  expect_error(
    plot_argmax_histogram(plot_res_multi_pmax(), metric = "bogus"),
    "metric must be a single value"
  )
})

test_that("plot_argmax_histogram errors informatively for an unmatched p_max filter", {
  expect_error(
    plot_argmax_histogram(plot_res_multi_pmax(), metric = "AE", p_maxs = 0.4),
    "No rows match the requested p_max"
  )
})

test_that("plot_argmax_histogram facets with alpha in columns and p_max in rows", {
  res <- run_simulation_experiment(
    alpha = c(2, 3, 4), K = 10L, n = 50L, B = 5L, taus = c(0.05, 0.10),
    proportion_method = "fixed_max_beta", p_max = c(0.6, 0.8), seed = 17L
  )
  layout <- ggplot2::ggplot_build(plot_argmax_histogram(res, metric = "AE"))$layout$layout
  expect_equal(length(unique(layout$COL)), length(unique(res$replicate_summaries$alpha)))
  expect_equal(length(unique(layout$ROW)), length(unique(res$replicate_summaries$p_max)))
})

# Success-plot target lines: consistent default (0.95) and dotted linetype ---------------------------

hline_layers <- function(p) {
  Filter(function(l) inherits(l$geom, "GeomHline"), p$layers)
}

expect_dotted_target <- function(p, yintercept) {
  layers <- hline_layers(p)
  expect_length(layers, 1L)
  expect_equal(layers[[1]]$data$yintercept, yintercept)
  expect_equal(layers[[1]]$aes_params$linetype, "dotted")
}

curves_vs_n_df <- function() {
  data.frame(alpha = 2, n = c(10, 20, 30), success_rate = c(0.5, 0.8, 0.96))
}

curves_tau_df <- function() {
  data.frame(
    alpha = 2, n_people = 1L, metric = "AE",
    tau = c(0.05, 0.1, 0.2), success_rate = c(0.5, 0.8, 0.96)
  )
}

curves_n_df <- function() {
  data.frame(
    alpha = 2, n_people = 1L, n_per_person = c(10, 100, 1000), metric = "AE",
    tau = 0.1, success_rate = c(0.5, 0.8, 0.96)
  )
}

test_that("plot_success_rate_curve draws a dotted 0.95 target by default and hides it for NULL", {
  res <- plot_res_beta()
  expect_dotted_target(plot_success_rate_curve(res, metric = "AE"), 0.95)
  expect_length(hline_layers(plot_success_rate_curve(res, metric = "AE", target = NULL)), 0L)
})

test_that("plot_success_rate_vs_n draws a dotted 0.95 target by default and hides it for NULL", {
  expect_dotted_target(plot_success_rate_vs_n(curves_vs_n_df()), 0.95)
  expect_length(hline_layers(plot_success_rate_vs_n(curves_vs_n_df(), target = NULL)), 0L)
})

test_that("plot_success_rate_vs_n uses inputs$success_rate_target only when target is not supplied", {
  res <- list(curves = curves_vs_n_df(), inputs = list(success_rate_target = 0.9))
  expect_dotted_target(plot_success_rate_vs_n(res), 0.9)
  expect_dotted_target(plot_success_rate_vs_n(res, target = 0.8), 0.8)
  expect_length(hline_layers(plot_success_rate_vs_n(res, target = NULL)), 0L)
})

test_that("plot_success_vs_tau draws a dotted 0.95 target by default and hides it for NULL", {
  expect_dotted_target(plot_success_vs_tau(curves_tau_df(), "AE"), 0.95)
  expect_length(hline_layers(plot_success_vs_tau(curves_tau_df(), "AE", target = NULL)), 0L)
})

test_that("plot_success_vs_n draws a dotted 0.95 target by default and hides it for NULL", {
  expect_dotted_target(plot_success_vs_n(curves_n_df(), "AE"), 0.95)
  expect_length(hline_layers(plot_success_vs_n(curves_n_df(), "AE", target = NULL)), 0L)
})


sample_size_df <- function() {
  data.frame(
    alpha = rep(c(1, 5), each = 6),
    n_people = rep(rep(c(10, 20, 40), 2), 4),
    metric = rep(c("AE", "ARE"), each = 3, times = 2),
    sample_size = c(
      100L, 60L, 40L, 300L, 200L, 150L,
      NA, 500L, 400L, 800L, NA, 700L
    ),
    stopping_reason = c(
      "tolerance", "tolerance", "max_iterations", "tolerance", "tolerance", "tolerance",
      "infeasible", "tolerance", "tolerance", "n_max", "infeasible", "tolerance"
    ),
    stringsAsFactors = FALSE
  )
}

test_that("plot_sample_size_by_people builds a log10, alpha-faceted ggplot without warnings", {
  p <- plot_sample_size_by_people(sample_size_df(), "AE", subtitle = "sub")
  expect_s3_class(p, "ggplot")
  expect_no_warning(b <- ggplot2::ggplot_build(p))
  expect_s3_class(p$facet, "FacetWrap")
  y_scale <- p$scales$get_scales("y")
  trans <- if (!is.null(y_scale$trans)) y_scale$trans else y_scale$transformation
  expect_equal(trans$name, "log-10")
  expect_equal(p$labels$subtitle, "sub")
  expect_equal(length(unique(b$layout$layout$PANEL)), 2L)
})

test_that("plot_sample_size_by_people works when a metric has no feasible values", {
  d <- sample_size_df()
  d <- d[d$metric == "ARE", ]
  d$sample_size <- NA_integer_
  d$stopping_reason <- "infeasible"
  p <- plot_sample_size_by_people(d, "ARE")
  expect_s3_class(p, "ggplot")
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("plot_sample_size_by_people validates its inputs", {
  d <- sample_size_df()
  expect_error(plot_sample_size_by_people(list(), "AE"), "data.frame")
  expect_error(plot_sample_size_by_people(d[, c("alpha", "metric")], "AE"), "missing required columns")
  expect_error(plot_sample_size_by_people(d, c("AE", "ARE")), "single character")
  expect_error(plot_sample_size_by_people(d, "nope"), "No rows")
})
