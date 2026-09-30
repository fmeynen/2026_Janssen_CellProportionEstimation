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

test_that("fixed_max_beta proportions plot keeps the fixed maximum off the curve", {
  b <- ggplot2::ggplot_build(plot_proportions_curve(plot_res_fixed()))
  curve_x <- b$data[[1]]$x
  point_x <- b$data[[2]]$x
  expect_true(any(abs(point_x - 1) < 1e-12))
  expect_false(any(abs(curve_x - 1) < 1e-12))
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
