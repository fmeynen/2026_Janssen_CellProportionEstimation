# simulation_success_surface.R
#
# Driver: compute and visualise success-probability contours over an
# AE x ARE threshold grid for a given hybrid cutoff.
#
# The workflow:
#   1. Generate true proportions for the given alpha / K.
#   2. Simulate B replicates.
#   3. Call estimate_success_surface() for a fixed cutoff c.
#   4. Plot contour lines via plot_success_contours().
#
# Usage (interactive):
#   source(here::here("scripts", "simulations", "simulation_success_surface.R"))
#   config <- success_surface_defaults()
#   result <- run_success_surface(config)
#   print(result$contour_plot)

simulation_helper_files <- list.files(here::here("scripts", "simulation_layers"))
lapply(simulation_helper_files, function(f) {
  source(here::here("scripts", "simulation_layers", f))
})


#' Default configuration for the success-surface workflow.
#'
#' @return Named list of default parameters.
success_surface_defaults <- function() {
  list(
    alpha    = 2,
    K        = 10L,
    n        = 1000L,
    B        = 500L,
    cutoff   = 0.05,
    ae_grid  = seq(0.001, 0.10,  length.out = 60L),
    are_grid = seq(0.01,  2.0,   length.out = 60L),
    seed     = 260926L
  )
}


#' Run the success-surface estimation and produce a contour plot.
#'
#' Reuses [run_replicates()] for simulation and [estimate_success_surface()]
#' for the grid evaluation.  Results are optionally cached as an `.rds` file.
#'
#' @param config          Named list; see [success_surface_defaults()].
#' @param cache           Logical; save / load result from `cache_dir`.
#' @param force_recompute Logical; ignore any existing cached result.
#' @param cache_dir       Directory for cached `.rds` files.
#' @param contour_levels  Numeric vector of probability levels for contour
#'   lines (default `seq(0.1, 0.9, 0.1)`).
#' @param x_limits        AE axis display window passed to
#'   [plot_success_contours()] (default `c(0.001, 0.05)`).
#' @param y_limits        ARE axis display window passed to
#'   [plot_success_contours()] (default `c(0.01, 1)`).
#'
#' @return List with elements:
#'   \describe{
#'     \item{config}{Copy of the configuration used.}
#'     \item{p}{True proportion vector of length K.}
#'     \item{surface_df}{Data.frame from [estimate_success_surface()].}
#'     \item{contour_plot}{ggplot contour-line plot.}
#'   }
run_success_surface <- function(config = success_surface_defaults(),
                                cache = TRUE,
                                force_recompute = FALSE,
                                cache_dir = here::here("results", "simresults"),
                                contour_levels = seq(0.1, 0.9, 0.1),
                                x_limits = c(0.001, 0.05),
                                y_limits = c(0.01, 1)) {
  result_file <- simulation_result_path(
    config = config,
    dir    = cache_dir,
    name   = "success_surface"
  )

  if (cache && !force_recompute && file.exists(result_file)) {
    return(readRDS(result_file))
  }

  p <- generate_proportions(alpha = config$alpha, K = config$K)

  rep_out <- run_replicates(
    p       = p,
    n       = config$n,
    B       = config$B,
    metrics = c("AE", "ARE"),
    seed    = config$seed
  )

  surface_df <- estimate_success_surface(
    phat_mat = rep_out$phat,
    p        = p,
    cutoff   = config$cutoff,
    ae_grid  = config$ae_grid,
    are_grid = config$are_grid
  )

  contour_plot <- plot_success_contours(
    surface_df     = surface_df,
    contour_levels = contour_levels,
    x_limits       = x_limits,
    y_limits       = y_limits
  )

  result <- list(
    config       = config,
    p            = p,
    surface_df   = surface_df,
    contour_plot = contour_plot
  )

  if (cache) {
    saveRDS(result, result_file)
  }

  result
}
