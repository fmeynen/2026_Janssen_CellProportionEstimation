source(here::here("scripts", "simulations", "simulation_isoband_strict_hybrid.R"))

config <- simulation_isoband_strict_defaults()
# Optional manual override:
# config$p0 <- 0.60

result <- run_simulation_isoband_strict(config)



if (!is.null(result$calibration)) {
  print(result$calibration)
}
print(head(result$design_history))
print(head(result$final_band_points))

plot_isoband_point_cloud(result)


result$design_history
# Scratchpad ------------------------------------------------------------------------------------------------------

result$final_band_points_unrefined



ranges <- validate_isoband_ranges(config$ranges)
midpoint <- compute_isoband_midpoint(ranges)
calibration_design <- data.frame(
  tau_AE = 0.005,
  tau_ARE = 0.05,
  cutoff = 0.01,
  stringsAsFactors = FALSE
)
calibration_result <- simulate_isoband_design_points(
  design_df = calibration_design,
  alpha = config$alpha,
  K = config$K,
  n = config$n,
  B_point = as.integer(config$calibration_B),
  proportion_method = config$proportion_method,
  p_max = config$p_max,
  model = config$model,
  seed = config$calibration_seed
)
calibration_result
