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
