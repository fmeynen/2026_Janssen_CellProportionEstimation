source(here::here("scripts", "simulations", "simulation_isoband_strict_hybrid.R"))

config <- simulation_isoband_strict_defaults()
# Optional manual override:
# config$p0 <- 0.60
# If config$p0 remains NULL, provide a calibration point inside config$ranges:
config$calibration_point <- list(tau_AE = 0.04, tau_ARE = 0.60, cutoff = 0.10)
result <- run_simulation_isoband_strict(config)

if (!is.null(result$calibration)) {
  print(result$calibration)
}
print(head(result$design_history))
print(head(result$final_band_points))
