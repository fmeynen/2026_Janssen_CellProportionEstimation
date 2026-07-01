source(here::here("scripts", "simulations", "simulation_isoband_strict_hybrid.R"))

config <- simulation_isoband_strict_defaults()
result <- run_simulation_isoband_strict(config)

print(head(result$design_history))
print(head(result$final_band_points))
