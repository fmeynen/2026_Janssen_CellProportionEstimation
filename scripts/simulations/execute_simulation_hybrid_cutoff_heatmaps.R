source(here::here("scripts", "simulations", "simulation_hybrid_cutoff_heatmaps.R"))


config <- simulation_hybrid_heatmap_defaults()
result <- run_simulation_hybrid_heatmaps(config)
  
print(result$p_table)
