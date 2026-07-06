source(here::here("scripts", "simulations", "simulation_cutoff_threshold_contours.R"))

config <- simulation_contour_defaults()
t_start <- proc.time()
result <- run_simulation_cutoff_threshold_contours(config)
elapsed <- (proc.time() - t_start)[["elapsed"]]

n_grid_points <- nrow(result$grid_results)
cat(sprintf(
  "Run complete in %.1f s | scenarios: %d | cutoffs: %d | grid points: %d\n",
  elapsed,
  nrow(result$p_table),
  length(config$cutoffs),
  n_grid_points
))
cat(sprintf("p_success range: [%.3f, %.3f]\n",
            min(result$grid_results$p_success),
            max(result$grid_results$p_success)))
cat("\nGrid results (first 10 rows):\n")
print(head(result$grid_results, 10))

cat("\nContour plot names:\n")
print(names(result$contour_plots))
