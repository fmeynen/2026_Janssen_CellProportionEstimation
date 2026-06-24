source(here::here("scripts", "simulations", "simulation_errorchoice.R"))

config <- simulation_errorchoice_defaults()
result <- run_simulation_errorchoice(config)
  
cat("True proportions (p):\n")
print(round(result$p_table, 6))
plot_proportions_curve(result = result)
  
cat("\nSuccess-rate curves (first 10 rows):\n")
print(head(result$curves, 10))

cat("\nArgmax summary (first 20 rows):\n")
print(head(result$argmax_summary, 20))
  
plot_success_rate_curve(result, metric = "AE")
plot_success_rate_curve(result, metric = "ARE")
plot_success_rate_curve(result, metric = "LAE")
plot_success_rate_curve(result, metric = "TSE")
plot_success_rate_curve(result, metric = NULL)
plot_success_rate_curve(result, metric = "AE", alphas = c(2, 5))
  
plot_argmax_histogram(result, metric = "AE")
print(result$p_table)
print(head(result$curves))

