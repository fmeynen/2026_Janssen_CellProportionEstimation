source(here::here("scripts", "simulations", "simulation_dm_errorchoice.R"))

config <- simulation_dm_errorchoice_defaults()
result <- run_simulation_dm_errorchoice(config)

cat("True proportions (p_table):\n")
print(round(result$p_table, 6))

cat("\ncurves_tau (first 10 rows):\n")
print(head(result$curves_tau, 10))

cat("\ncurves_n (first 10 rows):\n")
print(head(result$curves_n, 10))

for (m in c("AE", "ARE")) {
  print(plot_success_vs_tau(
    result$curves_tau,
    m,
    subtitle = sprintf(
      "n_per_person = %s, concentration = %s",
      config$n_per_person_fixed,
      config$concentration
    )
  ))
  print(plot_success_vs_n(
    result$curves_n,
    m,
    subtitle = sprintf(
      "tau_%s = %s, concentration = %s",
      m,
      config$taus_fixed[[m]],
      config$concentration
    ),
    target = config$target
  ))
}
