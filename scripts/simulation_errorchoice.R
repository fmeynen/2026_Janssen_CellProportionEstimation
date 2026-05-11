# simulation_errorchoice.R
#
# MVP demonstration: multinomial sampling error for K cell-type proportions.
#
# Research question: how often does the *maximum* absolute (or relative) error
# across all K cell types exceed a threshold tau?
#
# Workflow
#   1. Generate true proportions from a monotone Beta curve (deterministic).
#   2. Simulate B multinomial draws; extract max AE and max ARE per replicate.
#   3. Sweep tau grid to obtain success-rate curves.
#   4. Summarise which cell type most often drives the max error.
#
# Possible future changes:
#   * Allow overdispersion  -> model = "dirichlet_multinomial"
#   * Allow correlations    -> model = "logistic_normal_multinomial"
#   * Other monotone curves (exponential, power, ...)
#   * Additional error metrics
#   * Distribution of max errors (not just success rates)
# ---------------------------------------------------------------------------
source(file.path("scripts", "simulation_functions.R"))

# ---- Parameters ------------------------------------------------------------
simulation_errorchoice_defaults <- function() {
  taus_AE <- exp(seq(0, log(1 + 0.05), by = 0.0005)) - 1
  taus_ARE <- 10^seq(0, log10(1 + 3), by = 0.0005) - 1
  list(
    alpha = c(2, 2.5, 3, 4, 5),
    K = 10L,
    n = 1000L,
    B = 5000L,
    taus = list(AE = taus_AE, ARE = taus_ARE),
    metrics = c("AE", "ARE"),
    model = "multinomial",
    tie_method = "random",
    seed = 260926L
  )
}

# ---- Run experiment --------------------------------------------------------
run_simulation_errorchoice <- function(config = simulation_errorchoice_defaults()) {
  run_simulation_experiment(
    alpha = config$alpha,
    K = config$K,
    n = config$n,
    B = config$B,
    taus = config$taus,
    metrics = config$metrics,
    model = config$model,
    tie_method = config$tie_method,
    seed = config$seed
  )
}

# ---- Script execution ------------------------------------------------------
is_simulation_errorchoice_main <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  any(grepl("--file=.*simulation_errorchoice\\.R$", args))
}

if (is_simulation_errorchoice_main()) {
  result <- run_simulation_errorchoice()

  saveRDS(result, file = "results/simresults/simulation_errorchoice.RData")

  cat("True proportions (p):\n")
  print(round(result$p_table, 6))
  generate_proportions_curve(result$inputs$alpha)

  cat("\nSuccess-rate curves (first 10 rows):\n")
  print(head(result$curves, 10))

  cat("\nArgmax summary (first 20 rows):\n")
  print(head(result$argmax_summary, 20))

  plot_success_rate_curve(result, metric = "AE")
  plot_success_rate_curve(result, metric = "ARE")
  plot_success_rate_curve(result, metric = NULL)
  plot_success_rate_curve(result, metric = "AE", alphas = c(2, 5))
}
