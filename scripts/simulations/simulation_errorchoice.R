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
simulation_helper_files <- list.files(here::here("scripts", "simulation_layers"))
lapply(simulation_helper_files, function(f){
  source(here::here("scripts", "simulation_layers", f))
})

# ---- Parameters ------------------------------------------------------------
simulation_errorchoice_defaults <- function() {
  taus_AE <- exp(seq(0, log(1 + 0.02), by = 0.0005)) - 1
  taus_ARE <- 10^seq(0, log10(1 + 1.5), by = 0.0005) - 1
  taus_TSE <- 10^seq(0, log10(1 + 7.5), by = 0.005) - 1
  taus_LAE <- 10^seq(0, log10(1 + 0.001), by = 0.000005) - 1
  list(
    alpha = c(2, 2.5, 3, 4, 5),
    K = 10L,
    n = 1000L,
    B = 5000L,
    taus = list(AE = taus_AE, ARE = taus_ARE, TSE = taus_TSE, LAE = taus_LAE),
    metrics = c("AE", "ARE", "LAE", "TSE"),
    model = "multinomial",
    tie_method = "random",
    proportion_method = "beta",
    #p_max = c(0.2, 0.3,0.4, 0.5),
    seed = 260925L
  )
}

# ---- Run experiment --------------------------------------------------------
run_simulation_errorchoice <- function(config = simulation_errorchoice_defaults(),
                                       cache = TRUE,
                                       force_recompute = FALSE,
                                       cache_dir = here::here("results", "simresults")) {
  result_file <- simulation_result_path(
    config = config,
    dir    = cache_dir,
    name   = "errorchoice"
  )

  if (cache && !force_recompute && file.exists(result_file)) {
    return(readRDS(result_file))
  }

  result <- run_simulation_experiment(
    alpha = config$alpha,
    K = config$K,
    n = config$n,
    B = config$B,
    taus = config$taus,
    metrics = config$metrics,
    model = config$model,
    tie_method = config$tie_method,
    proportion_method = config$proportion_method,
    p_max = config$p_max,
    seed = config$seed
  )

  if (cache) {
    saveRDS(result, result_file)
  }

  result
}
