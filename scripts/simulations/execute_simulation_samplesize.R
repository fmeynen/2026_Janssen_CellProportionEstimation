#+ echo=TRUE, results='hide'


# simulation_samplesize.R
#
# Iterative sample-size estimation for each alpha.
#
# Workflow:
#   1. For each alpha, run 3 pilot sample sizes (95%, 100%, 105% of n_init).
#   2. Fit glm(success ~ n, family = binomial) to observed success counts.
#   3. Invert the curve to estimate the sample size hitting the target success rate.
#   4. Repeat until convergence or max iterations; always round up.
#   5. Propagate the final estimate as the next alpha's n_init.
# ---------------------------------------------------------------------------

#+ echo=TRUE, results='hide'
simulation_helper_files <- list.files(here::here("scripts", "simulation_layers"))
lapply(simulation_helper_files, function(f) {
  source(here::here("scripts", "simulation_layers", f))
})
library(ggplot2)

# Setup config file -----------------------------------------------------------------------------------------------

simulation_sample_size_defaults <- function(n_init = 200000) {
  tau_AE <- 0.002
  tau_ARE <- 0.05
  list(
    alpha                = seq(from = 2, to = 5, by = 0.05),
    K                    = 10L,
    n                    = n_init,
    B                    = 500L,
    taus                 = list(AE = tau_AE, ARE = tau_ARE),
    metrics              = c("AE", "ARE"),
    model                = "multinomial",
    tie_method           = "random",
    proportion_method    = "beta",
    seed                 = 260925L,
    success_rate_target  = 0.95,
    sample_size_tolerance = 100L,
    max_iterations       = 20L
  )
}

config <- simulation_sample_size_defaults()


# Calculate sample sizes ------------------------------------------------------------------------------------------


run_simulation_samplesize <- function(config = simulation_sample_size_defaults(),
                                      cache = TRUE,
                                      force_recompute = FALSE,
                                      cache_dir = here::here("results", "simresults")) {
  result_file <- simulation_result_path(
    config = config,
    dir    = cache_dir,
    name   = "sample_size"
  )
  if (cache && !force_recompute && file.exists(result_file)) {
    return(readRDS(result_file))
  }

  alphas        <- config$alpha
  n_samples     <- numeric(length(alphas))
  diag_list     <- vector("list", length(alphas))
  n_init        <- config$n

  for (a in seq_along(alphas)) {
    iter_result        <- iterate_sample_size_for_alpha(alphas[a], n_init, config)
    n_samples[a]       <- iter_result$final_n
    diag_list[[a]]     <- iter_result$diagnostics
    n_init             <- iter_result$final_n
  }

  result <- list(
    sample_size  = as.integer(n_samples),
    diagnostics  = do.call(rbind, diag_list)
  )

  if (cache) {
    saveRDS(result, result_file)
  }
  result
}



res <- run_simulation_samplesize()

#+ echo=TRUE, results='markup'
res$sample_size
ggplot(data.frame(alpha = config$alpha, sample_size = res$sample_size),
       aes(x = alpha, y = sample_size)) +
  geom_line() +
  labs(x = "Alpha", y = "Sample Size (log scale)") +
  theme_minimal() +
  scale_y_log10()

ggplot(data.frame(alpha = config$alpha, sample_size = res$sample_size),
       aes(x = alpha, y = sample_size)) +
  geom_line() +
  labs(x = "Alpha", y = "Sample Size") +
  theme_minimal()
