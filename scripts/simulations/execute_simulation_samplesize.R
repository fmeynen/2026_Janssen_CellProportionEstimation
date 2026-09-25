#+ echo=TRUE, results='hide'


# simulation_samplesize.R
#
# Sample-size estimation for each alpha, using the iterative solver on log(n) (`estimate_sample_size()` in
# scripts/simulation_layers/calculation.R), orchestrated across the alpha grid by
# `run_sample_size_experiment()` (scripts/simulation_layers/orchestration.R).
#
# Workflow:
#   1. For each alpha, in order, run the solver: fit pilots around a centre, invert the logistic success curve on
#      log(n), and iterate until the relative change in n is within tolerance or max_iterations is reached.
#   2. Warm start: the first alpha starts from config$n_init; every later alpha starts from the previous alpha's
#      final sample size.
#   3. Each alpha's result is cached to its own file (per-alpha cache, keyed on a config that includes that
#      alpha's warm-start n_init), so re-running only recomputes alphas whose inputs actually changed.
# ---------------------------------------------------------------------------

#+ echo=TRUE, results='hide'
simulation_helper_files <- list.files(here::here("scripts", "simulation_layers"))
lapply(simulation_helper_files, function(f) {
  source(here::here("scripts", "simulation_layers", f))
})
library(ggplot2)

# Setup config file -----------------------------------------------------------------------------------------------

simulation_sample_size_defaults <- function(n_init = 200000) {
  list(
    alpha                = seq(from = 2, to = 5, by = 0.5),
    K                    = 10L,
    n_init               = n_init,
    B                    = 500L,
    taus                 = list(AE = 0.002, ARE = 0.05),
    metrics              = c("AE", "ARE"),
    n_people             = 2,
    concentration        = 50,
    model                = "dirichlet_multinomial",
    # `n_init` (and, at solve time, `n`) then means cells sampled per person, not total cells.
    n_people             = NULL,
    concentration        = NULL,
    tie_method           = "random",
    proportion_method    = "beta",
    seed                 = 260925L,
    success_rate_target  = 0.95,
    rel_tol              = 0.01,
    max_iterations       = 20L,
    f0                   = 2,
    f_floor              = 1.1,
    n_max                = 1e9
  )
}

config <- simulation_sample_size_defaults()


# Calculate sample sizes ------------------------------------------------------------------------------------------

res <- run_sample_size_experiment(config)

#+ echo=TRUE, results='markup'
print(res$sample_size)

ggplot(res$sample_size, aes(x = alpha, y = sample_size)) +
  geom_line() +
  labs(x = "Alpha", y = "Sample Size (log scale)") +
  theme_minimal() +
  scale_y_log10()

ggplot(res$sample_size, aes(x = alpha, y = sample_size)) +
  geom_line() +
  labs(x = "Alpha", y = "Sample Size") +
  theme_minimal()
