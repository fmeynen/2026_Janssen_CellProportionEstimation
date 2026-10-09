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
#   2. Warm start: the first alpha starts from config$n_init (or config$concentration when n_init is NULL); later
#      alphas start from the last feasible alpha's final sample size.
#      Alphas whose success rate at n_max cannot reach the target are reported as infeasible (NA sample size) and
#      skipped.
#   3. Each alpha's result is cached to its own file (per-alpha cache, keyed on a config that includes that
#      alpha's warm-start n_init), so re-running only recomputes alphas whose inputs actually changed.
# ---------------------------------------------------------------------------

#+ echo=TRUE, results='hide'
source(here::here("scripts", "load_layers.R"))
library(ggplot2)

# Setup config file -----------------------------------------------------------------------------------------------

simulation_sample_size_defaults <- function() {
  list(
    alpha = 2,
    K = 10L,
    n_init = NULL, # NULL means "defaults to concentration"
    B = 500L,
    taus = list(AE = 0.02, ARE = 2),
    metrics = c("AE", "ARE"),
    n_people = 2,
    concentration = 50,
    model = "dirichlet_multinomial",
    tie_method = "random",
    proportion_method = "beta",
    # For a fixed smallest proportion: proportion_method = "fixed_min_beta", p_min = 0.01
    p_min = NULL,
    p_max = NULL,
    seed = 260925L,
    success_rate_target = 0.95,
    rel_tol = 0.01,
    max_iterations = 20L,
    f0 = 2,
    f_floor = 1.1,
    n_max = 1e9
  )
}

config <- simulation_sample_size_defaults()


# Calculate sample sizes ------------------------------------------------------------------------------------------

res <- run_sample_size_experiment(config)

#+ echo=TRUE, results='markup'
print(res$sample_size)

feasible <- subset(res$sample_size, !is.na(sample_size))

ggplot(feasible, aes(x = alpha, y = sample_size)) +
  geom_line() +
  geom_point() +
  labs(x = "Alpha", y = "Sample Size (log scale)") +
  theme_minimal() +
  scale_y_log10()

ggplot(feasible, aes(x = alpha, y = sample_size)) +
  geom_line() +
  geom_point() +
  labs(x = "Alpha", y = "Sample Size") +
  theme_minimal()
