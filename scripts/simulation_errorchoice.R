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
alpha     <- 2                          # Beta(1, 2): moderately skewed proportions
K         <- 10L                        # number of cell types
n         <- 1000L                      # reads / total count per sample
B         <- 500L                       # simulation replicates
taus      <- seq(0, 0.2, by = 0.005)   # threshold grid for success-rate curves
seed      <- 42L

# ---- Run experiment --------------------------------------------------------
result <- run_simulation_experiment(
  alpha      = alpha,
  K          = K,
  n          = n,
  B          = B,
  taus       = taus,
  metrics    = c("AE", "ARE"),
  model      = "multinomial",
  tie_method = "random",
  seed       = seed
)

# ---- Output ----------------------------------------------------------------
cat("True proportions (p):\n")
print(round(result$p, 4))

cat("\nSuccess-rate curves (first 10 rows):\n")
print(head(result$curves, 10))

cat("\nArgmax summary (first 10 rows):\n")
print(head(result$argmax_summary, 10))
