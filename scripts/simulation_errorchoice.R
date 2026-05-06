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
rm(list = ls())

source(file.path("scripts", "simulation_functions.R"))

# ---- Parameters ------------------------------------------------------------
alpha     <- 3                           # Beta(1, alpha)
K         <- 10L                         # number of cell types
n         <- 1000L                       # reads / total count per sample
B         <- 5000L                       # simulation replicates
# Separate threshold grids for each error metric.
# AE  is bounded in [0, 1]; a range of 0 to 0.5 is sufficient.
# ARE can exceed 1; a range of 0 to 2 covers most practical cases.
taus_AE   <- 10^seq(0, log10(1.5), by = 0.0005) - 1  # threshold grid for AE  (0 to 0.5)
taus_ARE  <- 10^seq(0, log10(3),   by = 0.0005) - 1  # threshold grid for ARE (0 to 2)
taus      <- list(AE = taus_AE, ARE = taus_ARE)
seed      <- 260926L

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

cat("\nArgmax summary (first 20 rows):\n")
print(head(result$argmax_summary, 20))

ggplot(data = result$curves, aes(x = tau, y = success_rate, color = metric)) + geom_point()
ggplot(data = result$argmax, aes(x = AE)) + geom_histogram()
