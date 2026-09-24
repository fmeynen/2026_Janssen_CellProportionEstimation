# scripts/deprecated/test_sample_size_deprecated.R
#
# Base-R test script for the deprecated iterative-GLM sample-size helpers
# (fit_success_glm, solve_sample_size_from_glm, iterate_sample_size_for_alpha).
# These functions were moved out of the active simulation layers into
# scripts/deprecated/deprecated.R and are being replaced by a new solver.
# This test script is not part of the active test suite; it must source
# scripts/deprecated/deprecated.R (and the simulation layer it depends on,
# e.g. simulate_success_at_n from scripts/simulation_layers/simulation.R)
# before running.
#
# Run from the repository root:
#   Rscript scripts/deprecated/test_sample_size_deprecated.R
# ---------------------------------------------------------------------------

source(file.path("scripts", "simulation_layers", "validation_utils.R"))
source(file.path("scripts", "simulation_layers", "calculation.R"))
source(file.path("scripts", "simulation_layers", "extraction.R"))
source(file.path("scripts", "simulation_layers", "simulation.R"))
source(file.path("scripts", "simulation_layers", "orchestration.R"))
source(file.path("scripts", "deprecated", "deprecated.R"))

pass <- function(name) cat(sprintf("  PASS  %s\n", name))

# ---- Sample-size estimation helpers ----------------------------------------
cat("\nSample-size estimation helpers\n")

ss_config <- list(
  K                    = 5L,
  B                    = 20L,
  taus                 = list(AE = 0.05, ARE = 2.0),
  metrics              = c("AE", "ARE"),
  model                = "multinomial",
  tie_method           = "random",
  proportion_method    = "beta",
  seed                 = 42L,
  success_rate_target  = 0.8,
  sample_size_tolerance = 10L,
  max_iterations       = 3L
)

# simulate_success_at_n: basic shape checks
sim_succ <- simulate_success_at_n(alpha = 2, n = 500L, config = ss_config)
if (!is.logical(sim_succ$success)) stop("success must be a logical vector")
if (length(sim_succ$success) != ss_config$B) stop("success must have length B")
if (!is.numeric(sim_succ$success_rate) || sim_succ$success_rate < 0 || sim_succ$success_rate > 1) {
  stop("success_rate must be in [0, 1]")
}
if (sim_succ$success_count != sum(sim_succ$success)) stop("success_count must equal sum(success)")
pass("simulate_success_at_n returns correct structure")

# Success requires both AE and ARE thresholds when both are in metrics
# Use extreme thresholds to force known outcomes
cfg_all_pass <- ss_config
cfg_all_pass$taus <- list(AE = 100, ARE = 100)  # everything passes
sim_all <- simulate_success_at_n(alpha = 2, n = 200L, config = cfg_all_pass)
if (!all(sim_all$success)) stop("all replicates should succeed with very loose thresholds")
pass("all replicates succeed with very loose thresholds")

cfg_none_pass <- ss_config
cfg_none_pass$taus <- list(AE = 0, ARE = 0)  # nothing passes
sim_none <- simulate_success_at_n(alpha = 2, n = 200L, config = cfg_none_pass)
if (any(sim_none$success)) stop("no replicates should succeed with zero thresholds")
pass("no replicates succeed with zero thresholds")

# Only AE metric: success based only on AE threshold
cfg_ae_only <- ss_config
cfg_ae_only$metrics <- "AE"
cfg_ae_only$taus    <- list(AE = 0.05)
sim_ae <- simulate_success_at_n(alpha = 2, n = 500L, config = cfg_ae_only)
if (length(sim_ae$success) != cfg_ae_only$B) stop("wrong length for AE-only success")
pass("simulate_success_at_n works with single metric")

# fit_success_glm
n_vals <- c(100L, 200L, 300L)
s_cnt  <- c(5L,   12L,  18L)
glm_fit <- fit_success_glm(n_vals, s_cnt, B = 20L)
if (!inherits(glm_fit, "glm")) stop("fit_success_glm must return a glm object")
if (length(stats::coef(glm_fit)) != 2L) stop("glm must have intercept and slope")
pass("fit_success_glm returns a binomial glm with two coefficients")

# solve_sample_size_from_glm: result is always rounded up
solved <- solve_sample_size_from_glm(glm_fit, target_success_rate = 0.8)
if (!is.list(solved) || is.null(solved$n_raw) || is.null(solved$n_rounded)) {
  stop("solve_sample_size_from_glm must return a list with n_raw and n_rounded")
}
if (solved$n_rounded != as.integer(ceiling(solved$n_raw))) {
  stop("n_rounded must equal ceiling(n_raw)")
}
pass("solve_sample_size_from_glm rounds up")

# iterate_sample_size_for_alpha: diagnostics structure
iter_result <- iterate_sample_size_for_alpha(alpha = 2, n_init = 200L, config = ss_config)
if (!is.list(iter_result)) stop("iterate_sample_size_for_alpha must return a list")
if (!is.integer(iter_result$final_n) || iter_result$final_n < 1L) {
  stop("final_n must be a positive integer")
}
if (!iter_result$stopping_reason %in% c("tolerance", "max_iterations")) {
  stop("stopping_reason must be 'tolerance' or 'max_iterations'")
}
req_cols <- c("alpha", "iteration", "pilot_index", "n", "success_count",
              "success_rate", "target_success_rate", "glm_intercept",
              "glm_slope", "n_raw", "n_rounded", "n_clamped", "stopping_reason")
missing_cols <- setdiff(req_cols, names(iter_result$diagnostics))
if (length(missing_cols) > 0L) {
  stop(paste("diagnostics missing columns:", paste(missing_cols, collapse = ", ")))
}
if (nrow(iter_result$diagnostics) %% 3L != 0L) {
  stop("diagnostics must have a multiple of 3 rows (3 pilot points per iteration)")
}
pass("iterate_sample_size_for_alpha returns correct structure")

# Clamping: when success rate is near 0%, new_n should not exceed 2x the largest pilot
cfg_clamp_low <- ss_config
cfg_clamp_low$taus           <- list(AE = 0, ARE = 0)  # force 0% success
cfg_clamp_low$max_iterations <- 1L
cfg_clamp_low$sample_size_tolerance <- Inf
n_init_clamp <- 200L
iter_clamp_low <- iterate_sample_size_for_alpha(alpha = 2, n_init = n_init_clamp, config = cfg_clamp_low)
upper_bound <- as.integer(ceiling(2.0 * as.integer(ceiling(1.05 * n_init_clamp))))
if (iter_clamp_low$final_n > upper_bound) {
  stop("final_n must not exceed 2x the largest pilot when success rate is near 0%")
}
pass("clamping prevents new_n from exceeding 2x largest pilot when success is ~0%")

# Clamping: when success rate is near 100%, new_n should not be less than 0.5x the smallest pilot
cfg_clamp_high <- ss_config
cfg_clamp_high$taus           <- list(AE = 100, ARE = 100)  # force 100% success
cfg_clamp_high$max_iterations <- 1L
cfg_clamp_high$sample_size_tolerance <- Inf
iter_clamp_high <- iterate_sample_size_for_alpha(alpha = 2, n_init = n_init_clamp, config = cfg_clamp_high)
lower_bound <- as.integer(ceiling(0.5 * as.integer(ceiling(0.95 * n_init_clamp))))
if (iter_clamp_high$final_n < lower_bound) {
  stop("final_n must not go below 0.5x the smallest pilot when success rate is ~100%")
}
pass("clamping prevents new_n from going below 0.5x smallest pilot when success is ~100%")

# Stopping at tolerance: use tolerance = Inf so it stops in 1 iteration
cfg_tol_inf <- ss_config
cfg_tol_inf$sample_size_tolerance <- Inf
iter_inf <- iterate_sample_size_for_alpha(alpha = 2, n_init = 200L, config = cfg_tol_inf)
if (iter_inf$stopping_reason != "tolerance") {
  stop("stopping_reason should be 'tolerance' when tolerance = Inf")
}
if (iter_inf$iterations_used != 1L) stop("should stop after 1 iteration when tolerance = Inf")
pass("stops after first iteration when tolerance is Inf")

# max_iterations stopping: use tolerance = 0 and 1 iteration
cfg_max1 <- ss_config
cfg_max1$sample_size_tolerance <- 0L
cfg_max1$max_iterations        <- 1L
iter_max1 <- iterate_sample_size_for_alpha(alpha = 2, n_init = 500L, config = cfg_max1)
if (iter_max1$iterations_used > 1L) stop("should use at most 1 iteration")
pass("max_iterations respected")

# Alpha propagation: final_n of first alpha feeds next alpha's n_init
cfg_prop <- ss_config
cfg_prop$sample_size_tolerance <- Inf  # always stop after 1 iteration

iter_a1 <- iterate_sample_size_for_alpha(alpha = 2,   n_init = 300L, config = cfg_prop)
iter_a2 <- iterate_sample_size_for_alpha(alpha = 2.5, n_init = iter_a1$final_n, config = cfg_prop)

# The second alpha's diagnostics$n values should be based on iter_a1$final_n
expected_pilots <- as.integer(ceiling(c(0.95, 1.00, 1.05) * iter_a1$final_n))
actual_pilots   <- iter_a2$diagnostics$n[iter_a2$diagnostics$iteration == 1L]
if (!identical(sort(actual_pilots), sort(expected_pilots))) {
  stop("alpha propagation: second alpha must use first alpha's final_n as n_init")
}
pass("final_n from one alpha is correctly propagated to the next alpha's n_init")

# ---- Done ------------------------------------------------------------------
cat("\nAll deprecated sample-size tests passed.\n")
