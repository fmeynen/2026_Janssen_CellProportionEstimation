# Shared cases for config-level bound validation (validate_proportion_bounds() and its entry points).
# Each case: proportion_method, p_min, p_max, and a regex the error must match (K = 10).
bound_error_cases <- list(
  min_missing    = list(method = "fixed_min_beta", p_min = NULL, p_max = NULL, regex = "p_min must be provided"),
  max_missing    = list(method = "fixed_max_beta", p_min = NULL, p_max = NULL, regex = "p_max must be provided"),
  min_invalid    = list(method = "fixed_min_beta", p_min = 1.5, p_max = NULL, regex = "p_min must contain"),
  max_invalid    = list(method = "fixed_max_beta", p_min = NULL, p_max = 0, regex = "p_max must contain"),
  min_vector     = list(method = "fixed_min_beta", p_min = c(0.01, 0.02), p_max = NULL, regex = "p_min must be a single"),
  max_vector     = list(method = "fixed_max_beta", p_min = NULL, p_max = c(0.3, 0.4), regex = "p_max must be a single"),
  min_unused_beta = list(method = "beta", p_min = 0.01, p_max = NULL, regex = "p_min is not used by method = 'beta'"),
  max_unused_beta = list(method = "beta", p_min = NULL, p_max = 0.4, regex = "p_max is not used by method = 'beta'"),
  max_unused_min = list(method = "fixed_min_beta", p_min = 0.01, p_max = 0.4,
                        regex = "p_max is not used by method = 'fixed_min_beta'"),
  min_unused_max = list(method = "fixed_max_beta", p_min = 0.01, p_max = 0.4,
                        regex = "p_min is not used by method = 'fixed_max_beta'"),
  min_infeasible = list(method = "fixed_min_beta", p_min = 0.2, p_max = NULL, regex = "Impossible fixed_min_beta"),
  max_infeasible = list(method = "fixed_max_beta", p_min = NULL, p_max = 0.05, regex = "Impossible fixed_max_beta")
)
