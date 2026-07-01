# simulation_isoband_strict_hybrid.R
#
# Driver script for the strict-hybrid isoband automation pipeline.

simulation_helper_files <- list.files(here::here("scripts", "simulation_layers"))
lapply(simulation_helper_files, function(f) {
  source(here::here("scripts", "simulation_layers", f))
})

simulation_isoband_strict_defaults <- function() {
  list(
    alpha = 2,
    p0 = NULL,
    ranges = list(
      tau_AE = c(0.01, 0.08),
      tau_ARE = c(0.10, 1.20),
      cutoff = c(0.02, 0.20)
    ),
    eps = 0.10,
    seed = 260925L,
    n = 200L,
    K = 5L,
    proportion_method = "beta",
    p_max = NULL,
    model = "multinomial",
    n_init = 12L,
    R_init = 20L,
    n_rounds_max = 1L,
    n_add = 8L,
    R_final = 30L,
    n_candidates = 500L,
    calibration_B = 5000L,
    calibration_seed = 270001L
  )
}

isoband_midpoint_point <- function(ranges) {
  ranges <- validate_isoband_ranges(ranges)
  list(
    tau_AE = mean(ranges$tau_AE),
    tau_ARE = mean(ranges$tau_ARE),
    cutoff = mean(ranges$cutoff)
  )
}

compute_isoband_calibration_scenario_hash <- function(config) {
  midpoint <- isoband_midpoint_point(config$ranges)
  scenario <- list(
    alpha = config$alpha,
    n = config$n,
    K = config$K,
    proportion_method = config$proportion_method,
    p_max = config$p_max,
    model = config$model,
    ranges = validate_isoband_ranges(config$ranges),
    calibration_point = midpoint
  )
  hash_config(scenario)
}

isoband_calibrated_p0_path <- function(config, dir = here::here("results", "calibrated_p")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  scenario_hash <- compute_isoband_calibration_scenario_hash(config)
  file.path(dir, paste0("p0_", scenario_hash, ".rds"))
}

calibrate_isoband_p0 <- function(config,
                                 cache = TRUE,
                                 force_recompute = FALSE,
                                 cache_dir = here::here("results", "calibrated_p")) {
  ranges <- validate_isoband_ranges(config$ranges)
  midpoint <- isoband_midpoint_point(ranges)
  cache_file <- isoband_calibrated_p0_path(config = config, dir = cache_dir)

  if (cache && !force_recompute && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  calibration_design <- data.frame(
    tau_AE = midpoint$tau_AE,
    tau_ARE = midpoint$tau_ARE,
    cutoff = midpoint$cutoff,
    stringsAsFactors = FALSE
  )
  calibration_result <- simulate_isoband_design_points(
    design_df = calibration_design,
    alpha = config$alpha,
    K = config$K,
    n = config$n,
    B_point = as.integer(config$calibration_B),
    proportion_method = config$proportion_method,
    p_max = config$p_max,
    model = config$model,
    seed = config$calibration_seed
  )

  p0_calibrated <- calibration_result$success_rate[[1L]]
  out <- list(
    p0 = p0_calibrated,
    n_success = calibration_result$n_success[[1L]],
    n_total = calibration_result$n_total[[1L]],
    success_rate = calibration_result$success_rate[[1L]],
    calibration_point = midpoint,
    scenario_hash = compute_isoband_calibration_scenario_hash(config)
  )

  if (cache) {
    saveRDS(out, cache_file)
  }

  out
}

run_simulation_isoband_strict <- function(config = simulation_isoband_strict_defaults(),
                                          cache = TRUE,
                                          force_recompute = FALSE,
                                          cache_dir = here::here("results", "simresults")) {
  config$ranges <- validate_isoband_ranges(config$ranges)

  calibration_info <- NULL
  if (is.null(config$p0)) {
    calibration_info <- calibrate_isoband_p0(config = config, cache = cache, force_recompute = force_recompute)
    config$p0 <- calibration_info$p0
  }

  result_file <- isoband_result_path(
    config = config,
    dir = cache_dir,
    name = "simulation_isoband_strict_hybrid"
  )

  if (cache && !force_recompute && file.exists(result_file)) {
    return(readRDS(result_file))
  }

  result <- run_isoband_pipeline(
    alpha = config$alpha,
    p0 = config$p0,
    ranges = config$ranges,
    eps = config$eps,
    seed = config$seed,
    n = config$n,
    K = config$K,
    proportion_method = config$proportion_method,
    p_max = config$p_max,
    model = config$model,
    n_init = config$n_init,
    R_init = config$R_init,
    n_rounds_max = config$n_rounds_max,
    n_add = config$n_add,
    R_final = config$R_final,
    n_candidates = config$n_candidates
  )

  if (cache) {
    save_isoband_result(
      result = result,
      config = config,
      dir = cache_dir,
      name = "simulation_isoband_strict_hybrid"
    )
  }

  result$calibration <- calibration_info

  result
}
