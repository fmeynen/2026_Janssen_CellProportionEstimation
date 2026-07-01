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
    calibration_point = NULL,
    # High replicate count for stable single-point p0 calibration (still user-configurable).
    calibration_B = 5000L,
    # Fixed default seed for reproducible calibration draws (arbitrary, user-configurable).
    calibration_seed = 270001L
  )
}

validate_isoband_calibration_point <- function(point, ranges) {
  ranges <- validate_isoband_ranges(ranges)
  required <- c("tau_AE", "tau_ARE", "cutoff")
  if (!is.list(point) || !all(required %in% names(point))) {
    stop(
      sprintf("calibration_point must be a named list containing: %s", paste(required, collapse = ", ")),
      call. = FALSE
    )
  }

  validated_point <- lapply(required, function(name) {
    value <- point[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop(
        sprintf("calibration_point$%s must be numeric, length 1, and finite.", name),
        call. = FALSE
      )
    }
    as.numeric(value)
  })
  names(validated_point) <- required

  in_range <- function(value, bounds) {
    value >= bounds[[1L]] && value <= bounds[[2L]]
  }
  if (!in_range(validated_point$tau_AE, ranges$tau_AE)) {
    stop("calibration_point$tau_AE must lie within ranges$tau_AE.", call. = FALSE)
  }
  if (!in_range(validated_point$tau_ARE, ranges$tau_ARE)) {
    stop("calibration_point$tau_ARE must lie within ranges$tau_ARE.", call. = FALSE)
  }
  if (!in_range(validated_point$cutoff, ranges$cutoff)) {
    stop("calibration_point$cutoff must lie within ranges$cutoff.", call. = FALSE)
  }

  validated_point
}

compute_isoband_calibration_scenario_hash <- function(config) {
  calibration_point <- validate_isoband_calibration_point(config$calibration_point, config$ranges)
  scenario <- list(
    alpha = config$alpha,
    n = config$n,
    K = config$K,
    proportion_method = config$proportion_method,
    p_max = config$p_max,
    model = config$model,
    ranges = validate_isoband_ranges(config$ranges),
    calibration_point = calibration_point
  )
  hash_config(scenario)
}

isoband_calibrated_p0_path <- function(config, dir = here::here("results", "calibrated_p")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  scenario_hash <- compute_isoband_calibration_scenario_hash(config)
  file.path(dir, paste0("p0_", scenario_hash, ".rds"))
}

validate_open_probability <- function(value, label = "p0") {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0 || value >= 1) {
    stop(sprintf("%s must be a finite scalar strictly between 0 and 1.", label), call. = FALSE)
  }
}

calibrate_isoband_p0 <- function(config,
                                 cache = TRUE,
                                 force_recompute = FALSE,
                                 cache_dir = here::here("results", "calibrated_p")) {
  ranges <- validate_isoband_ranges(config$ranges)
  calibration_point <- validate_isoband_calibration_point(config$calibration_point, ranges)
  cache_file <- isoband_calibrated_p0_path(config = config, dir = cache_dir)

  if (cache && !force_recompute && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  calibration_design <- data.frame(
    tau_AE = calibration_point$tau_AE,
    tau_ARE = calibration_point$tau_ARE,
    cutoff = calibration_point$cutoff,
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
  if (!is.data.frame(calibration_result) ||
      !("success_rate" %in% names(calibration_result)) ||
      nrow(calibration_result) < 1L) {
    stop("Calibration failed: expected non-empty calibration result with success_rate.", call. = FALSE)
  }

  p0_calibrated <- calibration_result$success_rate[[1L]]
  validate_open_probability(p0_calibrated, label = "Calibrated p0")
  out <- list(
    p0 = p0_calibrated,
    n_success = calibration_result$n_success[[1L]],
    n_total = calibration_result$n_total[[1L]],
    calibration_point = calibration_point,
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
  config_local <- config
  config_local$ranges <- validate_isoband_ranges(config_local$ranges)

  calibration_info <- NULL
  if (is.null(config_local$p0)) {
    if (is.null(config_local$calibration_point)) {
      stop("config$calibration_point must be provided when config$p0 is NULL.", call. = FALSE)
    }
    calibration_info <- calibrate_isoband_p0(config = config_local, cache = cache, force_recompute = force_recompute)
    config_local$p0 <- calibration_info$p0
  }
  validate_open_probability(config_local$p0, label = "Resolved p0")

  result_file <- isoband_result_path(
    config = config_local,
    dir = cache_dir,
    name = "simulation_isoband_strict_hybrid"
  )

  if (cache && !force_recompute && file.exists(result_file)) {
    return(readRDS(result_file))
  }

  result <- run_isoband_pipeline(
    alpha = config_local$alpha,
    p0 = config_local$p0,
    ranges = config_local$ranges,
    eps = config_local$eps,
    seed = config_local$seed,
    n = config_local$n,
    K = config_local$K,
    proportion_method = config_local$proportion_method,
    p_max = config_local$p_max,
    model = config_local$model,
    n_init = config_local$n_init,
    R_init = config_local$R_init,
    n_rounds_max = config_local$n_rounds_max,
    n_add = config_local$n_add,
    R_final = config_local$R_final,
    n_candidates = config_local$n_candidates
  )

  if (cache) {
    save_isoband_result(
      result = result,
      config = config_local,
      dir = cache_dir,
      name = "simulation_isoband_strict_hybrid"
    )
  }

  result$calibration <- calibration_info

  result
}
