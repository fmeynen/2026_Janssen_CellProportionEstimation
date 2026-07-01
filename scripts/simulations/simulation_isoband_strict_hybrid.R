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
    p0 = 0.60,
    ranges = list(
      tau_AE = c(0.01, 0.08),
      tau_ARE = c(0.10, 1.20),
      cutoff = c(0.02, 0.20)
    ),
    eps = 0.10,
    seed = 260926L,
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
    n_candidates = 500L
  )
}

run_simulation_isoband_strict <- function(config = simulation_isoband_strict_defaults(),
                                          cache = TRUE,
                                          force_recompute = FALSE,
                                          cache_dir = here::here("results", "simresults")) {
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

  result
}
