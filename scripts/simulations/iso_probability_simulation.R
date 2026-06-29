# iso_probability_simulation.R
#
# Showcase script for the iso-probability pipeline.
#
# Research question: which (tau_AE, tau_ARE, cutoff) triples produce a joint
# replicate-success probability within `delta` of a fixed baseline triple?
#
# Workflow
#   1. Define a baseline (tau_AE, tau_ARE, cutoff) and estimate its success
#      probability with B_baseline replicates.
#   2. Build a Cartesian candidate design over user-supplied grids of each
#      parameter using make_iso_parameter_design().
#   3. Evaluate every candidate with B_screen replicates and flag those whose
#      estimated probability is within delta of the baseline (point rule).
#   4. Summarise passing candidates (counts, ranges, closest near-misses) with
#      summarize_iso_candidates().
#   5. Visualise results with plot_iso_slice() and
#      plot_iso_probability_comparison().
# ---------------------------------------------------------------------------

simulation_helper_files <- list.files(here::here("scripts", "simulation_layers"))
lapply(simulation_helper_files, function(f) {
  source(here::here("scripts", "simulation_layers", f))
})

# ---- Default configuration -------------------------------------------------

iso_probability_simulation_defaults <- function() {
  list(
    # Proportion-generation scenario
    alpha             = 2,
    K                 = 10L,
    n                 = 500L,
    proportion_method = "beta",
    p_max             = NULL,

    # Baseline parameter triple
    baseline_tau_AE  = 0.10,
    baseline_tau_ARE = 0.20,
    baseline_cutoff  = 0.10,

    # Monte Carlo budget
    B_baseline = 5000L,
    B_screen   = 1000L,

    # Equivalence tolerance
    delta = 0.02,

    # Candidate grids
    tau_AE_grid  = seq(0.05, 0.20, by = 0.05),
    tau_ARE_grid = seq(0.10, 0.40, by = 0.10),
    cutoff_grid  = seq(0.05, 0.20, by = 0.05),

    # Summarisation
    top_k = 10L,

    # Reproducibility
    seed = 260926L
  )
}

# ---- Run iso-probability simulation ----------------------------------------

run_iso_probability_simulation <- function(
    config = iso_probability_simulation_defaults(),
    cache = TRUE,
    force_recompute = FALSE,
    cache_dir = here::here("results", "simresults")) {

  result_file <- simulation_result_path(
    config = config,
    dir    = cache_dir,
    name   = "iso_probability"
  )

  if (cache && !force_recompute && file.exists(result_file)) {
    return(readRDS(result_file))
  }

  # Step 1: build the candidate design
  design <- make_iso_parameter_design(
    tau_AE_grid  = config$tau_AE_grid,
    tau_ARE_grid = config$tau_ARE_grid,
    cutoff_grid  = config$cutoff_grid
  )

  # Step 2: run the search
  iso_result <- run_iso_success_search(
    alpha             = config$alpha,
    K                 = config$K,
    n                 = config$n,
    baseline_tau_AE   = config$baseline_tau_AE,
    baseline_tau_ARE  = config$baseline_tau_ARE,
    baseline_cutoff   = config$baseline_cutoff,
    B_baseline        = config$B_baseline,
    design            = design,
    B_screen          = config$B_screen,
    delta             = config$delta,
    proportion_method = config$proportion_method,
    p_max             = config$p_max,
    seed              = config$seed
  )

  # Step 3: summarise results
  summary <- summarize_iso_candidates(iso_result, top_k = config$top_k)

  # Step 4: produce plots
  plots <- list(
    comparison = plot_iso_probability_comparison(iso_result),
    slice_tau_AE  = plot_iso_slice(
      iso_result, fix = "tau_AE",
      fix_value = config$baseline_tau_AE
    ),
    slice_tau_ARE = plot_iso_slice(
      iso_result, fix = "tau_ARE",
      fix_value = config$baseline_tau_ARE
    ),
    slice_cutoff  = plot_iso_slice(
      iso_result, fix = "cutoff",
      fix_value = config$baseline_cutoff
    )
  )

  result <- list(
    config     = config,
    design     = design,
    iso_result = iso_result,
    summary    = summary,
    plots      = plots
  )

  if (cache) {
    saveRDS(result, result_file)
  }

  result
}
