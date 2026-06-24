# simulation_hybrid_cutoff_heatmaps.R
#
# Small driver script to simulate hybrid AE/ARE cutoff heatmaps
# for multiple alpha values.

simulation_helper_files <- list.files(here::here("scripts", "simulation_layers"))
lapply(simulation_helper_files, function(f){
  source(here::here("scripts", "simulation_layers", f))
})

simulation_hybrid_heatmap_defaults <- function() {
  list(
    alpha = c(2, 2.5, 3, 4, 5),
    K = 10L,
    n = 1000L,
    B = 500L,
    cutoffs = seq(0, 0.5, by = 0.002),
    tau_AE_values = c(0.01, 0.02, 0.03, 0.04, 0.05),
    tau_ARE_values = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2),
    maximize = "replicate",
    seed = 260926L
  )
}

run_simulation_hybrid_heatmaps <- function(config = simulation_hybrid_heatmap_defaults(),
                                           cache = TRUE,
                                           force_recompute = FALSE,
                                           cache_dir = here::here("results", "simresults")) {
  result_file <- simulation_result_path(
    config = config,
    dir    = cache_dir,
    name   = "hybrid_heatmaps"
  )

  if (cache && !force_recompute && file.exists(result_file)) {
    return(readRDS(result_file))
  }

  sim_out <- run_simulation_experiment(
    alpha = config$alpha,
    K = config$K,
    n = config$n,
    B = config$B,
    taus = list(
      AE = config$tau_AE_values,
      ARE = config$tau_ARE_values
    ),
    metrics = c("AE", "ARE"),
    seed = config$seed
  )

  p_table <- sim_out$p_table
  index_cols <- grep("^index_", names(p_table), value = TRUE)

  plots <- lapply(config$alpha, function(alpha_i) {
    p_row <- p_table[p_table$alpha == alpha_i, , drop = FALSE]
    p_i <- as.numeric(p_row[1, index_cols, drop = FALSE])

    phat_subset <- sim_out$phat_long[
      sim_out$phat_long$alpha == alpha_i,
      c("replicate", "index", "phat"),
      drop = FALSE
    ]
    phat_subset <- phat_subset[order(phat_subset$index, phat_subset$replicate), , drop = FALSE]
    phat_mat_i <- matrix(phat_subset$phat, nrow = config$B, ncol = config$K)

    plot_hybrid_best_cutoff_heatmap(
      phat_mat = phat_mat_i,
      p = p_i,
      cutoffs = config$cutoffs,
      tau_AE_values = config$tau_AE_values,
      tau_ARE_values = config$tau_ARE_values,
      maximize = config$maximize
    ) +
      ggplot2::labs(title = paste0("Best hybrid cutoff heatmap (alpha = ", alpha_i, ")"))
  })
  names(plots) <- paste0("alpha_", config$alpha)

  result <- list(
    config = config,
    p_table = p_table,
    phat_long = sim_out$phat_long,
    heatmaps = plots
  )

  if (cache) {
    saveRDS(result, result_file)
  }

  result
}
