# Updated Simulation Code

# Refactored code as per the requirements:
# 1. 'c' is now a vector of thresholds and not part of the simulation grid.

params <- list(
  n_cell = seq(100, 1000, by = 100),
  p = seq(0.1, 1, by = 0.1),
  B = 10,
  c = c(1)  # Vector of thresholds
)

# Create sim_grid without 'c'
sim_grid <- expand.grid(n_cell = params$n_cell, p = params$p, B = params$B)


## Helper functions ------------------------------------------------------------------------------------------------
### Simulation Functions ----
run_single_simulation <- function(n_cell, p,  B, c, seed) {
  set.seed(seed)
  #Simulation
  counts <- rbinom(B, n_cell, p)
  p_hat <- counts / n_cell
  #Output
  list(
    counts = counts,
    scopit = mean(counts >= c),
    ARE = abs(p_hat - p) / p,
    ALOR = abs(log(p_hat / p)),
    ATE = abs(asin(sqrt(p_hat)) - asin(sqrt(p))),
    CV = sd(p_hat) / mean(p_hat)
  )
}


run_all_simulations <- function(param_grid) {
  results <- pmap(param_grid, function(n_cell, p,c, B) {
    seed <- 2609 * which(
      param_grid$n_cell == n_cell & param_grid$p == p &
        param_grid$c == c & param_grid$B == B
    ) #the equivalent of choosing i in a for loop over the grid.
    run_single_simulation(n_cell, p, c, B, seed)
  })
  names(results) <- glue::glue(
    "n_cell{param_grid$n_cell}_p{param_grid$p}_c{param_grid$c}_B{param_grid$B}"
  )
  results


# Modified run_all_simulations
run_all_simulations <- function() {
  results <- list()

  for (i in 1:nrow(sim_grid)) {
    params <- sim_grid[i, ]
    set.seed(i)  # Generate seed using row index
    # Store results named like n_cell{n_cell}_p{p}_B{B}
    result_name <- paste0("n_cell", params$n_cell, "_p", params$p, "_B", params$B)
    #... [run simulations and save results] ...
  }
  # Compute scopit for each threshold
  # Produce long-format summary table sum_stats
}

# Update plotting functions accordingly
plot_are <- function(results, c = NULL) {
  result_name <- paste0("n_cell", results$n_cell, "_p", results$p, "_B", results$B)
  # Optionally accept 'c' only for labeling
}

plot_alor <- function(results, param_df, ncol = NULL) {
  
  # param_df is a dataframe with columns n, p, and optionally c, B
  # Example: data.frame(n_cell = c(100, 1000), p = c(0.05, 0.05))
  
  # Add default values for c and B if not present
  if (!"c" %in% names(param_df)) param_df$c <- 1
  if (!"B" %in% names(param_df)) param_df$B <- 1000
  
  alor_data <- pmap_dfr(param_df, function(n_cell, p, c, B) {
    
    result_name <- glue::glue("n_cell{n_cell}_p{p}_c{c}_B{B}")
    
    if (!result_name %in% names(results)) {
      warning(glue::glue("No simulation found for {result_name}"))
      return(NULL)
    }
    
    alor_values <- results[[result_name]]$ALOR
    data.frame(
      n_cell = rep(n_cell, length(alor_values)),
      p = rep(p, length(alor_values)),
      c = rep(c, length(alor_values)),
      B = rep(B, length(alor_values)),
      ALOR = alor_values,
      label = rep(glue::glue("n_cell = {n_cell}, p = {p}"), length(alor_values))
    )
  })
  if (nrow(alor_data) == 0) {
    stop("No valid simulations found for the provided parameters")
  }
  
  # Determine grid layout
  n_plots <- length(unique(alor_data$label))
  ncol <- ncol %||% min(n_plots, 3)
  
  ggplot(alor_data, aes(x = ALOR)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
    facet_wrap(~ label, scales = "fixed", ncol = ncol) +
    labs(
      title = "Comparison of Absolute Log Error Distributions",
      x = "Absolute Log Error (ALOR)",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# Simulation ------------------------------------------------------------------------------------------------------

res_sim <- run_all_simulations(param_grid)

# Summary stats ---------------------------------------------------------------------------------------------------

param_grid  |>
  mutate(
    scopit       = map_dbl(res_sim, "scopit"),
    mean_ARE     = map_dbl(res_sim, ~ mean(.x$ARE)),
    quant90_ARE  = map_dbl(res_sim, ~quantile(.x$ARE, probs = c(0.9))),
    quant95_ARE  = map_dbl(res_sim, ~quantile(.x$ARE, probs = c(0.95))),
    mean_ATE     = map_dbl(res_sim, ~ mean(.x$ATE)),
    quant90_ATE  = map_dbl(res_sim, ~quantile(.x$ATE, probs = c(0.9))),
    quant95_ATE  = map_dbl(res_sim, ~quantile(.x$ATE, probs = c(0.95))),
    mean_ALOR    = map_dbl(res_sim, ~ mean(.x$ALOR)),
    quant90_ALOR = map_dbl(res_sim, ~quantile(.x$ALOR, probs = c(0.9))),
    quant95_ALOR = map_dbl(res_sim, ~quantile(.x$ALOR, probs = c(0.95))),
    CV           = map_dbl(res_sim, "CV")
  ) -> sum_stats

# Visualization ---------------------------------------------------------------------------------------------------

vis_df <- data.frame(
  n_cell = c(10, 100, 1000, 10000),
  p = rep(0.05, 4)
)

plot_are(res_sim, param_grid[param_grid$n_cell == 1000, ])
plot_alor(res_sim, param_grid[param_grid$n_cell == 1000, ])
plot_ate(res_sim, param_grid[param_grid$n_cell == 1000, ])
sum_stats[sum_stats$n_cel == 1000, ]


# Scratchpad ------------------------------------------------------------------------------------------------------

p <- seq(1e-5,1-1e-5,1e-4)
a <- asin(sqrt(p))
l <- log(p/(1-p))

plot(p,a)
plot(p,l)


# Fix typo
sum_stats$n_cell <- sum_stats$n_cell

# Kept behavior similar otherwise.
