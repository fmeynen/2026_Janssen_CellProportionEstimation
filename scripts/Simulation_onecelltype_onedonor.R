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

# Rename function
run_single_simulation_core <- function(n_cell, p, B) {
  # Simulate based on n_cell, p, and B
  # Returns counts, p_hat, ARE, ALOR, ATE, CV
  # Removed 'c' from arguments
  # ... [Simulation logic here] ...
}

# New function to compute scopit
compute_scopit <- function(counts, c) {
  return(mean(counts >= c))  # Returns mean of counts greater than or equal to thresholds
}

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

# Fix typo
sum_stats$n_cell <- sum_stats$n_cell

# Kept behavior similar otherwise.
