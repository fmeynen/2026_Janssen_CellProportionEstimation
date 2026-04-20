
# Idea ------------------------------------------------------------------------------------------------------------

#Start small: one cell type in one donor (later: increase cell types and/or increase number of donors)
#Use binomial distribution to model cell type (later: beta binomial or more difficult frameworks)
#Go for precision instead of hypothesis testing (no hypothesis testing possible in one cell type on donor)
#Compare based on:
  #* Scopit-like where c = 1 (probability of at least observing 1 count)
  #* ARE: absolute relative error
  #* Coefficient of variation (simple since easy formula for binomial distribution)

#Idea of sampling:
  #* Choose sampling parameters
  #* Construct parameter grid
  #* Choose number of samples B
  #* For each combination of sampling parameters (=row in parameter grid) sample B times
  #* Calculate the comparison statistics

#Presentation
  #* Present summary stats
  #* Graphs of summary stats
  #* Distribution of ARE and ALOR

# Setup -----------------------------------------------------------------------------------------------------------
rm(list = ls())


## Libraries -------------------------------------------------------------------------------------------------------
library(tidyverse)  # For cleaner data manipulation and plotting

## Parameters -------------------------------------------------------------------------------------------------------
params <- list(
  n_cell = c(10, 50, 100, 1000, 10000),  # number of cells
  p = c(0.01, 0.05, 0.1, 0.3, 0.5),           # proportions
  c = 1,                             # Scopit limit
  B = 10000                           # repetitions
)

param_grid <- expand.grid(
  n_cell = params$n_cell,
  p = params$p,
  c = params$c,
  B = params$B
)


## Helper functions ------------------------------------------------------------------------------------------------
### Simulation Functions ----
run_single_simulation <- function(n_cell, p, c, B, seed) {
  set.seed(seed)
  #Simulation
  counts <- rbinom(B, n_cell, p)
  p_hat <- counts / n_cell
  #Output
  list(
    scopit = mean(counts >= c),
    ARE = abs(p_hat - p) / p,
    ALOR = abs(log(p_hat / p)),
    ATE = abs(asin(sqrt(p_hat)) - asin(sqrt(p))),
    CV = sd(p_hat) / mean(p_hat)
  )
}

run_all_simulations <- function(param_grid) {
  results <- pmap(param_grid, function(n_cell, p, c, B) {
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
}

### Visualization Functions -----------------------------------------------------------------------------------------

plot_are <- function(results, param_df, ncol = NULL) {
  # param_df is a dataframe with columns n, p, and optionally c, B
  
  # Add default values for c and B if not present
  if (!"c" %in% names(param_df)) param_df$c <- 1
  if (!"B" %in% names(param_df)) param_df$B <- 1000
  
  are_data <- pmap_dfr(param_df, function(n_cell, p, c, B) {
    
    result_name <- glue::glue("n_cell{n_cell}_p{p}_c{c}_B{B}")
    
    if (!result_name %in% names(results)) {
      warning(glue::glue("No simulation found for {result_name}"))
      return(NULL)
    }
    
    are_values <- results[[result_name]]$ARE
    
    data.frame(
      n_cell = rep(n_cell, length(are_values)),
      p = rep(p, length(are_values)),
      c = rep(c, length(are_values)),
      B = rep(B, length(are_values)),
      ARE = are_values,
      label = rep(glue::glue("n_cell = {n_cell}, p = {p}"), length(are_values))
    )
  })
  if (nrow(are_data) == 0) {
    stop("No valid simulations found for the provided parameters")
  }
  
  # Determine grid layout
  n_plots <- length(unique(are_data$label))
  ncol <- ncol %||% min(n_plots, 3)
  
  ggplot(are_data, aes(x = ARE)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
    facet_wrap(~ label, scales = "fixed", ncol = ncol) +
    labs(
      title = "Comparison of Absolute Relative Error Distributions",
      x = "Absolute Relative Error (ARE)",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

plot_ate <- function(results, param_df, ncol = NULL) {
  
  # param_df is a dataframe with columns n, p, and optionally c, B
  # Add default values for c and B if not present
  if (!"c" %in% names(param_df)) param_df$c <- 1
  if (!"B" %in% names(param_df)) param_df$B <- 1000
  
  ate_data <- pmap_dfr(param_df, function(n_cell, p, c, B) {
    
    result_name <- glue::glue("n_cell{n_cell}_p{p}_c{c}_B{B}")
    
    if (!result_name %in% names(results)) {
      warning(glue::glue("No simulation found for {result_name}"))
      return(NULL)
    }
    
    ate_values <- results[[result_name]]$ATE
    data.frame(
      n_cell = rep(n_cell, length(ate_values)),
      p = rep(p, length(ate_values)),
      c = rep(c, length(ate_values)),
      B = rep(B, length(ate_values)),
      ATE = ate_values,
      label = rep(glue::glue("n_cell = {n_cell}, p = {p}"), length(ate_values))
    )
  })
  if (nrow(ate_data) == 0) {
    stop("No valid simulations found for the provided parameters")
  }
  
  # Determine grid layout
  n_plots <- length(unique(ate_data$label))
  ncol <- ncol %||% min(n_plots, 3)
  
  ggplot(ate_data, aes(x = ATE)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
    facet_wrap(~ label, scales = "fixed", ncol = ncol) +
    labs(
      title = "Comparison of Absolute Transformed Error Distributions",
      x = "Absolute Transformed Error (ATE)",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
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
    scopit      = map_dbl(res_sim, "scopit"),
    mean_ARE    = map_dbl(res_sim, ~ mean(.x$ARE)),
    mean_ATE    = map_dbl(res_sim, ~ mean(.x$ATE)),
    mean_ALOR   = map_dbl(res_sim, ~ mean(.x$ALOR)),
    CV          = map_dbl(res_sim, "CV")
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
