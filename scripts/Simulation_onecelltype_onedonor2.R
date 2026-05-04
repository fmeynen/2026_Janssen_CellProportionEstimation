
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
library(VGAM)       # For betabinomial sampling
library(inflection) #For finding the inflection point of mAE/mARE plots

## Parameters -------------------------------------------------------------------------------------------------------
params <- list(
  n_cell = c(10, 50, 100, 1000, 10000, 10007),  # number of cells
  p = c(0.01, 0.03, 0.05, 0.1, 0.3, 0.5),           # proportions
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
  p_hat  <- counts / n_cell
  diff   <- p - p_hat
  diff   <- ifelse(diff != 0, diff, diff+1/(2*n_cell))
  #Output
  list(
    scopit = mean(counts >= c),
    AE     = abs(diff),                              #Absolute Error
    ARE    = abs(diff) / p,                          #Absolute Relative Error
    ASE    = abs(diff) / sqrt((p*(1-p))),            #Absolute Standardised error
    ALOR   = abs(log(p_hat / p)),                    #absolute log observed to expected ratio
    ATE    = abs(asin(sqrt(p_hat)) - asin(sqrt(p))), #absolute transformed error
    LAE    = log(abs(diff)),                         #log absolute error
    TSE    = asinh(sqrt(n_cell*(2*n_cell - 6)*(diff)^2)), #Transformed Square Error, VST for squared error
    CV     = sd(p_hat) / mean(p_hat)
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

plot_metric <- function(results,
                        param_df,
                        metric = c("AE", "ARE", "ASE", "ATE", "LAE", "ALOR", "TSE"),
                        ncol = NULL,
                        bins = 30,
                        fill = "steelblue",
                        alpha = 0.8,
                        fixed_scales = TRUE) {
  
  metric <- match.arg(metric)
  
  # Ensure param_df has the expected columns
  if (!"n_cell" %in% names(param_df)) stop("param_df must contain column `n_cell`.")
  if (!"p" %in% names(param_df)) stop("param_df must contain column `p`.")
  
  # Add default values for c and B if not present
  if (!"c" %in% names(param_df)) param_df$c <- 1
  if (!"B" %in% names(param_df)) param_df$B <- 1000
  
  plot_data <- purrr::pmap_dfr(param_df, function(n_cell, p, c, B) {
    
    result_name <- glue::glue("n_cell{n_cell}_p{p}_c{c}_B{B}")
    
    if (!result_name %in% names(results)) {
      warning(glue::glue("No simulation found for {result_name}"))
      return(NULL)
    }
    
    metric_values <- results[[result_name]][[metric]]
    
    if (is.null(metric_values)) {
      warning(glue::glue("Metric `{metric}` not found for {result_name}"))
      return(NULL)
    }
    
    data.frame(
      n_cell = rep(n_cell, length(metric_values)),
      p      = rep(p, length(metric_values)),
      c      = rep(c, length(metric_values)),
      B      = rep(B, length(metric_values)),
      value  = metric_values,
      label  = rep(glue::glue("n_cell = {n_cell}, p = {p}"), length(metric_values))
    )
  })
  
  if (nrow(plot_data) == 0) {
    stop("No valid simulations found for the provided parameters")
  }
  
  # Determine grid layout
  n_plots <- length(unique(plot_data$label))
  ncol <- ncol %||% min(n_plots, 3)
  
  scales <- if (fixed_scales) "fixed" else "free"
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = bins, fill = fill, color = "white", alpha = alpha) +
    ggplot2::facet_wrap(~ label, scales = scales, ncol = ncol) +
    ggplot2::labs(
      title = glue::glue("Comparison of {metric} Distributions"),
      x = metric,
      y = "Frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
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

param_grid  |>
  mutate(
    mean_AE   = map_dbl(res_sim, ~ mean(.x$AE)),
    var_AE    = map_dbl(res_sim, ~ var(.x$AE)),
    q95_AE    = map_dbl(res_sim, ~ quantile(.x$AE, probs = 0.95)),
    mean_ARE  = map_dbl(res_sim, ~ mean(.x$ARE)),
    var_ARE   = map_dbl(res_sim, ~ var(.x$ARE)),
    q95_ARE   = map_dbl(res_sim, ~ quantile(.x$ARE, probs = 0.95)),
    mean_ASE  = map_dbl(res_sim, ~ mean(.x$ASE)),
    var_ASE   = map_dbl(res_sim, ~ var(.x$ASE)),
    q95_ASE   = map_dbl(res_sim, ~ quantile(.x$ASE, probs = 0.95)),
    mean_LAE  = map_dbl(res_sim, ~ mean(.x$LAE)),
    var_LAE   = map_dbl(res_sim, ~ var(.x$LAE)),
    q95_LAE   = map_dbl(res_sim, ~ quantile(.x$LAE, probs = 0.95)),
    mean_TSE  = map_dbl(res_sim, ~ mean(.x$TSE)),
    var_TSE   = map_dbl(res_sim, ~ var(.x$TSE)),
    q95_TSE   = map_dbl(res_sim, ~ quantile(.x$TSE, probs = 0.95)),
  ) |> 
  select(-c, -B) |> 
  filter(n_cell == 10000) |>
  signif(digits = 6) |> 
  arrange(n_cell, p) -> sum_stats_mean

sum_stats_mean
write.csv(sum_stats_mean, file = "sum_stats.csv", quote = FALSE, sep = ",", row.names = FALSE)
writexl::write_xlsx(sum_stats_mean, "sum_stats.xlsx")
sum_stats_mean

# Visualization ---------------------------------------------------------------------------------------------------

plot_metric(res_sim, param_grid[param_grid$n_cell == 10000, ], metric = "AE")
plot_metric(res_sim, param_grid[param_grid$n_cell == 10000, ], metric = "ARE")
plot_metric(res_sim, param_grid[param_grid$n_cell == 10000, ], metric = "ASE")
plot_metric(res_sim, param_grid[param_grid$n_cell == 10000, ], metric = "LAE")
plot_metric(res_sim, param_grid[param_grid$n_cell == 10000, ], metric = "TSE")
plot_metric(res_sim, param_grid[param_grid$n_cell == 10000, ], metric = "ATE")


# plot mAE vs mARE ------------------------------------------------------------------------------------------------

props <- 10^seq(from = -5, to = log10(0.5), by = 0.005)

mARE <- mAE <- c()
qARE <- qAE <- c()
mLAE <- mLARE <- c()

for (p in props) {
  set.seed(p*2609)
  counts <- rbinom(n = 1000, size = 10000, prob = p)
  p_hat  <- counts/10000
  AE     <- abs(p_hat - p)
  ARE    <- AE/p
  LAE    <- log(AE)
  mAE    <- c(mAE, mean(AE))
  mARE   <- c(mARE, mean(ARE))
  mLAE   <- c(mLAE, mean(LAE))
  qAE    <- c(qAE, quantile(AE, 0.95))
  qARE   <- c(qARE, quantile(ARE, 0.95))
}

m_data <- data.frame(props, mAE, mARE)
ggplot(m_data, aes(x = mAE, y = mARE, colour = props)) +
  geom_point() +
  ggtitle("mean Absolute Error vs mean Absolute Relative Error for different proportions") +
  labs(x = "mean Absolute Error", y = "mean Absolute Relative Error")

q_data <- data.frame(props, qAE, qARE)
ggplot(q_data, aes(x = qAE, y = qARE, colour = props)) + geom_point()


# plot AE vs ARE (differing p hat) --------------------------------------------------------------------------------

p_hat         <- c(0.01, 0.03, 0.05, 0.1, 0.3, 0.5)
# p <- seq(from = 0.001, to = 0.5, by = 0.001)
p             <- 10^seq(from = -2.5, to = log10(0.5), by = 0.01)
p_error       <- expand.grid(p_hat = p_hat, p = p)
p_error |> 
  mutate(AE  = abs(p_hat- p),
         ARE = abs(p_hat- p)/p) |> 
  mutate(p_hat = as.factor(p_hat))-> p_error


ggplot(p_error, aes(x = AE, y = ARE, colour = p)) +
  facet_wrap(p_hat, scales = "free_y") +
  geom_point() +
  ggtitle("AE vs ARE for different values of p and p_hat")


# Scratchpad ------------------------------------------------------------------------------------------------------


plot(mAE, mARE)
plot(qAE, qARE, fit)
plot(fit)


kneedle <- ede(scale(mAE), scale(mARE), index = 0)
index <- kneedle[1]
props[index]

kneedle <- ede(scale(qAE), scale(qARE), index = 0)
index <- kneedle[1]
props[index]

fit <- smooth.spline(qAE, qARE)
d2 <- predict(fit, qAE, deriv = 2)
props[which.max(abs(d2$y))]