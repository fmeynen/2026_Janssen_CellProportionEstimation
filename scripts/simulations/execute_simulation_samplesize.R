# Set up a config file
# Go through all alpha's
## For the first alpha, estimate sample size
## Use this as an init for the second sample size

# Sample size calculation for a single alpha:
# Take 3 initial sample sizes with small repetition (500)
# fit logistic regression glm(success ~ n, family=binomial)
# estimate sample size needed for success rate from logistic regression
# again 500 repetitions at this sample size
# repeat until sample size estimate is within 10 of most recent sample size estimate.



# Setup config file -----------------------------------------------------------------------------------------------

simulation_sample_size_defaults <- function(n_init=200000) {
  tau_AE <- 0.02
  tau_ARE <- 0.5
  list(
    alpha = c(2, 2.5, 3, 4, 5),
    K = 10L,
    n = n_init,
    B = 500L,
    taus = list(AE = tau_AE, ARE = tau_ARE),
    metrics = c("AE", "ARE"),
    model = "multinomial",
    tie_method = "random",
    proportion_method = "beta",
    seed = 260926L
  )
}

run_simulation_samplesize <- function(config = simulation_sample_size_defaults(),
                                      cache = TRUE,
                                      force_recompute = FALSE,
                                      cache_dir = here::here("results", "simresults")){
  result_file <- simulation_result_path(
    config = config,
    dir    = cache_dir,
    name   = "sample_size"
  )
  if (cache && !force_recompute && file.exists(result_file)) {
    return(readRDS(result_file))
  }
  alphas      <- config$alpha
  sample_size <- alphas
  n_init      <- config$n_init
  for (a in 1:length(alphas)) {
    sample_size[a] <- calculate_sample_size(alphas[a], n_init, config)
    n_init <- sample_size[a]
  }
  sample_size
}

calculate_sample_size <- function(alpha, n_init, config){
  
}