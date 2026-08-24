

``` r
# simulation_samplesize.R
#
# Iterative sample-size estimation for each alpha.
#
# Workflow:
#   1. For each alpha, run 3 pilot sample sizes (95%, 100%, 105% of n_init).
#   2. Fit glm(success ~ n, family = binomial) to observed success counts.
#   3. Invert the curve to estimate the sample size hitting the target success rate.
#   4. Repeat until convergence or max iterations; always round up.
#   5. Propagate the final estimate as the next alpha's n_init.
```

``` r
simulation_helper_files <- list.files(here::here("scripts", "simulation_layers"))
lapply(simulation_helper_files, function(f) {
  source(here::here("scripts", "simulation_layers", f))
})
```

```
## [[1]]
## [[1]]$value
## function(alpha, n_init, config) {
##   stopifnot(
##     is.numeric(alpha), length(alpha) == 1L, alpha > 0,
##     is.numeric(n_init), length(n_init) == 1L, n_init >= 1L
##   )
##   target    <- config$success_rate_target
##   tolerance <- config$sample_size_tolerance
##   max_iter  <- config$max_iterations
##   B         <- config$B
##   stopifnot(
##     is.numeric(target),    length(target)    == 1L, target > 0,    target < 1,
##     is.numeric(tolerance), length(tolerance) == 1L, tolerance >= 0,
##     is.numeric(max_iter) || is.integer(max_iter),
##     length(max_iter) == 1L, max_iter >= 1L
##   )
##   max_iter <- as.integer(max_iter)
## 
##   current_n      <- as.integer(ceiling(n_init))
##   diag_rows      <- vector("list", max_iter * 3L)
##   diag_idx       <- 0L
##   stopping_reason <- "max_iterations"
## 
##   for (iter in seq_len(max_iter)) {
##     pilot_ns <- as.integer(ceiling(c(0.95, 1.00, 1.05) * current_n))
##     pilot_ns <- pmax(pilot_ns, 1L)   # guard against n < 1
## 
##     success_counts <- integer(3L)
##     success_rates  <- numeric(3L)
## 
##     for (j in seq_along(pilot_ns)) {
##       sim_j            <- simulate_success_at_n(alpha, pilot_ns[j], config)
##       success_counts[j] <- sim_j$success_count
##       success_rates[j]  <- sim_j$success_rate
##     }
## 
##     glm_fit  <- fit_success_glm(pilot_ns, success_counts, B)
##     solved   <- tryCatch(
##       solve_sample_size_from_glm(glm_fit, target),
##       error = function(e) list(n_raw = current_n, n_rounded = current_n)
##     )
##     new_n <- as.integer(ceiling(solved$n_rounded))
## 
##     # Clamp the new estimate to prevent extreme jumps when the pilot success
##     # rates are all near 0% or all near 100%.
##     # - If mean success rate < target (need more n): cap at 2x the largest pilot.
##     # - If mean success rate >= target (need less n): floor at 0.5x the smallest pilot.
##     mean_success_rate <- mean(success_rates)
##     if (mean_success_rate < target) {
##       upper_bound <- as.integer(ceiling(2.0 * max(pilot_ns)))
##       new_n <- min(new_n, upper_bound)
##       if(is.na(new_n)){new_n <- upper_bound}
##     } else {
##       lower_bound <- as.integer(ceiling(0.5 * min(pilot_ns)))
##       new_n <- max(new_n, lower_bound)
##       if(is.na(new_n)){new_n <- lower_bound}
##     }
##     new_n <- pmax(new_n, 1L)
## 
##     coefs <- stats::coef(glm_fit)
##     n_unclamped <- as.integer(ceiling(solved$n_rounded))
##     clamped     <- (new_n != n_unclamped)
## 
##     for (j in seq_along(pilot_ns)) {
##       diag_idx <- diag_idx + 1L
##       diag_rows[[diag_idx]] <- data.frame(
##         alpha               = alpha,
##         iteration           = iter,
##         pilot_index         = j,
##         n                   = pilot_ns[j],
##         success_count       = success_counts[j],
##         success_rate        = success_rates[j],
##         target_success_rate = target,
##         glm_intercept       = coefs[[1L]],
##         glm_slope           = coefs[[2L]],
##         n_raw               = solved$n_raw,
##         n_rounded           = new_n,
##         n_clamped           = clamped,
##         stopping_reason     = NA_character_,
##         stringsAsFactors    = FALSE
##       )
##     }
## 
##     if (abs(new_n - current_n) <= tolerance) {
##       stopping_reason <- "tolerance"
##       current_n       <- new_n
##       break
##     }
##     current_n <- new_n
##   }
## 
##   diagnostics <- do.call(rbind, diag_rows[seq_len(diag_idx)])
##   diagnostics$stopping_reason[nrow(diagnostics)] <- stopping_reason
## 
##   list(
##     final_n         = current_n,
##     stopping_reason = stopping_reason,
##     iterations_used = as.integer(diag_idx / 3L),
##     diagnostics     = diagnostics
##   )
## }
## 
## [[1]]$visible
## [1] FALSE
## 
## 
## [[2]]
## [[2]]$value
## function(rep_out, alpha_i, p_max_i, B) {
##   data.frame(
##     alpha = alpha_i,
##     p_max = p_max_i,
##     replicate = rep(seq_len(B), times = ncol(rep_out$phat)),
##     index = rep(seq_len(ncol(rep_out$phat)), each = B),
##     phat = as.vector(rep_out$phat),
##     stringsAsFactors = FALSE
##   )
## }
## 
## [[2]]$visible
## [1] FALSE
## 
## 
## [[3]]
## [[3]]$value
## function(alpha, K, n, B, cutoffs,
##                                                      tau_AE, tau_ARE,
##                                                      proportion_method = "beta",
##                                                      p_max = NULL,
##                                                      model = "multinomial",
##                                                      tie_method = "random",
##                                                      maximize = c("cell", "replicate"),
##                                                      seed = NULL, ...) {
##   maximize <- match.arg(maximize)
##   sim_out <- run_simulation_experiment(
##     alpha = alpha,
##     K = K,
##     n = n,
##     B = B,
##     taus = list(AE = tau_AE, ARE = tau_ARE),
##     metrics = c("AE", "ARE"),
##     proportion_method = proportion_method,
##     p_max = p_max,
##     model = model,
##     tie_method = tie_method,
##     seed = seed,
##     ...
##   )
## 
##   p_table <- sim_out$p_table
##   index_cols <- grep("^index_", names(p_table), value = TRUE)
##   cutoff_curves_list <- vector("list", nrow(p_table))
##   best_summary_list <- vector("list", nrow(p_table))
## 
##   for (i in seq_len(nrow(p_table))) {
##     alpha_i <- p_table$alpha[[i]]
##     p_max_i <- if ("p_max" %in% names(p_table)) p_table$p_max[[i]] else NA_real_
##     p_i <- as.numeric(p_table[i, index_cols, drop = FALSE])
## 
##     same_p_max <- if (is.na(p_max_i)) {
##       is.na(sim_out$phat_long$p_max)
##     } else {
##       sim_out$phat_long$p_max == p_max_i
##     }
##     phat_subset <- sim_out$phat_long[
##       sim_out$phat_long$alpha == alpha_i & same_p_max,
##       c("replicate", "index", "phat"),
##       drop = FALSE
##     ]
##     phat_subset <- phat_subset[order(phat_subset$index, phat_subset$replicate), , drop = FALSE]
##     phat_mat_i <- matrix(phat_subset$phat, nrow = B, ncol = length(p_i))
## 
##     hybrid_i <- run_hybrid_cutoff_analysis(
##       phat_mat = phat_mat_i,
##       p = p_i,
##       cutoffs = cutoffs,
##       tau_AE = tau_AE,
##       tau_ARE = tau_ARE,
##       maximize = maximize,
##       alpha = alpha_i,
##       p_max = p_max_i
##     )
## 
##     curve_i <- hybrid_i$cutoff_curve
##     curve_i$alpha <- alpha_i
##     curve_i$p_max <- p_max_i
##     curve_i$tau_AE <- tau_AE
##     curve_i$tau_ARE <- tau_ARE
##     cutoff_curves_list[[i]] <- curve_i[, c(
##       "alpha", "p_max", "cutoff", "tau_AE", "tau_ARE",
##       "success_rate_cell", "success_rate_replicate",
##       "prop_using_AE", "prop_using_ARE"
##     )]
## 
##     best_summary_list[[i]] <- hybrid_i$best_summary
##   }
## 
##   list(
##     inputs = list(
##       alpha = alpha, K = K, n = n, B = B, cutoffs = cutoffs,
##       tau_AE = tau_AE, tau_ARE = tau_ARE,
##       proportion_method = proportion_method, p_max = p_max,
##       model = model, tie_method = tie_method, maximize = maximize, seed = seed
##     ),
##     p_table = sim_out$p_table,
##     phat_long = sim_out$phat_long,
##     cutoff_curves = do.call(rbind, cutoff_curves_list),
##     best_cutoff_summary = do.call(rbind, best_summary_list)
##   )
## }
## 
## [[3]]$visible
## [1] FALSE
## 
## 
## [[4]]
## [[4]]$value
## function(alpha, n, config) {
##   p <- generate_proportions(
##     alpha  = alpha,
##     K      = config$K,
##     method = config$proportion_method
##   )
##   rep_out <- run_replicates(
##     p          = p,
##     n          = n,
##     B          = config$B,
##     metrics    = config$metrics,
##     model      = config$model,
##     tie_method = config$tie_method,
##     seed       = config$seed
##   )
##   max_errors <- rep_out$max_errors
##   metrics    <- config$metrics
##   taus       <- config$taus
## 
##   # Build a B-length logical success vector: replicate passes iff every metric
##   # with a threshold has its max error <= that threshold.
##   success <- rep(TRUE, config$B)
##   for (m in metrics) {
##     if (is.null(taus[[m]])) {
##       warning(sprintf(
##         "simulate_success_at_n: metric '%s' has no threshold in config$taus; it will not contribute to the success criterion.",
##         m
##       ))
##     } else if (length(taus[[m]]) == 1L) {
##       success <- success & (max_errors[, m] <= taus[[m]])
##     }
##   }
## 
##   list(
##     success       = success,
##     success_count = sum(success),
##     success_rate  = mean(success),
##     rep_out       = rep_out
##   )
## }
## 
## [[4]]$visible
## [1] FALSE
## 
## 
## [[5]]
## [[5]]$value
## function(config, dir, name) {
##   dir.create(dir, recursive = TRUE, showWarnings = FALSE)
##   file.path(dir, paste0(name, "_", hash_config(config), ".rds"))
## }
## 
## [[5]]$visible
## [1] FALSE
## 
## 
## [[6]]
## [[6]]$value
## function(non_max, alpha, K, p_max) {
##   warning(
##     sprintf(
##       paste(
##         "Impossible fixed_max_beta combination for alpha=%s, K=%s, p_max=%s:",
##         "after scaling the Beta-shaped remainder to sum to 1 - p_max,",
##         "at least one non-max component is >= p_max (max non-max = %s)."
##       ),
##       format(alpha, trim = TRUE),
##       K,
##       format(p_max, trim = TRUE),
##       format(max(non_max), trim = TRUE)
##     ),
##     call. = FALSE
##   )
##   stop(structure(
##     list(message = fixed_max_beta_impossible_error, call = NULL),
##     class = c("impossible_fixed_max_error", "error", "condition")
##   ))
## }
## 
## [[6]]$visible
## [1] FALSE
## 
## 
## [[7]]
## [[7]]$value
## function(phat_mat, p, cutoffs,
##                                             tau_AE_values, tau_ARE_values,
##                                             maximize = c("cell", "replicate"),
##                                             tie_break = c("smallest", "largest", "median"),
##                                             label_digits = 3L) {
##   maximize <- match.arg(maximize)
##   tie_break <- match.arg(tie_break)
## 
##   stopifnot(
##     is.matrix(phat_mat),
##     is.numeric(p),
##     length(p) == ncol(phat_mat),
##     is.numeric(cutoffs),
##     length(cutoffs) >= 1L,
##     all(is.finite(cutoffs)),
##     is.numeric(tau_AE_values),
##     length(tau_AE_values) >= 1L,
##     all(is.finite(tau_AE_values)),
##     is.numeric(tau_ARE_values),
##     length(tau_ARE_values) >= 1L,
##     all(is.finite(tau_ARE_values)),
##     is.numeric(label_digits),
##     length(label_digits) == 1L,
##     is.finite(label_digits),
##     label_digits %% 1 == 0,
##     label_digits >= 0
##   )
## 
##   threshold_grid <- expand.grid(
##     tau_AE = tau_AE_values,
##     tau_ARE = tau_ARE_values,
##     KEEP.OUT.ATTRS = FALSE,
##     stringsAsFactors = FALSE
##   )
## 
##   rows <- lapply(seq_len(nrow(threshold_grid)), function(i) {
##     tau_AE_i <- threshold_grid$tau_AE[[i]]
##     tau_ARE_i <- threshold_grid$tau_ARE[[i]]
##     curve_i <- sweep_hybrid_cutoffs_cell_level(
##       phat_mat = phat_mat,
##       p = p,
##       cutoffs = cutoffs,
##       tau_AE = tau_AE_i,
##       tau_ARE = tau_ARE_i
##     )
##     best_i <- find_best_hybrid_cutoff(
##       curve_df = curve_i,
##       maximize = maximize,
##       tie_break = tie_break
##     )
##     data.frame(
##       tau_AE = tau_AE_i,
##       tau_ARE = tau_ARE_i,
##       best_cutoff = best_i$best_cutoff,
##       stringsAsFactors = FALSE
##     )
##   })
## 
##   heatmap_df <- do.call(rbind, rows)
##   heatmap_df$tau_AE <- factor(heatmap_df$tau_AE, levels = sort(tau_AE_values, decreasing = TRUE))
##   heatmap_df$tau_ARE <- factor(heatmap_df$tau_ARE, levels = sort(tau_ARE_values))
##   heatmap_df$label <- formatC(heatmap_df$best_cutoff, digits = label_digits, format = "f")
## 
##   ggplot2::ggplot(
##     heatmap_df,
##     ggplot2::aes(x = tau_ARE, y = tau_AE, fill = best_cutoff)
##   ) +
##     ggplot2::geom_tile(color = "white", linewidth = 0.4) +
##     ggplot2::geom_text(ggplot2::aes(label = label), size = 3) +
##     ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08306B") +
##     ggplot2::labs(
##       x = "ARE threshold (tau_ARE)",
##       y = "AE threshold (tau_AE)",
##       fill = "Best cutoff",
##       title = "Best hybrid cutoff heatmap"
##     ) +
##     ggplot2::theme_bw()
## }
## 
## [[7]]$visible
## [1] FALSE
```

``` r
library(ggplot2)

# Setup config file -----------------------------------------------------------------------------------------------

simulation_sample_size_defaults <- function(n_init = 200000) {
  tau_AE <- 0.002
  tau_ARE <- 0.05
  list(
    alpha                = seq(from = 2, to = 5, by = 0.05),
    K                    = 10L,
    n                    = n_init,
    B                    = 500L,
    taus                 = list(AE = tau_AE, ARE = tau_ARE),
    metrics              = c("AE", "ARE"),
    model                = "multinomial",
    tie_method           = "random",
    proportion_method    = "beta",
    seed                 = 260925L,
    success_rate_target  = 0.95,
    sample_size_tolerance = 100L,
    max_iterations       = 20L
  )
}

config <- simulation_sample_size_defaults()

run_simulation_samplesize <- function(config = simulation_sample_size_defaults(),
                                      cache = TRUE,
                                      force_recompute = FALSE,
                                      cache_dir = here::here("results", "simresults")) {
  result_file <- simulation_result_path(
    config = config,
    dir    = cache_dir,
    name   = "sample_size"
  )
  if (cache && !force_recompute && file.exists(result_file)) {
    return(readRDS(result_file))
  }

  alphas        <- config$alpha
  n_samples     <- numeric(length(alphas))
  diag_list     <- vector("list", length(alphas))
  n_init        <- config$n

  for (a in seq_along(alphas)) {
    iter_result        <- iterate_sample_size_for_alpha(alphas[a], n_init, config)
    n_samples[a]       <- iter_result$final_n
    diag_list[[a]]     <- iter_result$diagnostics
    n_init             <- iter_result$final_n
  }

  result <- list(
    sample_size  = as.integer(n_samples),
    diagnostics  = do.call(rbind, diag_list)
  )

  if (cache) {
    saveRDS(result, result_file)
  }
  result
}



res <- run_simulation_samplesize()
res$sample_size
```

```
##  [1]    217237    276935    263890    259694    285077         4        10        22        48       102         4        10        22        48
## [15]       102         4        10        22        48       102   1935103   2339113   2707466   2992213   3564779   4036635   4467108   5295578
## [29]   6092764   7114089   7816152   9142184  10574225  12074466  13885018  15819065  18397874  20770297  23542756  27294534  31101753  35906961
## [43]  40696203  46931388  54264861  61791707  70283847  81787611  93167542 108996453 123998516 141018845 165998223 188668194 218539394 246777261
## [57] 279938651 330093035 377737075 438885388 497584720
```

``` r
plot(config$alpha, log10(res$sample_size))
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

``` r
ggplot(data.frame(alpha = config$alpha, sample_size = res$sample_size),
       aes(x = alpha, y = sample_size)) +
  geom_line() +
  labs(x = "Alpha", y = "log10(Sample Size)") +
  theme_minimal() +
  scale_y_log10()
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-2.png)

``` r
ggplot(data.frame(alpha = config$alpha, sample_size = res$sample_size),
       aes(x = alpha, y = sample_size)) +
  geom_line() +
  labs(x = "Alpha", y = "Sample Size") +
  theme_minimal()
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-3.png)

