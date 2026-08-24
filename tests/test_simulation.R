# tests/test_simulation.R
#
# Base-R test script for the layered simulation modules under R/.
#
# Run from the repository root:
#   Rscript tests/test_simulation.R
#
# Each test calls `stop()` with a descriptive message on failure so that CI
# (or a local developer) can immediately identify what broke.
# ---------------------------------------------------------------------------

source(file.path("R", "validation_utils.R"))
source(file.path("R", "calculation.R"))
source(file.path("R", "extraction.R"))
source(file.path("R", "simulation.R"))
source(file.path("R", "orchestration.R"))
source(file.path("R", "visualisation.R"))

pass <- function(name) cat(sprintf("  PASS  %s\n", name))

# ---- 1. generate_proportions_beta ------------------------------------------
cat("1. generate_proportions_beta\n")

# p sums to 1
p <- generate_proportions_beta(alpha = 2, K = 10)
if (abs(sum(p) - 1) > 1e-12) stop("p does not sum to 1")
pass("sum(p) == 1")

# length K
if (length(p) != 10L) stop("length(p) != 10")
pass("length(p) == K")

# all(p > 0) — strictly positive
if (any(p <= 0)) stop("p contains zero or negative values")
pass("all(p > 0)")

# deterministic given alpha
p2 <- generate_proportions_beta(alpha = 2, K = 10)
if (!identical(p, p2)) stop("generate_proportions_beta is not deterministic")
pass("deterministic given alpha")

# different alpha gives different p
p3 <- generate_proportions_beta(alpha = 5, K = 10)
if (identical(p, p3)) stop("different alpha should yield different p")
pass("different alpha -> different p")

# ---- 2. generate_proportions dispatcher ------------------------------------
cat("\n2. generate_proportions dispatcher\n")

p_dispatch <- generate_proportions(alpha = 2, K = 10, method = "beta")
if (!identical(p, p_dispatch)) stop("beta dispatcher route should match generate_proportions_beta")
pass("beta dispatcher route works and is backward-compatible")

fixed_K <- 10L
err_missing_pmax <- tryCatch(
  generate_proportions(alpha = 2, K = fixed_K, method = "fixed_max_beta"),
  error = function(e) e$message
)
if (!grepl("p_max must be provided", err_missing_pmax)) {
  stop("fixed_max_beta should require p_max")
}
pass("fixed_max_beta requires p_max")

p_fixed <- generate_proportions(alpha = 2, K = fixed_K, method = "fixed_max_beta", p_max = 0.4)
if (abs(sum(p_fixed) - 1) > 1e-12) stop("fixed_max_beta proportions do not sum to 1")
pass("fixed_max_beta proportions sum to 1")

if (length(p_fixed) != fixed_K) stop("fixed_max_beta returned wrong length")
pass("fixed_max_beta returns length K")

if (abs(p_fixed[fixed_K] - 0.4) > 1e-12) stop("fixed_max_beta should place p_max at highest index")
pass("fixed_max_beta fixes largest proportion at highest index")

if (which.max(p_fixed) != fixed_K || sum(p_fixed == max(p_fixed)) != 1L) {
  stop("fixed_max_beta should produce a strictly unique maximum at the highest index")
}
pass("fixed_max_beta maximum is strictly unique")

p_fixed_multi <- generate_proportions_fixed_max_beta(alpha = 2, K = fixed_K, p_max = c(0.3, 0.4))
if (!is.matrix(p_fixed_multi) || nrow(p_fixed_multi) != 2L || ncol(p_fixed_multi) != fixed_K) {
  stop("fixed_max_beta should return a matrix for vector p_max")
}
pass("fixed_max_beta returns a matrix for vector p_max")
if (any(abs(rowSums(p_fixed_multi) - 1) > 1e-12)) {
  stop("fixed_max_beta vector p_max rows should sum to 1")
}
pass("fixed_max_beta vector p_max rows sum to 1")
if (any(abs(p_fixed_multi[, fixed_K] - c(0.3, 0.4)) > 1e-12)) {
  stop("fixed_max_beta vector p_max should place each p_max at highest index")
}
if (!all(apply(p_fixed_multi, 1, function(row) which.max(row) == fixed_K && sum(row == max(row)) == 1L))) {
  stop("fixed_max_beta vector p_max should preserve strict unique maximum at highest index")
}
pass("fixed_max_beta vector p_max preserves strict unique maximum")

warn_fixed_max <- NULL
err_fixed_max <- tryCatch(
  withCallingHandlers(
    generate_proportions(alpha = 2, K = 2, method = "fixed_max_beta", p_max = 0.4),
    warning = function(w) {
      warn_fixed_max <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e) e$message
)
if (is.null(warn_fixed_max) || !grepl("Impossible fixed_max_beta combination", warn_fixed_max)) {
  stop("fixed_max_beta should warn on impossible alpha/K/p_max combinations")
}
pass("fixed_max_beta warns on impossible combinations")

if (!grepl("not strictly unique", err_fixed_max)) {
  stop("fixed_max_beta should fail after warning for impossible combinations")
}
pass("fixed_max_beta fails after warning on impossible combinations")

# ---- 3. simulate_counts_multinomial ----------------------------------------
cat("\n3. simulate_counts_multinomial\n")

set.seed(1L)
y <- simulate_counts_multinomial(p, n = 200L)

# correct length
if (length(y) != 10L) stop("y has wrong length")
pass("length(y) == K")

# sums to n
if (sum(y) != 200L) stop("sum(y) != n")
pass("sum(y) == n")

# nonnegative integers
if (any(y < 0L)) stop("y contains negative counts")
if (!is.integer(y)) stop("y is not integer")
pass("nonneg integer counts")

# ---- 4. simulate_counts dispatcher -----------------------------------------
cat("\n4. simulate_counts dispatcher\n")

set.seed(1L)
y2 <- simulate_counts(p, n = 200L, model = "multinomial")
if (length(y2) != 10L || sum(y2) != 200L) stop("dispatcher output wrong")
pass("multinomial route works")

# unimplemented models throw errors
err_dm <- tryCatch(
  simulate_counts(p, n = 200L, model = "dirichlet_multinomial"),
  error = function(e) e$message
)
if (!grepl("not yet implemented", err_dm)) stop("dirichlet_multinomial should error")
pass("dirichlet_multinomial errors with clear message")

err_ln <- tryCatch(
  simulate_counts(p, n = 200L, model = "logistic_normal_multinomial"),
  error = function(e) e$message
)
if (!grepl("not yet implemented", err_ln)) stop("logistic_normal_multinomial should error")
pass("logistic_normal_multinomial errors with clear message")

# ---- 5. compute_errors ------------------------------------------------------
cat("\n5. compute_errors\n")

# known example
p_ex   <- c(0.5, 0.3, 0.2)
phat_ex <- c(0.6, 0.25, 0.15)
err <- compute_errors(phat_ex, p_ex, metrics = c("AE", "ARE"))

expected_ae  <- abs(phat_ex - p_ex)          # 0.10, 0.05, 0.05
expected_are <- abs(phat_ex - p_ex) / p_ex   # 0.20, 0.1667, 0.25

if (!isTRUE(all.equal(err$AE,  expected_ae)))  stop("AE mismatch")
pass("AE correct")
if (!isTRUE(all.equal(err$ARE, expected_are))) stop("ARE mismatch")
pass("ARE correct")

# validate_proportions rejects zero proportions
err_zero_prop <- tryCatch(
  validate_proportions(c(0.5, 0.5, 0.0)),
  error = function(e) e$message
)
if (!is.character(err_zero_prop)) stop("validate_proportions should reject zero proportions")
pass("validate_proportions rejects zero proportions")

# ---- 6. max_error_summary ---------------------------------------------------
cat("\n6. max_error_summary\n")

ev <- c(0.1, 0.5, 0.5, 0.3)

# first
s_first <- max_error_summary(ev, tie_method = "first")
if (s_first$max_error_value != 0.5) stop("wrong max value (first)")
if (s_first$argmax_index != 2L)     stop("tie_method='first' should return index 2")
pass("tie_method='first' returns smallest tied index")

# last
s_last <- max_error_summary(ev, tie_method = "last")
if (s_last$max_error_value != 0.5) stop("wrong max value (last)")
if (s_last$argmax_index != 3L)     stop("tie_method='last' should return index 3")
pass("tie_method='last' returns largest tied index")

# random: over many draws, both tied indices should be selected
set.seed(99L)
random_indices <- replicate(200L, max_error_summary(ev, tie_method = "random")$argmax_index)
if (!all(c(2L, 3L) %in% random_indices)) stop("tie_method='random' should pick both tied indices")
pass("tie_method='random' samples among tied indices")

# unique max (no tie)
ev2 <- c(0.1, 0.9, 0.4)
s_unique <- max_error_summary(ev2, tie_method = "first")
if (s_unique$argmax_index != 2L) stop("unique max should return index 2")
pass("unique max handled correctly")

# ---- 7. evaluate_thresholds -------------------------------------------------
cat("\n7. evaluate_thresholds\n")

# construct toy max_errors matrix
me <- matrix(c(0.05, 0.15, 0.25,    # AE column
               0.10, 0.20, 0.30),   # ARE column
             nrow = 3L, ncol = 2L,
             dimnames = list(NULL, c("AE", "ARE")))

curves <- evaluate_thresholds(me, taus = c(0.10, 0.20))

# AE: 1/3 <= 0.10, 2/3 <= 0.20
ae_10 <- curves$success_rate[curves$metric == "AE" & curves$tau == 0.10]
ae_20 <- curves$success_rate[curves$metric == "AE" & curves$tau == 0.20]
if (!isTRUE(all.equal(ae_10, 1/3))) stop("AE success_rate at tau=0.10 wrong")
if (!isTRUE(all.equal(ae_20, 2/3))) stop("AE success_rate at tau=0.20 wrong")
pass("AE success rates correct")

# ARE: 1/3 <= 0.10, 2/3 <= 0.20
are_10 <- curves$success_rate[curves$metric == "ARE" & curves$tau == 0.10]
are_20 <- curves$success_rate[curves$metric == "ARE" & curves$tau == 0.20]
if (!isTRUE(all.equal(are_10, 1/3))) stop("ARE success_rate at tau=0.10 wrong")
if (!isTRUE(all.equal(are_20, 2/3))) stop("ARE success_rate at tau=0.20 wrong")
pass("ARE success rates correct")

# correct columns
expected_cols <- c("metric", "tau", "success_rate")
if (!all(expected_cols %in% names(curves))) stop("curves missing expected columns")
pass("curves has correct columns")

# named-list taus: AE and ARE use different threshold grids
curves_list <- evaluate_thresholds(
  me,
  taus = list(AE = c(0.10, 0.20), ARE = c(0.15, 0.25))
)

ae_list_10 <- curves_list$success_rate[curves_list$metric == "AE" & curves_list$tau == 0.10]
ae_list_20 <- curves_list$success_rate[curves_list$metric == "AE" & curves_list$tau == 0.20]
if (!isTRUE(all.equal(ae_list_10, 1/3))) stop("named-list: AE success_rate at tau=0.10 wrong")
if (!isTRUE(all.equal(ae_list_20, 2/3))) stop("named-list: AE success_rate at tau=0.20 wrong")
pass("named-list taus: AE success rates correct")

are_list_15 <- curves_list$success_rate[curves_list$metric == "ARE" & curves_list$tau == 0.15]
are_list_25 <- curves_list$success_rate[curves_list$metric == "ARE" & curves_list$tau == 0.25]
if (!isTRUE(all.equal(are_list_15, 1/3))) stop("named-list: ARE success_rate at tau=0.15 wrong")
if (!isTRUE(all.equal(are_list_25, 2/3))) stop("named-list: ARE success_rate at tau=0.25 wrong")
pass("named-list taus: ARE success rates correct")

# named-list with unequal grid lengths produces correct row counts
# 2 thresholds for AE + 2 thresholds for ARE = 4 rows total
if (nrow(curves_list) != 4L) stop("named-list curves: wrong row count")
pass("named-list taus: row count correct with different grid lengths")

# missing metric in named-list should error
err_missing <- tryCatch(
  evaluate_thresholds(me, taus = list(AE = c(0.10))),
  error = function(e) e$message
)
if (!grepl("ARE", err_missing)) stop("missing metric in taus list should error with metric name")
pass("named-list taus: missing metric produces informative error")

# ---- 8. run_simulation_experiment (smoke test) ------------------------------
cat("\n8. run_simulation_experiment smoke test\n")

set.seed(7L)
res <- run_simulation_experiment(
  alpha = 2, K = 10L, n = 100L, B = 50L,
  taus = c(0.05, 0.10, 0.20),
  metrics = c("AE", "ARE"),
  seed = 7L
)

expected_fields <- c("inputs", "p_table", "replicate_summaries", "errors_long", "phat_long", "curves", "argmax_summary")
if (!all(expected_fields %in% names(res))) stop("result missing expected fields")
pass("result has all expected fields")

if (!identical(res$inputs$proportion_method, "beta")) stop("default proportion_method should be beta")
pass("default proportion_method is beta")

if (abs(sum(as.numeric(res$p_table[1, grep("^index_", names(res$p_table))])) - 1) > 1e-12) {
  stop("p_table first row does not sum to 1")
}
pass("p_table rows encode valid proportions")

if (nrow(res$replicate_summaries) != 50L * 2L) stop("replicate_summaries wrong row count")
pass("replicate_summaries row count correct")

if (nrow(res$phat_long) != 50L * 10L) stop("phat_long wrong row count")
if (!all(c("alpha", "p_max", "replicate", "index", "phat") %in% names(res$phat_long))) {
  stop("phat_long missing expected columns")
}
pass("phat_long has expected shape and columns")

if (nrow(res$curves) != 3L * 2L) stop("curves has wrong row count")
pass("curves row count correct")
if (!("p_max" %in% names(res$curves))) stop("curves should include p_max metadata")
pass("curves includes p_max metadata")

res_fixed_method <- run_simulation_experiment(
  alpha = 2, K = 10L, n = 100L, B = 5L,
  taus = c(0.05, 0.10),
  proportion_method = "fixed_max_beta",
  p_max = 0.4,
  seed = 9L
)
if (!identical(res_fixed_method$inputs$proportion_method, "fixed_max_beta")) {
  stop("fixed_max_beta should be selectable through run_simulation_experiment")
}
if (!identical(res_fixed_method$inputs$p_max, 0.4)) {
  stop("run_simulation_experiment should record p_max in inputs")
}
if (!all(res_fixed_method$p_table$p_max == 0.4)) {
  stop("run_simulation_experiment should record p_max in p_table rows")
}
# Extract the proportions from the first simulation's index_1 ... index_K columns in p_table.
proportion_columns <- paste0("index_", seq_len(res_fixed_method$inputs$K))
simulated_proportions <- as.numeric(res_fixed_method$p_table[1, proportion_columns])
if (abs(simulated_proportions[length(simulated_proportions)] - 0.4) > 1e-12) {
  stop("run_simulation_experiment should use fixed_max_beta proportions")
}
pass("run_simulation_experiment routes fixed_max_beta through dispatcher")

# smoke test with named-list taus (different grid lengths per metric)
set.seed(7L)
res_list <- run_simulation_experiment(
  alpha = 2, K = 10L, n = 100L, B = 50L,
  taus  = list(AE = c(0.05, 0.10), ARE = c(0.05, 0.10, 0.20)),
  metrics = c("AE", "ARE"),
  seed = 7L
)

if (nrow(res_list$curves[res_list$curves$metric == "AE",  ]) != 2L) stop("named-list: AE curve row count wrong")
if (nrow(res_list$curves[res_list$curves$metric == "ARE", ]) != 3L) stop("named-list: ARE curve row count wrong")
pass("run_simulation_experiment: named-list taus produces correct per-metric row counts")

warn_multi_pmax <- character()
res_multi_pmax <- withCallingHandlers(
  run_simulation_experiment(
    alpha = c(2, 3),
    K = 2L, n = 100L, B = 5L,
    taus = c(0.05, 0.10),
    proportion_method = "fixed_max_beta",
    p_max = c(0.4, 0.8),
    seed = 11L
  ),
  warning = function(w) {
    warn_multi_pmax <<- c(warn_multi_pmax, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
if (!any(grepl("Impossible fixed_max_beta combination", warn_multi_pmax))) {
  stop("multi-p_max run should warn for impossible combinations")
}
pass("multi-p_max run warns for impossible combinations")
if (!setequal(unique(res_multi_pmax$curves$p_max), 0.8)) {
  stop("multi-p_max run should skip impossible combinations and keep feasible p_max only")
}
pass("multi-p_max run skips impossible combinations and keeps feasible ones")

# ---- 9. hybrid cutoff helpers ------------------------------------------------
cat("\n9. hybrid cutoff helpers\n")

phat_toy <- matrix(
  c(
    0.10, 0.90,
    0.02, 0.98
  ),
  nrow = 2L,
  byrow = TRUE
)
p_toy <- c(0.01, 0.99)

hybrid_eval <- evaluate_hybrid_success_cell_level(
  phat_mat = phat_toy,
  p = p_toy,
  cutoff = 0.05,
  tau_AE = 0.05,
  tau_ARE = 0.20
)
if (!isTRUE(all.equal(hybrid_eval$success_rate_cell, 0.75))) {
  stop("hybrid cell-level success rate mismatch")
}
if (!isTRUE(all.equal(hybrid_eval$success_rate_replicate, 0.5))) {
  stop("hybrid replicate-level success rate mismatch")
}
pass("evaluate_hybrid_success_cell_level computes expected rates")

hybrid_curve <- sweep_hybrid_cutoffs_cell_level(
  phat_mat = phat_toy,
  p = p_toy,
  cutoffs = c(0.02, 0.05),
  tau_AE = 0.05,
  tau_ARE = 0.20
)
if (nrow(hybrid_curve) != 2L) stop("sweep_hybrid_cutoffs_cell_level should return one row per cutoff")
pass("sweep_hybrid_cutoffs_cell_level row count correct")

best_hybrid <- find_best_hybrid_cutoff(
  curve_df = hybrid_curve,
  maximize = "cell",
  tie_break = "smallest"
)
if (!isTRUE(all.equal(best_hybrid$best_cutoff, 0.02))) {
  stop("find_best_hybrid_cutoff should select smallest tied cutoff")
}
if (best_hybrid$n_tied_maximizers != 2L) {
  stop("find_best_hybrid_cutoff should report tied maximizer count")
}
pass("find_best_hybrid_cutoff tie handling works")

hybrid_run <- run_hybrid_cutoff_analysis(
  phat_mat = phat_toy,
  p = p_toy,
  cutoffs = c(0.02, 0.05),
  tau_AE = 0.05,
  tau_ARE = 0.20,
  maximize = "cell",
  alpha = 2,
  p_max = 0.5
)
if (nrow(hybrid_run$best_summary) != 1L || !("best_cutoff" %in% names(hybrid_run$best_summary))) {
  stop("run_hybrid_cutoff_analysis should return one-row best_summary")
}
pass("run_hybrid_cutoff_analysis returns best_summary")

hybrid_exp <- run_simulation_hybrid_cutoff_experiment(
  alpha = 2,
  K = 5L,
  n = 100L,
  B = 10L,
  cutoffs = c(0.02, 0.05),
  tau_AE = 0.05,
  tau_ARE = 0.20,
  seed = 21L
)
expected_hybrid_fields <- c("inputs", "p_table", "phat_long", "cutoff_curves", "best_cutoff_summary")
if (!all(expected_hybrid_fields %in% names(hybrid_exp))) {
  stop("run_simulation_hybrid_cutoff_experiment missing expected fields")
}
pass("run_simulation_hybrid_cutoff_experiment returns expected fields")

# heatmap over AE/ARE threshold pairs should expose best cutoffs per cell
phat_heatmap <- matrix(
  c(
    0.01, 0.99,
    0.06, 0.94,
    0.20, 0.80
  ),
  nrow = 3L,
  byrow = TRUE
)
heatmap_plot <- plot_hybrid_best_cutoff_heatmap(
  phat_mat = phat_heatmap,
  p = c(0.02, 0.98),
  cutoffs = c(0.02, 0.05, 0.10, 0.20),
  tau_AE_values = c(0.01, 0.05),
  tau_ARE_values = c(0.05, 0.20),
  maximize = "cell",
  label_digits = 3L
)
if (!inherits(heatmap_plot, "ggplot")) stop("plot_hybrid_best_cutoff_heatmap should return a ggplot object")
heatmap_build <- ggplot2::ggplot_build(heatmap_plot)
if (nrow(heatmap_build$data[[1]]) != 4L) stop("heatmap should include one tile per AE/ARE threshold pair")
if (!setequal(unique(heatmap_build$data[[2]]$label), c("0.020", "0.100"))) {
  stop("heatmap labels should include expected best-cutoff values")
}
pass("plot_hybrid_best_cutoff_heatmap returns expected tiles and labels")

heatmap_plot_k3 <- plot_hybrid_best_cutoff_heatmap(
  phat_mat = matrix(c(0.10, 0.20, 0.70, 0.15, 0.25, 0.60), nrow = 2L, byrow = TRUE),
  p = c(0.10, 0.20, 0.70),
  cutoffs = c(0.05, 0.10, 0.20),
  tau_AE_values = c(0.02, 0.05),
  tau_ARE_values = c(0.10, 0.20),
  maximize = "cell"
)
if (!inherits(heatmap_plot_k3, "ggplot")) {
  stop("plot_hybrid_best_cutoff_heatmap should support K > 2 inputs")
}
pass("plot_hybrid_best_cutoff_heatmap supports K > 2 inputs")

# ---- 10. plotting helpers: p_max filters ------------------------------------
cat("\n10. plotting helper filters\n")

proportions_plot_beta <- plot_proportions_curve(res)
if (!inherits(proportions_plot_beta, "ggplot")) stop("plot_proportions_curve should return a ggplot object for beta results")
pass("plot_proportions_curve supports beta simulation results")

proportions_plot_fixed <- plot_proportions_curve(res_fixed_method)
if (!inherits(proportions_plot_fixed, "ggplot")) stop("plot_proportions_curve should return a ggplot object for fixed_max_beta results")
pass("plot_proportions_curve supports fixed_max_beta simulation results")

fixed_plot_build <- ggplot2::ggplot_build(proportions_plot_fixed)
fixed_curve_x <- fixed_plot_build$data[[1]]$x
fixed_point_x <- fixed_plot_build$data[[2]]$x
comparison_tol <- 1e-12
if (!any(abs(fixed_point_x - 1) < comparison_tol)) {
  stop("fixed_max_beta plot should include the fixed-maximum point at x = 1")
}
if (any(abs(fixed_curve_x - 1) < comparison_tol)) {
  stop("fixed_max_beta plot curve should exclude x = 1 so the fixed maximum is not on the curve")
}
pass("fixed_max_beta plot excludes the fixed maximum point from the curve")

success_plot <- plot_success_rate_curve(
  res_multi_pmax,
  metric = "AE",
  alphas = c(2, 3),
  p_maxs = 0.8
)
if (!inherits(success_plot, "ggplot")) stop("plot_success_rate_curve should return a ggplot object")
pass("plot_success_rate_curve supports p_max filtering")

argmax_plot <- plot_argmax_histogram(
  res_multi_pmax,
  metric = "AE",
  alphas = c(2, 3),
  p_maxs = 0.8
)
if (!inherits(argmax_plot, "ggplot")) stop("plot_argmax_histogram should return a ggplot object")
pass("plot_argmax_histogram supports p_max filtering")

err_plot_missing_metric <- tryCatch(
  plot_argmax_histogram(res_multi_pmax),
  error = function(e) e$message
)
if (!identical(err_plot_missing_metric, "metric must be a single value: 'AE', 'ARE', 'TSE', or 'LAE'.")) {
  stop("plot_argmax_histogram should require a valid metric")
}
pass("plot_argmax_histogram requires metric")

err_plot_pmax <- tryCatch(
  plot_argmax_histogram(res_multi_pmax, metric = "AE", p_maxs = 0.4),
  error = function(e) e$message
)
if (!grepl("No rows match the requested p_max", err_plot_pmax)) {
  stop("plot_argmax_histogram should error for unmatched p_max filter")
}
pass("plot_argmax_histogram gives informative error for unmatched p_max filter")

res_hist_grid <- run_simulation_experiment(
  alpha = c(2, 3, 4),
  K = 10L, n = 50L, B = 5L,
  taus = c(0.05, 0.10),
  proportion_method = "fixed_max_beta",
  p_max = c(0.6, 0.8),
  seed = 17L
)
argmax_plot_grid <- plot_argmax_histogram(res_hist_grid, metric = "AE")
layout_grid <- ggplot2::ggplot_build(argmax_plot_grid)$layout$layout
expected_alpha_n <- length(unique(res_hist_grid$replicate_summaries$alpha))
expected_p_max_n <- length(unique(res_hist_grid$replicate_summaries$p_max))
if (length(unique(layout_grid$COL)) != expected_alpha_n || length(unique(layout_grid$ROW)) != expected_p_max_n) {
  stop("plot_argmax_histogram should facet with alpha in columns and p_max in rows")
}
pass("plot_argmax_histogram facets with alpha columns and p_max rows")

# ---- Sample-size estimation helpers ----------------------------------------
cat("\nSample-size estimation helpers\n")

ss_config <- list(
  K                    = 5L,
  B                    = 20L,
  taus                 = list(AE = 0.05, ARE = 2.0),
  metrics              = c("AE", "ARE"),
  model                = "multinomial",
  tie_method           = "random",
  proportion_method    = "beta",
  seed                 = 42L,
  success_rate_target  = 0.8,
  sample_size_tolerance = 10L,
  max_iterations       = 3L
)

# simulate_success_at_n: basic shape checks
sim_succ <- simulate_success_at_n(alpha = 2, n = 500L, config = ss_config)
if (!is.logical(sim_succ$success)) stop("success must be a logical vector")
if (length(sim_succ$success) != ss_config$B) stop("success must have length B")
if (!is.numeric(sim_succ$success_rate) || sim_succ$success_rate < 0 || sim_succ$success_rate > 1) {
  stop("success_rate must be in [0, 1]")
}
if (sim_succ$success_count != sum(sim_succ$success)) stop("success_count must equal sum(success)")
pass("simulate_success_at_n returns correct structure")

# Success requires both AE and ARE thresholds when both are in metrics
# Use extreme thresholds to force known outcomes
cfg_all_pass <- ss_config
cfg_all_pass$taus <- list(AE = 100, ARE = 100)  # everything passes
sim_all <- simulate_success_at_n(alpha = 2, n = 200L, config = cfg_all_pass)
if (!all(sim_all$success)) stop("all replicates should succeed with very loose thresholds")
pass("all replicates succeed with very loose thresholds")

cfg_none_pass <- ss_config
cfg_none_pass$taus <- list(AE = 0, ARE = 0)  # nothing passes
sim_none <- simulate_success_at_n(alpha = 2, n = 200L, config = cfg_none_pass)
if (any(sim_none$success)) stop("no replicates should succeed with zero thresholds")
pass("no replicates succeed with zero thresholds")

# Only AE metric: success based only on AE threshold
cfg_ae_only <- ss_config
cfg_ae_only$metrics <- "AE"
cfg_ae_only$taus    <- list(AE = 0.05)
sim_ae <- simulate_success_at_n(alpha = 2, n = 500L, config = cfg_ae_only)
if (length(sim_ae$success) != cfg_ae_only$B) stop("wrong length for AE-only success")
pass("simulate_success_at_n works with single metric")

# fit_success_glm
n_vals <- c(100L, 200L, 300L)
s_cnt  <- c(5L,   12L,  18L)
glm_fit <- fit_success_glm(n_vals, s_cnt, B = 20L)
if (!inherits(glm_fit, "glm")) stop("fit_success_glm must return a glm object")
if (length(stats::coef(glm_fit)) != 2L) stop("glm must have intercept and slope")
pass("fit_success_glm returns a binomial glm with two coefficients")

# solve_sample_size_from_glm: result is always rounded up
solved <- solve_sample_size_from_glm(glm_fit, target_success_rate = 0.8)
if (!is.list(solved) || is.null(solved$n_raw) || is.null(solved$n_rounded)) {
  stop("solve_sample_size_from_glm must return a list with n_raw and n_rounded")
}
if (solved$n_rounded != as.integer(ceiling(solved$n_raw))) {
  stop("n_rounded must equal ceiling(n_raw)")
}
pass("solve_sample_size_from_glm rounds up")

# iterate_sample_size_for_alpha: diagnostics structure
iter_result <- iterate_sample_size_for_alpha(alpha = 2, n_init = 200L, config = ss_config)
if (!is.list(iter_result)) stop("iterate_sample_size_for_alpha must return a list")
if (!is.integer(iter_result$final_n) || iter_result$final_n < 1L) {
  stop("final_n must be a positive integer")
}
if (!iter_result$stopping_reason %in% c("tolerance", "max_iterations")) {
  stop("stopping_reason must be 'tolerance' or 'max_iterations'")
}
req_cols <- c("alpha", "iteration", "pilot_index", "n", "success_count",
              "success_rate", "target_success_rate", "glm_intercept",
              "glm_slope", "n_raw", "n_rounded", "n_clamped", "stopping_reason")
missing_cols <- setdiff(req_cols, names(iter_result$diagnostics))
if (length(missing_cols) > 0L) {
  stop(paste("diagnostics missing columns:", paste(missing_cols, collapse = ", ")))
}
if (nrow(iter_result$diagnostics) %% 3L != 0L) {
  stop("diagnostics must have a multiple of 3 rows (3 pilot points per iteration)")
}
pass("iterate_sample_size_for_alpha returns correct structure")

# Clamping: when success rate is near 0%, new_n should not exceed 2x the largest pilot
cfg_clamp_low <- ss_config
cfg_clamp_low$taus           <- list(AE = 0, ARE = 0)  # force 0% success
cfg_clamp_low$max_iterations <- 1L
cfg_clamp_low$sample_size_tolerance <- Inf
n_init_clamp <- 200L
iter_clamp_low <- iterate_sample_size_for_alpha(alpha = 2, n_init = n_init_clamp, config = cfg_clamp_low)
upper_bound <- as.integer(ceiling(2.0 * as.integer(ceiling(1.05 * n_init_clamp))))
if (iter_clamp_low$final_n > upper_bound) {
  stop("final_n must not exceed 2x the largest pilot when success rate is near 0%")
}
pass("clamping prevents new_n from exceeding 2x largest pilot when success is ~0%")

# Clamping: when success rate is near 100%, new_n should not be less than 0.5x the smallest pilot
cfg_clamp_high <- ss_config
cfg_clamp_high$taus           <- list(AE = 100, ARE = 100)  # force 100% success
cfg_clamp_high$max_iterations <- 1L
cfg_clamp_high$sample_size_tolerance <- Inf
iter_clamp_high <- iterate_sample_size_for_alpha(alpha = 2, n_init = n_init_clamp, config = cfg_clamp_high)
lower_bound <- as.integer(ceiling(0.5 * as.integer(ceiling(0.95 * n_init_clamp))))
if (iter_clamp_high$final_n < lower_bound) {
  stop("final_n must not go below 0.5x the smallest pilot when success rate is ~100%")
}
pass("clamping prevents new_n from going below 0.5x smallest pilot when success is ~100%")

# Stopping at tolerance: use tolerance = Inf so it stops in 1 iteration
cfg_tol_inf <- ss_config
cfg_tol_inf$sample_size_tolerance <- Inf
iter_inf <- iterate_sample_size_for_alpha(alpha = 2, n_init = 200L, config = cfg_tol_inf)
if (iter_inf$stopping_reason != "tolerance") {
  stop("stopping_reason should be 'tolerance' when tolerance = Inf")
}
if (iter_inf$iterations_used != 1L) stop("should stop after 1 iteration when tolerance = Inf")
pass("stops after first iteration when tolerance is Inf")

# max_iterations stopping: use tolerance = 0 and 1 iteration
cfg_max1 <- ss_config
cfg_max1$sample_size_tolerance <- 0L
cfg_max1$max_iterations        <- 1L
iter_max1 <- iterate_sample_size_for_alpha(alpha = 2, n_init = 500L, config = cfg_max1)
if (iter_max1$iterations_used > 1L) stop("should use at most 1 iteration")
pass("max_iterations respected")

# Alpha propagation: final_n of first alpha feeds next alpha's n_init
cfg_prop <- ss_config
cfg_prop$sample_size_tolerance <- Inf  # always stop after 1 iteration

iter_a1 <- iterate_sample_size_for_alpha(alpha = 2,   n_init = 300L, config = cfg_prop)
iter_a2 <- iterate_sample_size_for_alpha(alpha = 2.5, n_init = iter_a1$final_n, config = cfg_prop)

# The second alpha's diagnostics$n values should be based on iter_a1$final_n
expected_pilots <- as.integer(ceiling(c(0.95, 1.00, 1.05) * iter_a1$final_n))
actual_pilots   <- iter_a2$diagnostics$n[iter_a2$diagnostics$iteration == 1L]
if (!identical(sort(actual_pilots), sort(expected_pilots))) {
  stop("alpha propagation: second alpha must use first alpha's final_n as n_init")
}
pass("final_n from one alpha is correctly propagated to the next alpha's n_init")

# ---- Done ------------------------------------------------------------------
cat("\nAll tests passed.\n")
