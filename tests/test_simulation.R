# tests/test_simulation.R
#
# Base-R test script for simulation_functions.R.
#
# Run from the repository root:
#   Rscript tests/test_simulation.R
#
# Each test calls `stop()` with a descriptive message on failure so that CI
# (or a local developer) can immediately identify what broke.
# ---------------------------------------------------------------------------

source(file.path("scripts", "simulation_functions.R"))

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

# ---- 2. simulate_counts_multinomial ----------------------------------------
cat("\n2. simulate_counts_multinomial\n")

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

# ---- 3. simulate_counts dispatcher -----------------------------------------
cat("\n3. simulate_counts dispatcher\n")

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

# ---- 4. compute_errors ------------------------------------------------------
cat("\n4. compute_errors\n")

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

# ---- 5. max_error_summary ---------------------------------------------------
cat("\n5. max_error_summary\n")

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

# ---- 6. evaluate_thresholds -------------------------------------------------
cat("\n6. evaluate_thresholds\n")

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
if (nrow(curves_list) != 2L + 2L) stop("named-list curves: wrong row count")
pass("named-list taus: row count correct with different grid lengths")

# missing metric in named-list should error
err_missing <- tryCatch(
  evaluate_thresholds(me, taus = list(AE = c(0.10))),
  error = function(e) e$message
)
if (!grepl("ARE", err_missing)) stop("missing metric in taus list should error with metric name")
pass("named-list taus: missing metric produces informative error")

# ---- 7. run_simulation_experiment (smoke test) ------------------------------
cat("\n7. run_simulation_experiment smoke test\n")

set.seed(7L)
res <- run_simulation_experiment(
  alpha = 2, K = 10L, n = 100L, B = 50L,
  taus = c(0.05, 0.10, 0.20),
  metrics = c("AE", "ARE"),
  seed = 7L
)

expected_fields <- c("inputs", "p", "max_errors", "argmax", "curves", "argmax_summary")
if (!all(expected_fields %in% names(res))) stop("result missing expected fields")
pass("result has all expected fields")

if (abs(sum(res$p) - 1) > 1e-12) stop("res$p does not sum to 1")
pass("res$p sums to 1")

if (!all(dim(res$max_errors) == c(50L, 2L))) stop("max_errors wrong dimensions")
pass("max_errors dimensions correct")

if (nrow(res$curves) != 3L * 2L) stop("curves has wrong row count")
pass("curves row count correct")

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

# ---- Done ------------------------------------------------------------------
cat("\nAll tests passed.\n")
