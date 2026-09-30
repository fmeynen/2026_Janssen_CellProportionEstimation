# Small smoke tests for run_simulation_experiment (multinomial model).

test_that("run_simulation_experiment returns the documented fields and shapes", {
  res <- run_simulation_experiment(
    alpha = 2, K = 10L, n = 100L, B = 10L,
    taus = c(0.05, 0.10, 0.20), metrics = c("AE", "ARE"), seed = 7L
  )
  expect_true(all(
    c("inputs", "p_table", "replicate_summaries", "errors_long", "phat_long", "curves", "argmax_summary") %in%
      names(res)
  ))
  expect_identical(res$inputs$proportion_method, "beta")
  expect_equal(sum(as.numeric(res$p_table[1, grep("^index_", names(res$p_table))])), 1, tolerance = 1e-12)
  expect_equal(nrow(res$replicate_summaries), 10L * 2L)
  expect_equal(nrow(res$phat_long), 10L * 10L)
  expect_true(all(c("alpha", "p_max", "replicate", "index", "phat") %in% names(res$phat_long)))
  expect_equal(nrow(res$curves), 3L * 2L)
  expect_true("p_max" %in% names(res$curves))
})

test_that("run_simulation_experiment routes fixed_max_beta through the dispatcher and records p_max", {
  res <- run_simulation_experiment(
    alpha = 2, K = 10L, n = 100L, B = 5L, taus = c(0.05, 0.10),
    proportion_method = "fixed_max_beta", p_max = 0.4, seed = 9L
  )
  expect_identical(res$inputs$proportion_method, "fixed_max_beta")
  expect_identical(res$inputs$p_max, 0.4)
  expect_true(all(res$p_table$p_max == 0.4))
  props <- as.numeric(res$p_table[1, paste0("index_", seq_len(res$inputs$K))])
  expect_equal(props[length(props)], 0.4, tolerance = 1e-12)
})

test_that("run_simulation_experiment with named-list taus gives per-metric curve row counts", {
  res <- run_simulation_experiment(
    alpha = 2, K = 10L, n = 100L, B = 10L,
    taus = list(AE = c(0.05, 0.10), ARE = c(0.05, 0.10, 0.20)),
    metrics = c("AE", "ARE"), seed = 7L
  )
  expect_equal(sum(res$curves$metric == "AE"), 2L)
  expect_equal(sum(res$curves$metric == "ARE"), 3L)
})

test_that("a multi-p_max run warns on impossible combinations and keeps only feasible p_max", {
  expect_warning(
    res <- run_simulation_experiment(
      alpha = c(2, 3), K = 2L, n = 100L, B = 5L, taus = c(0.05, 0.10),
      proportion_method = "fixed_max_beta", p_max = c(0.4, 0.8), seed = 11L
    ),
    "Impossible fixed_max_beta combination"
  )
  expect_setequal(unique(res$curves$p_max), 0.8)
})
