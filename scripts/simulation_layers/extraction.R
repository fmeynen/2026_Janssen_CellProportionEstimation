
# Extraction Layer ------------------------------------------------------------------------------------------------
# Extraction layer: pull structured outputs from simulation results.


#' Extract the maximum error value and its index from one error vector.
#'
#' @param error_vec  Numeric error vector (length K).
#' @param tie_method How to break ties among indices sharing the maximum:
#'   \describe{
#'     \item{"random"}{sample uniformly among tied indices (default)}
#'     \item{"first"}{smallest tied index}
#'     \item{"last"}{largest tied index}
#'   }
#'
#' @return Named list: `max_error_value` (numeric scalar) and `argmax_index` (integer scalar).
max_error_summary <- function(error_vec, tie_method = c("random", "first", "last")) {
  tie_method <- match.arg(tie_method)
  max_val <- max(error_vec, na.rm = TRUE)
  tied <- which(error_vec == max_val)          # always a proper integer vector
  argmax_idx <- switch(tie_method,
    random = tied[sample.int(length(tied), 1L)],   # safe even when length(tied)==1
    first  = tied[[1L]],
    last   = tied[[length(tied)]]
  )
  list(max_error_value = max_val, argmax_index = argmax_idx)
}


#' Tabulate which cell-type index most often achieved the max error.
#'
#' @param argmax  B x M integer matrix of argmax indices (from `run_replicates()`).
#' @param p       True proportion vector (length K).
#'
#' @return Tidy data.frame with columns: metric, index, count, fraction, p_value.
summarize_argmax <- function(argmax, p) {
  validate_named_matrix(argmax, "argmax")
  K <- length(p)
  metrics <- colnames(argmax)
  rows <- vector("list", length(metrics))
  for (i in seq_along(metrics)) {
    m <- metrics[[i]]
    counts <- tabulate(argmax[, m], nbins = K)
    rows[[i]] <- data.frame(
      metric   = m,
      index    = seq_len(K),
      count    = counts,
      fraction = counts / sum(counts),
      p_value  = p,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


#' Build one row of the p_table for a single alpha/p_max scenario.
#'
#' @param alpha_i  Numeric scalar; the alpha value for this scenario.
#' @param p_max_i  Numeric scalar (or NA); the p_max value for this scenario.
#' @param p        Numeric vector of length K; true proportions for this scenario.
#' @param K        Integer; number of cell types.
#' @param before   Optional named list of extra columns placed before `alpha`.
#' @param after    Optional named list of extra columns placed between `p_max` and `index_1`.
#'
#' @return A single-row data.frame with columns: [before], alpha, p_max, [after], index_1, ..., index_K.
extract_p_table_row <- function(alpha_i, p_max_i, p, K, before = NULL, after = NULL) {
  do.call(
    data.frame,
    c(
      before,
      list(alpha = alpha_i, p_max = p_max_i),
      after,
      as.list(stats::setNames(as.numeric(p), paste0("index_", seq_len(K)))),
      list(stringsAsFactors = FALSE, check.names = FALSE)
    )
  )
}


#' Build the replicate_summaries data.frame for one alpha/p_max scenario.
#'
#' @param rep_out  Output list from `run_replicates()`.
#' @param alpha_i  Numeric scalar; the alpha value for this scenario.
#' @param p_max_i  Numeric scalar (or NA); the p_max value for this scenario.
#' @param B        Integer; number of replicates.
#' @param metrics  Character vector of metric names.
#'
#' @return Tidy data.frame with columns:
#'   alpha, p_max, replicate, metric, max_error, argmax_index.
extract_replicate_summaries <- function(rep_out, alpha_i, p_max_i, B, metrics) {
  data.frame(
    alpha = alpha_i,
    p_max = p_max_i,
    replicate = rep(seq_len(B), times = length(metrics)),
    metric = rep(metrics, each = B),
    max_error = as.vector(rep_out$max_errors[, metrics, drop = FALSE]),
    argmax_index = as.vector(rep_out$argmax[, metrics, drop = FALSE]),
    stringsAsFactors = FALSE
  )
}


#' Build the errors_long data.frame for one alpha/p_max scenario.
#'
#' @param rep_out  Output list from `run_replicates()`.
#' @param alpha_i  Numeric scalar; the alpha value for this scenario.
#' @param p_max_i  Numeric scalar (or NA); the p_max value for this scenario.
#' @param B        Integer; number of replicates.
#' @param metrics  Character vector of metric names.
#'
#' @return Tidy data.frame with columns:
#'   alpha, p_max, replicate, metric, index, error.
extract_errors_long <- function(rep_out, alpha_i, p_max_i, B, metrics) {
  errors_m_list <- vector("list", length(metrics))
  for (j in seq_along(metrics)) {
    m <- metrics[[j]]
    errors_m <- rep_out$errors[, , m, drop = TRUE]
    if (is.null(dim(errors_m))) errors_m <- matrix(errors_m, nrow = B)
    errors_m_list[[j]] <- data.frame(
      alpha = alpha_i,
      p_max = p_max_i,
      replicate = rep(seq_len(B), times = ncol(errors_m)),
      metric = m,
      index = rep(seq_len(ncol(errors_m)), each = B),
      error = as.vector(errors_m),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, errors_m_list)
}


#' Build the phat_long data.frame for one alpha/p_max scenario.
#'
#' @param rep_out  Output list from `run_replicates()`.
#' @param alpha_i  Numeric scalar; the alpha value for this scenario.
#' @param p_max_i  Numeric scalar (or NA); the p_max value for this scenario.
#' @param B        Integer; number of replicates.
#'
#' @return Tidy data.frame with columns:
#'   alpha, p_max, replicate, index, phat.
extract_phat_long <- function(rep_out, alpha_i, p_max_i, B) {
  data.frame(
    alpha = alpha_i,
    p_max = p_max_i,
    replicate = rep(seq_len(B), times = ncol(rep_out$phat)),
    index = rep(seq_len(ncol(rep_out$phat)), each = B),
    phat = as.vector(rep_out$phat),
    stringsAsFactors = FALSE
  )
}


# Dirichlet-multinomial success rate ------------------------------------------------------------------------------


#' Identifier of the current replicate success rule.
#'
#' Included in the cache keys of experiments whose cached results depend on the rule (`run_sample_size_experiment()`,
#' `run_simulation_dm_errorchoice()`), so that changing the rule never silently reuses results computed under an
#' earlier one. Change this string whenever `replicate_pooled_error()` changes what it computes.
#'
#' @return Character scalar.
success_rule_id <- function() {
  "pooled_proportion_v1"
}


#' Compute the per-replicate pooled-proportion error statistic for one metric.
#'
#' Single source of truth for the success statistic used by `replicate_success()` (and by any caller that needs the
#' raw per-replicate statistic, e.g. to evaluate many tau thresholds via `mean(stat <= tau)` without
#' recomputation). Per (scenario, replicate, cell type j), the estimated proportions are first averaged over persons,
#' giving the pooled estimate `pbar_j`, which is compared with the population-level true proportion `p_j`:
#'   * AE:  `|pbar_j - p_j|`
#'   * ARE: `|pbar_j - p_j| / p_j`
#' The largest of these cell-type errors is the per-replicate "stat"; a replicate passes the metric iff
#' `stat <= tau`. `NaN` (ARE with `pbar_j` and `p_j` both 0) counts as 0; `Inf` (ARE with `p_j = 0` and
#' `pbar_j > 0`) is kept. With a single person per replicate this reduces to the per-draw error against `p`.
#'
#' The computation is fully vectorised: integer group keys, `rowsum()` for the per-cell-type means, and a single
#' `order()` + `duplicated()` pass for the max-by-group step. There is no `stats::aggregate()` call, no
#' `tapply()`, and no per-row loop.
#'
#' @param person_results Data.frame with at least the columns `replicate`, `cell_type`, `metric`,
#'   `observed_proportion` and `population_mean_proportion`, one row per (replicate, person, cell type, metric). A
#'   `scenario_id` column is optional: when it is absent, or present but entirely `NA`, all rows are treated as a
#'   single scenario and the output's `scenario_id` is `NA_character_`.
#' @param metric Single metric name, `"AE"` or `"ARE"`. Rows of `person_results` are filtered to this metric (so
#'   each person's proportions are counted once); an error is raised if none match.
#'
#' @return Data.frame with one row per (scenario_id, replicate), sorted by scenario_id then replicate, with
#'   columns scenario_id, replicate, stat (max over cell types of the pooled-proportion error).
replicate_pooled_error <- function(person_results, metric) {
  validate_required_columns(person_results, person_results_required_cols)
  if (!is.character(metric) || length(metric) != 1L || is.na(metric)) {
    stop("metric must be a single character string.", call. = FALSE)
  }
  if (!metric %in% c("AE", "ARE")) {
    stop(sprintf("The pooled success rule is defined for AE and ARE only, not '%s'.", metric), call. = FALSE)
  }

  idx <- person_results$metric == metric
  if (!any(idx)) {
    stop(sprintf("metric '%s' is not present in person_results.", metric), call. = FALSE)
  }

  # A missing/all-NA scenario_id is replaced by a non-NA sentinel for grouping purposes (rowsum() drops NA
  # groups), and mapped back to NA_character_ in the output.
  has_scenario <- "scenario_id" %in% names(person_results)
  scenario_raw <- if (has_scenario) {
    as.character(person_results$scenario_id)
  } else {
    rep(NA_character_, nrow(person_results))
  }
  no_scenario <- all(is.na(scenario_raw))
  scenario_key <- if (no_scenario) rep("__single_scenario__", length(scenario_raw)) else scenario_raw

  scenario_i   <- scenario_key[idx]
  replicate_i  <- person_results$replicate[idx]
  cell_type_i  <- person_results$cell_type[idx]
  observed_i   <- person_results$observed_proportion[idx]
  population_i <- person_results$population_mean_proportion[idx]

  scenario_lv  <- sort(unique(scenario_i))
  replicate_lv <- sort(unique(replicate_i))
  cell_type_lv <- sort(unique(cell_type_i))
  nS <- length(scenario_lv)
  nR <- length(replicate_lv)
  nC <- length(cell_type_lv)

  s_idx <- match(scenario_i, scenario_lv)
  r_idx <- match(replicate_i, replicate_lv)
  c_idx <- match(cell_type_i, cell_type_lv)

  # Stage 1: pooled estimate (mean over persons) and population proportion per (scenario, replicate, cell_type) via
  # one integer group key and rowsum(); then the per-cell-type error of the pooled estimate.
  key1 <- (s_idx - 1L) * nR * nC + (r_idx - 1L) * nC + c_idx
  counts1     <- rowsum(rep(1L, length(observed_i)), key1)[, 1L]
  pooled1     <- rowsum(observed_i, key1)[, 1L] / counts1
  population1 <- rowsum(population_i, key1)[, 1L] / counts1
  error1 <- abs(pooled1 - population1)
  if (identical(metric, "ARE")) {
    error1 <- error1 / population1
  }
  error1[is.nan(error1)] <- 0
  key1_sorted <- as.numeric(names(counts1))

  tmp1 <- (key1_sorted - 1L) %/% nC
  r1 <- (tmp1 %% nR) + 1L
  s1 <- (tmp1 %/% nR) + 1L

  # Stage 2: max over cell types per (scenario, replicate). Sort by group then by value descending and keep the
  # first row of each group -- fully vectorised, no tapply()/aggregate() and no per-row loop.
  key2 <- (s1 - 1L) * nR + r1
  ord  <- order(key2, -error1)
  keep <- !duplicated(key2[ord])
  max_key2 <- key2[ord][keep]
  stat     <- error1[ord][keep]

  r2 <- ((max_key2 - 1L) %% nR) + 1L
  s2 <- ((max_key2 - 1L) %/% nR) + 1L

  out <- data.frame(
    scenario_id = scenario_lv[s2],
    replicate   = replicate_lv[r2],
    stat        = stat,
    stringsAsFactors = FALSE
  )
  if (no_scenario) out$scenario_id <- NA_character_

  out <- out[order(out$scenario_id, out$replicate), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Determine per-replicate success against a set of metric thresholds.
#'
#' Replicate success rule shared by `extract_success_rate()` and `simulate_success_at_n()` (both models). For each
#' metric, the per-replicate statistic comes from `replicate_pooled_error()`, the single source of truth for the rule:
#' the estimated proportions are averaged over persons per cell type, compared with the population-level true
#' proportion, and the largest cell-type error is taken. A replicate passes that metric if this value is `<= tau`,
#' and passes overall if it passes every metric in `taus`. This function only turns each metric's `stat` into a
#' `pass_<metric>` flag and ANDs them together.
#'
#' @param person_results Data.frame with at least the columns `replicate`, `cell_type`, `metric`,
#'   `observed_proportion`, `population_mean_proportion`. A `scenario_id` column is optional: when it is absent, or present but entirely `NA`, all rows are treated as a
#'   single scenario and the output's `scenario_id` is `NA_character_`.
#' @param taus   Named list with one scalar threshold per metric (e.g. `list(AE = 0.02, ARE = 0.5)`). Metrics missing
#'   from `person_results$metric` are skipped with a warning; an error is raised if none of the metrics remain.
#'
#' @return Data.frame with one row per (scenario_id, replicate), sorted by scenario_id then replicate, with columns
#'   scenario_id, replicate, `pass_<metric>` for each metric used, and `pass` (logical AND across all `pass_<metric>`
#'   columns).
replicate_success <- function(person_results, taus) {
  validate_required_columns(person_results, person_results_required_cols)
  if (!is.list(taus) || is.null(names(taus)) || any(names(taus) == "")) {
    stop("taus must be a named list with one scalar threshold per metric.", call. = FALSE)
  }
  for (m in names(taus)) {
    if (!is.numeric(taus[[m]]) || length(taus[[m]]) != 1L || is.na(taus[[m]])) {
      stop(sprintf("taus$%s must be a single numeric threshold.", m), call. = FALSE)
    }
  }

  metrics <- names(taus)
  missing_metrics <- setdiff(metrics, unique(person_results$metric))
  for (m in missing_metrics) {
    warning(sprintf("replicate_success: metric '%s' is not in person_results; it is skipped.", m), call. = FALSE)
  }
  metrics <- setdiff(metrics, missing_metrics)
  if (length(metrics) == 0L) {
    stop("None of the metrics in taus are present in person_results.", call. = FALSE)
  }

  replicate_pass <- NULL
  for (m in metrics) {
    stat_df <- replicate_pooled_error(person_results, m)
    stat_df[[paste0("pass_", m)]] <- stat_df$stat <= taus[[m]]
    stat_df$stat <- NULL
    replicate_pass <- if (is.null(replicate_pass)) {
      stat_df
    } else {
      merge(replicate_pass, stat_df, by = c("scenario_id", "replicate"))
    }
  }
  pass_cols <- paste0("pass_", metrics)
  replicate_pass$pass <- Reduce(`&`, replicate_pass[pass_cols])

  replicate_pass <- replicate_pass[order(replicate_pass$scenario_id, replicate_pass$replicate), , drop = FALSE]
  rownames(replicate_pass) <- NULL
  replicate_pass[, c("scenario_id", "replicate", pass_cols, "pass")]
}


#' Extract the success rate per scenario from a Dirichlet-multinomial experiment.
#'
#' A replicate succeeds for a metric if the largest cell-type error of the person-pooled proportion estimate against
#' the population-level true proportion is `<= tau` (see `replicate_pooled_error()`). A
#' replicate succeeds jointly if it succeeds for every metric in `taus`. The success rate is the fraction of
#' successful replicates. Success is determined by `replicate_success()`, the single source of truth for this rule.
#'
#' @param result Output of `run_dirichlet_multinomial_experiment()` (or `run_simulation_experiment()` with
#'   `model = "dirichlet_multinomial"`), in memory or read back with `readRDS()`.
#' @param taus   Named list with one scalar threshold per metric (e.g. `list(AE = 0.02, ARE = 0.5)`). Metrics missing
#'   from `result$person_results` are skipped with a warning.
#'
#' @return Data.frame with one row per scenario and columns: scenario_id, alpha, p_max, n_people, concentration, B,
#'   success_count, success_rate (joint over all metrics), and `success_rate_<metric>` per metric.
extract_success_rate <- function(result, taus) {
  person_results <- result$person_results
  if (!is.data.frame(person_results)) {
    stop("result must contain a person_results data.frame (Dirichlet-multinomial experiment output).", call. = FALSE)
  }

  replicate_pass <- replicate_success(person_results, taus)
  metrics <- sub("^pass_", "", setdiff(names(replicate_pass), c("scenario_id", "replicate", "pass")))

  scenario_ids <- unique(replicate_pass$scenario_id)
  summary_rows <- lapply(scenario_ids, function(id) {
    rows <- replicate_pass[replicate_pass$scenario_id == id, , drop = FALSE]
    out <- data.frame(
      scenario_id = id,
      B = nrow(rows),
      success_count = sum(rows$pass),
      success_rate = mean(rows$pass),
      stringsAsFactors = FALSE
    )
    for (m in metrics) {
      out[[paste0("success_rate_", m)]] <- mean(rows[[paste0("pass_", m)]])
    }
    out
  })
  summary <- do.call(rbind, summary_rows)

  scenario_cols <- c("scenario_id", "alpha", "p_max", "n_people", "concentration")
  out <- merge(result$p_table[, scenario_cols], summary, by = "scenario_id", sort = FALSE)
  out[order(as.integer(sub("^scenario_", "", out$scenario_id))), , drop = FALSE]
}
