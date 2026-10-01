
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
#' @return Tidy data.frame with columns: metric, cell_type, count, fraction, true_proportion.
summarize_argmax <- function(argmax, p) {
  validate_named_matrix(argmax, "argmax")
  K <- length(p)
  metrics <- colnames(argmax)
  rows <- vector("list", length(metrics))
  for (i in seq_along(metrics)) {
    m <- metrics[[i]]
    counts <- tabulate(argmax[, m], nbins = K)
    rows[[i]] <- data.frame(
      metric          = m,
      cell_type       = seq_len(K),
      count           = counts,
      fraction        = counts / sum(counts),
      true_proportion = p,
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
#' @param after    Optional named list of extra columns placed between `p_max` and `cell_type_1`.
#'
#' @return A single-row data.frame with columns: [before], alpha, p_max, [after], cell_type_1, ..., cell_type_K.
extract_p_table_row <- function(alpha_i, p_max_i, p, K, before = NULL, after = NULL) {
  do.call(
    data.frame,
    c(
      before,
      list(alpha = alpha_i, p_max = p_max_i),
      after,
      as.list(stats::setNames(as.numeric(p), paste0("cell_type_", seq_len(K)))),
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
#'   alpha, p_max, replicate, metric, cell_type, error.
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
      cell_type = rep(seq_len(ncol(errors_m)), each = B),
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
#'   alpha, p_max, replicate, cell_type, phat.
extract_phat_long <- function(rep_out, alpha_i, p_max_i, B) {
  data.frame(
    alpha = alpha_i,
    p_max = p_max_i,
    replicate = rep(seq_len(B), times = ncol(rep_out$phat)),
    cell_type = rep(seq_len(ncol(rep_out$phat)), each = B),
    phat = as.vector(rep_out$phat),
    stringsAsFactors = FALSE
  )
}


# Dirichlet-multinomial success rate ------------------------------------------------------------------------------


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

  # One grouped pass: per-scenario sums of the pass columns (first-appearance order) divided by the replicate counts.
  g <- factor(replicate_pass$scenario_id, levels = unique(replicate_pass$scenario_id))
  pass_cols <- c("pass", paste0("pass_", metrics))
  sums <- rowsum(as.matrix(replicate_pass[pass_cols]) * 1, g, reorder = FALSE)
  B <- tabulate(g, nlevels(g))
  summary <- data.frame(
    scenario_id = levels(g),
    B = B,
    success_count = as.integer(sums[, "pass"]),
    success_rate = sums[, "pass"] / B,
    stringsAsFactors = FALSE
  )
  for (m in metrics) {
    summary[[paste0("success_rate_", m)]] <- sums[, paste0("pass_", m)] / B
  }

  scenario_cols <- c("scenario_id", "alpha", "p_max", "n_people", "concentration")
  out <- merge(result$p_table[, scenario_cols], summary, by = "scenario_id", sort = FALSE)
  out[order(as.integer(sub("^scenario_", "", out$scenario_id))), , drop = FALSE]
}
