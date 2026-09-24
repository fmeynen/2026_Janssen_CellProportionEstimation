
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
  stopifnot(is.matrix(argmax), !is.null(colnames(argmax)))
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
#'
#' @return A single-row data.frame with columns: alpha, p_max, index_1, ..., index_K.
extract_p_table_row <- function(alpha_i, p_max_i, p, K) {
  data.frame(
    alpha = alpha_i,
    p_max = p_max_i,
    as.list(stats::setNames(as.numeric(p), paste0("index_", seq_len(K)))),
    stringsAsFactors = FALSE,
    check.names = FALSE
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


#' Compute the per-replicate max-mean error for one metric of a Dirichlet-multinomial experiment.
#'
#' For each scenario and replicate, the error is averaged over all persons per cell type, and the largest of these
#' cell-type means is returned. `NaN` errors (ARE with true and observed proportion both 0) count as 0; `Inf` errors are
#' kept, so the corresponding mean is `Inf`.
#'
#' @param person_results `person_results` data.frame from `run_dirichlet_multinomial_experiment()`.
#' @param metric         Metric name present in `person_results$metric`.
#'
#' @return Data.frame with columns: scenario_id, replicate, max_mean_error.
dm_max_mean_errors <- function(person_results, metric) {
  rows <- person_results[person_results$metric == metric, c("scenario_id", "replicate", "cell_type", "error")]
  rows$error[is.nan(rows$error)] <- 0
  cell_means <- stats::aggregate(error ~ scenario_id + replicate + cell_type, data = rows, FUN = mean,
                                 na.action = stats::na.pass)
  out <- stats::aggregate(error ~ scenario_id + replicate, data = cell_means, FUN = max,
                          na.action = stats::na.pass)
  names(out)[names(out) == "error"] <- "max_mean_error"
  out
}


#' Extract the success rate per scenario from a Dirichlet-multinomial experiment.
#'
#' A replicate succeeds for a metric if the largest cell-type error, averaged over all persons, is `<= tau`. A replicate
#' succeeds jointly if it succeeds for every metric in `taus`. The success rate is the fraction of successful replicates.
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
    warning(sprintf("extract_success_rate: metric '%s' is not in person_results; it is skipped.", m), call. = FALSE)
  }
  metrics <- setdiff(metrics, missing_metrics)
  if (length(metrics) == 0L) {
    stop("None of the metrics in taus are present in person_results.", call. = FALSE)
  }

  replicate_pass <- NULL
  for (m in metrics) {
    max_means <- dm_max_mean_errors(person_results, m)
    max_means[[paste0("pass_", m)]] <- max_means$max_mean_error <= taus[[m]]
    max_means$max_mean_error <- NULL
    replicate_pass <- if (is.null(replicate_pass)) {
      max_means
    } else {
      merge(replicate_pass, max_means, by = c("scenario_id", "replicate"))
    }
  }
  pass_cols <- paste0("pass_", metrics)
  replicate_pass$pass <- Reduce(`&`, replicate_pass[pass_cols])

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
