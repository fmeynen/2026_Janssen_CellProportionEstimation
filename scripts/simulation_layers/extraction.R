# R/extraction.R
#
# Extraction layer: pull structured outputs from simulation results.
#
# Depends on: (none from this package; uses base R only)
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Iso-probability pipeline
# ---------------------------------------------------------------------------


#' Summarize iso-search yield, ranges, and closest candidates.
#'
#' @param iso_result Output list from `run_iso_success_search()`.
#' @param top_k      Integer; number of closest-to-baseline candidates to return
#'   (ranked by smallest `abs_diff_p0`; default 10L).
#'
#' @return Named list with elements:
#'   \describe{
#'     \item{counts}{One-row data.frame: n_candidates_total, n_pass_screen.}
#'     \item{ranges}{One-row data.frame: tau_AE_min, tau_AE_max, tau_ARE_min,
#'       tau_ARE_max, cutoff_min, cutoff_max among passing candidates.
#'       All columns are NA when no candidates pass.}
#'     \item{closest}{data.frame of up to top_k rows from screen_results ranked
#'       by ascending abs_diff_p0: candidate_id, tau_AE, tau_ARE, cutoff,
#'       p_hat, abs_diff_p0, ci_low, ci_high.}
#'   }
summarize_iso_candidates <- function(iso_result, top_k = 10L) {
  stopifnot(
    is.list(iso_result),
    all(c("screen_results", "diagnostics") %in% names(iso_result)),
    is.data.frame(iso_result$screen_results),
    is.numeric(top_k), length(top_k) == 1L, top_k >= 1L
  )

  sr  <- iso_result$screen_results
  top_k <- as.integer(top_k)

  counts <- data.frame(
    n_candidates_total = iso_result$diagnostics$n_candidates_total,
    n_pass_screen      = iso_result$diagnostics$n_pass_screen,
    stringsAsFactors   = FALSE
  )

  pass_rows <- sr[sr$pass, , drop = FALSE]
  if (nrow(pass_rows) == 0L) {
    ranges <- data.frame(
      tau_AE_min  = NA_real_,
      tau_AE_max  = NA_real_,
      tau_ARE_min = NA_real_,
      tau_ARE_max = NA_real_,
      cutoff_min  = NA_real_,
      cutoff_max  = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    ranges <- data.frame(
      tau_AE_min  = min(pass_rows$tau_AE),
      tau_AE_max  = max(pass_rows$tau_AE),
      tau_ARE_min = min(pass_rows$tau_ARE),
      tau_ARE_max = max(pass_rows$tau_ARE),
      cutoff_min  = min(pass_rows$cutoff),
      cutoff_max  = max(pass_rows$cutoff),
      stringsAsFactors = FALSE
    )
  }

  sr_sorted <- sr[order(sr$abs_diff_p0), , drop = FALSE]
  k         <- min(top_k, nrow(sr_sorted))
  closest   <- sr_sorted[seq_len(k),
                         c("candidate_id", "tau_AE", "tau_ARE", "cutoff",
                           "p_hat", "abs_diff_p0", "ci_low", "ci_high"),
                         drop = FALSE]
  rownames(closest) <- NULL

  list(
    counts  = counts,
    ranges  = ranges,
    closest = closest
  )
}
