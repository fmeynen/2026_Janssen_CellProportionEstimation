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

#' Validate and normalize isoband triad ranges.
#'
#' @param ranges Named list with numeric length-2 vectors:
#'   `tau_AE`, `tau_ARE`, and `cutoff`.
#'
#' @return Normalized list with elements `tau_AE`, `tau_ARE`, and `cutoff`.
validate_isoband_ranges <- function(ranges) {
  required <- c("tau_AE", "tau_ARE", "cutoff")
  if (!is.list(ranges) || !all(required %in% names(ranges))) {
    stop(
      sprintf("ranges must be a named list containing: %s", paste(required, collapse = ", ")),
      call. = FALSE
    )
  }

  normalize_pair <- function(x, name) {
    if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x))) {
      stop(sprintf("ranges$%s must be a numeric length-2 finite vector.", name), call. = FALSE)
    }
    x <- as.numeric(x)
    if (!(x[[1L]] < x[[2L]])) {
      stop(sprintf("ranges$%s must be strictly increasing.", name), call. = FALSE)
    }
    x
  }

  tau_AE <- normalize_pair(ranges$tau_AE, "tau_AE")
  tau_ARE <- normalize_pair(ranges$tau_ARE, "tau_ARE")
  cutoff <- normalize_pair(ranges$cutoff, "cutoff")

  if (!(tau_AE[[1L]] > 0 && tau_AE[[2L]] < 1)) {
    stop("tau_AE range must satisfy 0 < min < max < 1.", call. = FALSE)
  }
  if (!(tau_ARE[[1L]] > 0 && is.finite(tau_ARE[[2L]]) && tau_ARE[[2L]] < Inf)) {
    stop("tau_ARE range must satisfy 0 < min < max < Inf with finite bounds.", call. = FALSE)
  }
  if (!(cutoff[[1L]] > 0 && cutoff[[2L]] < 1)) {
    stop("cutoff range must satisfy 0 < min < max < 1.", call. = FALSE)
  }

  list(tau_AE = tau_AE, tau_ARE = tau_ARE, cutoff = cutoff)
}

#' Create an empty isoband result container.
#'
#' @param inputs Optional list of pipeline inputs.
#'
#' @return Standardized list skeleton for isoband outputs.
new_isoband_result_container <- function(inputs = NULL) {
  list(
    design_history = data.frame(),
    round_summaries = data.frame(),
    final_band_points = data.frame(),
    final_model = NULL,
    inputs = inputs
  )
}
