
# Orchestration ---------------------------------------------------------------------------------------------------

## Hybrid Cutoff Simulation ----------------------------------------------------------------------------------------

#' Run simulation and hybrid cutoff analysis across all alpha/p_max scenarios.
#'
#' @param alpha             Numeric vector of positive Beta shape values.
#' @param K                 Number of cell types.
#' @param n                 Total sample size per replicate.
#' @param B                 Number of replicates.
#' @param cutoffs           Numeric vector of candidate cutoffs.
#' @param tau_AE            Numeric scalar AE success threshold.
#' @param tau_ARE           Numeric scalar ARE success threshold.
#' @param proportion_method Proportion-generation method (`"beta"` or
#'   `"fixed_max_beta"`).
#' @param p_max             Optional fixed largest true proportion(s).
#' @param model             Sampling model passed to `run_simulation_experiment()`.
#' @param tie_method        Tie-breaking rule for max-error argmax.
#' @param maximize          Which success rate to maximize (`"cell"` or
#'   `"replicate"`).
#' @param seed              Optional integer seed for reproducibility.
#' @param ...               Additional arguments forwarded to
#'   `run_simulation_experiment()`.
#'
#' @return List with elements:
#'   \describe{
#'     \item{inputs}{Copy of all function inputs.}
#'     \item{p_table}{Scenario-level true proportions table.}
#'     \item{phat_long}{Tidy replicate-level observed proportions.}
#'     \item{cutoff_curves}{Combined hybrid cutoff curves across scenarios.}
#'     \item{best_cutoff_summary}{One-row-per-scenario best-cutoff summary.}
#'   }
run_simulation_hybrid_cutoff_experiment <- function(alpha, K, n, B, cutoffs,
                                                    tau_AE, tau_ARE,
                                                    proportion_method = "beta",
                                                    p_max = NULL,
                                                    model = "multinomial",
                                                    tie_method = "random",
                                                    maximize = c("cell", "replicate"),
                                                    seed = NULL, ...) {
  maximize <- match.arg(maximize)
  sim_out <- run_simulation_experiment(
    alpha = alpha,
    K = K,
    n = n,
    B = B,
    taus = list(AE = tau_AE, ARE = tau_ARE),
    metrics = c("AE", "ARE"),
    proportion_method = proportion_method,
    p_max = p_max,
    model = model,
    tie_method = tie_method,
    seed = seed,
    ...
  )
  
  p_table <- sim_out$p_table
  index_cols <- grep("^index_", names(p_table), value = TRUE)
  cutoff_curves_list <- vector("list", nrow(p_table))
  best_summary_list <- vector("list", nrow(p_table))
  
  for (i in seq_len(nrow(p_table))) {
    alpha_i <- p_table$alpha[[i]]
    p_max_i <- if ("p_max" %in% names(p_table)) p_table$p_max[[i]] else NA_real_
    p_i <- as.numeric(p_table[i, index_cols, drop = FALSE])
    
    same_p_max <- if (is.na(p_max_i)) {
      is.na(sim_out$phat_long$p_max)
    } else {
      sim_out$phat_long$p_max == p_max_i
    }
    phat_subset <- sim_out$phat_long[
      sim_out$phat_long$alpha == alpha_i & same_p_max,
      c("replicate", "index", "phat"),
      drop = FALSE
    ]
    phat_subset <- phat_subset[order(phat_subset$index, phat_subset$replicate), , drop = FALSE]
    phat_mat_i <- matrix(phat_subset$phat, nrow = B, ncol = length(p_i))
    
    hybrid_i <- run_hybrid_cutoff_analysis(
      phat_mat = phat_mat_i,
      p = p_i,
      cutoffs = cutoffs,
      tau_AE = tau_AE,
      tau_ARE = tau_ARE,
      maximize = maximize,
      alpha = alpha_i,
      p_max = p_max_i
    )
    
    curve_i <- hybrid_i$cutoff_curve
    curve_i$alpha <- alpha_i
    curve_i$p_max <- p_max_i
    curve_i$tau_AE <- tau_AE
    curve_i$tau_ARE <- tau_ARE
    cutoff_curves_list[[i]] <- curve_i[, c(
      "alpha", "p_max", "cutoff", "tau_AE", "tau_ARE",
      "success_rate_cell", "success_rate_replicate",
      "prop_using_AE", "prop_using_ARE"
    )]
    
    best_summary_list[[i]] <- hybrid_i$best_summary
  }
  
  list(
    inputs = list(
      alpha = alpha, K = K, n = n, B = B, cutoffs = cutoffs,
      tau_AE = tau_AE, tau_ARE = tau_ARE,
      proportion_method = proportion_method, p_max = p_max,
      model = model, tie_method = tie_method, maximize = maximize, seed = seed
    ),
    p_table = sim_out$p_table,
    phat_long = sim_out$phat_long,
    cutoff_curves = do.call(rbind, cutoff_curves_list),
    best_cutoff_summary = do.call(rbind, best_summary_list)
  )
}



# Calculation -----------------------------------------------------------------------------------------------------

## Hybrid Thresholds -----------------------------------------------------------------------------------------------

#' Evaluate hybrid AE/ARE success at the cell-type level for one cutoff.
#'
#' @param phat_mat B x K numeric matrix of observed proportions.
#' @param p        True proportion vector of length K.
#' @param cutoff   Numeric scalar cutoff on observed proportions.
#' @param tau_AE   Numeric scalar AE success threshold.
#' @param tau_ARE  Numeric scalar ARE success threshold.
#'
#' @return List with elements:
#'   \describe{
#'     \item{success_matrix}{B x K logical matrix of hybrid success indicators.}
#'     \item{use_AE_matrix}{B x K logical matrix; TRUE when AE is selected.}
#'     \item{success_rate_cell}{Mean success across all replicate-index pairs.}
#'     \item{success_rate_replicate}{Fraction of replicates where all K indices succeed.}
#'   }
evaluate_hybrid_success_cell_level <- function(phat_mat, p, cutoff, tau_AE, tau_ARE) {
  stopifnot(
    is.matrix(phat_mat),
    is.numeric(phat_mat),
    is.numeric(p),
    length(p) == ncol(phat_mat),
    all(is.finite(p)),
    all(p > 0),
    is.numeric(cutoff),
    length(cutoff) == 1L,
    is.finite(cutoff),
    is.numeric(tau_AE),
    length(tau_AE) == 1L,
    is.finite(tau_AE),
    is.numeric(tau_ARE),
    length(tau_ARE) == 1L,
    is.finite(tau_ARE)
  )
  
  p_mat <- matrix(p, nrow = nrow(phat_mat), ncol = ncol(phat_mat), byrow = TRUE)
  ae_mat <- abs(phat_mat - p_mat)
  are_mat <- ae_mat / p_mat
  use_AE_matrix <- phat_mat <= cutoff
  success_matrix <- ifelse(use_AE_matrix, ae_mat <= tau_AE, are_mat <= tau_ARE)
  
  list(
    success_matrix = success_matrix,
    use_AE_matrix = use_AE_matrix,
    success_rate_cell = mean(success_matrix, na.rm = TRUE),
    success_rate_replicate = mean(rowSums(success_matrix, na.rm = TRUE) == ncol(success_matrix))
  )
}

#' Sweep hybrid AE/ARE cutoffs and summarize success rates.
#'
#' @param phat_mat B x K numeric matrix of observed proportions.
#' @param p        True proportion vector of length K.
#' @param cutoffs  Numeric vector of candidate cutoffs.
#' @param tau_AE   Numeric scalar AE success threshold.
#' @param tau_ARE  Numeric scalar ARE success threshold.
#'
#' @return Tidy data.frame with one row per cutoff and columns:
#'   cutoff, success_rate_cell, success_rate_replicate,
#'   prop_using_AE, prop_using_ARE.
sweep_hybrid_cutoffs_cell_level <- function(phat_mat, p, cutoffs, tau_AE, tau_ARE) {
  stopifnot(is.numeric(cutoffs), length(cutoffs) >= 1L, all(is.finite(cutoffs)))
  
  rows <- lapply(cutoffs, function(cutoff_i) {
    eval_i <- evaluate_hybrid_success_cell_level(
      phat_mat = phat_mat,
      p = p,
      cutoff = cutoff_i,
      tau_AE = tau_AE,
      tau_ARE = tau_ARE
    )
    prop_using_AE_i <- mean(eval_i$use_AE_matrix, na.rm = TRUE)
    data.frame(
      cutoff = cutoff_i,
      success_rate_cell = eval_i$success_rate_cell,
      success_rate_replicate = eval_i$success_rate_replicate,
      prop_using_AE = prop_using_AE_i,
      prop_using_ARE = 1 - prop_using_AE_i,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, rows)
}

#' Find the best hybrid cutoff from a sweep curve.
#'
#' @param curve_df  Data.frame returned by `sweep_hybrid_cutoffs_cell_level()`.
#' @param maximize  Which success rate to maximize: `"cell"` or `"replicate"`.
#' @param tie_break Tie-breaker among maximizing cutoffs:
#'   `"smallest"`, `"largest"`, or `"median"`.
#'
#' @return List with elements:
#'   best_cutoff, best_success_rate_cell, best_success_rate_replicate,
#'   best_prop_using_AE, n_tied_maximizers.
find_best_hybrid_cutoff <- function(curve_df,
                                    maximize = c("cell", "replicate"),
                                    tie_break = c("smallest", "largest", "median")) {
  maximize <- match.arg(maximize)
  tie_break <- match.arg(tie_break)
  
  required_cols <- c(
    "cutoff", "success_rate_cell", "success_rate_replicate",
    "prop_using_AE", "prop_using_ARE"
  )
  if (!all(required_cols %in% names(curve_df))) {
    stop(
      sprintf(
        "curve_df must contain columns: %s",
        paste(required_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (nrow(curve_df) == 0L) {
    stop("curve_df must contain at least one row.", call. = FALSE)
  }
  
  target_col <- if (identical(maximize, "cell")) "success_rate_cell" else "success_rate_replicate"
  target_values <- curve_df[[target_col]]
  best_value <- max(target_values, na.rm = TRUE)
  tied_idx <- which(target_values == best_value)
  tied_cutoffs <- curve_df$cutoff[tied_idx]
  
  best_cutoff <- switch(
    tie_break,
    smallest = min(tied_cutoffs, na.rm = TRUE),
    largest = max(tied_cutoffs, na.rm = TRUE),
    median = stats::median(tied_cutoffs, na.rm = TRUE)
  )
  
  chosen_row <- curve_df[curve_df$cutoff == best_cutoff, , drop = FALSE][1L, , drop = FALSE]
  
  list(
    best_cutoff = best_cutoff,
    best_success_rate_cell = chosen_row$success_rate_cell[[1L]],
    best_success_rate_replicate = chosen_row$success_rate_replicate[[1L]],
    best_prop_using_AE = chosen_row$prop_using_AE[[1L]],
    n_tied_maximizers = length(tied_idx)
  )
}

#' Run hybrid cutoff analysis for one simulation scenario.
#'
#' @param phat_mat  B x K numeric matrix of observed proportions.
#' @param p         True proportion vector of length K.
#' @param cutoffs   Numeric vector of candidate cutoffs.
#' @param tau_AE    Numeric scalar AE success threshold.
#' @param tau_ARE   Numeric scalar ARE success threshold.
#' @param maximize  Which success rate to maximize: `"cell"` or `"replicate"`.
#' @param alpha     Optional scalar metadata carried into outputs.
#' @param p_max     Optional scalar metadata carried into outputs.
#'
#' @return List with two elements:
#'   \describe{
#'     \item{cutoff_curve}{Data.frame from `sweep_hybrid_cutoffs_cell_level()`.}
#'     \item{best_summary}{One-row data.frame with the selected best cutoff summary.}
#'   }
run_hybrid_cutoff_analysis <- function(phat_mat, p, cutoffs, tau_AE, tau_ARE,
                                       maximize = c("cell", "replicate"),
                                       alpha = NULL, p_max = NULL) {
  maximize <- match.arg(maximize)
  cutoff_curve <- sweep_hybrid_cutoffs_cell_level(
    phat_mat = phat_mat,
    p = p,
    cutoffs = cutoffs,
    tau_AE = tau_AE,
    tau_ARE = tau_ARE
  )
  best <- find_best_hybrid_cutoff(
    curve_df = cutoff_curve,
    maximize = maximize,
    tie_break = "smallest"
  )
  
  best_summary <- data.frame(
    alpha = if (is.null(alpha)) NA_real_ else alpha,
    p_max = if (is.null(p_max)) NA_real_ else p_max,
    maximize = maximize,
    tau_AE = tau_AE,
    tau_ARE = tau_ARE,
    best_cutoff = best$best_cutoff,
    best_success_rate_cell = best$best_success_rate_cell,
    best_success_rate_replicate = best$best_success_rate_replicate,
    best_prop_using_AE = best$best_prop_using_AE,
    n_tied_maximizers = best$n_tied_maximizers,
    stringsAsFactors = FALSE
  )
  
  list(
    cutoff_curve = cutoff_curve,
    best_summary = best_summary
  )
}



