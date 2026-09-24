# Orchestration ---------------------------------------------------------------------------------------------------
# Orchestration layer: high-level experiment runners that coordinate all other layers.
#
# Depends on: simulation.R, calculation.R, extraction.R

# Simulation - I/O -------------------------------------------------------------------------------------------------

# I/O helpers for persisting simulation results.
#
# Each unique combination of simulation parameters is saved as a separate .rds file whose name contains a short MD5 hash
# of those parameters.

#' Compute an MD5 hash of an R object.
#'
#' Serialises `x` to a temporary file, computes its MD5 checksum, and returns the hex string.
#' Only `tools` (part of base R) is required.
#'
#' @param x Any R object.
#' @return A 32-character hexadecimal string.
hash_config <- function(x) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  saveRDS(x, tmp)
  unname(tools::md5sum(tmp))
}


#' Build the file path for a simulation result.
#'
#' The file name is `<name>_<md5>.rds`, where the MD5 is derived from the serialised `config`.
#' Identical parameters always resolve to the same path; any change produces a different path.
#'
#' @param config  Named list of simulation parameters.
#' @param dir     Directory that will hold the result files.  Created automatically if it does not yet exist.
#' @param name    Short label prefixed to the file name
#'   (e.g. `"simulation_errorchoice"`).
#'
#' @return Absolute path string.
simulation_result_path <- function(config, dir, name) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(name, "_", hash_config(config), ".rds"))
}

# Run individual experiments --------------------------------------------------------------------------------------
run_dirichlet_multinomial_experiment <- function(
  alpha,
  K,
  B,
  metrics,
  proportion_method = "beta",
  p_max = NULL,
  n_people,
  n_per_person,
  concentration,
  seed
) {
  n_people <- validate_positive_integer(
    n_people,
    "n_people",
    allow_vector = TRUE
  )
  n_per_person <- validate_positive_integer(n_per_person, "n_per_person")
  concentration <- validate_positive_numeric(
    concentration,
    "concentration",
    allow_vector = TRUE
  )

  if (identical(proportion_method, "fixed_max_beta")) {
    if (is.null(p_max)) {
      stop(
        "p_max must be provided when proportion_method = 'fixed_max_beta'.",
        call. = FALSE
      )
    }
    if (
      !is.numeric(p_max) ||
        any(!is.finite(p_max)) ||
        any(p_max <= 0) ||
        any(p_max >= 1)
    ) {
      stop(
        "p_max must contain numbers strictly between 0 and 1.",
        call. = FALSE
      )
    }
    p_max_values <- as.numeric(p_max)
  } else {
    p_max_values <- NA_real_
  }

  population_scenarios <- expand.grid(
    alpha = alpha,
    p_max = p_max_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  person_scenarios <- expand.grid(
    n_people = n_people,
    concentration = concentration,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  scenario_grid <- merge(
    population_scenarios,
    person_scenarios,
    by = NULL,
    sort = FALSE
  )
  scenario_grid$scenario_id <- paste0("scenario_", seq_len(nrow(scenario_grid)))

  p_table_list <- vector("list", nrow(scenario_grid))
  person_results_list <- vector("list", nrow(scenario_grid))
  keep <- logical(nrow(scenario_grid))
  population_compositions <- list()
  impossible_population_keys <- character()
  skip_impossible <- identical(proportion_method, "fixed_max_beta") &&
    length(p_max_values) > 1L

  for (i in seq_len(nrow(scenario_grid))) {
    scenario <- scenario_grid[i, , drop = FALSE]
    population_key <- paste(
      format(scenario$alpha[[1L]], scientific = FALSE, trim = TRUE),
      if (is.na(scenario$p_max[[1L]])) {
        "NA"
      } else {
        format(scenario$p_max[[1L]], scientific = FALSE, trim = TRUE)
      },
      sep = "__"
    )
    if (population_key %in% impossible_population_keys) {
      next
    }
    p <- population_compositions[[population_key]]
    if (is.null(p)) {
      p <- tryCatch(
        generate_proportions(
          alpha = scenario$alpha[[1L]],
          K = K,
          method = proportion_method,
          p_max = if (is.na(scenario$p_max[[1L]])) {
            NULL
          } else {
            scenario$p_max[[1L]]
          }
        ),
        error = function(e) {
          if (skip_impossible && inherits(e, "impossible_fixed_max_error")) {
            return(NULL)
          }
          stop(e)
        }
      )
      if (!is.null(p)) {
        population_compositions[[population_key]] <- p
      }
    }
    if (is.null(p)) {
      impossible_population_keys <- c(
        impossible_population_keys,
        population_key
      )
    }
    if (is.null(p)) {
      next
    }

    rep_out <- run_replicates(
      p = p,
      B = B,
      metrics = metrics,
      model = "dirichlet_multinomial",
      seed = if (is.null(seed)) NULL else seed + i - 1L,
      n_people = scenario$n_people[[1L]],
      n_per_person = n_per_person,
      concentration = scenario$concentration[[1L]],
      scenario_id = scenario$scenario_id[[1L]]
    )
    person_results <- rep_out$person_results
    person_results$alpha <- scenario$alpha[[1L]]
    person_results$p_max <- scenario$p_max[[1L]]
    person_results_list[[i]] <- person_results[, c(
      "scenario_id",
      "alpha",
      "p_max",
      "n_people",
      "concentration",
      "replicate",
      "person_id",
      "cell_type",
      "metric",
      "count",
      "observed_proportion",
      "person_true_proportion",
      "population_mean_proportion",
      "error"
    )]
    p_table_list[[i]] <- data.frame(
      scenario_id = scenario$scenario_id[[1L]],
      alpha = scenario$alpha[[1L]],
      p_max = scenario$p_max[[1L]],
      n_people = scenario$n_people[[1L]],
      concentration = scenario$concentration[[1L]],
      as.list(stats::setNames(as.numeric(p), paste0("index_", seq_len(K)))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    keep[[i]] <- TRUE
  }

  if (!any(keep)) {
    stop(
      "No feasible alpha/p_max combinations produced simulation output.",
      call. = FALSE
    )
  }

  list(
    inputs = list(
      alpha = alpha,
      K = K,
      B = B,
      metrics = metrics,
      proportion_method = proportion_method,
      p_max = p_max,
      model = "dirichlet_multinomial",
      n_people = n_people,
      n_per_person = n_per_person,
      concentration = concentration,
      seed = seed
    ),
    p_table = do.call(rbind, p_table_list[keep]),
    person_results = do.call(rbind, person_results_list[keep])
  )
}


# Original Simulation ---------------------------------------------------------------------------------------------

#' Run the original full simulation experiment end-to-end.
#'
#' @param alpha      Numeric vector; one or more positive shape values used by the selected proportion-generation method
#'      (default method is Beta-based).
#' @param K          Number of cell types (default 10).
#' @param n          Total sample size per replicate for the multinomial
#'   model.
#' @param B          Number of replicates.
#' @param taus       Numeric vector of thresholds (same for all metrics) or a named list with one numeric vector per
#'      metric (e.g. `list(AE = c(...), ARE = c(...))`).
#' @param metrics    Error metrics; any subset of `c("AE", "ARE", "TSE", "LAE")`.
#' @param proportion_method Proportion-generation method (`"beta"` or `"fixed_max_beta"`).
#'    The fixed-max Beta method places `p_max` at the highest index and warns then fails for impossible combinations.
#' @param p_max      Fixed largest true proportion(s) used by `proportion_method = "fixed_max_beta"`.
#'      Can be a numeric vector.
#'      When multiple values are provided, all alpha × p_max combinations are attempted; impossible fixed-max
#'      combinations are warned and skipped.
#' @param model      Sampling model (`"multinomial"` or
#'   `"dirichlet_multinomial"`).
#' @param n_people   Positive integer vector of people per replicate for the
#'   Dirichlet-multinomial model.
#' @param n_per_person Positive integer sampled-cell count for every person
#'   in the Dirichlet-multinomial model.
#' @param concentration Positive numeric vector of Dirichlet concentration
#'   values for the Dirichlet-multinomial model.
#' @param tie_method Tie-breaking rule for max-error argmax.
#' @param seed       Optional integer seed for reproducibility.
#' @param ...        Additional arguments forwarded to `simulate_counts()`.
#'
#' @return List with elements:
#'   \describe{
#'     \item{inputs}{All input arguments.}
#'     \item{p_table}{Data.frame with one row per simulated alpha/p_max combination,
#'       an `alpha` column, a `p_max` column, and one column per index
#'       (`index_1`, ..., `index_K`) containing the corresponding p values.}
#'     \item{replicate_summaries}{Tidy data.frame:
#'       alpha, p_max, replicate, metric, max_error, argmax_index.}
#'     \item{errors_long}{Tidy data.frame:
#'       alpha, p_max, replicate, metric, index, error.}
#'     \item{phat_long}{Tidy data.frame:
#'       alpha, p_max, replicate, index, phat.}
#'     \item{curves}{Tidy data.frame:
#'       alpha, p_max, metric, tau, success_rate, mean_n_above.}
#'     \item{argmax_summary}{Tidy data.frame:
#'       alpha, p_max, metric, index, count, fraction, p_value.}
#'   }
run_simulation_experiment <- function(
  alpha,
  K = 10,
  n = NULL,
  B,
  taus,
  metrics = c("AE", "ARE"),
  proportion_method = "beta",
  p_max = NULL,
  model = "multinomial",
  tie_method = "random",
  seed = NULL,
  n_people = NULL,
  n_per_person = NULL,
  concentration = NULL,
  ...
) {
  stopifnot(is.numeric(alpha), length(alpha) >= 1L, all(alpha > 0))
  model <- match.arg(model, c("multinomial", "dirichlet_multinomial"))

  if (identical(model, "dirichlet_multinomial")) {
    return(run_dirichlet_multinomial_experiment(
      alpha = alpha,
      K = K,
      B = B,
      metrics = metrics,
      proportion_method = proportion_method,
      p_max = p_max,
      n_people = n_people,
      n_per_person = n_per_person,
      concentration = concentration,
      seed = seed
    ))
  }
  if (is.null(n)) {
    stop("n must be provided for model = 'multinomial'.", call. = FALSE)
  }

  if (identical(proportion_method, "fixed_max_beta")) {
    if (is.null(p_max)) {
      stop(
        "p_max must be provided when proportion_method = 'fixed_max_beta'.",
        call. = FALSE
      )
    }
    if (
      !is.numeric(p_max) ||
        any(!is.finite(p_max)) ||
        any(p_max <= 0) ||
        any(p_max >= 1)
    ) {
      stop(
        "p_max must contain numbers strictly between 0 and 1.",
        call. = FALSE
      )
    }
    p_max_values <- as.numeric(p_max)
  } else {
    p_max_values <- NA_real_
  }

  combinations <- expand.grid(
    alpha = alpha,
    p_max = p_max_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  n_combinations <- nrow(combinations)
  p_table_list <- vector("list", n_combinations)
  replicate_summaries_list <- vector("list", n_combinations)
  errors_long_list <- vector("list", n_combinations)
  phat_long_list <- vector("list", n_combinations)
  curves_list <- vector("list", n_combinations)
  argmax_summary_list <- vector("list", n_combinations)
  keep <- logical(n_combinations)
  should_skip_impossible_combinations <- identical(
    proportion_method,
    "fixed_max_beta"
  ) &&
    length(p_max_values) > 1L

  for (i in seq_len(n_combinations)) {
    alpha_i <- combinations$alpha[[i]]
    p_max_i <- combinations$p_max[[i]]
    seed_i <- if (is.null(seed)) NULL else seed + i - 1L
    p <- tryCatch(
      generate_proportions(
        alpha = alpha_i,
        K = K,
        method = proportion_method,
        p_max = if (is.na(p_max_i)) NULL else p_max_i
      ),
      error = function(e) {
        if (
          should_skip_impossible_combinations &&
            inherits(e, "impossible_fixed_max_error")
        ) {
          return(NULL)
        }
        stop(e)
      }
    )
    if (is.null(p)) {
      next
    } # Skip impossible alpha/p_max combinations caught by the error handler.

    rep_out <- run_replicates(
      p,
      n,
      B,
      metrics = metrics,
      model = model,
      tie_method = tie_method,
      seed = seed_i,
      ...
    )
    keep[[i]] <- TRUE

    p_table_list[[i]] <- extract_p_table_row(alpha_i, p_max_i, p, K)

    replicate_summaries_list[[i]] <- extract_replicate_summaries(
      rep_out,
      alpha_i,
      p_max_i,
      B,
      metrics
    )

    errors_long_list[[i]] <- extract_errors_long(
      rep_out,
      alpha_i,
      p_max_i,
      B,
      metrics
    )

    phat_long_list[[i]] <- extract_phat_long(rep_out, alpha_i, p_max_i, B)

    curves_i <- evaluate_thresholds(
      rep_out$max_errors,
      taus,
      errors = rep_out$errors
    )
    curves_i$alpha <- alpha_i
    curves_i$p_max <- p_max_i
    curves_list[[i]] <- curves_i[, c(
      "alpha",
      "p_max",
      "metric",
      "tau",
      "success_rate",
      "mean_n_above"
    )]

    argmax_i <- summarize_argmax(rep_out$argmax, p)
    argmax_i$alpha <- alpha_i
    argmax_i$p_max <- p_max_i
    argmax_summary_list[[i]] <- argmax_i[, c(
      "alpha",
      "p_max",
      "metric",
      "index",
      "count",
      "fraction",
      "p_value"
    )]
  }

  if (!any(keep)) {
    stop(
      "No feasible alpha/p_max combinations produced simulation output.",
      call. = FALSE
    )
  }

  p_table_list <- p_table_list[keep]
  replicate_summaries_list <- replicate_summaries_list[keep]
  errors_long_list <- errors_long_list[keep]
  phat_long_list <- phat_long_list[keep]
  curves_list <- curves_list[keep]
  argmax_summary_list <- argmax_summary_list[keep]

  list(
    inputs = list(
      alpha = alpha,
      K = K,
      n = n,
      B = B,
      taus = taus,
      metrics = metrics,
      proportion_method = proportion_method,
      p_max = p_max,
      model = model,
      tie_method = tie_method,
      seed = seed
    ),
    p_table = do.call(rbind, p_table_list),
    replicate_summaries = do.call(rbind, replicate_summaries_list),
    errors_long = do.call(rbind, errors_long_list),
    phat_long = do.call(rbind, phat_long_list),
    curves = do.call(rbind, curves_list),
    argmax_summary = do.call(rbind, argmax_summary_list)
  )
}
