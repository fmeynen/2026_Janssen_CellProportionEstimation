# Backlog

Code-review findings (clarity, efficiency, consistency) from 2026-09-30, on branch `feat/dm-errorchoice`.
They come from reading the code and searching for definitions; the scripts and tests were not run.
Items are grouped by priority. Line numbers refer to the code as of commit `df54645`.
P1 pass done on branch `fix/p1-broken-code` (2026-10-01); all P1 items resolved, later items pruned for the deprecations.

## P1: Broken code

- [x] **Hybrid-heatmap path fails.** `plot_hybrid_best_cutoff_heatmap()`
  ([visualisation.R:482-489](scripts/simulation_layers/visualisation.R#L482-L489)) calls
  `sweep_hybrid_cutoffs_cell_level()` and `find_best_hybrid_cutoff()`, which now exist only in
  [deprecated.R](scripts/deprecated/deprecated.R) (never sourced). This breaks
  [simulation_hybrid_cutoff_heatmaps.R](scripts/deprecated/simulation_hybrid_cutoff_heatmaps.R), its execute script and
  both hybrid reports. Decide: restore the helpers into a layer, or deprecate the driver and the reports.
  Resolved: deprecated, not restored. `plot_hybrid_best_cutoff_heatmap()` moved to `scripts/deprecated/deprecated.R`;
  the driver and execute script moved to `scripts/deprecated/`; both hybrid reports (and PDFs) moved to `reports/deprecated/`.
- [x] **Wrong source path in report.** [simulation_cutoff_heatmaps_report.Rmd:29](reports/simulation_cutoff_heatmaps_report.Rmd#L29)
  sources `scripts/simulation_hybrid_cutoff_heatmaps.R`; the file is in `scripts/simulations/`.
  Resolved: fixed by the move to `reports/deprecated/` (source path corrected there).
- [x] **`p_max` not passed through.** `simulate_success_at_n()`
  ([simulation.R:727-731](scripts/simulation_layers/simulation.R#L727-L731)) doesn't pass `config$p_max` to
  `generate_proportions()`, so `proportion_method = "fixed_max_beta"` always stops with "p_max must be provided".
  Resolved: fixed; `simulate_success_at_n()` now passes `config$p_max`, with a regression test.
- [x] **Duplicate config keys.** [execute_simulation_samplesize.R:38-43](scripts/simulations/execute_simulation_samplesize.R#L38-L43)
  lists `n_people` and `concentration` twice (`2`/`50`, then `NULL`). `$` returns the first value, but both copies go into
  the cache hash. Remove the `NULL` copies and the leftover comment.
  Resolved: fixed; `NULL` copies removed, `2`/`50` kept.
- [x] **Undefined function.** `simulation_success_curve.R:105` calls
  `simulate_success_curve_for_alpha()`, which isn't defined anywhere; the `overwrite` argument (L91) is ignored.
  Fix or deprecate.
  Resolved: deprecated, not fixed; moved to `scripts/deprecated/`.
- [x] **Legacy test script is broken.** `tests/test_simulation.R` calls
  `generate_proportions_fixed_max_beta()` (L80; renamed to `generate_props_fixed_max_beta`) and the deprecated hybrid
  helpers (L400+). Move anything still useful into `tests/testthat/`, then delete the script.
  Resolved: fixed; non-hybrid, uncovered tests ported to `tests/testthat/`, script deleted.
- [x] **Script doesn't parse.** `Simulation_onecelltype_onedonor.R:36-48`
  is missing a `}` and uses `param_grid` and `plot_ate`, which are never defined.
  Resolved: deprecated; this script and `...onedonor2.R` moved unchanged to `scripts/deprecated/`.

## P2: Efficiency

- [ ] **Slim the Dirichlet-multinomial replicate output.**
  - Currently each replicate ([simulation.R:563-577](scripts/simulation_layers/simulation.R#L563-L577)) builds a
    14-column data frame with one row per person × cell type × metric (count and proportion columns repeated per
    metric), and the B frames are combined with `do.call(rbind, …)`.
  - Downstream, only the pooled estimate `colMeans(observed_p)` and `p` are used; the per-person `error` column is no
    longer read.
  - Change: each replicate returns the pooled K-vector; build a B × K matrix; compute each replicate's max error with
    matrix operations; make `person_results` optional.
  - Likely the main cost of `run_dm_errorchoice_experiment` (about 520 scenarios at B = 1000). Profile before and after.
- [ ] **Stop rebuilding data in the multinomial success path.**
  [simulation.R:772-794](scripts/simulation_layers/simulation.R#L772-L794) builds a B × K × M long data frame just to
  call `replicate_success()`. With one person the pooled rule equals `rep_out$max_errors[, m] <= tau`, which is already
  computed.
- [ ] **Parallelism on Windows.** [`replicate_cores()`](scripts/simulation_layers/simulation.R#L328) returns 1 on
  Windows, so all runs are serial. The streams don't depend on the worker, so a PSOCK/`parLapply` backend would keep
  results reproducible. The layer functions need to be sourced on each worker.
- [ ] **Vectorise the Dirichlet draws.** One `rgamma` call over `n_people × K` values followed by `rowSums`, instead of a
  `vapply` over people that re-validates its input on every call.
- [ ] **Replace the per-scenario `lapply` filter** in `extract_success_rate()`
  ([extraction.R:365](scripts/simulation_layers/extraction.R#L365)) with one `rowsum`/`aggregate` pass.

## P3: Structure and consistency

- [ ] **One way to load the layers.** There are currently three:
  - `list.files(pattern = "[.]R$")` in the testthat helper
  - `list.files()` without a pattern, copied into three drivers (would also source non-`.R` files)
  - the loader in the `.qmd`

  Preferred: make `simulation_layers/` a small package (a `DESCRIPTION` file plus an `R/` folder, loaded with
  `devtools::load_all()`). `R CMD check` would then flag undefined functions like those in P1.
- [ ] **One cache helper.** The cache-or-compute pattern is written out three times, with different behaviour:
  - `errorchoice` hashes the whole config, so plot-only changes such as `cutoffs` trigger a re-simulation.
  - `dm_errorchoice` hashes only the simulation fields plus `success_rule_id()`.

  Introduce e.g. `cached_result(key, name, compute, cache, force_recompute, dir)` and pass only simulation-relevant fields
  in the key.
- [ ] **Remove duplicated code:**
  - The RNG save/restore block in `replicate_streams` and `replicate_apply`
    ([simulation.R:369-381](scripts/simulation_layers/simulation.R#L369-L381),
    [:453-465](scripts/simulation_layers/simulation.R#L453-L465)).
  - The `p_max` validation in the generator and both orchestrators
    ([orchestration.R:158](scripts/simulation_layers/orchestration.R#L158),
    [:563](scripts/simulation_layers/orchestration.R#L563)).
  - Inline p_table row construction in two orchestrators; use `extract_p_table_row()` instead.
  - The `required_cols` check in `replicate_pooled_error` and `replicate_success`.
  - The scenario loop that skips impossible combinations, in `run_simulation_experiment` and
    `run_dirichlet_multinomial_experiment`.
- [ ] **One validation style.** Currently mixed: bare `stopifnot`, named `stopifnot`, the `validate_*` helpers, inline
  `if`/`stop`, and a local `is_scalar`. `evaluate_thresholds` is the only place that calls `stop()` without
  `call. = FALSE`.
- [ ] **Align naming:**
  - Cell types are `index`/`index_k` in the multinomial outputs but `cell_type`/`cell_type_k` in the
    Dirichlet-multinomial outputs.
  - The header of `visualisation.R` says "Visualization Layer".
  - Folder names differ: `Simulations_Alemu/` vs `simulations/`, and `simResults_alemu` vs `simresults`.
- [ ] **Align plot defaults.** `plot_success_vs_tau` defaults to `target = NULL` and draws a dashed line;
  `plot_success_vs_n` defaults to `0.95`; `plot_success_rate_curve` draws a dotted line.

## P4: Clarity

- [ ] **Put functions in the right layer files:**
  - `compute_errors` is in `simulation.R`, but the header of `calculation.R` says it holds the error metrics.
  - The success rule (`replicate_pooled_error`, `replicate_success`, `success_rule_id`) is in `extraction.R`.
  - `validate_positive_integer` and `validate_positive_numeric` are in `simulation.R`, not `validation_utils.R`.
- [ ] **Remove the stale note** at [calculation.R:3-5](scripts/simulation_layers/calculation.R#L3-L5) (the solver
  replacement is done).
- [ ] **Move the misplaced doc block.** The block at
  [simulation.R:487-518](scripts/simulation_layers/simulation.R#L487-L518) describes the multinomial output of
  `run_replicates` but sits on `run_replicates_dirichlet_multinomial`. Split it, and document `run_replicates` itself.
- [ ] **Rename `p_value`** in `summarize_argmax` ([extraction.R:49](scripts/simulation_layers/extraction.R#L49)). It
  holds the true proportion, not a p-value; e.g. `p_true`.
- [ ] **Fix `generate_proportions()`**: it accepts `grid` but doesn't pass it on for `fixed_max_beta`, where the default
  length would be wrong anyway.
- [ ] **Remove the placeholder option** `"logistic_normal_multinomial"` from `simulate_counts()`'s `match.arg` choices
  (it only stops with "not implemented"); a comment is enough.
