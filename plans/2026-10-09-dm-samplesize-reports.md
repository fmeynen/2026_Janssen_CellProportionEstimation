# Plan: required cells per person in the DM error-choice reports

Decisions from a grilling session on 2026-10-09. Goal: use the sample-size solver (see
`plans/2026-10-09-samplesize-pmin-feasibility.md`) in `reports/simulation_dm_errorchoice.qmd` and
`reports/simulation_dm_errorchoice_pmin.qmd`. Each report gets a graph of the cells per person needed to reach a 0.95
success rate, against the number of people:

- x-axis: `n_people`;
- y-axis: cells per person;
- one figure per metric (AE, ARE);
- one panel per alpha.

Line numbers are as of commit `b40a331` (`feat/samplesize-pmin`).

**Status (2026-10-09): planned**, not started.

## Background

- **Report settings.** Both reports take their settings from `simulation_dm_errorchoice_defaults()`
  (`scripts/simulations/simulation_dm_errorchoice.R:55-77`):
  - alpha {2, 3, 4, 5}, K = 10, B = 1000, c = 1e4, `n_people` {1, 2, 3, 5, 10};
  - `taus_fixed` AE = 0.005, ARE = 0.5;
  - target 0.95, seed 260926.
  - The p_min report uses `simulation_dm_errorchoice_pmin_defaults()` (`:82-87`), which adds `fixed_min_beta` and
    `p_min = 0.01`.
- **The solver checks every metric at once.** `run_sample_size_experiment()` (`scripts/simulation_layers/orchestration.R`)
  counts a replicate as a success only if it passes every metric in `config$taus`, and takes a single `n_people`.
  - So keeping AE and ARE separate, as the reports do, needs one solver run per (`n_people`, metric).
  - That is 10 runs per report, each solving 4 alphas.
- **Solver output.** `res$sample_size` has the columns `alpha`, `sample_size`, `stopping_reason`, `iterations_used` and
  `success_ceiling`.
  - `stopping_reason` is `"tolerance"` when the solver converged and `"infeasible"` when the target cannot be reached
    (with `sample_size` NA).
  - Other reasons, such as hitting `n_max` or `max_iterations`, mean the solver did not converge.
  - Check `estimate_sample_size()` in `scripts/simulation_layers/calculation.R` for the exact strings.

## Decisions

- **Placement.**
  - Add a new section, "Required cells per person", to each of the two reports. Each report shows its own AE and ARE
    figures.
  - `simulation_dm_errorchoice.qmd`: put it after "Success rate versus cells per person" (ends before line 250, "Why
    the curves level off").
  - `simulation_dm_errorchoice_pmin.qmd`: put it after "Success rate versus cells per person" (from line 183), before
    the `result-values` chunk if that chunk only feeds earlier text; otherwise put it at the end of the results.
- **Settings.**
  - Same as the reports: alpha, K, B = 1000, c, the `n_people` grid {1, 2, 3, 5, 10}, the `taus_fixed` thresholds,
    `target`, `seed`, and for the p_min report `fixed_min_beta` with `p_min`.
  - Add these solver fields to `simulation_dm_errorchoice_defaults()`, so the p_min defaults inherit them: `rel_tol =
    0.01`, `max_iterations = 20L`, `f0 = 2`, `f_floor = 1.1`, `n_max = 1e9`, `n_init = NULL`.
  - `n_init = NULL` makes the solver start from c.
  - The solver needs a `tie_method` field; use `"random"`, as in the sample-size driver.
- **Seed.** Use the same `config$seed` for every (`n_people`, metric) run. This gives common random numbers across N
  and across metrics.
- **Grid runner.**
  - New `run_dm_errorchoice_samplesize(config, cache = TRUE, force_recompute = FALSE, cache_dir = ...,
    simulate = simulate_success_at_n)` in `scripts/simulations/simulation_dm_errorchoice.R`.
  - It loops over `n_people` × `metrics`. For each pair it builds a solver config:
    - `taus = config$taus_fixed[metric]`, `metrics = metric`;
    - `model = "dirichlet_multinomial"`, `n_people`, `concentration = config$concentration`;
    - `success_rate_target = config$target`;
    - the solver fields, the alpha grid, `proportion_method`, `p_min` and `p_max`.
  - It calls `run_sample_size_experiment()` with that config, `cache`, `force_recompute`, `cache_dir` and `simulate`.
  - It returns one tidy data.frame with `alpha`, `n_people`, `metric`, `sample_size`, `stopping_reason`,
    `iterations_used` and `success_ceiling`.
  - It keeps infeasible rows.
  - `run_sample_size_experiment()` already caches per alpha.
  - Unreachable alphas emit warnings. The reports run with `warning: false`; the runner should not suppress them
    itself.
- **Plot.**
  - New `plot_sample_size_by_people(df, metric, subtitle = NULL)` in `scripts/simulation_layers/visualisation.R`,
    following the style of `plot_success_vs_n()` (`:402`).
  - Axes: x is `n_people`. y is cells per person on a **log10 scale, shared across panels**. Facet by alpha.
  - Converged points: a line plus filled points, drawn only through converged points.
  - **Not converged** (feasible, but `stopping_reason` is neither `"tolerance"` nor `"infeasible"`): an open marker at
    the solver's final n, labelled "not converged" in the legend.
  - **Unreachable** (`"infeasible"`): a × marker pinned to the top of the panel (`y = Inf`, or the panel's maximum),
    labelled "target unreachable" in the legend.
  - It must work when a panel has no converged points.
- **Report text.**
  - Figures with captions. One short sentence says that n is the smallest number of cells per person that reaches the
    target, found with the sample-size solver.
  - Captions mention that combinations whose success rate cannot reach the target even at very many cells per person
    (a check at `n_max`) are marked with × at the top. They do not describe the method beyond that.
  - No inline-number paragraph and no table.
- **Tests.**
  - `run_dm_errorchoice_samplesize()` with a fake `simulate` and a temporary `cache_dir`:
    - one row per (alpha, `n_people`, metric);
    - each solver run sees only its own metric in `taus` and `metrics` (the fake records the configs it receives);
    - infeasible rows are kept with NA `sample_size`;
    - `n_people` is forwarded.
  - `plot_sample_size_by_people()` on a small hand-made data.frame with converged, not-converged and infeasible rows
    across two alphas:
    - it returns a ggplot that builds (`ggplot2::ggplot_build()`) without error;
    - the y scale is log10;
    - it is faceted by alpha.
  - Tests go in `tests/testthat/test-dm-errorchoice.R` and `tests/testthat/test-plots.R`.
- **Delivery.**
  - Continue on `feat/samplesize-pmin` with focused commits: config fields + runner, plot, tests, reports, PDFs.
  - Before running the real simulations, **estimate the runtime and ask the user** (see "Runtime estimate").
  - Then render both PDFs and commit them, as with the earlier reports.

## Runtime estimate (to confirm with the user before running)

- **Reference:** the full DM error-choice simulation took about 52 s, which works out to about 0.1 s per scenario
  point. A scenario point is one (alpha, N, n) combination with B = 1000.
- **Per alpha:** one solver run costs one simulation at `n_max` plus 3 pilots per iteration.
  - Typically about 5 iterations, so about 16 points, or about 1.6 s.
  - Worst case 20 iterations, so about 61 points, or about 6 s.
- **Per report:** 4 alphas × 5 N × 2 metrics = 40 alpha-solves.
  - Typically about 1–2 min; worst case about 4 min.
  - The cost grows with N, because each replicate draws N people.
- **Both reports:** typically about 2–4 min, worst case about 8 min, plus about 25 s per PDF render.
- Re-renders after that are fast, because the results are cached.
- Before quoting a number to the user, time one (alpha, N = 10, metric) solve.

## Steps

1. **Config fields.** Add the solver fields and `tie_method` to `simulation_dm_errorchoice_defaults()` and document
   them in its roxygen block.
   - Check that the DM error-choice cache key is unaffected: `run_simulation_dm_errorchoice()` only keys on
     simulation-relevant fields.
2. **Grid runner.** Add `run_dm_errorchoice_samplesize()` with a roxygen block, and update the file header's workflow
   list.
3. **Plot.** Add `plot_sample_size_by_people()` with a roxygen block.
4. **Tests.** Add the runner and plot tests described above. The full suite must pass with 0 warnings; the baseline is
   893 passing.
5. **Reports.**
   - In each report's setup chunk, call `samplesize <- run_dm_errorchoice_samplesize(config)`.
   - Add the new section with two figure chunks (`fig-samplesize-ae`, `fig-samplesize-are`), each with `fig-cap` and
     `fig-alt`.
   - The subtitle shows c and τ.
6. **Run + render.**
   - Time one solve, give the user the runtime estimate and get their approval.
   - Run both reports' sample-size grids, which populates the cache.
   - Render both reports to PDF and check that the figures look right (markers, log axis, facets).
7. **Wrap-up.**
   - Commit the PDFs.
   - Update the status line of this plan.
   - Push the branch and open or update the PR.

## Out of scope

- Solving for the number of people.
- Combining AE and ARE in one success criterion.
- Showing both proportion settings in one figure.
- Tables of solved n.
- Changes to the existing report sections.
