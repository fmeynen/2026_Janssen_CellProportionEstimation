# Plan: DM error-choice report with p_min = 0.01

Decisions from a grilling session on 2026-10-09. Goal: a report like `reports/simulation_dm_errorchoice.qmd`, but
with the smallest population proportion fixed at `p_min = 0.01`, and without the concentration and level-off
explainers. Line numbers are as of commit `74856ea` (`main`).

**Status (2026-10-09): done** on `feat/dm-errorchoice-pmin`.
- Steps 1–5 and 7 are done.
- Step 6 was a no-op: `results/simresults/` held no cache files.
- The full DM simulation takes about 52 s, and the PDF render about 22 s with the simulation cached.
- The test suite has 726 passing tests.

## Decisions

- **Generator semantics (replace in place).** `fixed_min_beta` and `fixed_max_beta` keep their names, but their
  construction changes. The old construction pinned the bound at index 1 / K, built a Beta remainder over a K - 1 grid,
  and failed when the combination was impossible. That construction is removed.
  - Start from plain `generate_proportions_beta(alpha, K)` (K-point default grid), the same composition as the
    `"beta"` method.
  - **Clip + always pin the extreme.** The most extreme type (smallest for `p_min`, largest for `p_max`) is always set
    to the bound, together with every type that crosses it. The unpinned types are rescaled proportionally to fill
    `1 - n_pinned * bound`. Repeat until no unpinned type crosses the bound.
  - So the bound is **always attained exactly**. Ties at the bound are allowed. Pinned values are exactly the bound,
    not approximately ("exact" renormalisation, not `pmax(p, p_min) / sum(...)`).
  - If nothing crosses, only the single extreme type moves to the bound and the others scale proportionally. For
    example, `p_max = 0.3` at alpha = 2 sets type 10 to 0.30 and scales types 1..9 by 0.70 / 0.81.
  - Error only when infeasible: `K * p_min > 1` or `K * p_max < 1`.
  - A vector `p_min` / `p_max` still returns a matrix with one row per value.
  - Expected result at `p_min = 0.01`, K = 10:
    - alpha = 2: unchanged (beta min is already 0.01).
    - alpha = 3: 0.01, 0.01, 0.019, ...
    - alpha = 4: 0.01 x3, 0.017, ...
    - alpha = 5: 0.01 x4, 0.020, ...
- **Other simulation settings** are the same as `simulation_dm_errorchoice_defaults()`:
  - alpha 2-5, K = 10, B = 1000, c = 1e4, N in {1, 2, 3, 5, 10}
  - n grid 10^2-10^8, n_fixed = 2e5
  - tau_AE = 0.005, tau_ARE = 0.5, target 0.95, same seed
- **Caches.**
  - Add `p_min` and `p_max` to the DM cache key in `run_simulation_dm_errorchoice()`.
  - Delete the old `simulation_errorchoice` cache files built with `fixed_max_beta`, since their key no longer
    describes their proportions. List them and confirm with the user before deleting. No generator-version id is added
    to the keys.
- **Plot.** `plot_true_proportions()` supports both `fixed_max_beta` and `fixed_min_beta`: the final proportions as
  points over the scaled beta curve, plus a dashed horizontal line at the bound, with facet rows by `p_max` / `p_min`.
- **Untouched:** `scripts/deprecated/deprecated.R` and `reports/simulation_errorchoice_report.Rmd` (no wording change,
  no re-render).
- **Delivery.**
  - New branch off `main`.
  - Tests pass.
  - Render the new report to PDF. Estimate the render time and ask the user first.
  - Open a PR.

## Steps

1. **Generators** (`scripts/simulation_layers/simulation.R:29-154`).
   - Rewrite `generate_props_fixed_max_beta()` and `generate_props_fixed_min_beta()` with the pin-and-rescale loop.
     They could share one helper, e.g. `pin_to_bound(p, bound, side = c("min", "max"))`.
   - Their `grid` argument now defaults to the K-point beta grid.
   - Update the roxygen docs and the `generate_proportions()` docs (`grid` length is now K for every method).
   - In `validation_utils.R`, replace `fail_fixed_max_beta_impossible()` / `fail_fixed_min_beta_impossible()` with
     the feasibility check (`K * p_min > 1`, `K * p_max < 1`).
   - Check `feasible_scenarios()` (`orchestration.R:194`): the skip-impossible logic now only applies to the
     feasibility check.
2. **Tests** (`tests/testthat/test-proportions.R` and the others that reference `fixed_*_beta`). Rewrite them for
   the new semantics:
   - Bound attained exactly, with no value crossing it.
   - Sums to 1.
   - Ties allowed.
   - The no-crossing case pins only the extreme type.
   - Matches the plain beta when the bound already equals the extreme.
   - Infeasible bounds error.
   - Vector input returns a matrix.
3. **DM plumbing.**
   - Pass `p_min` and `p_max` through `run_dm_errorchoice_experiment()` (`orchestration.R:412`, the
     `generate_proportions()` call at line 433) and `run_simulation_dm_errorchoice()`
     (`scripts/simulations/simulation_dm_errorchoice.R`), including the cache key and the returned `inputs`.
   - Add `simulation_dm_errorchoice_pmin_defaults()` to the same script: the current defaults plus
     `proportion_method = "fixed_min_beta"` and `p_min = 0.01`.
4. **Plot.** Extend `plot_true_proportions()` (`visualisation.R:14-120`) as described above, and add a test.
5. **Report** `reports/simulation_dm_errorchoice_pmin.qmd`, a copy of `simulation_dm_errorchoice.qmd` with these
   changes:
   - Title: "Dirichlet-Multinomial Error-Choice Simulation (p_min = 0.01)".
   - Use `simulation_dm_errorchoice_pmin_defaults()`.
   - Add "Proportion method" and "p_min" rows to `tbl-setup`. Update the overview and setup text to describe how the
     proportions are generated with the fixed minimum.
   - Remove `## Concentration c` (`tbl-theta`, the `theta-example` chunk and its paragraph).
   - Remove `# Why the curves level off`, including `## Illustration with Additional Simulation` and the overlay
     chunks.
   - Remove helpers that are no longer used (`var_people`, `var_count`, `sd_total`, `technical_share`, `n_c`).
   - Trim the remaining c / level-off wording:
     - the `fig-n-ae` caption ("levels off from roughly n ≈ c") and alt text
     - the AE paragraph (the n ≈ c value, "Beyond roughly 10c ...", the `@sec-why-level` reference)
     - the ARE paragraph's `c p_j << 1` explanation, if it still holds; rewrite it against the new p table
   - Keep the one-line definition of c in the setup.
6. **Caches.** List the `results/simresults/*.rds` files from `fixed_max_beta` runs, confirm with the user, then
   delete them.
7. **Verify and deliver.**
   - Run the full test suite.
   - Estimate the render time (a new DM simulation run is needed) and ask before rendering.
   - Render the PDF, commit, and open a PR.
