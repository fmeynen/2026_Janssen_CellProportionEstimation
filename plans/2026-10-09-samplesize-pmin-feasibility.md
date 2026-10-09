# Plan: sample-size solver with p_min support and a feasibility check

Decisions from a grilling session on 2026-10-09. Goals:

- Support `fixed_min_beta` / `p_min` in the sample-size calculator.
- Tighten bound validation across the layers.
- Detect (c, n_people) combinations whose success-rate plateau never reaches the target, so the solver doesn't run on
  them.
- Add tests.

Line numbers are as of commit `c56463d` (`main`).

**Status (2026-10-09): planned**, not started.

## Background

- **p_min is dropped.** `simulate_success_at_n()` (`scripts/simulation_layers/simulation.R:820-826`) passes only
  `p_max` to `generate_proportions()`, so `proportion_method = "fixed_min_beta"` fails with "p_min must be provided".
  `p_min` is also missing from the sample-size cache key (`scripts/simulation_layers/orchestration.R:132-152`).
- **Vector bounds are not caught.** A vector bound makes `generate_proportions()` return a matrix, and nothing on the
  sample-size path stops that.
- **The success rate can level off below the target.** Under the Dirichlet-multinomial model the success rate levels
  off as cells per person grow (`reports/simulation_dm_errorchoice.qmd:246-273`). The pooled estimate converges to the
  average of the N people's θ_i ~ Dirichlet(c·p), not to p. For some (c, n_people) the plateau is below the target, so
  no n per person is enough, and today the solver then runs until `n_max` or `max_iterations`.
- **The existing tests use a fake simulation.** The sample-size tests inject a fake `simulate`
  (`tests/testthat/test-sample-size-experiment.R:28-36`). No sample-size test touches the real DM simulation,
  `p_min` or `p_max`.

## Decisions

- **Bounds are single values on the sample-size path.** Each run takes one `p_min` or one `p_max`. To sweep values,
  call `run_sample_size_experiment()` once per value. `p_max` gets the same treatment as `p_min`.
- **Validation has two levels.**
  - **Low level:** generalise `validate_p_max()` (`scripts/simulation_layers/validation_utils.R:118-136`) into
    `validate_p_bound(x, bound = c("p_min", "p_max"), method_arg = "proportion_method")`.
    - It checks that the bound is not NULL and every value is finite and strictly in (0, 1). It still accepts vectors.
    - It replaces `validate_p_max()` in `generate_props_fixed_max_beta()` (`simulation.R:90`) and
      `feasible_scenarios()` (`orchestration.R:197`).
    - It replaces the inline p_min checks in `generate_props_fixed_min_beta()` (`simulation.R:142-147`).
    - Remove `validate_p_max()` and port its tests (`tests/testthat/test-validation.R:18-30`).
  - **High level:** new `validate_proportion_bounds(proportion_method, K, p_min, p_max)`, for config-level callers.
    - The bound required by the method is present, via `validate_p_bound()`.
    - A bound the method does not use is an error.
    - Each bound is a single value.
    - The combination is feasible for K: `K * p_min <= 1` and `K * p_max >= 1`, with the generators' 1e-12 tolerance.
      Otherwise it stops with a plain `stop(call. = FALSE)`, worded like the generators' messages
      (`validation_utils.R:49-50, 77-78`), with no warning and no classed condition.
  - **Callers:** call `validate_proportion_bounds()` at the top of `run_sample_size_experiment()`,
    `simulate_success_at_n()` and `run_dm_errorchoice_experiment()` (`orchestration.R:415-438`).
- **An unused bound is an error everywhere.** The `generate_proportions()` dispatcher (`simulation.R:191-213`) also
  errors when a bound is set that the method does not use: `p_min` with `"beta"` or `"fixed_max_beta"`, and `p_max`
  with `"beta"` or `"fixed_min_beta"`. That covers `feasible_scenarios()` and the DM error-choice runner as well. A
  sweep of the scripts, reports and tests found no current call site that passes an unused bound.
- **p_min plumbing.**
  - `simulate_success_at_n()` passes `p_min = config$p_min`.
  - `run_sample_size_experiment()` adds `p_min = config$p_min` to its cache key.
  - There are no cached sample-size results in `results/simresults/`, so the key change invalidates nothing.
- **Feasibility check (per alpha, before the solver).**
  - **Where:** in `run_sample_size_experiment()`, inside the `cached_result()` compute, so it is cached together with
    the solver result.
  - **How:** call `simulate(alpha, config$n_max, config, config$seed)` through the same injectable `simulate` hook as
    the solver pilots. It runs for both models; the multinomial ceiling is about 1, so that model effectively always
    passes.
    - Using the same `config$seed` gives common random numbers with the solver's pilots.
  - **Rule:** with `x` successes out of `B`, the alpha is infeasible only if the one-sided 95% Clopper–Pearson upper
    bound is below `success_rate_target`. That bound is `qbeta(0.95, x + 1, B - x)`, or 1 when `x == B`. The
    confidence level is fixed at 0.95, with no config field.
  - **Infeasible:**
    - Skip the solver and return `sample_size = NA`, `stopping_reason = "infeasible"` and `iterations_used = 0`.
    - Return empty or `NULL` diagnostics; the implementer picks, and the `rbind` must still work.
    - `warning()` naming alpha, `concentration`, `n_people`, the success rate at `n_max` and the upper bound.
    - Continue with the next alpha.
  - **Borderline** (point estimate < target but upper bound >= target): run the solver as usual. It may still end at
    `n_max` with its existing warning.
  - **New output column:** `success_ceiling`, the success rate at `n_max`, filled for every alpha whether or not it is
    feasible.
  - **Warm start:** the next alpha starts from the last *feasible* alpha's `final_n`. If no alpha has been solved yet,
    it uses the resolved initial n.
- **n_init default.**
  - `config$n_init` becomes optional. When it is `NULL`, `run_sample_size_experiment()` uses `config$concentration`.
  - When it is `NULL` and the model is `"multinomial"`, there is no concentration to fall back on, so it errors and
    asks for `config$n_init`.
  - An explicit `n_init` wins.
  - The cache key records the resolved value (it already keys on `n_init`).
- **Driver** (`scripts/simulations/execute_simulation_samplesize.R`).
  - `simulation_sample_size_defaults()` gets `n_init = NULL`, `p_min = NULL` and `p_max = NULL`.
  - A comment shows `fixed_min_beta` usage (`proportion_method = "fixed_min_beta", p_min = 0.01`).
  - The default run stays plain `"beta"`.
  - The function's `n_init` argument goes away or defaults to `NULL`.
  - Plots must tolerate NA `sample_size` rows.
- **Delivery.**
  - Branch `feat/samplesize-pmin` off `main`.
  - Focused commits.
  - Full test suite green.
  - Open a PR.

## Steps

1. **Bound validation refactor** (`validation_utils.R`, `simulation.R`, `orchestration.R`).
   - Add `validate_p_bound()` and switch every `validate_p_max()` caller and the inline p_min checks to it.
   - Remove `validate_p_max()`.
   - Add `validate_proportion_bounds()`.
   - Port `validate_p_max()` tests to `validate_p_bound()`, for both bounds and the `method_arg` wording.
   - Keep the existing error-message substrings asserted in `tests/testthat/test-proportions.R:27,34` ("p_max must be
     provided", "p_min must be provided").
2. **Unused-bound error in `generate_proportions()`.**
   - Add the check to the dispatcher.
   - Tests: each method × each unused bound errors, and the existing valid calls still pass.
3. **Validate at config entry points.**
   - Call `validate_proportion_bounds()` at the top of `run_sample_size_experiment()`, `simulate_success_at_n()` and
     `run_dm_errorchoice_experiment()`.
   - Tests: missing, invalid, vector, unused and infeasible-for-K bounds error at each entry point before any
     simulation happens. Use a fake `simulate` that fails if called.
4. **p_min plumbing.** Pass `p_min` in `simulate_success_at_n()` and add it to the sample-size cache key.
5. **n_init default** in `run_sample_size_experiment()`.
   - Tests:
     - NULL resolves to `concentration` for DM, checked via the first solver centre or the cache key.
     - NULL errors for multinomial.
     - An explicit value wins.
6. **Feasibility check** in `run_sample_size_experiment()`.
   - Add the `success_ceiling` column.
   - Return an `"infeasible"` row with a warning when the alpha is infeasible.
   - Warm start from the last feasible alpha.
   - Update the roxygen docs, including the `simulate` hook now also being called at `n_max`.
   - Update existing call-count assertions in `test-sample-size-experiment.R` (+1 call per alpha) and the expected
     `res$sample_size` column names (`test-sample-size-experiment.R:50`).
7. **Driver.** Update `simulation_sample_size_defaults()` and make the plots NA-safe.
8. **New tests** (in `tests/testthat/test-sample-size-experiment.R`, or a new `test-sample-size-feasibility.R`).
   - **Plumbing**, with the real `simulate_success_at_n()` at tiny B and n:
     - with `fixed_min_beta` / `fixed_max_beta`, the generated p (from `rep_out` or `generate_proportions()`) attains
       the bound;
     - `p_min` and `p_max` each change the sample-size cache key.
   - **Model coverage:** `simulate_success_at_n()` runs for `model = "dirichlet_multinomial"` and `"multinomial"`,
     returning a logical `success` of length B and a consistent `success_count` / `success_rate`.
   - **Feasibility with a fake `simulate`:**
     - **Clearly infeasible:** the fake levels off at, e.g., 0.5 × B. Expect an NA row, `"infeasible"` and the
       warning.
     - **Borderline:** x at `n_max` gives a point estimate < target but a Clopper–Pearson upper bound >= target. The
       solver still runs.
     - **Warm start:** with alphas `c(feasible, infeasible, feasible)`, the third alpha starts from the first alpha's
       `final_n`.
     - `success_ceiling` is filled on every row.
   - **Feasibility with the real DM simulation:** small `concentration`, `n_people = 1` and a tight `tau_AE`, so the
     plateau is clearly below 0.95. Expect an NA row and the warning. The solver is skipped, so this is fast.
   - **End to end** with the real DM simulation, small B and K and a loose `rel_tol`:
     - `stopping_reason == "converged"`;
     - re-simulating at the solved n with a different seed gives a success rate within ±0.1 of the target.
     - Keep it to a few seconds; it runs in the normal suite with no skip.
9. **Wrap-up.**
   - Run the full suite.
   - Open a PR (`🤖 Generated with Claude Code` footer).
   - Update the status line of this plan.

## Out of scope

- Grids over p_min / p_max inside `run_sample_size_experiment()`.
- Solving for the number of people instead of cells per person.
- An analytic plateau formula.
- Packaging `simulation_layers/` and Windows parallelism (still in `backlog.md`).
