# Plan: P2 Efficiency

Decisions from a grilling session on 2026-10-01, covering the `## P2: Efficiency` items in [backlog.md](../backlog.md).
Line numbers are as of commit `c112f1c` (the tip of `refactor/p4-clarity`).

**Status (2026-10-01): done.** Steps 1–5 are committed on `perf/p2-efficiency`, one commit each.
- Benchmark: 145.7 s → 26.5 s (5.5×). Old and new success rates are statistically equivalent: 270 comparisons, max
  |z| = 2.73.
- The old code came from `git archive` instead of a worktree, because the worktree checkout hit Windows path-length
  limits.

## Decisions

- **Scope.** Four of the five P2 items: slim the Dirichlet-multinomial (DM) replicate output, stop rebuilding data in
  the multinomial success path, vectorise the Dirichlet draws, and replace the per-scenario filter in
  `extract_success_rate`. **Parallelism on Windows is out of scope.** It stays open in the backlog, to be decided in a
  later PR once the timings from this pass are known.
- **Reproducibility.** Results need only be **statistically equivalent**: the same seed may give different numbers
  than before. As a consequence:
  - Bump `CACHE_SCHEMA` from 2 to 3 (`orchestration.R:47`), so old cached results are never reused.
  - Keep `success_rule_id()` at `"pooled_proportion_v1"`; the rule's definition doesn't change.
  - Don't re-render `reports/simulation_dm_errorchoice.pdf` in this pass.
- **Verification.**
  - Exact unit tests for deterministic code: the new core must give exactly the same values as the old data.frame
    implementation on fixed inputs.
  - A one-off benchmark outside the test suite: time old against new code and compare their success rates (see
    step 6).
- **Delivery.**
  - Branch `perf/p2-efficiency` off `refactor/p4-clarity`.
  - One commit per step below.
  - Tick the backlog items, including the timings.
  - Write the PR description to `plans/p2-efficiency-pr.md`. Pushing and opening the PR are done by hand.

## Step 1: Matrix-based success-rule core

In `scripts/simulation_layers/calculation.R` (success-rule section):

- Add `pooled_error_stat(pbar, p, metric)`:
  - `pbar`: B × K numeric matrix of pooled estimates.
  - `p`: true proportions, length K.
  - `metric`: `"AE"` or `"ARE"`. Anything else errors with the same "AE and ARE only" message as
    `replicate_pooled_error`.
  - Returns a length-B numeric vector, the max over cell types of `|pbar - p|` (AE) or `|pbar - p| / p` (ARE).
  - NaN counts as 0 and Inf is kept, exactly as the current rule does.
  - Use matrix operations, e.g. `sweep` / `abs` / `matrixStats`-free `apply(…, 1, max)` or `do.call(pmax, …)` over
    columns. Avoid per-row loops; `do.call(pmax, as.data.frame(x))` is fast for small K.
- Rewrite `replicate_pooled_error(person_results, metric)` as a thin wrapper:
  - Keep the validation, the `scenario_id` handling and the output shape.
  - Build `pbar` per (scenario, replicate) with `rowsum`, then call `pooled_error_stat`.
- Leave `replicate_success` as is; it already calls `replicate_pooled_error`.
- Tests (`tests/testthat/test-success.R`):
  - Add `pooled_error_stat` cases: AE, ARE, the NaN → 0 case, the Inf case, a non-AE/ARE metric erroring, and B = 1.
  - The existing `replicate_pooled_error` / `replicate_success` tests must pass unchanged. That is the
    exact-equivalence check against the old implementation.
- Document `pooled_error_stat` as the single source of truth for the rule. Update the `replicate_pooled_error` doc,
  whose text says it is the single source of truth.

## Step 2: Slim the DM replicate output

Change `run_replicates_dirichlet_multinomial` (`simulation.R:526`) and its doc block:

- **Metric check.** Validate `metrics` up front: only `"AE"` and `"ARE"` are allowed, otherwise
  `stop(..., call. = FALSE)`. No current DM caller uses TSE or LAE. Mirror the check in `run_replicates`'s doc.
- **New argument.** Add `keep_person_results = FALSE`.
- **Per replicate:**
  - Compute `observed_p` (n_people × K).
  - Return the pooled estimate `colMeans(observed_p)`.
  - Only when `keep_person_results` is `TRUE`, also return the per-person data.frame, with the same columns as today
    **minus `error`**.
  - Drop the per-person `compute_errors` call entirely.
- **Assembly.** After `replicate_apply`, build `phat` (B × K, `rbind` of the pooled vectors). Then build
  `max_errors` (B × M, one column per metric, from `pooled_error_stat(phat, p, m)`, column names = metrics).
- **Return** `list(phat, max_errors, inputs)`. When `keep_person_results` is `TRUE`, also return `person_results`
  (`do.call(rbind, …)` of the per-replicate frames). Put `keep_person_results` in `inputs`.
- `run_replicates` (`simulation.R:635`): add and forward `keep_person_results`, and document it as DM only.

Update the callers:

- **`run_dm_errorchoice_experiment`** (`orchestration.R`, ~L440-470): read `rep_out$max_errors[, m]` instead of
  calling `replicate_pooled_error(rep_out$person_results, m)`. Remove the `rm(rep_out)` memory comment if it no
  longer applies.
- **`simulate_success_at_n`**, DM branch (`simulation.R:786-804`): compute pass/fail from `rep_out$max_errors`.
  - A replicate passes when, for every metric in `config$taus` that is present, `max_errors[, m] <= taus[[m]]`.
  - Match `replicate_success`'s semantics for metrics missing from `taus` (they're skipped; the warning is already
    raised above).
  - Factor this into a small helper (e.g. `pass_from_max_errors(max_errors, taus)` in `calculation.R`), so step 3
    can reuse it.
  - Check that it matches `replicate_success`'s handling of missing metrics and of NaN/Inf.
- **`run_dirichlet_multinomial_experiment`** (`orchestration.R:246`, ~L300-325): call with
  `keep_person_results = TRUE`, and remove `"error"` from the column selection. Update its doc (and
  `run_simulation_experiment`'s, if it lists the columns).
- **`CACHE_SCHEMA`:** bump to `3L` (`orchestration.R:47`). Check the doc and tests that mention the value
  (`tests/testthat/test-cache.R`).

Tests:

- **`test-replicates.R`:**
  - The DM shape tests (~L60-111) now cover:
    - the default output (`phat` B × K, `max_errors` B × M with metric column names, no `person_results`)
    - `keep_person_results = TRUE` (`person_results` columns without `error`, at ~L75)
  - Keep the reproducibility tests for the same seed (`expect_identical` on `phat`/`max_errors`).
  - Keep the prefix-stability test (8 vs 4 replicates, ~L107) on `phat`.
  - Add a test that TSE/LAE are rejected for the DM model.
  - Add a test that `max_errors` equals `pooled_error_stat(phat, p, m)`.
- **`test-simulate-success.R`:**
  - L44-47: compare against a pass computed from `rep_out$max_errors` / `pooled_error_stat`, instead of
    `replicate_success` on `person_results`.
  - L68: compare `phat` instead.
  - L71: compare `max_errors` instead of `person_results$error`.
  - L149: `population_mean_proportion` is gone from the default output; use `max(rep_out$inputs$p)` (or the
    returned `p`) instead.
- **Experiment tests** (`test-experiment.R`, `test-success.R` `extract_success_rate`) must still pass via the
  opt-in path.

## Step 3: Multinomial success path without the long data.frame

In `simulate_success_at_n`'s multinomial branch (`simulation.R`, ~L806-850):

- Delete the `expand.grid` / `person_results` construction.
- Compute `max_errors` for the success rule as `pooled_error_stat(rep_out$phat, p, m)` for each metric in
  `config$metrics`. Use the core on `phat`, **not** `rep_out$max_errors`: that one comes from `compute_errors` +
  `max`, which returns NaN instead of 0 for ARE 0/0.
- Apply the step-2 helper (`pass_from_max_errors`).
- Keep the comment explaining why the multinomial model reduces to a single synthetic person; shorten it to fit.
- Tests: the existing `test-simulate-success.R` multinomial tests must pass. Add one test that the result equals
  `replicate_success` on the old long-format construction, built inside the test from `rep_out$phat`. That keeps
  an exact-equivalence check.

## Step 4: Vectorise the Dirichlet draws

In `simulate_counts_dirichlet_multinomial` (`simulation.R:245`):

- Replace the per-person `vapply(sample_dirichlet(...))` with one draw:
  - `g <- matrix(stats::rgamma(n_people * K, shape = rep(concentration * p, each = n_people), rate = 1), n_people, K)`
    (or another layout that is clearly correct)
  - then `person_true_proportions <- g / rowSums(g)`
- Keep a single validation up front.
- The draw order may change (statistically equivalent is enough).
- Keep the invalid-total check: any `rowSums(g)` that is not finite or not `> 0` errors with the same message.
- Leave `sample_dirichlet` in place if anything else or any test uses it; otherwise remove it with its tests.
- Keep the per-person multinomial draws (`rmultinom` can't vectorise over different probabilities).
- Tests (`test-counts.R` / `test-replicates.R`):
  - Shapes, rows summing to 1, and same-seed reproducibility.
  - A loose moment check: the mean of `person_true_proportions` over many people is about `p`, within a generous
    tolerance and with a fixed seed.

## Step 5: `extract_success_rate` in one pass

In `extract_success_rate` (`extraction.R:172`), replace the per-scenario `lapply` filter:

- Use one `rowsum` (or `tapply`) over `replicate_pass` grouped by `scenario_id` to get `B`, `success_count`,
  `success_rate` and each `success_rate_<metric>`.
- Keep the scenario order and the merge with `p_table` as they are.
- The existing test "extract_success_rate() output format is unchanged" (`test-success.R:210`) must pass unchanged.

## Step 6: Benchmark and statistical check (one-off, not committed)

Before running anything:

- Estimate the run time.
- Ask the user whether that is acceptable, as required by their global instructions.

Setup:

- **Old code:** a `git worktree` at commit `c112f1c`, created in the session scratchpad.
- **New code:** the branch tip.
- **Script:** a small R script in the session scratchpad. Give its full absolute path when reporting.

The script:

- Sources each tree's `scripts/load_layers.R`.
- Runs one reduced `run_dm_errorchoice_experiment` config, e.g. 2 alphas × 2 `n_people` × 3 `n_per_person`,
  `B = 1000`, AE and ARE, fixed seed.
- Records the time with `system.time` for each tree. Run each once, or best of 2 if it's quick.
- Compares success rates at a few taus per scenario: old vs new should agree within Monte Carlo error. Report the
  max absolute difference next to the binomial standard error `sqrt(r(1-r)/B)`, and flag any difference above about
  3 standard errors.

Report the timings and the comparison table in the PR description and in the backlog's `Resolved:` note.

## Step 7: Wrap-up

- Run the full testthat suite: `Rscript -e "testthat::test_dir('tests/testthat', reporter='summary')"`.
- In `backlog.md`:
  - Tick the four P2 items, each with a `Resolved:` note; the "slim" note includes the before/after timings.
  - Leave "Parallelism on Windows" open and note that it's deferred until these timings are known.
- Write `plans/p2-efficiency-pr.md`: a title line plus a description that covers the behaviour changes. Don't commit
  it. The behaviour changes are:
  - the default DM output shape
  - DM metrics restricted to AE/ARE
  - `person_results` is opt-in and has no `error` column
  - the `CACHE_SCHEMA` bump
  - the same seed no longer reproduces old numbers
- Add a status line to this plan.
