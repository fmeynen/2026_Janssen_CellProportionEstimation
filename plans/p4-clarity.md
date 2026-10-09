# Plan: P4 Clarity

Decisions from a grilling session on 2026-10-01, covering the `## P4: Clarity` items in [backlog.md](../backlog.md).
Line numbers are as of commit `18df9b8`.

**Status (2026-10-01): done.** All six sections are implemented on `refactor/p4-clarity`, one commit each.
`results/simresults/` had no cache files, so section 4 deleted nothing.

## Delivery

- Branch: `refactor/p4-clarity` (from `main`).
- One commit per section below. The two extra fixes are folded into the related sections.
- Tick each P4 item in `backlog.md` with a short `Resolved:` note, as was done for P3.
- Run the testthat suite (`tests/testthat/`) before opening the PR.
- Open a PR against `main`.

## 1. Put functions in the right layer files

- Move `compute_errors` from `scripts/simulation_layers/simulation.R` (L340-377, doc block included) to
  `calculation.R`. It has two callers in `simulation.R` (L616, L704). They need no change, because every layer is
  sourced before anything runs.
- Move the success rule from `extraction.R` to `calculation.R`:
  - `success_rule_id`
  - `replicate_pooled_error`
  - `replicate_success`
  - their doc blocks and any helpers or comments that belong only to them

  Put them next to `evaluate_thresholds` / `success_rate_from_stat`.
- Leave `extract_success_rate` in `extraction.R`. It is the thin wrapper that calls `replicate_success()`.
- Leave `person_results_required_cols` in `validation_utils.R`.
- Extra fix: the `success_rule_id` doc mentions `run_simulation_dm_errorchoice()`, which doesn't exist. The correct
  name is `run_dm_errorchoice_experiment()`.
- Update comments that point to the old location. For example, `tests/testthat/test-simulate-success.R:2` says the
  success rule is in "(extraction.R)", and `tests/testthat/test-success.R:1` probably names the file too. Grep for
  `extraction.R` in tests and scripts.
- Note: P2's "slim the DM replicate output" will later rewrite `replicate_pooled_error`/`replicate_success`. This move
  only decides where they belong.

## 2. Remove the stale note and rewrite the headers

- Delete the NOTE at `calculation.R:3-5`. It says the GLM helpers "are being replaced by a new solver below", but that
  solver (`sample_size_pilots` … `estimate_sample_size`) now exists.
- Rewrite the `calculation.R` header to list what the file holds after section 1:
  - error metrics (`compute_errors`)
  - the success rule
  - threshold evaluation
  - the sample-size solver
- Update the `extraction.R` header if it needs it after the success rule moves out.

## 3. Split the misplaced `run_replicates` doc block

The block at `simulation.R:543-578` describes `run_replicates` (the multinomial return value) but sits on
`run_replicates_dirichlet_multinomial`. `run_replicates` (L656) has no doc block. It is the public dispatcher: it
delegates to `run_replicates_dirichlet_multinomial()` for `model = "dirichlet_multinomial"` and runs the multinomial
path itself.

- Give `run_replicates` the full doc block:
  - every `@param`
  - both return shapes: the multinomial list (`max_errors`, `argmax`, `errors`, `phat`, `inputs`) and, for the DM
    model, a pointer to `run_replicates_dirichlet_multinomial()`'s return value
- Give `run_replicates_dirichlet_multinomial` a short block of its own:
  - its params (`p`, `B`, `metrics`, `n_people`, `n_per_person`, `concentration`, `scenario_id`, `seed`)
  - its return value: `person_results`, with one row per replicate × person × cell type × metric, plus `inputs`

## 4. Rename `p_value` to `true_proportion`

- `summarize_argmax` (`extraction.R:49`): rename the column and its `@return` doc (L34).
- `run_simulation_experiment` (`orchestration.R:675`): update the column selection and the `argmax_summary` doc
  (around L556).
- No `CACHE_SCHEMA` bump. It would invalidate the expensive dm_errorchoice cache too.
- Instead, delete only the errorchoice cache `.rds` files under `results/simresults/`, which come from
  `scripts/simulations/simulation_errorchoice.R:72`. They regenerate with the new column name on the next run.
  - Check the file names before deleting, and leave the dm_errorchoice files alone.
- Grep for any test that checks the column name.

## 5. Fix `grid` in `generate_proportions()`

`generate_proportions()` (`simulation.R:171`) defaults `grid` to `default_beta_grid(K)`, but passes it only to the
`"beta"` method. The fixed-max and fixed-min generators expect length `K - 1`.

- Change the default to `grid = NULL`.
- If `grid` is `NULL`, call each generator without it, so each one uses its own default.
- Otherwise, pass `grid` to all three methods; each generator's own length check applies.
- Update the `@param grid` doc to say this.
- Add tests:
  - a custom grid reaches `fixed_max_beta` (and `fixed_min_beta`)
  - a wrong length is rejected with the generator's error

## 6. Remove the `logistic_normal_multinomial` placeholder

- `simulate_counts()` (`simulation.R:287-315`): remove `"logistic_normal_multinomial"` from the `match.arg` choices
  and remove its `switch` branch. Check the `@param model` doc as well.
- Keep the commented-out stub at `simulation.R:267` as the note for a future model.
- Delete the test at `tests/testthat/test-counts.R:25-31`.
- Extra fix: reword the comment at `scripts/simulations/simulation_errorchoice.R:15` ("Allow correlations -> model =
  'logistic_normal_multinomial'") as a future idea rather than an existing option.
