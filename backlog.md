# Backlog

Open items left over from the 2026-09-30 code review (clarity, efficiency, consistency). The P1–P4 passes are
done; resolved items were removed (see git history for their resolutions). Items are listed in priority order.

## Efficiency

- [ ] **Parallelism on Windows.** [`replicate_cores()`](scripts/simulation_layers/simulation.R#L337) returns 1 on
  Windows, so all runs are serial. The streams don't depend on the worker, so a PSOCK/`parLapply` backend would keep
  results reproducible. The layer functions need to be sourced on each worker.
  Deferred: to be decided in a separate PR now that the P2 timings are known.

## Structure

- [ ] **Consider packaging `simulation_layers/`.** Turning it into a package (`DESCRIPTION` plus `R/`, loaded with
  `devtools::load_all()`) would let `R CMD check` flag undefined functions. Not done; `scripts/load_layers.R` covers loading.
