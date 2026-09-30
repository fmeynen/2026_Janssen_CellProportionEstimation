# Load simulation layers -------------------------------------------------------------------------------------------
# Shared loader: sources every R file in scripts/simulation_layers/ (in sorted order) into the calling environment.
#
# Used by the simulation drivers in scripts/simulations/ and by tests/testthat/helper-layers.R.

invisible(lapply(
  sort(list.files(here::here("scripts", "simulation_layers"), pattern = "[.]R$", full.names = TRUE)),
  source
))
