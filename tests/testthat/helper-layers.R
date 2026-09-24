# Source every simulation layer so tests can call the functions directly.
layer_dir <- here::here("scripts", "simulation_layers")
for (f in list.files(layer_dir, pattern = "[.]R$", full.names = TRUE)) {
  source(f)
}
