# simulation_cache.R
#
# Cache helpers: avoid re-running expensive simulations when parameters
# have not changed.
#
# The cache key is an MD5 hash of the config list, so any change to any
# parameter automatically produces a new cache file and triggers a fresh run.
# ---------------------------------------------------------------------------


#' Compute an MD5 hash of an R object.
#'
#' Serialises `x` to a temporary file, computes its MD5 checksum, and returns
#' the hex string.  Only `tools` (part of base R) is required.
#'
#' @param x Any R object.
#' @return A 32-character hexadecimal string.
hash_config <- function(x) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  saveRDS(x, tmp)
  unname(tools::md5sum(tmp))
}


#' Load a cached simulation result, or run and cache it.
#'
#' Looks for an `.rds` file whose name encodes `name` and an MD5 hash of
#' `config` inside `cache_dir`.  If found, the file is loaded and returned
#' without re-running the simulation.  Otherwise `run_fn()` is called, the
#' result is saved to the cache, and then returned.
#'
#' @param run_fn    Zero-argument function that runs the simulation and returns
#'   a result.
#' @param config    Named list of all simulation parameters; used solely to
#'   derive the cache key.  Must be identical (bit-for-bit after serialisation)
#'   across runs to obtain a cache hit.
#' @param cache_dir Directory where `.rds` cache files are stored.  Created
#'   automatically if it does not yet exist.
#' @param name      Short label prepended to the cache file name
#'   (e.g. `"simulation_errorchoice"`).
#'
#' @return The simulation result, either loaded from cache or freshly computed.
cached_simulation <- function(run_fn, config, cache_dir, name) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  hash <- hash_config(config)
  cache_file <- file.path(cache_dir, paste0(name, "_", hash, ".rds"))

  if (file.exists(cache_file)) {
    message("Loading cached simulation result: ", basename(cache_file))
    return(readRDS(cache_file))
  }

  message("No cached result found for '", name, "'. Running simulation ...")
  result <- run_fn()
  saveRDS(result, cache_file)
  message("Simulation result cached: ", basename(cache_file))
  result
}
