# simulation_io.R
#
# I/O helpers for persisting simulation results.
#
# Each unique combination of simulation parameters is saved as a separate
# .rds file whose name contains a short MD5 hash of those parameters.
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


#' Build the file path for a simulation result.
#'
#' The file name is `<name>_<md5>.rds`, where the MD5 is derived from the
#' serialised `config`.  Identical parameters always resolve to the same path;
#' any change produces a different path.
#'
#' @param config  Named list of simulation parameters.
#' @param dir     Directory that will hold the result files.  Created
#'   automatically if it does not yet exist.
#' @param name    Short label prefixed to the file name
#'   (e.g. `"simulation_errorchoice"`).
#'
#' @return Absolute path string.
simulation_result_path <- function(config, dir, name) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(name, "_", hash_config(config), ".rds"))
}

#' Build deterministic path for isoband simulation output.
#'
#' @param config Named list of isoband parameters.
#' @param dir Directory for output files.
#' @param name Prefix for file names.
#'
#' @return Absolute path string.
isoband_result_path <- function(config, dir, name = "simulation_isoband") {
  simulation_result_path(config = config, dir = dir, name = name)
}

#' Save isoband result object to disk.
#'
#' @param result Isoband result object.
#' @param config Named list of isoband parameters.
#' @param dir Directory for output files.
#' @param name Prefix for file names.
#'
#' @return Invisibility returns the saved file path.
save_isoband_result <- function(result, config, dir, name = "simulation_isoband") {
  out_path <- isoband_result_path(config = config, dir = dir, name = name)
  saveRDS(result, out_path)
  invisible(out_path)
}
