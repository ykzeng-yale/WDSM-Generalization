#!/usr/bin/env Rscript
# Prepare an admitted current-source forward config. Historical execution pins
# remain unchanged; no observations, counts, fits, streams or jobs are created.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[1L]))
  stop("Usage: Rscript simulations/prepare_bounded_wdsm_replay.R NEW_CONFIG.json")
output <- args[1L]
if (!dir.exists(dirname(output)) || file.exists(output) || file.exists(paste0(output, ".partial")))
  stop("Output parent must exist and the config/partial paths must be new.")
had_rng <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
if (had_rng) rng_before <- get(".Random.seed", .GlobalEnv)
for (pkg in c("jsonlite", "digest", "wdsmatch"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing dependency: ", pkg)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Call as an Rscript file.")
script <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
historical <- file.path(root, "simulations/config/bounded_wdsm_calibration_sharded.json")
historical_sha <- "4e577017722c701844817f33e2c899c00c7e3b57bd8de58ede4097302812540f"
if (!identical(sha(historical), historical_sha)) stop("Historical config was changed.")
cfg <- jsonlite::read_json(historical, simplifyVector = TRUE)
admission_relative <- "simulations/config/bounded_wdsm_current_source_admission_20261001.json"
admission_path <- file.path(root, admission_relative)
admission_sha <- "c26c7aa53d28eab3bd03f1b55db4a632027c3bff8d7312cdea3a1922d428ae47"
if (!file.exists(admission_path) || !identical(sha(admission_path), admission_sha))
  stop("Current-source admission is missing or changed; independent review is required.")
admission <- jsonlite::read_json(admission_path, simplifyVector = TRUE)
if (!identical(admission$schema, "bounded_wdsm_current_source_admission_v1") ||
    !identical(admission$admission_id, "supplied-weight-current-source-20261001-v1") ||
    !identical(admission$historical_config$path,
               "simulations/config/bounded_wdsm_calibration_sharded.json") ||
    !identical(admission$historical_config$sha256, historical_sha) ||
    !identical(admission$package_version, cfg$package_version) ||
    !identical(admission$source_sha256, cfg$source_sha256))
  stop("Current-source admission identity or scientific-source pins differ.")
verify <- function(pins) {
  if (!is.list(pins) || !length(pins) || is.null(names(pins)) ||
      anyNA(names(pins)) || any(!nzchar(names(pins))) || anyDuplicated(names(pins)) ||
      any(!grepl("^[a-f0-9]{64}$", unlist(pins)))) stop("Malformed admitted source pins.")
  paths <- file.path(root, names(pins))
  if (!all(file.exists(paths)) ||
      !identical(unname(unlist(pins)), vapply(paths, sha, character(1L), USE.NAMES = FALSE)))
    stop("Pinned source mismatch; current-source admission requires independent review.")
}
verify(admission$source_sha256)
verify(admission$package_metadata_sha256)
verify(admission$package_source_sha256)
paths <- sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
relative <- paste0("R/", basename(paths))
if (!setequal(relative, names(admission$package_source_sha256)))
  stop("Unreviewed package R-source inventory.")

# Host-specific installed artifact hashes are resolved only after complete
# equality of the admitted source function bodies/formals and package metadata.
ns <- asNamespace("wdsmatch")
package_root <- find.package("wdsmatch")
if (!identical(as.character(utils::packageVersion("wdsmatch")), admission$package_version))
  stop("Installed package version mismatch.")
for (file in c("DESCRIPTION", "NAMESPACE"))
  if (!identical(sha(file.path(root, file)), sha(file.path(package_root, file)))) {
    if (file != "DESCRIPTION") stop("Installed/source NAMESPACE mismatch.")
    source_dcf <- read.dcf(file.path(root, file))
    installed_dcf <- read.dcf(file.path(package_root, file))
    normalize_dcf <- function(x) trimws(gsub("[[:space:]]+", " ", as.character(x)))
    if (!all(colnames(source_dcf) %in% colnames(installed_dcf)) ||
        !identical(normalize_dcf(source_dcf),
                   normalize_dcf(installed_dcf[, colnames(source_dcf), drop = FALSE])))
      stop("Installed/source DESCRIPTION mismatch.")
  }
env <- new.env(parent = ns)
for (path in paths) sys.source(path, env)
functions <- Filter(function(name) is.function(env[[name]]), ls(env, all.names = TRUE))
installed_functions <- Filter(function(name) is.function(get(name, ns)), ls(ns, all.names = TRUE))
if (!setequal(functions, installed_functions)) stop("Installed/source function name sets differ.")
for (name in functions)
  if (!exists(name, ns, inherits = FALSE) ||
      !identical(deparse(formals(env[[name]])), deparse(formals(get(name, ns)))) ||
      !identical(deparse(body(env[[name]])), deparse(body(get(name, ns)))))
    stop("Installed/source function mismatch: ", name)
package_files <- names(cfg$package_sha256)
cfg$package_sha256 <- as.list(setNames(vapply(file.path(package_root, package_files),
  sha, character(1L), USE.NAMES = FALSE), package_files))
cfg$package_source_sha256 <- admission$package_source_sha256[relative]
cfg$replay_provenance <- list(
  kind = "current_source_forward_run_not_historical_execution",
  historical_config_sha256 = historical_sha,
  current_source_admission = list(path = admission_relative, sha256 = admission_sha,
    admission_id = admission$admission_id),
  preparation_script_sha256 = sha(script),
  unchanged = c("scientific configuration", "seven simulation sources",
    "dataset and count seeds", "allocation", "M", "sample sizes", "B"),
  package_compatibility = admission$compatibility,
  scope = admission$scope,
  change = paste("The exact reviewed current 35-file inventory is admitted.",
    "Historical package pins/results remain unchanged. Source extensions and installed",
    "runtime hashes do not establish numerical parity or repeated-sampling calibration."),
  installed_functions_verified = functions,
  R_version = R.version.string, platform = R.version$platform,
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  simulations_run = FALSE, observations_generated = 0L, count_columns_generated = 0L,
  models_fitted = 0L, rng_streams_allocated = FALSE)
if (!identical(had_rng, exists(".Random.seed", .GlobalEnv, inherits = FALSE)) ||
    (had_rng && !identical(rng_before, get(".Random.seed", .GlobalEnv))))
  stop("Preparation unexpectedly changed RNG state.")
# Repeat the source gates before export; the output pins are the admitted values.
verify(admission$source_sha256)
verify(admission$package_metadata_sha256)
verify(admission$package_source_sha256)
if (!identical(sha(historical), historical_sha) ||
    !identical(sha(admission_path), admission_sha)) stop("Admission changed during preparation.")
temporary <- paste0(output, ".partial")
jsonlite::write_json(cfg, temporary, auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "null")
if (file.exists(output) || !file.rename(temporary, output))
  stop("Config finalization failed; partial file retained.")
cat("Created current-source forward config: ", normalizePath(output), "\n",
    "SHA256: ", sha(output), "\nNo simulation launched.\n", sep = "")
