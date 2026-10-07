# Load transferred statistical definitions and an explicitly supplied upstream.
# This loader does not generate observations/counts, fit models or run methods.
original_wdsm_output <- function(path, inputs) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || file.exists(path))
    stop("Output must be a new path.")
  link <- Sys.readlink(path)
  if (!is.na(link) && nzchar(link)) stop("Output must not be a symbolic link.")
  parent <- normalizePath(dirname(path), mustWork = TRUE)
  if (!dir.exists(parent) || file.access(parent, 2L) != 0L ||
      basename(path) %in% c(".", "..", "")) stop("Writable existing output parent required.")
  output <- file.path(parent, basename(path))
  for (input in inputs) {
    input <- normalizePath(input, mustWork = TRUE)
    if (output == input || startsWith(output, paste0(input, "/")) ||
        (dir.exists(input) && startsWith(input, paste0(output, "/"))))
      stop("Output must be outside the package and supplied input directories.")
  }
  output
}

original_wdsm_load <- function(directory, upstream) {
  directory <- normalizePath(directory, mustWork = TRUE)
  upstream <- normalizePath(upstream, mustWork = TRUE)
  for (package in c("digest", "jsonlite", "sampling", "WeightedMatching"))
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
  sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
  provenance <- jsonlite::read_json(file.path(directory, "provenance.json"))
  for (entry in provenance$statistical_modules)
    if (!identical(sha(file.path(directory, entry$file)), entry$sha256))
      stop("Transferred statistical source differs: ", entry$file)
  pins <- jsonlite::read_json(file.path(directory, "upstream_source_pins.json"))$files
  for (path in names(pins))
    if (!file.exists(file.path(upstream, path)) ||
        !identical(sha(file.path(upstream, path)), pins[[path]]))
      stop("Upstream source differs or is missing: ", path)
  engine <- new.env(parent = globalenv())
  for (name in c("weighted_batch", "native_adapter", "analytic_adapter",
                 "reporting", "lossless_storage", "execution_io", "sample_engine"))
    sys.source(file.path(directory, paste0("original_survey_", name, ".R")), engine)
  sys.source(file.path(directory, "..", "paired_comparison_metrics.R"), engine)
  original <- new.env(parent = globalenv())
  sys.source(file.path(upstream, "R", "wdsm_core_final.R"), original)
  sys.source(file.path(upstream, "R", "comparators_final.R"), original)
  original$comparator_final_check_dependencies()
  wm <- asNamespace("WeightedMatching")
  expected_version <- unname(read.dcf(
    file.path(directory, "..", "..", "DESCRIPTION"), fields = "Version")[1L, 1L])
  if (!identical(as.character(utils::packageVersion("WeightedMatching")), expected_version))
    stop("Install WeightedMatching version ", expected_version, " from this source checkout.")
  required <- c("wm_match", ".wm_wdsm_controls", ".wm_wdsm_fit_stack",
                ".wm_wdsm_prediction_stack", ".wm_wdsm_nuisance_stack",
                ".wm_wdsm_fitted_pipeline")
  if (!all(vapply(required, exists, logical(1), envir = wm, inherits = FALSE)))
    stop("Installed package lacks required original-study functions.")
  list(engine = engine, original = original, wm = wm, upstream = upstream)
}
