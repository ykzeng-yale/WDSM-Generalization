# Lenis-specific reporting contract and saved-record provenance adapter.
# No population generation, sampling, estimation, or count generation occurs here.
lenis_reporting_plan <- function(config, requested = config$formal_requested_R) {
  native <- as.vector(outer(config$native_methods$PS_and_weight_transfer,
                            config$native_methods$outcome_adjustments, paste, sep = "_"))
  methods <- c(paste0("Lenis_", native), config$adapted_matchit_methods$names,
    as.vector(outer(config$WM_methods$families, config$M,
      function(f, M) paste0("WM_", f, "_M", M))))
  stopifnot(length(methods) == 25L, !anyDuplicated(methods), length(requested) == 1L,
    is.finite(requested), requested >= 2L, requested <= config$formal_requested_R,
    requested == floor(requested))
  plan <- expand.grid(population_scenario = config$population_scenarios,
    multiplier = config$effect_multipliers, response = config$response_models,
    method = methods, stringsAsFactors = FALSE)
  plan$estimand <- "PATT"; plan$scenario <- "source_specification"
  plan$requested <- as.integer(requested); plan$interval_required <- TRUE
  plan$coverage_rule <- ifelse(startsWith(plan$method, "Lenis_"), "strict", "inclusive")
  plan
}

# The exact 22-column record() schema from the frozen Lenis shard driver.
# Explicit character classes retain "" as an empty string; the literal CSV
# missing-value token NA remains missing. Numeric fields use strict parsers.
lenis_saved_record_schema <- function() c(
  population_scenario = "integer", multiplier = "integer", response = "character",
  replicate = "integer", scenario = "character", estimand = "character",
  selected_n = "integer", respondent_n = "integer", method = "character",
  target = "numeric", estimate = "numeric", raw_estimate = "numeric",
  variance = "numeric", lower = "numeric", upper = "numeric", source_coverage = "numeric",
  status = "character", point_status = "character", interval_status = "character",
  warning_count = "integer", error = "character", inference_scope = "character")

lenis_read_saved_records <- function(path) {
  schema <- lenis_saved_record_schema()
  header <- names(read.csv(path, nrows = 0L, check.names = FALSE,
    colClasses = "character", stringsAsFactors = FALSE, row.names = NULL))
  stopifnot(identical(header, names(schema)))
  records <- read.csv(path, check.names = FALSE, colClasses = unname(schema),
    stringsAsFactors = FALSE, na.strings = "NA", row.names = NULL, fill = FALSE)
  stopifnot(identical(names(records), names(schema)))
  records
}

lenis_saved_reporting_records <- function(run, source_plan, accepted_source_pins) {
  sha <- function(f) digest::digest(file = f, algo = "sha256", serialize = FALSE)
  input_pins <- character()
  read_artifact <- function(path, type) {
    stopifnot(file.exists(path)); input_pins[[path]] <<- sha(path)
    switch(type, json = jsonlite::read_json(path, simplifyVector = TRUE),
      rds = readRDS(path), csv = lenis_read_saved_records(path))
  }
  receipt <- read_artifact(file.path(run, "receipt.json"), "json")
  manifest <- read_artifact(file.path(run, "manifest.json"), "json")
  stopifnot(receipt$status %in% c("PILOT_COMPLETE", "PILOT_COMPLETE_WITH_FAILURES_RETAINED",
    "SHARD_COMPLETE", "SHARD_COMPLETE_WITH_FAILURES_RETAINED"),
    identical(manifest$plan, source_plan),
    all(names(accepted_source_pins) %in% names(manifest$source_pins)))
  for (f in names(accepted_source_pins)) stopifnot(
    manifest$source_pins[[f]] == accepted_source_pins[[f]])
  for (f in names(manifest$source_pins)) {
    p <- file.path(run, "frozen", f); input_pins[[p]] <- sha(p)
    stopifnot(input_pins[[p]] == manifest$source_pins[[f]])
  }
  if (startsWith(receipt$status, "PILOT")) {
    scenarios <- source_plan$population_scenarios
    reps <- source_plan$local_pilot_replicates
  } else {
    scenarios <- manifest$job$population_scenarios; reps <- manifest$job$replicates
    stopifnot(identical(receipt$kind, manifest$job$kind))
  }
  stopifnot(all(scenarios %in% source_plan$population_scenarios),
    all(reps %in% seq_len(source_plan$formal_requested_R)),
    !anyDuplicated(scenarios), !anyDuplicated(reps))
  method_plan <- lenis_reporting_plan(source_plan)
  expected <- expand.grid(population_scenario = scenarios, replicate = reps,
    response = source_plan$response_models, multiplier = source_plan$effect_multipliers,
    method = unique(method_plan$method), stringsAsFactors = FALSE)
  key <- function(x) do.call(paste, c(unname(x[c("population_scenario", "replicate", "response",
                                               "multiplier", "method")]), list(sep = "|")))
  records <- read_artifact(file.path(run, "records.csv"), "csv")
  stopifnot(nrow(records) == nrow(expected), !anyDuplicated(key(records)),
    setequal(key(records), key(expected)), receipt$records == nrow(records),
    receipt$completed_cases == length(scenarios) * length(reps) * length(source_plan$response_models),
    all(records$estimand == "PATT"), all(records$scenario == "source_specification"),
    all(records$selected_n == source_plan$selected_n),
    sum(records$point_status != "ok") == receipt$point_failures,
    sum(records$interval_status != "ok") == receipt$interval_failures,
    sum(records$warning_count > 0L) == receipt$warning_records)
  records$pairing_id <- NA_character_
  for (r in reps) {
    design <- read_artifact(file.path(run, paste0("selected_design_rep", r, ".rds")), "rds")
    stopifnot(design$seed == r, length(design$ids) == source_plan$selected_n, !anyDuplicated(design$ids))
    for (s in scenarios) for (mechanism in source_plan$response_models) {
      case <- paste0("s", s, "_rep", r, "_", mechanism)
      response <- read_artifact(file.path(run, case, "response.rds"), "rds")
      case_records <- read_artifact(file.path(run, case, "records.csv"), "csv")
      ix <- which(records$population_scenario == s & records$replicate == r & records$response == mechanism)
      stopifnot(nrow(case_records) == 150L, !anyDuplicated(key(case_records)),
        setequal(key(case_records), key(records[ix, ])))
      from_global <- records[match(key(case_records), key(records)), setdiff(names(records), "pairing_id"), drop = FALSE]
      rownames(from_global) <- rownames(case_records) <- NULL
      stopifnot(length(ix) == 150L, identical(case_records, from_global))
      cache_sha <- manifest$cache_receipt$cache_hashes[[paste0("scenario_", s, ".rds")]]
      stopifnot(is.character(cache_sha), length(cache_sha) == 1L, nchar(cache_sha) == 64L)
      # The pinned population contains all literal potential outcomes and response
      # indicators. Together with the actual design and response mechanism, it
      # determines the shared realized data even if a response fit failed.
      records$pairing_id[ix] <- digest::digest(list(population_sha256 = cache_sha,
        selected_design = design, response_model = mechanism), algo = "sha256")
      if (isTRUE(response$ok)) {
        target <- as.numeric(response$targets[records$multiplier[ix]])
        # write.csv rounds doubles to its decimal representation; RDS retains
        # their binary values. This bound covers decimal serialization only.
        stopifnot(identical(response$selected_ids, design$ids),
          identical(response$respondent_ids, design$ids[response$respondent_rows]),
          all(records$respondent_n[ix] == length(response$respondent_ids)),
          all(abs(records$target[ix] - target) <= 32 * .Machine$double.eps * pmax(1, abs(target))))
      }
    }
  }
  stopifnot(!anyNA(records$pairing_id))
  list(records = records, input_pins = as.list(input_pins), receipt = receipt,
       manifest = manifest, source = run)
}
