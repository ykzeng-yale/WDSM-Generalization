# Direct execution of the existing statistical loop; no controller or scheduler.
lenis_run_preflight <- function(source_dir, population_dir, job_path, output, export_dir) {
  setup <- lenis_export_load(export_dir, "run")
  sources <- lenis_export_sources(source_dir, setup$export_dir)
  job <- lenis_export_job(job_path, setup$plan)
  cache <- lenis_export_cache(population_dir, job$population_scenarios,
    sources[["biosts-17039-File015.R"]])
  out <- lenis_export_output(output, c(source_dir, population_dir, setup$code_root, job_path))
  list(setup = setup, sources = sources, job = job, cache = cache, output = out)
}

lenis_run_comparison <- function(source_dir, population_dir, job_path, output, export_dir) {
  checked <- lenis_run_preflight(source_dir, population_dir, job_path, output, export_dir)
  setup <- checked$setup; sources <- checked$sources; job <- checked$job
  out <- checked$output; plan <- setup$plan; wm <- asNamespace("WeightedMatching")
  list2env(as.list(setup$functions), envir = environment())
  sha <- lenis_export_sha
  cache_files <- checked$cache$paths; cache_receipt <- checked$cache$receipt
  runtime <- list(id = "user_installed_namespace", R = as.character(getRversion()))
  namespaces <- setup$packages
  source_names <- c("entrypoint.R", "generate.R", "run.R", "report.R",
    "lenis_comparison_engine.R", "lenis_comparison_batch.R", "lenis_reporting.R",
    "lenis_canonical_counts.R", "design.json", "external_sources.json")
  files <- setNames(file.path(setup$export_dir, source_names), source_names)
  files <- c(files, paired_comparison_metrics.R = file.path(setup$export_dir, "..", "paired_comparison_metrics.R"),
    job.json = normalizePath(job_path, mustWork = TRUE))
  pins <- setNames(vapply(files, sha, character(1)), names(files))
  stopifnot(dir.create(out, recursive = FALSE, showWarnings = FALSE))
  dir.create(file.path(out, "frozen"))
  for (name in names(files)) stopifnot(file.copy(files[[name]], file.path(out, "frozen", name),
    overwrite = FALSE), sha(file.path(out, "frozen", name)) == pins[[name]])
  manifest <- list(status = "STARTED", plan = plan, source_pins = as.list(pins),
    external_source_pins = as.list(setNames(vapply(sources, sha, character(1)), names(sources))),
    runtime = runtime, job = job, namespaces = namespaces, cache_receipt = cache_receipt,
    started = format(Sys.time(), tz = "UTC", usetz = TRUE), session = capture.output(sessionInfo()))
  jsonlite::write_json(manifest, file.path(out, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  caches <- setNames(lapply(cache_files, readRDS), as.character(job$population_scenarios))
  for (s in seq_along(caches)) stopifnot(identical(caches[[1L]]$common[c("Cluster", "Strata")],
    caches[[s]]$common[c("Cluster", "Strata")]))
  native_ids <- as.vector(outer(plan$native_methods$PS_and_weight_transfer,
    plan$native_methods$outcome_adjustments, paste, sep = "_"))
  method_ids <- c(paste0("Lenis_", native_ids), plan$adapted_matchit_methods$names,
    as.vector(outer(plan$WM_methods$families, plan$M, function(x,y) paste0("WM_",x,"_M",y))))
  stopifnot(length(method_ids) == 25L, !anyDuplicated(method_ids))
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  setTimeLimit(cpu = job$max_seconds, elapsed = job$max_seconds, transient = FALSE)
  record <- function(s, r, mechanism, multiplier, method, target, n) data.frame(
    population_scenario = s, multiplier = multiplier, response = mechanism, replicate = r,
    scenario = "source_specification", estimand = "PATT", selected_n = 5000L, respondent_n = n,
    method = method, target = target, estimate = NA_real_, raw_estimate = NA_real_,
    variance = NA_real_, lower = NA_real_, upper = NA_real_, source_coverage = NA_real_,
    status = "failed", point_status = "unavailable", interval_status = "unavailable",
    warning_count = 0L, error = "Not evaluated", inference_scope = "", stringsAsFactors = FALSE)
  save_plain <- function(object, path) {
    lenis_assert_plain(object)
    saveRDS(object, path, version = 2)
  }
  records <- list(); timings <- list(); started <- proc.time()[[3L]]
  for (replicate in job$replicates) {
    sample <- lenis_source_sample(caches[[1L]], replicate,
      source = sources[["biosts-17039-File013.R"]])
    if (job$kind == "compatibility") stopifnot(identical(sample, reference$selected))
    save_plain(sample, file.path(out, paste0("selected_design_rep", replicate, ".rds")))
    for (s in job$population_scenarios) {
      input <- lenis_selected_from_cache(caches[[as.character(s)]], sample)
      for (mechanism in plan$response_models) {
        began <- proc.time()[[3L]]
        key <- paste0("s", s, "_rep", replicate, "_", mechanism)
        case_dir <- file.path(out, key); dir.create(case_dir)
        response <- lenis_capture(function() lenis_response_data(input$selected, mechanism))
        rows <- do.call(rbind, lapply(1:6, function(q) do.call(rbind, lapply(method_ids,
          function(k) record(s, replicate, mechanism, q, k, input$targets[q],
            if (response$ok) nrow(response$value$data) else NA_integer_)))))
        checkpoint <- lenis_compact_status(response)
        if (response$ok) {
          resp <- response$value; d <- resp$data
          Y <- input$Y[resp$selected_rows, , drop = FALSE]
          checkpoint$selected_ids <- resp$selected_ids
          checkpoint$respondent_rows <- resp$selected_rows
          checkpoint$respondent_ids <- d$id
          checkpoint$supplied_weights <- d$s.wt
          checkpoint$response_probability <- resp$response_probability
          checkpoint$response_fit <- lenis_numeric_fit(resp$response_fit)
          checkpoint$targets <- input$targets
          checkpoint$outcomes_sha256 <- digest::digest(Y, algo = "sha256")
          if (job$kind == "compatibility") {
            expected <- reference$cases[[key]]
            stopifnot(!is.null(expected), identical(resp$selected_ids, expected$selected_ids),
              identical(resp$selected_rows, expected$respondent_rows),
              identical(d$id, expected$respondent_ids), identical(d$z, expected$z),
              identical(input$targets, expected$targets),
              identical(checkpoint$outcomes_sha256, expected$outcomes_sha256))
          }
        } else rows$error <- paste("response stage:", response$error)
        save_plain(checkpoint, file.path(case_dir, "response.rds"))
        if (response$ok) {
          native <- lenis_capture(function() lenis_native_case(resp, input$targets[1L]))
          if (native$ok) {
            save_plain(lenis_compact_native(native$value), file.path(case_dir, "native_matching.rds"))
            nb <- lenis_capture(function() lenis_native_outcome_batch(native$value, resp, Y, input$targets))
            save_plain(if (nb$ok) nb$value else lenis_compact_status(nb), file.path(case_dir, "native_outcomes.rds"))
          } else {
            nb <- native
            save_plain(lenis_compact_status(native), file.path(case_dir, "native_matching.rds"))
            save_plain(lenis_compact_status(native), file.path(case_dir, "native_outcomes.rds"))
          }
          for (j in seq_along(plan$native_methods$PS_and_weight_transfer)) for (q in 1:6) {
            result <- if (nb$ok) nb$value[[j]][[q]] else lenis_compact_status(nb)
            for (adjustment in c("A", "U")) {
              k <- paste0("Lenis_", plan$native_methods$PS_and_weight_transfer[j], "_", adjustment)
              i <- which(rows$method == k & rows$multiplier == q)
              ps_kind <- substr(plan$native_methods$PS_and_weight_transfer[j], 1L, 1L)
              ps_warnings <- if (native$ok) length(native$value$ps[[ps_kind]]$warnings) else 0L
              rows$warning_count[i] <- length(result$warnings) + length(response$warnings) + ps_warnings
              rows$error[i] <- result$error
              rows$inference_scope[i] <- "Native Lenis outcome-regression survey SE, source 1.96 normal interval"
              if (!result$ok && !is.null(result$partial_results$values)) {
                partial <- result$partial_results$values
                est <- unname(partial[startsWith(names(partial), paste0("ATT|", adjustment, "|"))])
                if (length(est) == 1L && is.finite(est)) {
                  rows$estimate[i] <- est
                  rows$point_status[i] <- "ok"
                }
              }
              if (result$ok) {
                values <- result$values
                extract <- function(prefix) unname(values[startsWith(names(values), paste0(prefix, "|", adjustment, "|"))])
                est <- extract("ATT"); se <- extract("SE")
                stopifnot(length(est) == 1L, length(se) == 1L)
                rows$estimate[i] <- est; rows$variance[i] <- se^2
                rows$lower[i] <- est - 1.96 * se; rows$upper[i] <- est + 1.96 * se
                rows$source_coverage[i] <- extract("Cov")
                rows$status[i] <- rows$point_status[i] <- rows$interval_status[i] <- "ok"
              }
            }
          }
          for (weighted_ps in c(TRUE, FALSE)) {
            name <- if (weighted_ps) "weighted_PS" else "unit_PS"
            quick <- lenis_capture(function() lenis_quick_case(resp, weighted_ps))
            save_plain(lenis_compact_quick(quick), file.path(case_dir, paste0("quick_", name, "_matching.rds")))
            qb <- lenis_capture(function() lenis_quick_outcome_batch(quick, resp, Y, weighted_ps))
            save_plain(if (qb$ok) qb$value else lenis_compact_status(qb), file.path(case_dir, paste0("quick_", name, "_outcomes.rds")))
            for (q in 1:6) {
              result <- if (qb$ok) qb$value[[q]] else lenis_compact_status(qb)
              i <- which(rows$method == paste0("MatchIt_quick_", name) & rows$multiplier == q)
              rows$warning_count[i] <- length(result$warnings) +
                (if (q == 1L) 0L else length(quick$warnings)) + length(response$warnings)
              rows$error[i] <- result$error
              rows$inference_scope[i] <- "Adapted MatchIt quick ATT; subclass HC1 sandwich"
              if (result$ok) {
                v <- result$values
                rows$estimate[i] <- v$estimate; rows$variance[i] <- v$variance
                rows$lower[i] <- v$conf.int[1L]; rows$upper[i] <- v$conf.int[2L]
                rows$status[i] <- rows$point_status[i] <- rows$interval_status[i] <- "ok"
              }
            }
          }
          count_seed <- 805000000L + 100000L * s + 10L * replicate + match(mechanism, plan$response_models)
          suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding")); set.seed(count_seed)
          rng_before <- .Random.seed
          counts <- rmultinom(plan$B, nrow(d), rep(1/nrow(d), nrow(d)))
          count_record <- list(seed = count_seed, rng_kind = RNGkind(), rng_before = rng_before,
            rng_after = .Random.seed, n = nrow(d), B = plan$B,
            sha256 = digest::digest(counts, algo = "sha256"),
            column_sha256 = vapply(seq_len(ncol(counts)), function(b)
              digest::digest(counts[, b], algo = "sha256"), character(1)))
          save_plain(count_record, file.path(case_dir, "counts.rds"))
          if (job$kind == "compatibility") stopifnot(identical(count_record, reference$cases[[key]]$counts))
          candidate <- lenis_capture(function() lenis_wm_case(resp, wm, plan$M, counts))
          save_plain(if (candidate$ok) lenis_compact_wm(candidate$value) else lenis_compact_status(candidate),
            file.path(case_dir, "wm_matching_and_refits.rds"))
          batches <- list()
          for (family in plan$WM_methods$families) for (M in plan$M) {
            id <- paste0(family, "_M", M)
            batch <- if (!candidate$ok) candidate else lenis_capture(function()
              lenis_patt_outcome_batch(candidate$value$fits[[id]], Y, counts,
                candidate$value$refit_means[[if (family == "DSM") "DSM" else "full_x"]]))
            batches[[id]] <- if (batch$ok) batch$value else lenis_compact_status(batch)
            for (q in 1:6) {
              i <- which(rows$method == paste0("WM_", id) & rows$multiplier == q)
              group <- if (family == "DSM") "DSM" else "full_x"
              refit_warnings <- if (candidate$ok) sum(vapply(candidate$value$refit_diagnostics,
                function(x) length(x[[group]]$warnings), integer(1))) else 0L
              rows$warning_count[i] <- length(response$warnings) + length(candidate$warnings) +
                length(batch$warnings) + refit_warnings
              rows$error[i] <- batch$error
              rows$inference_scope[i] <- "Fixed-reuse prediction refit; supplied response weights frozen; external empirical comparison"
              if (batch$ok) {
                v <- batch$value
                rows$estimate[i] <- v$estimate[q]; rows$raw_estimate[i] <- v$raw_estimate[q]
                rows$variance[i] <- v$variance[q]; rows$lower[i] <- v$lower[q]; rows$upper[i] <- v$upper[q]
                rows$point_status[i] <- "ok"; rows$interval_status[i] <- v$status
                rows$status[i] <- if (v$status == "ok") "ok" else "failed"
                if (v$status != "ok") rows$error[i] <- paste("Failed prescribed count columns:", paste(v$failed_columns, collapse = ","))
              }
            }
          }
          save_plain(batches, file.path(case_dir, "wm_outcomes.rds"))
        }
        stopifnot(nrow(rows) == 150L, !anyDuplicated(rows[c("method", "multiplier")]))
        write.csv(rows, file.path(case_dir, "records.csv"), row.names = FALSE)
        records[[key]] <- rows
        timings[[key]] <- data.frame(case = key, seconds = proc.time()[[3L]] - began,
          records = nrow(rows), failed = sum(rows$status != "ok"),
          bytes = sum(file.info(list.files(case_dir, full.names = TRUE))$size))
        write.csv(do.call(rbind, timings), file.path(out, "timings.csv"), row.names = FALSE)
        cat(key, "completed", nrow(rows), "records; failures", sum(rows$status != "ok"),
            "seconds", timings[[key]]$seconds, "\n"); flush.console()
      }
    }
  }
  all <- do.call(rbind, records)
  expected_cases <- length(job$replicates) * length(job$population_scenarios) * length(plan$response_models)
  expected_records <- expected_cases * 150L
  stopifnot(nrow(all) == expected_records,
    identical(pins, setNames(vapply(files, sha, character(1)), names(files))))
  write.csv(all, file.path(out, "records.csv"), row.names = FALSE)
  receipt <- list(status = if (all(all$status == "ok")) "SHARD_COMPLETE" else "SHARD_COMPLETE_WITH_FAILURES_RETAINED",
    kind = job$kind, runtime_id = runtime$id,
    completed_cases = length(records), expected_cases = expected_cases,
    records = nrow(all), expected_records = expected_records,
    point_failures = sum(all$point_status != "ok"), interval_failures = sum(all$interval_status != "ok"),
    warning_records = sum(all$warning_count > 0), elapsed_seconds = proc.time()[[3L]] - started,
    output_bytes = sum(file.info(list.files(out, full.names = TRUE, recursive = TRUE))$size),
    source_pins_unchanged = TRUE, session = capture.output(sessionInfo()),
    scope = if (job$kind == "compatibility") "Saved-input runtime compatibility; no new Monte Carlo performance conclusion" else
      "Predeclared formal source-specification comparison; full aggregation and independent review required",
    finished = format(Sys.time(), tz = "UTC", usetz = TRUE))
  jsonlite::write_json(receipt, file.path(out, "receipt.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(receipt$status, receipt$records, "records in", receipt$elapsed_seconds, "seconds\n")
}

if (sys.nframe() == 0L) {
  file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  stopifnot(length(file) == 1L)
  export_dir <- dirname(normalizePath(sub("^--file=", "", file), mustWork = TRUE))
  source(file.path(export_dir, "entrypoint.R"))
  args <- commandArgs(TRUE)
  stopifnot(length(args) == 5L, args[1L] %in% c("preflight", "run"))
  if (args[1L] == "preflight") {
    lenis_run_preflight(args[2L], args[3L], args[4L], args[5L], export_dir)
    cat("PASS preflight only: no output creation, RNG, populations, samples, counts or fits\n")
  } else lenis_run_comparison(args[2L], args[3L], args[4L], args[5L], export_dir)
}
