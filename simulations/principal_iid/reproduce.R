#!/usr/bin/env Rscript
# Portable single-process entry: truth calculation or one declared sample.
# It never dispatches jobs, installs packages, retries a failure or overwrites output.
principal_reproduce_main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (!length(argv) || !argv[1L] %in% c("truth", "case"))
    stop("Usage: reproduce.R truth NEW_DIR | case TRUTH_DIR OVERLAP N REPLICATE NEW_DIR")
  truth_mode <- identical(argv[1L], "truth")
  required_length <- if (truth_mode) 2L else 6L
  if (length(argv) != required_length) stop("Incorrect number of arguments.")
  if (!requireNamespace("digest", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE))
    stop("Install digest and jsonlite before running; this entry installs nothing.")
  if (!identical(as.character(getRversion()), "4.4.2"))
    stop("This exact source-prefix reproduction uses R 4.4.2.")
  threads <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS")
  if (any(Sys.getenv(threads) != "1")) stop("Set all documented numerical thread variables to 1 before R starts.")
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Invoke with Rscript --vanilla.")
  here <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
  package_root <- normalizePath(file.path(here, "../.."), mustWork = TRUE)
  sha <- function(p) digest::digest(file = p, algo = "sha256", serialize = FALSE)
  manifest_path <- file.path(here, "source_inventory.json")
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  manifest_sha <- sha(manifest_path)
  source_paths <- vapply(manifest$files, function(x) file.path(package_root, x$path), character(1L))
  source_hashes <- vapply(manifest$files, function(x) x$sha256, character(1L))
  verify_sources <- function() identical(sha(manifest_path), manifest_sha) &&
    identical(unname(vapply(source_paths, sha, character(1L))), unname(source_hashes))
  if (!verify_sources()) stop("The reproduction and common-package source inventory do not match.")
  model_env <- new.env(parent = globalenv())
  sys.source(file.path(here, "principal_iid_model.R"), model_env)
  model <- targets <- NULL
  input_paths <- character(); input_hashes <- character()
  if (!truth_mode) {
    truth_dir <- normalizePath(argv[2L], mustWork = TRUE)
    truth_receipt <- jsonlite::read_json(file.path(truth_dir, "result.json"), simplifyVector = FALSE)
    if (!identical(truth_receipt$schema, "principal_iid_public_truth_v1") ||
        !identical(truth_receipt$status, "COMPLETE_TRUTH") ||
        !identical(truth_receipt$source_inventory_sha256, manifest_sha) ||
        !isTRUE(truth_receipt$source_input_end_verified)) stop("Supply a complete truth result from the same source.")
    truth_names <- vapply(truth_receipt$files, function(x) x$path, character(1L))
    if (anyDuplicated(truth_names) || !setequal(truth_names, c("model.rds", "mixture_means.csv",
        "truth_orders.csv", "truth_components.csv", "truth_targets.csv", "truth.rds")))
      stop("Truth artifact inventory is incomplete.")
    for (item in truth_receipt$files) {
      p <- file.path(truth_dir, item$path)
      if (!identical(sha(p), item$sha256) || file.info(p)$size != item$bytes)
        stop("Truth artifact size/hash mismatch.")
      input_paths <- c(input_paths, p); input_hashes <- c(input_hashes, sha(p))
    }
    input_paths <- c(input_paths, file.path(truth_dir, "result.json"))
    input_hashes <- c(input_hashes, sha(tail(input_paths, 1L)))
    model <- readRDS(file.path(truth_dir, "model.rds"))
    targets <- utils::read.csv(file.path(truth_dir, "truth_targets.csv"), stringsAsFactors = FALSE)
    if (!grepl("^(1000|5000)$", argv[4L]) || !grepl("^[1-9][0-9]{0,3}$", argv[5L]))
      stop("Use n1000/5000 and replicate1:1000.")
  }
  output <- path.expand(tail(argv, 1L))
  if (file.exists(output) || !dir.exists(dirname(output))) stop("Supply a fresh output under an existing parent.")
  if (!dir.create(output, recursive = FALSE)) stop("Cannot reserve fresh output directory.")
  out <- normalizePath(output, mustWork = TRUE)
  write_json <- function(x, p) jsonlite::write_json(x, p, auto_unbox = TRUE,
    pretty = TRUE, null = "null", na = "null", digits = NA)
  began <- proc.time()[[3L]]; result <- NULL
  tryCatch({
    if (truth_mode) {
      model <- model_env$principal_iid_model()
      saveRDS(model, file.path(out, "model.rds"), version = 2L, compress = "xz")
      utils::write.csv(cbind(model$labels, model$means), file.path(out, "mixture_means.csv"), row.names = FALSE)
      truth <- model_env$principal_iid_truth(model, on_order = function(orders, components) {
        utils::write.csv(orders, file.path(out, "truth_orders.csv"), row.names = FALSE)
        utils::write.csv(components, file.path(out, "truth_components.csv"), row.names = FALSE)
      })
      utils::write.csv(truth$orders, file.path(out, "truth_orders.csv"), row.names = FALSE)
      utils::write.csv(truth$components, file.path(out, "truth_components.csv"), row.names = FALSE)
      utils::write.csv(truth$targets, file.path(out, "truth_targets.csv"), row.names = FALSE)
      saveRDS(truth, file.path(out, "truth.rds"), compress = "gzip", version = 3L)
      result <- list(schema = "principal_iid_public_truth_v1",
        status = if (truth$admitted) "COMPLETE_TRUTH" else "TRUTH_NOT_ADMITTED",
        admitted = truth$admitted, accuracy_scope = truth$accuracy_scope)
    } else {
      wm <- new.env(parent = globalenv())
      r_paths <- sort(list.files(file.path(package_root, "R"), "\\.R$", full.names = TRUE))
      declared_R <- sort(source_paths[grepl("/R/[^/]+\\.R$", source_paths)])
      if (!identical(r_paths, declared_R)) stop("The common-package R source inventory changed.")
      for (p in r_paths) sys.source(p, wm)
      original_helpers <- new.env(parent = globalenv())
      sys.source(file.path(package_root, "simulations/original_wdsm/original_survey_weighted_batch.R"), original_helpers)
      stack_helpers <- new.env(parent = globalenv())
      sys.source(file.path(package_root, "simulations/wm_saved_full_x_stack.R"), stack_helpers)
      recipes <- new.env(parent = globalenv())
      sys.source(file.path(here, "principal_iid_case.R"), recipes)
      result <- recipes$principal_iid_case(model, targets, argv[3L], as.integer(argv[4L]),
        as.integer(argv[5L]), out, wm, original_helpers, stack_helpers, model_env,
        progress = function(phase) { cat(phase, "\n"); flush.console() })
    }
  }, error = function(e) {
    result <<- list(schema = "principal_iid_public_failure_v1", status = "FAILED",
      error = conditionMessage(e), partial_artifacts_retained = TRUE)
  })
  result$source_input_end_verified <- tryCatch(verify_sources() &&
    identical(unname(vapply(input_paths, sha, character(1L))), unname(input_hashes)), error = function(e) FALSE)
  if (!result$source_input_end_verified) result$status <- "FAILED_INTEGRITY"
  result$source_inventory_sha256 <- manifest_sha
  result$elapsed_seconds <- proc.time()[[3L]] - began
  result$R_version <- R.version.string
  result$numerical_threads <- as.list(Sys.getenv(threads))
  result$inputs <- lapply(seq_along(input_paths), function(i) list(path = basename(input_paths[i]), sha256 = input_hashes[i]))
  members <- sort(list.files(out, recursive = TRUE, full.names = TRUE))
  members <- members[!file.info(members)$isdir]
  result$artifact_manifest <- lapply(members, function(p) list(
    path = substring(p, nchar(out) + 2L), sha256 = sha(p), bytes = unname(file.info(p)$size)))
  if (truth_mode) result$files <- result$artifact_manifest
  write_json(result, file.path(out, "result.json"))
  cat(result$status, "\n")
  if (!result$status %in% c("COMPLETE_TRUTH", "COMPLETE_CASE", "COMPLETE_WITH_STATISTICAL_FAILURES"))
    stop("Reproduction did not complete; retain its artifacts and error.", call. = FALSE)
  invisible(result)
}
if (sys.nframe() == 0L) principal_reproduce_main()
