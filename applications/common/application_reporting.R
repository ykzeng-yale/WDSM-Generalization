# SPDX-License-Identifier: GPL-3.0-only
# Descriptive reporting on existing point objects; no fits/counts/RNG/rematching.
wm_application_read_balance_schema <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required for schema JSON")
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

wm_application_balance <- function(result, schema) {
  a <- result$sample
  diagnostic_data <- a$data
  if (!is.data.frame(diagnostic_data) || nrow(diagnostic_data) != a$n)
    stop("Diagnostic data must retain exact prepared rows/order")
  diagnostics <- wb_diagnostics(diagnostic_data, schema)
  map <- diagnostics$mapping
  names(map)[names(map) == "kind"] <- "diagnostic_schema_kind"
  field <- function(name, default) vapply(schema, function(x)
    if (is.null(x[[name]])) default else x[[name]], character(1))
  map$diagnostic_family <- field("diagnostic_family", "declared_diagnostic")
  map$measurement_class <- field("measurement_class", "declared_diagnostic")
  map$declared_categories <- vapply(schema, function(x) if (is.null(x$allowed_codes)) "" else
    paste(unlist(x$allowed_codes, use.names = FALSE), collapse = ","), character(1))
  binary <- which(diagnostics$binary)
  map$unweighted_n_total <- map$unweighted_n_control <- map$unweighted_n_treated <- NA_integer_
  map$unweighted_n_total[binary] <- as.integer(colSums(diagnostics$X[, binary, drop = FALSE]))
  map$unweighted_n_control[binary] <- as.integer(colSums(diagnostics$X[a$Z == 0L, binary, drop = FALSE]))
  map$unweighted_n_treated[binary] <- as.integer(colSums(diagnostics$X[a$Z == 1L, binary, drop = FALSE]))
  diagnostics$mapping <- map
  reference <- setNames(lapply(c("PATE", "PATT"), function(e)
    wb_reference(diagnostics$X, a$Z, a$W, diagnostics$binary, e)), c("PATE", "PATT"))
  tables <- statuses <- values <- summaries <- list()
  requested <- result$requested
  for (i in seq_len(nrow(requested))) {
    item <- requested[i, , drop = FALSE]
    key <- paste0(item$family, "_M", item$M, "_", item$estimand)
    value <- app_capture(function() {
      point <- result$points[[key]]
      if (!isTRUE(point$ok)) stop("Saved point unavailable; no regeneration: ", point$error)
      wm_application_check_point(a, result$components, point, item$family, item$M, item$estimand)
      out <- wb_saved_graph(point$value, diagnostics$X, a$Z, a$W, reference[[item$estimand]])
      out$table$method_family <- item$family; out$table$M <- item$M; out$table$graph_key <- key
      out$table$sample_n <- a$n; out$table$sample_baseline <- a$baseline
      j <- match(out$table$variable, map$diagnostic); stopifnot(!anyNA(j))
      out$table <- cbind(out$table, map[j, setdiff(names(map), "diagnostic"), drop = FALSE])
      out
    })
    values[[key]] <- value
    if (value$ok) {
      tab <- value$value$table; tables[[key]] <- tab
      valid <- is.finite(tab$abs_SMD_before) & is.finite(tab$abs_SMD_after)
      summaries[[key]] <- data.frame(graph_key = key, n_diagnostics = nrow(tab),
        finite_SMD_pairs = sum(valid), undefined_SMD_after = sum(!is.finite(tab$SMD_after)),
        max_abs_SMD_before = if (any(valid)) max(tab$abs_SMD_before[valid]) else NA_real_,
        max_abs_SMD_after = if (any(valid)) max(tab$abs_SMD_after[valid]) else NA_real_,
        worsened = sum(tab$abs_SMD_after[valid] > tab$abs_SMD_before[valid]),
        worst_after = if (any(valid)) tab$variable[which(valid)[which.max(tab$abs_SMD_after[valid])]] else NA_character_)
    }
    statuses[[key]] <- data.frame(graph_key = key, family = item$family, M = item$M, estimand = item$estimand,
      status = if (value$ok) "saved_graph_balance_retained" else "unavailable_saved_graph_balance",
      error = if (value$ok) "" else value$error, stringsAsFactors = FALSE)
  }
  table <- if (length(tables)) do.call(rbind, tables) else data.frame()
  categories <- if (nrow(table)) table[table$diagnostic_family == "ordinal_category_proportion", , drop = FALSE] else table
  if (nrow(categories)) {
    categories$proportion_control_before <- categories$mean_control_before
    categories$proportion_treated_before <- categories$mean_treated_before
    categories$proportion_control_after <- categories$mean_control_after
    categories$proportion_treated_after <- categories$mean_treated_after
    categories$proportion_difference_before <- categories$difference_before
    categories$proportion_difference_after <- categories$difference_after
    categories$target_proportion <- categories$target_mean
  }
  list(status = do.call(rbind, statuses), tables = table, ordinal_category_proportions = categories,
    summary = if (length(summaries)) do.call(rbind, summaries) else data.frame(),
    graphs = values, reference = reference, diagnostic_mapping = map,
    scope = "Saved raw graph distributions; fixed target-specific pre-W SD; no causal/overlap/coverage certification",
    scientific_acceptance = FALSE)
}

wm_application_save <- function(result, new_directory, balance = NULL) {
  if (file.exists(new_directory) || !dir.exists(dirname(new_directory)) ||
      !dir.create(new_directory, recursive = FALSE)) stop("Fresh local output directory required; no overwrite/resume")
  files <- character()
  put <- function(name, value) {
    p <- file.path(new_directory, name); app_saved(p, value)
    files <<- c(files, p); invisible(p)
  }
  csv <- function(name, value) {
    p <- file.path(new_directory, name); app_csv(p, value)
    files <<- c(files, p); invisible(p)
  }
  tryCatch({
    put("prepared_input.rds", result$sample); put("shared_counts.rds", result$count_input)
    put("point_components.rds", result$components); put("shared_prediction_refits.rds", result$refits)
    put("requested_WM_rows.rds", result$requested)
    row_keys <- paste0(sub("^WM_", "", result$rows$method), "_M", result$rows$M, "_", result$rows$estimand)
    if (anyDuplicated(row_keys) || !setequal(row_keys, names(result$points))) stop("Saved row/point keys differ")
    for (key in names(result$points)) put(paste0(key, ".rds"),
      list(point = result$points[[key]], replication = result$replications[[key]],
        row = result$rows[match(key, row_keys), , drop = FALSE],
        count_validation = result$count_input$validation, counts_sha256 = result$count_input$sha256,
        scientific_acceptance = FALSE))
    csv("WM_comparisons.csv", result$rows)
    if (!is.null(balance)) {
      put("saved_graph_balance.rds", balance); csv("requested_graph_status.csv", balance$status)
      if (nrow(balance$tables)) csv("same_sample_weighted_balance.csv", balance$tables)
      if (nrow(balance$ordinal_category_proportions)) csv("ordinal_category_proportions.csv", balance$ordinal_category_proportions)
      if (nrow(balance$summary)) csv("balance_summary.csv", balance$summary)
    }
    put("phase_complete.rds", list(status = "APPLICATION_OUTPUTS_RETAINED_REVIEW_REQUIRED",
      WM_rows = nrow(result$rows), unavailable = result$rows[result$rows$status != "empirical_complete", , drop = FALSE],
      B_requested = 200L, counts_sha256 = result$count_input$sha256,
      sample_sha256 = result$sample$sample_sha256, reuse_binding = result$reuse_binding,
      nuisance_sha256 = result$sample$nuisance_sha256, nuisance_recipe = result$nuisance_recipe,
      output_pins = setNames(vapply(files, function(p) digest::digest(file = p, algo = "sha256", serialize = FALSE), character(1)), basename(files)),
      scientific_acceptance = FALSE))
  }, error = function(e) {
    app_saved(file.path(new_directory, "phase_failed.rds"), list(status = "FAILED_PARTIAL_OUTPUTS_RETAINED",
      error = conditionMessage(e), no_automatic_retry = TRUE, scientific_acceptance = FALSE))
    stop(e)
  })
  invisible(files)
}

# Unchanged local atomic no-overwrite IO helpers from the accepted engine.
app_saved <- function(path, value) {
  if (file.exists(path) || file.exists(paste0(path, ".tmp"))) stop("Refusing overwrite: ", path)
  saveRDS(value, paste0(path, ".tmp"), version = 3L)
  if (!file.rename(paste0(path, ".tmp"), path)) stop("Atomic save failed: ", path)
  invisible(path)
}
app_csv <- function(path, value) {
  if (file.exists(path) || file.exists(paste0(path, ".tmp"))) stop("Refusing overwrite: ", path)
  utils::write.csv(value, paste0(path, ".tmp"), row.names = FALSE)
  if (!file.rename(paste0(path, ".tmp"), path)) stop("Atomic CSV save failed: ", path)
}
