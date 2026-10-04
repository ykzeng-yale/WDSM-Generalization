# SPDX-License-Identifier: GPL-3.0-only
# Pure saved-graph balance binding; source defines functions only.
# Forty-nine diagnostic columns do not change the 35 original metric columns.

meps_diagnostic_matrix <- function(a, s) {
  columns <- list(); mapping <- list(); binary <- logical()
  for (field in s$fields) {
    name <- field$source_field
    if (identical(field$type, "continuous")) {
      stopifnot(identical(unlist(field$columns, use.names = FALSE), name))
      columns[[name]] <- a$X[, name]
      mapping[[name]] <- data.frame(variable = name, source_field = name,
        category = NA_integer_, diagnostic_kind = "continuous",
        original_matching_coordinate = TRUE, stringsAsFactors = FALSE)
      binary <- c(binary, FALSE)
    } else {
      stopifnot(identical(field$type, "categorical"), field$reference_level == 1L)
      levels <- as.integer(unlist(field$declared_source_levels, use.names = FALSE))
      declared <- paste0(name, "__", levels[levels != 1L])
      stopifnot(identical(unlist(field$columns, use.names = FALSE), declared))
      dummy <- a$X[, declared, drop = FALSE]
      stopifnot(all(dummy %in% c(0, 1)), all(rowSums(dummy) <= 1))
      for (level in levels) {
        label <- paste0(name, "__", level)
        columns[[label]] <- if (level == 1L) 1 - rowSums(dummy) else dummy[, label]
        mapping[[label]] <- data.frame(variable = label, source_field = name,
          category = level, diagnostic_kind = "category_proportion",
          original_matching_coordinate = level != 1L, stringsAsFactors = FALSE)
        binary <- c(binary, TRUE)
      }
    }
  }
  X <- do.call(cbind, columns)
  map <- do.call(rbind, mapping); rownames(map) <- NULL
  stopifnot(nrow(X) == a$n, ncol(X) == 49L, sum(binary) == 44L,
    !anyDuplicated(colnames(X)), all(is.finite(X)),
    sum(map$original_matching_coordinate) == 35L)
  list(X = X, binary = binary, mapping = map)
}

# Authenticate saved lambda and K without recomputing any neighbor or distance.
meps_graph_binding <- function(point, a, target, family) {
  stopifnot(inherits(point, "wm_match"), point$n == a$n, point$M == 3L,
    identical(point$method, "self_normalized"), identical(point$estimand, target),
    identical(point$graph$tie_rule, "Exact distance, then original row index"),
    identical(as.numeric(point$data$Y), as.numeric(a$Y)),
    identical(as.integer(point$data$Z), a$Z),
    identical(as.numeric(point$weights), as.numeric(a$W)),
    isTRUE(point$info$corrected))
  d <- switch(family, PS = 1L, DSM = 2L, full_X = 35L)
  for (z in if (target == "PATE") 0:1 else 0L)
    stopifnot(identical(dim(point$graph[[paste0("scores", z)]]), c(a$n, d)))
  queries <- if (target == "PATE") seq_len(a$n) else which(a$Z == 1L)
  edges <- point$graph$edges
  stopifnot(nrow(edges) == length(queries) * 3L,
    identical(as.integer(edges$query), rep(queries, each = 3L)),
    all(edges$donor >= 1L & edges$donor <= a$n),
    all(edges$donor == as.integer(edges$donor)),
    all(a$Z[edges$donor] == 1L - a$Z[edges$query]),
    all(edges$arm == a$Z[edges$donor]),
    all(lengths(point$graph$neighbors[queries]) == 3L),
    identical(as.integer(edges$donor), as.integer(unlist(point$graph$neighbors[queries], use.names = FALSE))),
    !anyDuplicated(paste(edges$query, edges$donor)),
    all(is.finite(edges$share)), all(edges$share > 0))
  w <- a$W / point$weight_scale
  donor_w <- matrix(w[edges$donor], ncol = 3L, byrow = TRUE)
  lambda <- as.vector(t(donor_w / rowSums(donor_w)))
  agree <- function(x, y) {
    stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)))
    error <- max(abs(x - y) / pmax(1, abs(x), abs(y)))
    if (error > 1e-10) stop("Saved graph fraction/load identity failed")
    error
  }
  fraction_error <- agree(edges$share, lambda)
  outcome_share_error <- agree(edges$outcome_share, edges$share)
  K <- matrix(0, a$n, 2L)
  for (z in if (target == "PATE") 0:1 else 0L) {
    take <- edges$arm == z
    summed <- rowsum(w[edges$query[take]] * edges$share[take],
                    edges$donor[take], reorder = FALSE)
    K[as.integer(rownames(summed)), z + 1L] <- summed[, 1L]
  }
  load_error <- agree(as.numeric(point$loads$incoming), as.numeric(K))
  list(fraction_error = fraction_error, outcome_share_error = outcome_share_error,
       incoming_load_error = load_error)
}

meps_saved_balance <- function(a, schema, points, balance) {
  required <- as.vector(outer(c("PS", "DSM", "full_X"), c("PATE", "PATT"), paste, sep = "_"))
  stopifnot(is.list(points), setequal(names(points), required), !anyDuplicated(names(points)))
  before_points <- digest::digest(points, algo = "sha256")
  had_seed <- exists(".Random.seed", globalenv(), inherits = FALSE)
  seed_before <- if (had_seed) get(".Random.seed", globalenv(), inherits = FALSE) else NULL
  kind_before <- RNGkind()
  diagnostics <- meps_diagnostic_matrix(a, schema)
  refs <- stats::setNames(lapply(c("PATE", "PATT"), function(e)
    balance$wb_reference(diagnostics$X, a$Z, a$W, diagnostics$binary, e)), c("PATE", "PATT"))
  records <- summaries <- statuses <- bindings <- list()
  pair <- a$baseline
  for (family in c("PS", "DSM", "full_X")) for (target in c("PATE", "PATT")) {
    name <- paste(family, target, sep = "_"); key <- paste(pair, name, sep = "_")
    saved <- points[[name]]
    if (!isTRUE(saved$ok)) {
      statuses[[key]] <- data.frame(key = key, status = "SAVED_POINT_UNAVAILABLE",
        error = saved$error, diagnostic_rows = 0L, stringsAsFactors = FALSE)
      next
    }
    point <- saved$value
    bindings[[key]] <- meps_graph_binding(point, a, target, family)
    result <- balance$wb_saved_graph(point, diagnostics$X, a$Z, a$W, refs[[target]])
    stopifnot(max(abs(result$mass$relative_residual_each_arm)) <= 1e-10)
    tab <- cbind(pair = pair, family = family, M = 3L, n = a$n, n_white = sum(a$Z),
      result$table, stringsAsFactors = FALSE)
    at <- match(tab$variable, diagnostics$mapping$variable); stopifnot(!anyNA(at))
    tab <- cbind(tab, diagnostics$mapping[at, -1L, drop = FALSE])
    tab$unweighted_n_total <- colSums(diagnostics$X == 1)[at]
    tab$unweighted_n_control <- colSums(diagnostics$X[a$Z == 0L, , drop = FALSE] == 1)[at]
    tab$unweighted_n_White <- colSums(diagnostics$X[a$Z == 1L, , drop = FALSE] == 1)[at]
    tab$unweighted_n_total[tab$diagnostic_kind == "continuous"] <- NA_integer_
    tab$unweighted_n_control[tab$diagnostic_kind == "continuous"] <- NA_integer_
    tab$unweighted_n_White[tab$diagnostic_kind == "continuous"] <- NA_integer_
    records[[key]] <- tab
    defined <- is.finite(tab$SMD_before) & is.finite(tab$SMD_after)
    summaries[[key]] <- data.frame(pair = pair, family = family, estimand = target,
      M = 3L, diagnostics = nrow(tab), undefined = sum(!defined),
      max_abs_SMD_before = if (all(defined)) max(tab$abs_SMD_before) else NA_real_,
      max_abs_SMD_after = if (all(defined)) max(tab$abs_SMD_after) else NA_real_,
      worsened_diagnostics = if (all(defined)) sum(tab$abs_SMD_after > tab$abs_SMD_before) else NA_integer_,
      worst_after = if (all(defined)) tab$variable[which.max(tab$abs_SMD_after)] else NA_character_,
      control_post_ESS = result$after$control$ESS, White_post_ESS = result$after$treated$ESS,
      maximum_relative_mass_error = max(abs(result$mass$relative_residual_each_arm)),
      stringsAsFactors = FALSE)
    statuses[[key]] <- data.frame(key = key, status = "SAVED_GRAPH_BALANCE_AVAILABLE",
      error = "", diagnostic_rows = nrow(tab), stringsAsFactors = FALSE)
  }
  status_table <- do.call(rbind, statuses); rownames(status_table) <- NULL
  stopifnot(identical(before_points, digest::digest(points, algo = "sha256")),
    identical(RNGkind(), kind_before),
    identical(exists(".Random.seed", globalenv(), inherits = FALSE), had_seed))
  if (had_seed) stopifnot(identical(get(".Random.seed", globalenv(), inherits = FALSE), seed_before))
  list(statuses = status_table, covariate_balance = if (length(records)) do.call(rbind, records) else NULL,
    balance_summary = if (length(summaries)) do.call(rbind, summaries) else NULL,
    bindings = bindings, diagnostic_mapping = diagnostics$mapping,
    status = if (all(status_table$status == "SAVED_GRAPH_BALANCE_AVAILABLE"))
      "COMPLETE_DESCRIPTIVE_BALANCE" else "INCOMPLETE_DESCRIPTIVE_BALANCE",
    original_metric_columns = 35L, diagnostic_columns = 49L,
    RNG_unchanged = TRUE, no_fit_rematch_bootstrap = TRUE, scientific_acceptance = FALSE)
}
