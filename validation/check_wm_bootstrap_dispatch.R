# Saved objects/counts only. No fitting, graph rebuilding, DGP or random draws.
# Args: package root, optional receipt RDS, optional isolated library path,
# optional existing original-survey integration-pilot results directory.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1L] else ".", mustWork = TRUE)
installed <- length(args) >= 3L && nzchar(args[3L])
if (installed) {
  lib <- normalizePath(args[3L], mustWork = TRUE)
  ns <- loadNamespace("WeightedMatching", lib.loc = lib)
  stopifnot(identical(normalizePath(getNamespaceInfo(ns, "path")),
    normalizePath(file.path(lib, "WeightedMatching"))))
  wm_bootstrap <- getExportedValue("WeightedMatching", "wm_bootstrap")
  wm_bootstrap_refit <- getExportedValue("WeightedMatching", "wm_bootstrap_refit")
} else {
  for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
}
labels <- character()
check <- function(ok, label) {
  if (!isTRUE(ok)) stop(label, call. = FALSE)
  labels <<- c(labels, label)
}
near <- function(x, y, label, tolerance = 2e-10) {
  check(identical(dim(x), dim(y)) && length(x) == length(y) &&
    all(is.finite(c(x, y))) &&
    all(abs(x - y) <= tolerance * pmax(1, abs(x), abs(y))), label)
}
fails <- function(expression, label, pattern = NULL) {
  result <- try(expression, silent = TRUE)
  check(inherits(result, "try-error") &&
    (is.null(pattern) || grepl(pattern, as.character(result), fixed = TRUE)), label)
}
had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
if (had_rng) rng_before <- get(".Random.seed", envir = .GlobalEnv)
start <- proc.time()[["elapsed"]]
forbidden <- function(...) stop("Forbidden fitting, graph rebuilding or random draw.")
wm_match <- wm_fit <- wm_wdsm_fit <- wdsm_fit_ps <- wm_scalar_logistic_match <- forbidden
runif <- sample <- rmultinom <- rnorm <- forbidden

check_fit <- function(fit, counts, label, refit = NULL, recorded = NULL) {
  before <- serialize(fit, NULL)
  counts_before <- serialize(counts, NULL)
  direct <- wm_bootstrap_refit(fit, counts, refit = refit, conf.level = .9)
  common <- wm_bootstrap(fit, method = "fixed_reuse", counts = counts,
    refit = refit, conf.level = .9)
  check(identical(common, direct), paste(label, "direct/common complete normal output"))
  check(identical(common$estimate, fit$estimate) && identical(class(common), class(direct)),
    paste(label, "same original point and result class"))
  explicit <- wm_bootstrap(fit, B = ncol(counts), method = "fixed_reuse",
    counts = counts, refit = refit, conf.level = .9)
  check(identical(explicit, direct), paste(label, "explicit B agrees with inferred B"))
  omitted <- wm_bootstrap(fit, method = "fixed_reuse", counts = counts,
    refit = refit, interval = "none", conf.level = .9)
  expected_none <- direct; expected_none$conf.int <- NULL
  check(identical(omitted, expected_none), paste(label, "none removes only conf.int"))
  near(common$variance, mean((common$draws - mean(common$draws))^2),
    paste(label, "original B-divisor variance"))
  near(common$root_n_draws, sqrt(fit$n) * (common$draws - fit$estimate),
    paste(label, "full observed n root scaling"))
  near(common$conf.int, fit$estimate + c(lower = -1, upper = 1) *
    stats::qnorm(.95) * sqrt(common$variance), paste(label, "source normal interval"))
  if (!is.null(recorded)) near(common$draws, recorded, paste(label, "existing saved refit draws"))
  row_labeled <- counts; rownames(row_labeled) <- as.character(seq_len(fit$n))
  check(identical(wm_bootstrap(fit, method = "fixed_reuse", counts = row_labeled,
    refit = refit, conf.level = .9), direct), paste(label, "literal original row labels"))
  check(identical(serialize(fit, NULL), before) &&
    identical(serialize(counts, NULL), counts_before), paste(label, "input objects unchanged"))
  list(n = fit$n, M = fit$M, estimand = fit$estimand,
    d = ncol(fit$graph$scores0), B = common$B, callback = !is.null(refit))
}

# Portable existing raw scalar fits, with the literal count transfers used by
# the preceding saved-array bootstrap validation. No model is refitted.
fixture <- readRDS(file.path(root, "validation", "fixtures", "scalar_replication_24.rds"))
cases <- list()
for (target in c("PATE", "PATT")) {
  fit <- fixture$fits[[target]]$fit
  n <- fit$n; counts <- matrix(1, n, 3L)
  i <- which(fit$data$Z == 0L)[1L]; j <- which(fit$data$Z == 1L)[1L]
  counts[i, 2L] <- 2; counts[j, 2L] <- 0
  counts[i, 3L] <- 0; counts[j, 3L] <- 2
  cases[[target]] <- check_fit(fit, counts, paste("portable", target))
  w <- fit$analysis_weights
  K <- fit$loads$incoming[cbind(seq_len(n), fit$data$Z + 1L)]
  expected <- vapply(seq_len(ncol(counts)), function(b) {
    m <- counts[, b]; z <- fit$data$Z; y <- fit$data$Y
    if (target == "PATE") sum(m * (2 * z - 1) * (w + K) * y) / sum(m * w)
    else (sum(m * z * w * y) - sum(m * (1 - z) * K * y)) / sum(m * z * w)
  }, numeric(1))
  near(wm_bootstrap(fit, method = "fixed_reuse", counts = counts)$draws,
    expected, paste(target, "independent raw random-denominator formula"))
  calls <- 0L
  callback <- function(m) {calls <<- calls + 1L; stop("must not execute")}
  fails(wm_bootstrap(fit, counts = counts, method = "fixed_reuse", refit = callback),
    paste(target, "raw fit rejects prediction refit"))
  check(calls == 0L, paste(target, "raw-fit rejection before callback"))
}

# Optional existing pilot contains genuine saved scalar PS, distinct-arm DSM,
# and full six-covariate graphs at M=1/3/5 for both estimands. Return predictions
# rebuilt from its stored coefficient arrays, never re-estimate a coefficient.
if (length(args) >= 4L && nzchar(args[4L])) {
  pilot <- normalizePath(args[4L], mustWork = TRUE)
  dirs <- c("GoodOverlap_retrospective_CorCor_PATE", "GoodOverlap_prospective_MisCor_PATT")
  for (directory in dirs) {
    path <- file.path(pilot, directory)
    weighted <- readRDS(file.path(path, "weighted.rds"))
    data <- readRDS(file.path(path, "data.rds"))
    saved_counts <- readRDS(file.path(path, "counts.rds"))
    counts <- saved_counts$counts[, 1:4, drop = FALSE]
    ps_terms <- c(paste0("X", 1:6),
      if (weighted$specification %in% c("CorCor", "CorMis")) "X1:X2")
    pg_terms <- c(paste0("X", 1:6),
      if (weighted$specification %in% c("CorCor", "MisCor")) "X1:X2")
    Dps <- stats::model.matrix(stats::reformulate(ps_terms), data)
    Dpg <- stats::model.matrix(stats::reformulate(pg_terms), data)
    x <- as.matrix(data[paste0("X", 1:6)])
    x <- sweep(x, 2L, colMeans(x), "-")
    x <- sweep(x, 2L, sqrt(colMeans(x^2)), "/")
    D6 <- cbind(1, x)
    for (j in 1:6) for (k in j:6) D6 <- cbind(D6, x[, j] * x[, k])
    for (family in c("PS", "DSM", "X6")) {
      predictions <- lapply(1:4, function(b) {
        saved <- weighted$refit_diagnostics[[b]]$fits[[family]]
        stopifnot(saved$ok)
        out <- list()
        theta <- saved$parameter
        ps <- if (family == "DSM") stats::plogis(Dps %*% theta[grep("^ps:", names(theta))]) else NULL
        for (z in if (weighted$estimand == "PATE") 0:1 else 0L) {
          if (family == "DSM") {
            label <- paste0("pg", z)
            pg <- as.vector(Dpg %*% theta[grep(paste0("^", label, ":"), names(theta))])
            a <- (as.numeric(ps) - theta[["center:ps"]]) / sqrt(theta[["variance:ps"]])
            g <- (pg - theta[[paste0("center:", label)]]) / sqrt(theta[[paste0("variance:", label)]])
            basis <- cbind(1, a, g, a^2, a * g, g^2)
            mean <- as.vector(basis %*% theta[grep(paste0("^bc", z, ":"), names(theta))])
          } else {
            design <- if (family == "PS") Dpg else D6
            mean <- as.vector(design %*% saved$coefficients[[paste0("arm", z)]])
          }
          out[[paste0("mean", z)]] <- mean
        }
        out
      })
      calls <- list()
      callback <- function(m) {
        b <- which(vapply(1:4, function(j) identical(unname(m), unname(counts[, j])), logical(1)))
        stopifnot(length(b) == 1L)
        calls[[length(calls) + 1L]] <<- m
        predictions[[b]]
      }
      for (M in c(1L, 3L, 5L)) {
        key <- paste0("WM_", family, "_M", M)
        point <- weighted$points[[key]]; stopifnot(point$ok)
        fit <- point$value
        d <- match(family, c("PS", "DSM", "X6")); d <- c(1L, 2L, 6L)[d]
        check(fit$M == M && ncol(fit$graph$scores0) == d,
          paste(directory, key, "saved M and actual graph dimension"))
        if (family == "DSM" && fit$estimand == "PATE")
          check(!identical(fit$graph$scores0, fit$graph$scores1), "saved DSM retains distinct arm maps")
        cases[[paste(directory, key)]] <- check_fit(fit, counts,
          paste(directory, key), callback, weighted$summaries[[key]]$draws[1:4])
        check_fit(fit, counts, paste(directory, key, "fixed predictions"))
        check(length(calls) == 20L && all(vapply(calls, function(m)
          any(vapply(1:4, function(b) identical(unname(m), unname(counts[, b])), logical(1))), logical(1))),
          paste(directory, key, "callback receives unchanged full-n columns"))
        calls <- list()
      }
    }
  }
}

fit <- fixture$fits$PATT$fit
counts <- matrix(1, fit$n, 3L)
calls <- 0L
callback <- function(m) {calls <<- calls + 1L; stop("must not execute")}
for (change in c("missing_counts", "B_mismatch", "B_one", "B_NULL", "seed", "basic",
                 "row_labels", "negative", "fractional", "nonfinite", "wrong_sum",
                 "wrong_rows", "empty_arm", "fitted_inference", "wrong_method", "contribution_refit")) {
  request <- list(object = fit, method = "fixed_reuse", counts = counts)
  if (change %in% c("missing_counts", "B_mismatch", "B_one", "B_NULL", "seed", "basic",
                    "row_labels", "fitted_inference", "wrong_method", "contribution_refit"))
    request$refit <- callback
  if (change == "missing_counts") request$counts <- NULL
  if (change == "B_mismatch") request$B <- 2L
  if (change == "B_one") {request$counts <- counts[, 1L, drop = FALSE]; request$B <- 1L}
  if (change == "B_NULL") request <- c(request, list(B = NULL))
  if (change == "seed") request$seed <- 7L
  if (change == "basic") request$interval <- "basic"
  if (change == "row_labels") rownames(request$counts) <- as.character(fit$n:1)
  if (change == "negative") request$counts[1, 1] <- -1
  if (change == "fractional") request$counts[1, 1] <- .5
  if (change == "nonfinite") request$counts[1, 1] <- NA_real_
  if (change == "wrong_sum") request$counts[1, 1] <- 2
  if (change == "wrong_rows") request$counts <- counts[-1, , drop = FALSE]
  if (change == "empty_arm") request$counts[, 1] <- fit$data$Z * fit$n / sum(fit$data$Z)
  if (change == "fitted_inference") class(request$object) <- c("wm_fitted_inference", class(fit))
  if (change == "wrong_method") request$method <- "unrecognized"
  if (change == "contribution_refit") request$method <- "contribution"
  pattern <- switch(change,
    missing_counts = "requires a supplied", B_mismatch = "B must be",
    B_one = "B must be", B_NULL = "B must be", seed = "Do not supply seed",
    basic = "basic is unsupported", row_labels = "Count row labels",
    fitted_inference = "ordinary self_normalized",
    contribution_refit = "refit is supported only", empty_arm = "Both arms",
    negative = "n-by-B", fractional = "n-by-B", nonfinite = "n-by-B",
    wrong_sum = "n-by-B", wrong_rows = "n-by-B", NULL)
  fails(do.call(wm_bootstrap, request), paste(change, "dispatch rejection"), pattern)
}
check(calls == 0L, "malformed requests rejected before callback")
fails(wm_bootstrap(fit, counts = counts), "ordinary wm_match default contribution does not silently choose refit")

# Same literal row/edge/training array as the existing contribution validator.
# Substitute a finite prescribed normal sequence temporarily, not an RNG draw.
with_literal_normals <- function(values, expression) {
  ns <- asNamespace("stats"); previous <- get("rnorm", ns)
  locked <- bindingIsLocked("rnorm", ns); cursor <- 0L
  literal <- function(n, mean = 0, sd = 1) {
    stopifnot(mean == 0, sd == 1)
    index <- seq.int(cursor + 1L, cursor + n); cursor <<- cursor + n
    stopifnot(cursor <= length(values)); values[index]
  }
  if (locked) unlockBinding("rnorm", ns)
  assign("rnorm", literal, ns); if (locked) lockBinding("rnorm", ns)
  on.exit({if (bindingIsLocked("rnorm", ns)) unlockBinding("rnorm", ns)
    assign("rnorm", previous, ns); if (locked) lockBinding("rnorm", ns)}, add = TRUE)
  result <- force(expression)
  check(cursor == length(values), "prescribed normal values consumed exactly")
  result
}
row <- c(1, 2, -1, 3); edge <- c(-2, .5); training <- c(-.4, .2, .8)
coefficients <- c(row, edge, training)
legacy <- structure(list(n = 4L, estimate = 1.25,
  root_n_variance = sum(coefficients^2) / 4,
  contributions = list(row = row, edge = edge, training = training)), class = c("wm_match", "list"))
values <- seq(-.9, 1.3, length.out = length(coefficients) * 2L)
default <- with_literal_normals(values, wm_bootstrap(legacy, B = 2L, chunk_size = 2L))
explicit <- with_literal_normals(values,
  wm_bootstrap(legacy, B = 2L, chunk_size = 2L, method = "contribution"))
check(identical(default, explicit), "default and explicit contribution outputs unchanged")
expected <- vapply(1:2, function(b) sum(coefficients *
  values[seq.int((b - 1L) * length(coefficients) + 1L, b * length(coefficients))]) / 2, numeric(1))
near(default$root_n_draws, expected, "unchanged Gaussian row/edge/training ordering")
check(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE) == had_rng &&
  (!had_rng || identical(rng_before, get(".Random.seed", envir = .GlobalEnv))), "RNG state unchanged")
check(length(tools::parse_Rd(file.path(root, "man", "wm_bootstrap.Rd"))) > 0, "Rd parses")
receipt <- list(status = "PASS", checks = length(labels), labels = labels, cases = cases,
  code_mode = if (installed) "isolated installed exports" else "workspace source",
  library_path = if (installed) lib else NULL,
  pilot_results = if (exists("pilot")) pilot else NULL,
  elapsed_seconds = proc.time()[["elapsed"]] - start, new_fits = 0L,
  graph_rebuilds = 0L, random_draws = 0L, new_datasets = 0L, RNG_changed = FALSE,
  scope = "Dispatch and original saved replicate algebra only; no sampling-validity or coverage claim")
if (length(args) >= 2L && nzchar(args[2L])) saveRDS(receipt, args[2L])
cat(length(labels), "checks PASS; existing objects and columns, no fits or RNG draws.\n")
