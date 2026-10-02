# Saved inference rows and literal columns only; no fits, DGP or random draws.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1L] else ".", mustWork = TRUE)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
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
fails <- function(expr, label) check(inherits(try(expr, silent = TRUE), "try-error"), label)
had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
if (had_rng) rng_before <- get(".Random.seed", envir = .GlobalEnv)
forbidden <- function(...) stop("Forbidden statistical fit, rematching or RNG call.")
wdsm_fit_ps <- wm_match <- .wm_scalar_nw <- forbidden
rnorm <- runif <- sample <- rmultinom <- forbidden
start <- proc.time()[["elapsed"]]
fixture <- readRDS(file.path(root, "validation", "fixtures", "scalar_replication_24.rds"))
inference <- if (length(args) >= 3L) readRDS(args[3L])$results else
  lapply(fixture$fits, wm_fitted_inference, covariance_scope = "scalar_psm")

# Exercise the existing Gaussian execution path with prescribed literal
# multiplier values in this process. No normal random generator is called.
with_literal_normals <- function(values, expression) {
  namespace <- asNamespace("stats")
  original <- get("rnorm", envir = namespace)
  original_lock <- bindingIsLocked("rnorm", namespace)
  cursor <- 0L
  literal <- function(n, mean = 0, sd = 1) {
    stopifnot(length(n) == 1L, n >= 0, n == floor(n), mean == 0, sd == 1)
    end <- cursor + n
    if (end > length(values)) stop("Prescribed normal sequence exhausted.")
    selected <- if (n) values[seq.int(cursor + 1L, end)] else numeric()
    cursor <<- end
    selected
  }
  if (original_lock) unlockBinding("rnorm", namespace)
  assign("rnorm", literal, envir = namespace)
  if (original_lock) lockBinding("rnorm", namespace)
  on.exit({
    if (bindingIsLocked("rnorm", namespace)) unlockBinding("rnorm", namespace)
    assign("rnorm", original, envir = namespace)
    if (original_lock) lockBinding("rnorm", namespace)
  }, add = TRUE)
  result <- force(expression)
  check(cursor == length(values), "literal Gaussian path consumes exact prescribed sequence")
  result
}

results <- list()
for (target in c("PATE", "PATT")) {
  x <- inference[[target]]
  check(x$available && identical(x$fit, fixture$fits[[target]]$fit),
    paste(target, "existing saved scalar inference bound to fixture"))
  before <- serialize(x, NULL)
  n <- x$n
  rows <- x$augmented_rows - mean(x$augmented_rows)
  i <- which(x$fit$data$Z == 0L)[1L]
  j <- which(x$fit$data$Z == 1L)[1L]
  counts <- matrix(1, n, 4L)
  counts[i, 2L] <- 2
  counts[j, 2L] <- 0
  counts[i, 3L] <- 0
  counts[j, 3L] <- 2
  counts[, 4L] <- 0
  counts[i, 4L] <- n
  expected <- vapply(seq_len(ncol(counts)), function(b)
    sum((counts[, b] - 1) * rows) / sqrt(n), numeric(1))
  result <- wm_bootstrap(x, counts = counts)
  check(result$B == 4L && result$n == n &&
    result$draw_input == "supplied_full_n_multinomial_counts" &&
    !result$supplied_draw_law_verified && !result$population_assumptions_verified &&
    !result$original_refit_bootstrap, paste(target, "full-n count scope and inferred B"))
  near(result$root_n_draws, expected, paste(target, "literal counts-minus-one root algebra"))
  near(result$draws, x$estimate + expected / sqrt(n), paste(target, "same raw point draw algebra"))
  near(result$root_n_draws[4L], sqrt(n) * rows[i],
    paste(target, "all counts on one control row admitted without arm denominator"))
  near(result$conditional_root_n_variance, mean(rows^2),
    paste(target, "analytic contribution variance"))
  near(result$conditional_variance, mean(rows^2) / n,
    paste(target, "full observed n for sampling variance"))
  near(result$monte_carlo_root_n_variance, stats::var(expected),
    paste(target, "separate empirical B-minus-one draw diagnostic"))
  near(result$conf.int, x$estimate + c(lower = -1, upper = 1) *
    stats::qnorm(.025, lower.tail = FALSE) * sqrt(mean(rows^2) / n),
    paste(target, "normal interval uses analytic variance"))
  check(identical(result$source_inference, x) && identical(serialize(x, NULL), before),
    paste(target, "source inference fully preserved"))
  one <- wm_bootstrap(x, counts = counts[, 1L, drop = FALSE])
  check(one$B == 1L && is.na(one$monte_carlo_root_n_variance) &&
    is.na(one$monte_carlo_variance), paste(target, "one-column diagnostic unavailable"))
  near(one$conditional_root_n_variance, result$conditional_root_n_variance,
    paste(target, "one-column retains analytic variance"))
  near(one$conf.int, result$conf.int, paste(target, "normal interval independent of count columns"))
  same <- wm_bootstrap(x, counts = matrix(1, n, 2L))
  check(same$monte_carlo_root_n_variance == 0 && same$conditional_root_n_variance > 0,
    paste(target, "identical supplied draws cannot replace positive analytic variance"))
  explicit <- wm_bootstrap(x, B = 4L, counts = counts, chunk_size = 1L)
  near(explicit$root_n_draws, result$root_n_draws, paste(target, "explicit B and small row chunks"))
  basic <- wm_bootstrap(x, counts = counts, interval = "basic", conf.level = .8)
  near(basic$conf.int, x$estimate - as.numeric(stats::quantile(expected,
    c(.9, .1), names = FALSE)) / sqrt(n), paste(target, "literal basic interval quantiles"))
  check(is.null(wm_bootstrap(x, counts = counts, interval = "none")$conf.int),
    paste(target, "no interval branch"))
  values <- seq(-1.2, .9, length.out = n * 3L)
  gaussian <- with_literal_normals(values,
    wm_bootstrap(x, B = 3L, interval = "normal", chunk_size = 5L))
  gaussian_expected <- vapply(seq_len(3L), function(b)
    sum(rows * values[seq.int((b - 1L) * n + 1L, b * n)]) / sqrt(n), numeric(1))
  near(gaussian$root_n_draws, gaussian_expected, paste(target, "literal Gaussian multipliers complete rows"))
  near(gaussian$conditional_root_n_variance, result$conditional_root_n_variance,
    paste(target, "Gaussian and multinomial conditional variance identity"))
  near(gaussian$conf.int, result$conf.int, paste(target, "Gaussian and count normal intervals agree"))
  fails(wm_bootstrap(x, counts = counts, B = 3L), paste(target, "B mismatch rejected"))
  fails(wm_bootstrap(x, counts = counts, seed = 9L), paste(target, "frozen counts reject seed"))
  for (change in c("negative", "fractional", "nonfinite", "wrong_sum", "wrong_n",
                   "complex", "zero_columns", "wrong_row_names")) {
    bad <- counts
    if (change == "negative") bad[1L, 1L] <- -1
    if (change == "fractional") bad[1L, 1L] <- .5
    if (change == "nonfinite") bad[1L, 1L] <- Inf
    if (change == "wrong_sum") bad[1L, 1L] <- 2
    if (change == "wrong_n") bad <- bad[-1L, , drop = FALSE]
    if (change == "complex") bad <- bad + 1i
    if (change == "zero_columns") bad <- matrix(numeric(), n, 0L)
    if (change == "wrong_row_names") rownames(bad) <- as.character(n:1L)
    fails(wm_bootstrap(x, counts = bad), paste(target, change, "counts rejected"))
  }
  for (change in c("unavailable", "different_scope", "point", "n", "source",
                   "variance", "row_shift", "row_nonfinite", "row_length")) {
    bad <- x
    if (change == "unavailable") bad$available <- FALSE
    if (change == "different_scope") bad$covariance_scope <- "full_x"
    if (change == "point") bad$estimate <- bad$estimate + 1
    if (change == "n") bad$n <- bad$n - 1L
    if (change == "source") bad$source_object$estimate <- bad$estimate + 1
    if (change == "variance") bad$root_n_variance <- 2 * bad$root_n_variance
    if (change == "row_shift") bad$augmented_rows <- bad$augmented_rows + 1
    if (change == "row_nonfinite") bad$augmented_rows[1L] <- Inf
    if (change == "row_length") bad$augmented_rows <- bad$augmented_rows[-1L]
    fails(wm_bootstrap(bad, counts = counts), paste(target, change, "inference rejected"))
  }
  fails(wm_bootstrap(x, counts = counts, chunk_size = 0), paste(target, "zero chunk rejected"))
  fails(wm_bootstrap(x, counts = counts, conf.level = 1), paste(target, "invalid confidence level rejected"))
  fails(wm_bootstrap(x, counts = counts, B = NULL), paste(target, "explicit NULL B rejected"))
  fails(wm_bootstrap(x$fit, counts = counts), paste(target, "counts unsupported for legacy wm_match"))
  check(identical(serialize(x, NULL), before), paste(target, "all calls preserve original inference"))
  results[[target]] <- result
}

# Legacy wm_match generator: exact original uncentered rows, edges and training
# order, tested with literal multiplier values rather than a random simulation.
row <- c(1, 2, -1, 3)
edge <- c(-2, .5)
training <- c(-.4, .2, .8)
coefficients <- c(row, edge, training)
legacy <- structure(list(n = 4L, estimate = 1.25,
  root_n_variance = sum(coefficients^2) / 4,
  contributions = list(row = row, edge = edge, training = training)),
  class = c("wm_match", "list"))
values <- seq(-.9, 1.3, length.out = length(coefficients) * 2L)
legacy_out <- with_literal_normals(values,
  wm_bootstrap(legacy, B = 2L, interval = "none", chunk_size = 2L))
expected <- vapply(1:2, function(b) sum(coefficients *
  values[seq.int((b - 1L) * length(coefficients) + 1L, b * length(coefficients))]) / 2,
  numeric(1))
near(legacy_out$root_n_draws, expected, "legacy Gaussian row-edge-training ordering preserved")
near(legacy_out$conditional_root_n_variance, sum(coefficients^2) / 4,
  "legacy uncentered row-edge-training variance preserved")
check(legacy_out$method == "Gaussian contribution multipliers" &&
  is.null(legacy_out$source_inference), "legacy result shape preserved")
check(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE) == had_rng &&
  (!had_rng || identical(rng_before, get(".Random.seed", envir = .GlobalEnv))),
  "RNG unchanged, including prior absence")
check(length(tools::parse_Rd(file.path(root, "man", "wm_bootstrap.Rd"))) > 0,
  "bootstrap documentation parses")
receipt <- list(status = "PASS", checks = length(labels), labels = labels,
  elapsed_seconds = proc.time()[["elapsed"]] - start,
  propensity_refits = 0, outcome_model_refits = 0, matching_reruns = 0,
  new_datasets = 0, random_draws = 0, RNG_changed = FALSE,
  inference_source = if (length(args) >= 3L) "existing saved inference results" else
    "deterministic variance nuisances on existing portable fixture",
  scope = "Finite-contribution algebra, analytic variance and guards; no new simulation or calibration claim")
if (length(args) >= 2L) saveRDS(list(receipt = receipt, results = results), args[2L])
cat(length(labels), "checks PASS; zero refits, rematching, random draws or datasets.\n")
