# Differential validation against the frozen, complete pre-optimization core.
# Run with base R: Rscript code-release/validation/check_core_reference.R
# Optional bounded timing: append --benchmark (n=400,1600, dimensions 2,6).
# Algorithm: candidate donor rows are cached once per exact cell/arm. Every
# eligible distance is still computed; for cell sizes q_c,d_c this costs
# O(sum_c q_c*d_c*d) arithmetic. The old radix ordering was already linear
# for fixed-width numbers: partial selection improves ranking/allocation
# constants, not the unrestricted quadratic distance-computation class.
# The new cutoff scan orders only M selected rows and retains every exact
# boundary tie long enough to choose original row indices. The cache removes
# the old O(n*sum_c q_c) eligibility rescans, useful for many small cells.
# Working storage is O(max_c d_c*d) for distances, O(n) for cached row indices,
# and O(M*sum_c q_c) for directed edges; no n-by-n distance matrix is formed.
# Norm row maxima are columnwise vectorized, with the original scaling and
# accumulation order retained. Edge vectors are allocated once after the
# first query's donor and distance validation. No spatial pruning, tolerance,
# approximate match, metric rescaling, or compiled/library dependency is used.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this validation with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
reference <- new.env(parent = globalenv())
candidate <- new.env(parent = globalenv())
sys.source(file.path(root, "validation", "wm_core_reference.R"), envir = reference)
sys.source(file.path(root, "R", "wm_core.R"), envir = candidate)
checked <- 0L
check <- function(args) {
  old <- do.call(reference$wm_match, args)
  new <- do.call(candidate$wm_match, args)
  # Graph identity is stricter than numerical equality: includes ordering,
  # integer row types, original labels, distances, shares, and unordered flags.
  stopifnot(identical(new$graph, old$graph))
  difference <- all.equal(new, old, tolerance = 1e-13, check.attributes = TRUE)
  if (!isTRUE(difference)) stop(paste(difference, collapse = "; "))
  checked <<- checked + 1L
  invisible(new)
}
make_args <- function(n, d, M) {
  scores <- matrix(runif(n * d, -2, 2), n, d)
  list(Y = rnorm(n), Z = rep(0:1, n / 2), weights = exp(runif(n, -2, 2)),
       scores0 = scores, M = M, mean0 = sin(scores[, 1]),
       mean1 = 1 + scores[, 1]^2)
}
set.seed(928314)
for (d in c(1L, 2L, 3L, 6L)) for (M in c(1L, 3L)) {
  base <- make_args(96L, d, M)
  base$fold_id <- rep(c("a.b", "a"), each = 48L)
  base$strata <- rep(rep(c("c", "b.c"), each = 24L), 2)
  for (scores in list(base$scores0, round(base$scores0), round(base$scores0) * 2^660)) {
    for (estimand in c("PATE", "PATT")) for (method in c("self_normalized", "stabilized")) {
      a <- base
      a$scores0 <- scores
      a$estimand <- estimand
      a$method <- method
      if (estimand == "PATT") a$mean1 <- NULL
      if (method == "stabilized") {
        a$rho0 <- 1.2 + runif(96)
        if (estimand == "PATE") a$rho1 <- 1.5 + runif(96)
      }
      check(a)
      if (method == "self_normalized") {
        a$nuisance_influence <- matrix(rnorm(192), 96, 2)
        a$sensitivity <- c(.3, -.2)
        if (estimand == "PATE") a$scores1 <- cbind(scores, 0.25 * scores[, 1])
        check(a)
      }
    }
  }
}
extra <- make_args(12L, 3L, 6L)
check(extra) # Every eligible donor selected; no partial-selection branch.
extra$fold_id <- rep(c(1, 1 + .Machine$double.eps), each = 6L)
extra$M <- 3L
check(extra) # Adjacent doubles remain distinct restriction labels.
extra$mean0 <- extra$mean1 <- NULL
extra$variance <- FALSE
check(extra)

same_error <- function(args) {
  message_for <- function(fun) tryCatch({ do.call(fun, args); NULL }, error = conditionMessage)
  a <- message_for(reference$wm_match)
  b <- message_for(candidate$wm_match)
  stopifnot(is.character(a), identical(a, b))
}
valid <- make_args(12L, 2L, 3L)
bad <- valid; bad$M <- .Machine$integer.max; same_error(bad)
bad <- valid; bad$M <- 1 + 1i; same_error(bad)
bad <- valid; bad$weights[1] <- 0; same_error(bad)
bad <- valid; bad$mean1 <- NULL; same_error(bad)
bad <- valid; bad$mean0[2] <- Inf; same_error(bad)
bad <- valid; bad$fold_id <- seq_along(bad$Y); same_error(bad)
bad <- valid; bad$fold_id <- c(rep(1, 6), 2:7); same_error(bad)
bad <- valid; bad$scores0[1, 1] <- -1e308; bad$scores0[2, 1] <- 1e308; same_error(bad)
bad <- valid; bad$scores0[1, ] <- 0; bad$scores0[2, ] <- 1.4e308; same_error(bad)
cat("PASS:", checked, "complete reference comparisons and 9 identical error messages.\n")

if ("--benchmark" %in% commandArgs(trailingOnly = TRUE)) {
  # Graph/estimator equality is checked before each timing. Three repetitions
  # are bounded, and no result is asserted from noisy wall-clock measurements.
  timing <- list()
  for (n in c(400L, 1600L)) for (d in c(2L, 6L)) {
    a <- make_args(n, d, 3L)
    check(a)
    for (version in c("reference", "candidate")) {
      fun <- if (version == "reference") reference$wm_match else candidate$wm_match
      elapsed <- replicate(3L, system.time(do.call(fun, a))[["elapsed"]])
      timing[[length(timing) + 1L]] <- data.frame(n = n, d = d, M = 3L,
        version = version, median_seconds = median(elapsed),
        minimum_seconds = min(elapsed), maximum_seconds = max(elapsed))
    }
  }
  print(do.call(rbind, timing), row.names = FALSE)
}
print(sessionInfo())
