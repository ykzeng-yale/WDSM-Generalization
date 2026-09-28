#!/usr/bin/env Rscript
# Reconstruct saved numerical geometry and repeat aggregate benchmark joins.
# No generated subjects, matching fits, numerical integration or RNG is used.

wm_results_geometry <- function(directory) {
  read <- function(name) utils::read.csv(file.path(directory, name), stringsAsFactors = FALSE)
  index <- read("geometry.csv")
  beta <- read("beta.csv")
  covariance <- read("beta_covariance.csv")
  if (anyDuplicated(index$key) || anyNA(index) ||
      !setequal(index$key, c("M1_d1", "M1_d2", "M1_d3", "M1_d5",
                            "M3_d1", "M3_d2", "M3_d3", "M3_d5")))
    stop("Unexpected saved geometry grid")
  result <- setNames(vector("list", nrow(index)), index$key)
  for (i in seq_len(nrow(index))) {
    row <- index[i, , drop = FALSE]
    M <- as.integer(row$M); d <- as.integer(row$d)
    if (!identical(row$key, sprintf("M%d_d%d", M, d))) stop("Geometry key differs")
    b <- beta[beta$key == row$key, , drop = FALSE]
    c <- covariance[covariance$key == row$key, , drop = FALSE]
    if (nrow(b) != M || anyDuplicated(b$overlap) ||
        !setequal(b$overlap, seq.int(0L, M - 1L)) || any(b$M != M | b$d != d) ||
        nrow(c) != M * M || anyDuplicated(c[c("overlap_row", "overlap_column")]) ||
        any(!c$overlap_row %in% b$overlap | !c$overlap_column %in% b$overlap))
      stop("Incomplete saved component/covariance grid")
    b <- b[order(b$overlap), , drop = FALSE]
    V <- matrix(NA_real_, M, M)
    V[cbind(c$overlap_row + 1L, c$overlap_column + 1L)] <- c$covariance
    close <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-11, check.attributes = FALSE))
    if (any(!is.finite(c(b$beta, b$mcse, V))) || any(b$beta <= 0 | b$mcse < 0) ||
        !close(V, t(V)) || min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) < -1e-12 ||
        !close(sqrt(pmax(0, diag(V))), b$mcse) || !close(sum(b$beta), row$alpha) ||
        !close(sqrt(max(0, sum(V))), row$alpha_mcse)) stop("Inconsistent saved geometry moments")
    if (d == 1L) {
      if (row$draws != 0 || any(V != 0) || any(b$precision != "exact_formula"))
        stop("Exact one-dimensional geometry has numerical uncertainty")
    } else if (row$draws <= 0 || anyNA(b$nonzero_draws) ||
               any(b$nonzero_draws <= 0 | b$nonzero_draws > row$draws) ||
               any(b$precision != "monte_carlo_estimate") || any(diag(V) <= 0))
      stop("Numerical geometry has unresolved precision")
    # Only fields consumed by the benchmark interfaces are reconstructed.
    result[[row$key]] <- structure(list(M = M, d = d, overlap = b$overlap,
      beta = b$beta, beta_mcse = b$mcse, covariance = V,
      alpha = row$alpha, alpha_mcse = row$alpha_mcse, method = row$method,
      precision_status = row$precision, component_precision = b$precision,
      draws = row$draws, requested_draws = row$requested_draws, seed = row$seed,
      diagnostics = list(nonzero_draws = if (d == 1L) NULL else b$nonzero_draws,
        below_jensen = row$alpha < M * M),
      representation = "Reconstructed saved CSV moments; not a new geometry computation"),
      class = c("wm_geometry", "list"))
  }
  result
}

wm_results_reproduce <- function(package_directory, output_directory = NULL) {
  root <- normalizePath(package_directory, mustWork = TRUE)
  saved <- file.path(root, "results")
  hashes <- utils::read.csv(file.path(saved, "artifact_hashes.csv"), stringsAsFactors = FALSE)
  if (!identical(names(hashes), c("path", "bytes", "sha256", "md5")) ||
      anyNA(hashes) || anyDuplicated(hashes$path) ||
      any(grepl("(^/|(^|/)\\.\\.(/|$))", hashes$path)) ||
      any(!grepl("^[0-9a-f]{32}$", hashes$md5))) stop("Invalid public artifact hash table")
  paths <- file.path(saved, hashes$path)
  if (any(!file.exists(paths)) || !identical(unname(tools::md5sum(paths)), hashes$md5) ||
      any(file.info(paths)$size != hashes$bytes)) stop("Public result artifact hash mismatch")
  environment <- new.env(parent = globalenv())
  for (file in c("benchmarks.R", "join_benchmarks.R", "first_stage.R",
                 "first_stage_benchmarks.R", "supplement.R", "supplement_join.R"))
    sys.source(file.path(root, "simulations", file), envir = environment)
  geometry <- wm_results_geometry(file.path(saved, "geometry"))
  primary <- utils::read.csv(file.path(saved, "primary", "summary.csv"), stringsAsFactors = FALSE)
  supplement <- utils::read.csv(file.path(saved, "supplement", "summary.csv"), stringsAsFactors = FALSE)
  primary_joined <- environment$wm_join_benchmarks(primary, geometry)
  supplement_joined <- environment$wm_supplement_join(supplement, geometry)
  check <- function(value, family, expected_rows) {
    reference <- utils::read.csv(file.path(saved, family, "summary_joined.csv"), stringsAsFactors = FALSE)
    if (nrow(reference) != expected_rows || !identical(names(value), names(reference)))
      stop("Saved aggregate schema differs: ", family)
    for (name in names(reference)) if (!isTRUE(all.equal(value[[name]], reference[[name]],
        tolerance = 2e-10, check.attributes = FALSE))) stop("Benchmark reproduction differs: ", family, "/", name)
  }
  check(primary_joined, "primary", 576L)
  check(supplement_joined, "supplement", 54L)
  if (!is.null(output_directory)) {
    if (file.exists(output_directory)) stop("Use a new output directory")
    dir.create(output_directory, recursive = TRUE)
    saveRDS(geometry, file.path(output_directory, "geometry.rds"))
    utils::write.csv(primary_joined, file.path(output_directory, "primary_joined.csv"), row.names = FALSE)
    utils::write.csv(supplement_joined, file.path(output_directory, "supplement_joined.csv"), row.names = FALSE)
  }
  cat("PASS: public artifact hashes and all 576 primary / 54 supplemental benchmark rows agree.\n")
  invisible(list(primary = primary_joined, supplement = supplement_joined, geometry = geometry))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("usage: Rscript results/reproduce_benchmarks.R [new_output_directory]")
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(script) != 1L) stop("Run this interface with Rscript")
  wm_results_reproduce(file.path(dirname(normalizePath(script)), ".."),
                       if (length(args)) args[1L] else NULL)
}
