#!/usr/bin/env Rscript
# Fixed independent production grid for the validated boundary-rank integral.
# Rscript geometry_boundary_constants.R config.csv new_output_directory [source_package]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) stop("usage: config.csv new_output_directory [source_package]")
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
root <- if (length(args) == 3L) normalizePath(args[3], mustWork = TRUE) else normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
paths <- file.path(root, c("R/wm_geometry.R", "R/wm_geometry_overlap.R", "R/wm_geometry_alpha.R", "simulations/benchmarks.R"))
for (path in paths) source(path)
config <- utils::read.csv(args[1], stringsAsFactors = FALSE)
if (!identical(names(config), c("M", "d", "draws", "seed")) || nrow(config) != 8L ||
    anyNA(config) || any(!is.finite(as.matrix(config))) || any(as.matrix(config) != floor(as.matrix(config))) ||
    any(!config$M %in% c(1L, 3L)) || any(!config$d %in% c(1L, 2L, 3L, 5L)) ||
    any(config$draws < 1 | config$draws > 600000) || any(config$seed < 1) ||
    anyDuplicated(config[c("M", "d")]) || anyDuplicated(config$seed)) stop("Invalid fixed eight-case bounded grid")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
output <- args[2]
if (file.exists(output)) stop("Use a new output directory")
if (!dir.create(output, recursive = TRUE)) stop("Cannot create output directory")
file.copy(args[1], file.path(output, "config.csv"), overwrite = FALSE)
metadata <- list(started = Sys.time(), config = config, code_md5 = tools::md5sum(c(script, paths)),
  config_md5 = tools::md5sum(args[1]), rng_kind = RNGkind(), session = utils::sessionInfo(),
  immutable_sha256_manifest = Sys.getenv("WDSM_RUN_MANIFEST", unset = NA_character_),
  design = paste("Seeds and draw counts were fixed before these integrations. Independent of estimator and numerical-pilot streams.",
    "600k M3 draws selected from pilot MCSE to target at most1% propagated benchmark MCSE; this is not guaranteed.",
    "Every result remains unprojected; unmet precision is reported, not masked."))
saveRDS(metadata, file.path(output, "metadata.rds"))
objects <- list(); precision <- list(); beta <- list(); alpha <- list()
for (i in seq_len(nrow(config))) {
  cfg <- config[i, ]; key <- sprintf("M%d_d%d", cfg$M, cfg$d)
  result <- wm_geometry_overlap(cfg$M, cfg$d, cfg$draws, cfg$seed, chunk_size = 4096L)
  objects[[key]] <- result
  saveRDS(result, file.path(output, paste0(key, ".rds")))
  beta[[i]] <- data.frame(key = key, M = cfg$M, d = cfg$d, overlap = result$overlap,
    beta = unname(result$beta), mcse = unname(result$beta_mcse), precision = result$component_precision)
  alpha[[i]] <- data.frame(key = key, M = cfg$M, d = cfg$d, alpha = result$alpha,
    mcse = result$alpha_mcse, below_jensen = result$diagnostics$below_jensen,
    method = result$method, precision = result$precision_status)
  for (design in c("strong", "weak")) {
    covariance <- if (cfg$d == 1L) matrix(0, cfg$M, cfg$M) else result$covariance
    b <- if (design == "strong") wm_benchmark_strong(cfg$M, result$beta, covariance) else
      wm_benchmark_weak(cfg$M, result$beta, covariance)
    z <- b$table
    z$relative_geometry_mcse <- z$geometry_mcse / z$root_n_variance
    z$goal_met <- is.finite(z$relative_geometry_mcse) & z$relative_geometry_mcse <= .01 &
      !any(grepl("^unresolved", result$component_precision))
    z$goal_met[!is.finite(z$root_n_variance)] <- NA
    precision[[length(precision) + 1L]] <- data.frame(key = key, M = cfg$M, d = cfg$d, design = design, z)
  }
  saveRDS(objects, file.path(output, "geometry.rds"))
  cat(key, "alpha", result$alpha, "MCSE", result$alpha_mcse, "\n")
}
utils::write.csv(do.call(rbind, beta), file.path(output, "beta.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, alpha), file.path(output, "alpha.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, precision), file.path(output, "benchmark_precision.csv"), row.names = FALSE)
metadata$completed <- Sys.time(); metadata$cases <- length(objects)
metadata$code_unchanged <- identical(metadata$code_md5, tools::md5sum(c(script, paths)))
saveRDS(metadata, file.path(output, "metadata.rds"))
stopifnot(metadata$code_unchanged, metadata$cases == 8L)
