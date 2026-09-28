#!/usr/bin/env Rscript
# Sourceable summaries; CLI: output.csv input1.csv [input2.csv ...]
.wm_fss_scalar_same <- function(x) length(unique(x)) == 1L
.wm_fss_key <- function(x, fields) do.call(paste, c(x[fields], sep = "\r"))
.wm_fss_hash <- function(x, label) {
  if (!is.character(x) || !length(x) || anyNA(x) || is.null(names(x)) ||
      any(!grepl("^[[:xdigit:]]{32}$", x))) stop("invalid ", label, " hashes")
  names(x) <- basename(names(x))
  if (anyDuplicated(names(x))) stop("ambiguous ", label, " source names")
  x[order(names(x))]
}

# Optional validation-only allowance for declared cross-host recovery datasets.
# It never changes records or statistical formulas; default pairing is exact.
wm_first_stage_validate_records <- function(r, paired_parameter_roundoff_keys = NULL) {
  fields <- c("schema_version", "scenario", "replication", "branch", "n", "m", "d", "M",
    "base_seed", "estimand", "comparison", "pairing_id", "target", "status", "inference_status",
    "estimate", "root_n_variance", "sampling_variance", "lower", "upper",
    "evaluation_rng_state", "training_rng_state", "first_stage_status", "parameter1", "parameter2",
    "sigma2", "r", "evaluation_training_ratio", "warning", "error")
  if (!is.data.frame(r) || !nrow(r) || !all(fields %in% names(r))) stop("invalid first-stage records")
  keys <- c("scenario", "replication", "estimand", "comparison")
  if (anyNA(r[keys]) || anyDuplicated(r[keys]) || any(!nzchar(r$scenario))) stop("duplicate or missing record keys")
  for (name in c("schema_version", "replication", "n", "m", "d", "M", "base_seed"))
    if (!is.numeric(r[[name]]) || any(!is.finite(r[[name]])) || any(r[[name]] != floor(r[[name]])) ||
        any(r[[name]] < if (name == "base_seed") 0 else 1)) stop("invalid integer record field: ", name)
  if (any(r$schema_version != 1L) || any(!r$status %in% c("ok", "error")) ||
      any(!r$first_stage_status %in% c("ok", "error")) ||
      any(!r$branch %in% c("same_sample", "independent_training", "estimated_weights", "gaussian_same_sample")) ||
      any(!r$estimand %in% c("PATE", "PATT")) ||
      any(!r$comparison %in% c("known_first_stage", "fitted_naive", "fitted_adjusted")))
    stop("unknown record schema or category")
  label <- ifelse(r$comparison == "fitted_naive", "first_stage_omission_diagnostic", "within_declared_parametric_scope")
  if (anyNA(r$inference_status) || any(r$inference_status != label)) stop("inference label mismatch")
  if (anyNA(r$pairing_id) || any(r$pairing_id != paste(r$scenario, r$replication, r$estimand, sep = "/")))
    stop("pairing identity mismatch")
  roundoff_ids <- character()
  if (!is.null(paired_parameter_roundoff_keys)) {
    k <- paired_parameter_roundoff_keys
    if (!is.data.frame(k) || !identical(names(k), c("scenario", "replication")) ||
        anyNA(k) || anyDuplicated(k) || !is.character(k$scenario) ||
        !is.numeric(k$replication) || any(k$replication != floor(k$replication)))
      stop("invalid paired-parameter roundoff whitelist")
    roundoff_ids <- paste(k$scenario, k$replication, sep = "\r")
    existing <- paste(r$scenario, r$replication, sep = "\r")
    if (any(!roundoff_ids %in% existing) ||
        any(r$branch[existing %in% roundoff_ids] != "gaussian_same_sample"))
      stop("paired-parameter roundoff whitelist has an unknown or non-Gaussian dataset")
  }
  config <- c("branch", "n", "m", "d", "M", "base_seed", "target", "sigma2", "r", "evaluation_training_ratio")
  for (s in unique(r$scenario)) {
    x <- r[r$scenario == s, , drop = FALSE]
    if (!all(vapply(x[config], .wm_fss_scalar_same, logical(1)))) stop("inconsistent scenario configuration")
    if (any(!is.finite(x$target)) || x$n[1] < 2 || x$m[1] < 2 ||
        (x$branch[1] != "independent_training" && x$n[1] != x$m[1]) ||
        (x$branch[1] == "same_sample" && x$d[1] == 1L)) stop("unsupported scenario scope")
    weight_branch <- x$branch[1] == "estimated_weights"
    if (x$target[1] != (if (weight_branch) 1.5 else 1) || x$sigma2[1] != 1/16 ||
        (weight_branch && !is.na(x$r[1])) || (!weight_branch && !isTRUE(x$r[1] == .6)) ||
        (x$branch[1] == "independent_training" && !isTRUE(all.equal(x$evaluation_training_ratio[1], x$n[1]/x$m[1]))) ||
        (x$branch[1] != "independent_training" && !is.na(x$evaluation_training_ratio[1])))
      stop("record DGP constants disagree with declared design")
    expected <- expand.grid(replication = unique(x$replication), estimand = c("PATE", "PATT"),
      comparison = c("known_first_stage", "fitted_naive", "fitted_adjusted"), stringsAsFactors = FALSE)
    f <- c("replication", "estimand", "comparison")
    if (!setequal(.wm_fss_key(x, f), .wm_fss_key(expected, f))) stop("incomplete six-record replication grid")
    for (rep in unique(x$replication)) {
      z <- x[x$replication == rep, , drop = FALSE]
      if (!all(vapply(z[c("evaluation_rng_state", "training_rng_state", "first_stage_status", "parameter2")],
                      .wm_fss_scalar_same, logical(1)))) stop("paired first-stage or RNG states disagree")
      if (paste(s, rep, sep = "\r") %in% roundoff_ids) {
        values <- z$parameter1
        if (any(!is.finite(values)) || max(values)-min(values) >
            64*.Machine$double.eps*max(1,abs(values)))
          stop("paired parameter1 exceeds declared cross-host roundoff allowance")
      } else if (!.wm_fss_scalar_same(z$parameter1)) stop("paired first-stage or RNG states disagree")
      if (anyNA(z$evaluation_rng_state) || any(!nzchar(z$evaluation_rng_state))) stop("missing evaluation RNG state")
      if (anyNA(z$training_rng_state) ||
          (z$branch[1] == "independent_training" && any(!nzchar(z$training_rng_state))) ||
          (z$branch[1] != "independent_training" && any(nzchar(z$training_rng_state)))) stop("training RNG scope mismatch")
    }
  }
  numeric_fields <- c("estimate", "root_n_variance", "sampling_variance", "lower", "upper")
  ok <- r$status == "ok"
  if (any(!is.finite(as.matrix(r[ok, numeric_fields, drop = FALSE]))) ||
      any(r$root_n_variance[ok] < 0) || any(r$sampling_variance[ok] < 0) ||
      any(r$lower[ok] > r$estimate[ok]) || any(r$upper[ok] < r$estimate[ok])) stop("invalid successful numeric record")
  if (any(abs(r$root_n_variance[ok] - r$n[ok] * r$sampling_variance[ok]) >
      1e-10 * pmax(1, abs(r$root_n_variance[ok])))) stop("variance scaling mismatch")
  half <- stats::qnorm(.975) * sqrt(r$sampling_variance[ok])
  if (any(abs(r$lower[ok] - (r$estimate[ok] - half)) > 1e-10 * pmax(1, abs(r$estimate[ok]))) ||
      any(abs(r$upper[ok] - (r$estimate[ok] + half)) > 1e-10 * pmax(1, abs(r$estimate[ok]))))
    stop("normal interval disagrees with reported variance")
  if (any(!is.na(as.matrix(r[!ok, numeric_fields, drop = FALSE])))) stop("failed record contains numerical results")
  if (any(!ok & (is.na(r$error) | !nzchar(r$error)))) stop("failed record lacks an error")
  if (any(ok & r$comparison != "known_first_stage" & r$first_stage_status != "ok")) stop("fitted success after first-stage failure")
  for (id in unique(r$pairing_id)) {
    x <- r[r$pairing_id == id, , drop = FALSE]
    a <- x[x$comparison == "fitted_adjusted", ]; b <- x[x$comparison == "fitted_naive", ]
    if (a$status == "ok" && b$status != "ok") stop("adjusted fit succeeded without naive fit")
    if (a$status == "ok" && abs(a$estimate - b$estimate) > 1e-12 * max(1, abs(b$estimate)))
      stop("adjusted and naive point identity violated")
  }
  invisible(r)
}

wm_first_stage_read_batches <- function(inputs, require_complete = TRUE) {
  if (!is.character(inputs) || !length(inputs) || anyDuplicated(normalizePath(inputs, mustWork = TRUE)))
    stop("missing or duplicate batch inputs")
  records <- manifests <- vector("list", length(inputs))
  common_code <- common_design <- NULL
  scenario_configs <- list()
  for (i in seq_along(inputs)) {
    # Empty quoted fields are valid for same-sample training streams and for
    # successful error/warning text. Automatic type conversion can turn an
    # entirely empty column into logical NAs; preserve their character type.
    file <- inputs[i]; r <- utils::read.csv(file, stringsAsFactors = FALSE,
      colClasses = c(training_rng_state = "character", warning = "character", error = "character"))
    wm_first_stage_validate_records(r)
    m <- readRDS(paste0(file, ".metadata.rds"))
    if (!identical(m$schema_version, "first_stage_v1") || is.null(m$completed) ||
        !is.data.frame(m$config) || nrow(m$config) != 1L) stop("incomplete or invalid batch manifest")
    cfg <- m$config
    f <- c("scenario", "branch", "n", "m", "d", "M", "seed", "replications")
    if (!all(f %in% names(cfg)) || anyNA(cfg[f])) stop("invalid manifest configuration")
    for (name in c("n", "m", "d", "M", "seed", "replications"))
      if (!is.numeric(cfg[[name]]) || !is.finite(cfg[[name]]) || cfg[[name]] != floor(cfg[[name]]) ||
          cfg[[name]] < if (name == "seed") 0 else 1) stop("invalid manifest integer")
    bounds <- c(m$first, m$last)
    if (length(bounds) != 2L || any(!is.finite(bounds)) || any(bounds != floor(bounds)) ||
        bounds[1] < 1 || bounds[2] < bounds[1] || bounds[2] > cfg$replications) stop("invalid batch bounds")
    requested <- 6 * (m$last - m$first + 1)
    if (!identical(names(r), m$columns) || !isTRUE(m$requested_records == requested) ||
        !isTRUE(m$records_written == requested) || nrow(r) != requested ||
        !isTRUE(m$last_completed_replication == m$last) ||
        !setequal(unique(r$replication), seq.int(m$first, m$last))) stop("incomplete manifest record grid")
    for (name in setdiff(f, "replications")) {
      field <- if (name == "seed") "base_seed" else name
      if (anyNA(r[[field]]) || any(r[[field]] != cfg[[name]])) stop("record/manifest configuration mismatch")
    }
    code <- .wm_fss_hash(m$code_md5, "code"); design <- .wm_fss_hash(m$design_md5, "design")
    config_hash <- .wm_fss_hash(m$config_file_md5, "configuration")
    if (is.null(common_code)) { common_code <- code; common_design <- design }
    if (!identical(code, common_code) || !identical(design, common_design)) stop("source or design hashes differ across batches")
    previous <- scenario_configs[[cfg$scenario]]
    declaration <- list(config = cfg[f], configuration_hash = unname(config_hash))
    if (!is.null(previous) && !isTRUE(all.equal(previous, declaration, check.attributes = FALSE)))
      stop("scenario manifest declarations differ")
    scenario_configs[[cfg$scenario]] <- declaration
    records[[i]] <- r; manifests[[i]] <- m
  }
  r <- do.call(rbind, records); rownames(r) <- NULL
  wm_first_stage_validate_records(r)
  if (require_complete) for (s in names(scenario_configs))
    if (!setequal(unique(r$replication[r$scenario == s]), seq_len(scenario_configs[[s]]$config$replications)))
      stop("scenario lacks planned replications; summary refuses incomplete study")
  attr(r, "batch_provenance") <- list(inputs = normalizePath(inputs),
    input_md5 = tools::md5sum(inputs), metadata_md5 = tools::md5sum(paste0(inputs, ".metadata.rds")),
    declarations = scenario_configs, code_md5 = common_code, design_md5 = common_design,
    manifests = manifests, complete_scenarios_required = require_complete)
  r
}

wm_first_stage_summarize <- function(records, paired_parameter_roundoff_keys = NULL) {
  wm_first_stage_validate_records(records, paired_parameter_roundoff_keys)
  mcse <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
  avg <- function(x) if (length(x)) mean(x) else NA_real_
  groups <- unique(records[c("scenario", "estimand", "comparison")]); rows <- list()
  for (i in seq_len(nrow(groups))) {
    use <- .wm_fss_key(records, names(groups)) == .wm_fss_key(groups[i, , drop = FALSE], names(groups))
    r <- records[use, , drop = FALSE]; x <- r[r$status == "ok", , drop = FALSE]; R <- nrow(x)
    v <- if (R > 1) stats::var(x$estimate) else NA_real_
    sq <- if (R > 1) (x$estimate - mean(x$estimate))^2 * R / (R - 1) else numeric()
    v_se <- mcse(sq); mean_v <- avg(x$sampling_variance)
    ratio <- if (is.finite(v) && v > 0) mean_v / v else NA_real_
    ratio_if <- if (is.finite(ratio)) (x$sampling_variance - mean_v) / v - mean_v / v^2 * (sq - v) else numeric()
    hit <- x$lower <= x$target & x$upper >= x$target; coverage <- avg(hit)
    # Clopper-Pearson inversion remains informative with zero/all successes.
    cp <- if (R) stats::binom.test(sum(hit), R, conf.level = .95)$conf.int else c(NA_real_, NA_real_)
    pair <- records[records$scenario == r$scenario[1] & records$estimand == r$estimand[1], , drop = FALSE]
    a <- pair[pair$comparison == "fitted_adjusted", ]; b <- pair[pair$comparison == "fitted_naive", ]
    b <- b[match(a$replication, b$replication), ]; joint <- a$status == "ok" & b$status == "ok"
    delta <- a$root_n_variance[joint] - b$root_n_variance[joint]
    rows[[i]] <- cbind(groups[i, ], data.frame(branch = r$branch[1], n = r$n[1], m = r$m[1], d = r$d[1], M = r$M[1],
      base_seed = r$base_seed[1], target = r$target[1], inference_status = r$inference_status[1],
      requested = nrow(r), completed = R, failed = nrow(r) - R, failure_rate = mean(r$status != "ok"),
      independent_replications = nrow(r), records_per_replication = 6L,
      mean_estimate = avg(x$estimate), bias = avg(x$estimate - x$target), bias_mcse = mcse(x$estimate),
      empirical_variance = v, empirical_variance_mcse_asymptotic = v_se,
      empirical_root_n_variance = r$n[1] * v, empirical_root_n_variance_mcse_asymptotic = r$n[1] * v_se,
      mean_root_n_variance = avg(x$root_n_variance), mean_root_n_variance_mcse = mcse(x$root_n_variance),
      variance_ratio = ratio, variance_ratio_mcse_delta = mcse(ratio_if),
      intervals = R, coverage = coverage, coverage_mcse = if (R) sqrt(coverage * (1 - coverage) / R) else NA_real_,
      coverage_exact95_lower = cp[1], coverage_exact95_upper = cp[2],
      mean_interval_width = avg(x$upper - x$lower), mean_interval_width_mcse = mcse(x$upper - x$lower),
      paired_variance_requested = nrow(a), paired_variance_completed = sum(joint), paired_variance_failed = sum(!joint),
      adjusted_minus_naive_root_n_variance = avg(delta), adjusted_minus_naive_root_n_variance_mcse = mcse(delta),
      warning_records = sum(!is.na(r$warning) & nzchar(r$warning)), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  attr(out, "interpretation") <- paste("Six records share each replication; counts and MCSE use independent replication units.",
    "Numerical summaries condition on successful records. Variance and ratio MCSE are asymptotic/delta approximations.",
    "Coverage has both plug-in binomial MCSE and exact pointwise 95% Clopper-Pearson intervals.",
    "Paired adjusted-minus-naive summaries use joint successes and repeat across comparison rows; they are not independent.")
  out
}

.wm_first_stage_summarize_main <- function(args) {
  if (length(args) < 2L) stop("usage: output.csv input1.csv [input2.csv ...]")
  output <- args[1]; metadata <- paste0(output, ".metadata.rds")
  if (any(file.exists(c(output, metadata)))) stop("summary or companion metadata already exists")
  records <- wm_first_stage_read_batches(args[-1]); summary <- wm_first_stage_summarize(records)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary, output, row.names = FALSE, na = "NA")
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  saveRDS(list(created = Sys.time(), session = utils::sessionInfo(),
    summary_source_md5 = tools::md5sum(script), summary_md5 = tools::md5sum(output),
    provenance = attr(records, "batch_provenance"), interpretation = attr(summary, "interpretation")), metadata)
  cat("Summarized", nrow(records), "paired records in", nrow(summary), "comparison rows; failures retained.\n")
}
if (sys.nframe() == 0L) .wm_first_stage_summarize_main(commandArgs(trailingOnly = TRUE))
