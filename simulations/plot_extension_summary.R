#!/usr/bin/env Rscript
# Plot already-joined summaries only; no package, estimator or geometry runs.
# CLI: Rscript plot_extension_summary.R joined.csv new_directory [--no-png]
# Sourceable: wm_plot_extension_summary(input, output, source_file = this_file).

.wmep_numeric <- function(dat, names, nonnegative = character()) {
  for (name in names) {
    x <- dat[[name]]
    # read.csv guesses logical for an all-NA column.
    if (is.logical(x) && all(is.na(x))) x <- as.numeric(x)
    if (!is.numeric(x) || is.complex(x) || any(!is.na(x) & !is.finite(x)) ||
        (name %in% nonnegative && any(x < 0, na.rm = TRUE)))
      stop("Invalid numeric summary column: ", name, call. = FALSE)
    dat[[name]] <- x
  }
  dat
}

.wmep_prepare <- function(dat) {
  if (!is.data.frame(dat) || !nrow(dat) || anyDuplicated(names(dat)))
    stop("A nonempty joined summary with unique column names is required")
  first <- "comparison" %in% names(dat)
  supplement <- all(c("method", "correction") %in% names(dat))
  if (first == supplement) stop("Cannot uniquely identify first-stage or supplemental schema")
  common <- c("scenario", "branch", "estimand", "n", "d", "M", "requested",
    "completed", "failed", "bias", "bias_mcse", "coverage", "coverage_exact95_lower",
    "coverage_exact95_upper", "mean_root_n_variance", "mean_root_n_variance_mcse",
    "empirical_root_n_variance", "benchmark_root_n_variance", "benchmark_geometry_mcse",
    "geometry_status", "warning_records")
  extra <- if (first) c("comparison", "m", "inference_status",
    "empirical_root_n_variance_mcse_asymptotic", "reported_root_n_variance_limit",
    "benchmark_comparison_status", "paired_variance_requested", "paired_variance_completed",
    "paired_variance_failed", "adjusted_minus_naive_root_n_variance",
    "adjusted_minus_naive_root_n_variance_mcse", "adjusted_minus_naive_reported_limit") else
      c("method", "correction", "analysis_n", "fraction1", "clipped_records",
        "empirical_root_n_variance_mcse", "variance_ratio", "variance_ratio_mcse",
        "geometry_goal_met")
  if (!all(c(common, extra) %in% names(dat)))
    stop("Missing joined-summary columns: ", paste(setdiff(c(common, extra), names(dat)), collapse = ", "))
  text_fields <- c("scenario", "branch", "estimand", "geometry_status",
    if (first) c("comparison", "inference_status", "benchmark_comparison_status") else c("method", "correction"))
  for (name in text_fields) if (!is.character(dat[[name]]) || anyNA(dat[[name]]) || any(!nzchar(dat[[name]])))
    stop("Invalid text summary column: ", name)
  number_fields <- setdiff(c(common, extra), c(text_fields, "geometry_goal_met"))
  signed <- c("bias", "adjusted_minus_naive_root_n_variance", "adjusted_minus_naive_reported_limit")
  dat <- .wmep_numeric(dat, number_fields, setdiff(number_fields, signed))
  counts <- c("n", "d", "M", "requested", "completed", "failed", "warning_records",
    if (first) c("m", "paired_variance_requested", "paired_variance_completed", "paired_variance_failed") else
      c("analysis_n", "clipped_records"))
  for (name in counts) if (anyNA(dat[[name]]) || any(dat[[name]] != floor(dat[[name]])))
    stop("Invalid count column: ", name)
  if (any(dat$n < 2 | dat$d < 1 | dat$M < 1 | dat$requested < 1) ||
      any(dat$completed + dat$failed != dat$requested) ||
      any(dat$warning_records > dat$requested)) stop("Inconsistent counts")
  valid <- is.finite(dat$coverage)
  if (any(dat$completed > 0 & !valid) || any(dat$completed == 0 & valid) ||
      any(!is.finite(dat$coverage_exact95_lower[valid]) | !is.finite(dat$coverage_exact95_upper[valid])) ||
      any(dat$coverage_exact95_lower[valid] < 0 | dat$coverage_exact95_upper[valid] > 1 |
          dat$coverage_exact95_lower[valid] > dat$coverage[valid] + 1e-12 |
          dat$coverage_exact95_upper[valid] < dat$coverage[valid] - 1e-12))
    stop("Invalid exact coverage interval")
  for (name in c("bias", "mean_root_n_variance", "empirical_root_n_variance"))
    if (any(is.finite(dat[[name]][dat$completed == 0])))
      stop("All-failed groups must have missing estimator summaries")
  if (first) {
    if (any(!dat$branch %in% c("same_sample", "independent_training", "estimated_weights", "gaussian_same_sample")) ||
        any(!dat$estimand %in% c("PATE", "PATT")) ||
        any(!dat$comparison %in% c("known_first_stage", "fitted_naive", "fitted_adjusted")) ||
        any(dat$m < 2) || any(dat$paired_variance_completed + dat$paired_variance_failed != dat$paired_variance_requested) ||
        any(dat$paired_variance_requested != dat$requested)) stop("Unsupported first-stage grid or paired counts")
    dat$series <- match(dat$comparison, c("known_first_stage", "fitted_naive", "fitted_adjusted"))
    dat$analysis_n <- dat$n
    dat$training_ratio <- ifelse(dat$branch == "independent_training", dat$m/dat$n, 1)
    dat$clipped_records <- NA_real_
    dat$empirical_se <- dat$empirical_root_n_variance_mcse_asymptotic
    dat$geometry_resolved <- dat$benchmark_comparison_status == "precision_goal_met"
    dat$reported_limit <- dat$reported_root_n_variance_limit
    group_fields <- c("branch", "d", "M", "estimand", "training_ratio")
    # Paired summaries repeat across comparisons; show them once, not as independent curves.
    for (s in unique(dat$scenario)) for (e in unique(dat$estimand[dat$scenario == s])) {
      r <- dat[dat$scenario == s & dat$estimand == e, , drop = FALSE]
      if (nrow(r) != 3L || !setequal(r$series, 1:3)) stop("Incomplete first-stage comparison group")
      for (name in c("n", "m", "d", "M", "branch", "requested",
        "paired_variance_requested", "paired_variance_completed", "paired_variance_failed",
        "adjusted_minus_naive_root_n_variance", "adjusted_minus_naive_root_n_variance_mcse",
        "adjusted_minus_naive_reported_limit"))
        if (length(unique(r[[name]])) != 1L) stop("Inconsistent repeated first-stage field: ", name)
    }
  } else {
    if (any(!dat$branch %in% c("spline", "strata", "split")) ||
        any(!dat$estimand %in% c("PATE", "PATT", "potential0", "potential1")) ||
        any(!dat$method %in% c("self_normalized", "stabilized")) ||
        any(!dat$correction %in% c("oracle", "feasible")) ||
        any(dat$analysis_n < 2 | dat$analysis_n > dat$n) ||
        any(dat$clipped_records > dat$completed)) stop("Unsupported supplemental grid or analysis counts")
    is_split <- dat$branch == "split"
    if (any(!is.finite(dat$fraction1[is_split]) | dat$fraction1[is_split] <= 0 | dat$fraction1[is_split] >= 1) ||
        any(!is.na(dat$fraction1[!is_split]))) stop("Invalid split fractions")
    dat$series <- match(paste(dat$method, dat$correction),
      c("self_normalized oracle", "self_normalized feasible", "stabilized oracle", "stabilized feasible"))
    dat$empirical_se <- dat$empirical_root_n_variance_mcse
    if (!is.logical(dat$geometry_goal_met) || anyNA(dat$geometry_goal_met))
      stop("geometry_goal_met must be nonmissing logical values")
    dat$geometry_resolved <- dat$geometry_goal_met
    dat$reported_limit <- dat$benchmark_root_n_variance
    group_fields <- c("branch", "d", "M", "estimand", "fraction1")
    for (s in unique(dat$scenario)) for (e in unique(dat$estimand[dat$scenario == s])) {
      r <- dat[dat$scenario == s & dat$estimand == e, , drop = FALSE]
      if (r$branch[1] == "split") {
        if (nrow(r) != 1L || r$series != 3L) stop("Split comparisons must be stabilized oracle only")
      } else if (nrow(r) != 4L || !setequal(r$series, 1:4)) stop("Incomplete supplemental comparison group")
      for (name in c("n", "analysis_n", "d", "M", "branch", "requested", "fraction1"))
        if (length(unique(r[[name]])) != 1L) stop("Inconsistent repeated supplemental field: ", name)
    }
  }
  if (anyDuplicated(dat[c(group_fields, "n", "series")]))
    stop("Multiple scenarios occupy the same plotted group, sample size and series; do not pool them implicitly")
  panels <- unique(dat[group_fields])
  panels <- panels[do.call(order, c(unname(panels), list(na.last = TRUE))), , drop = FALSE]
  rownames(panels) <- NULL
  list(data = dat, panels = panels, fields = group_fields, first_stage = first,
    schema = if (first) "joined_first_stage" else "joined_supplement")
}

.wmep_rows <- function(dat, facet) {
  keep <- rep(TRUE, nrow(dat))
  for (field in names(facet)) {
    value <- facet[[field]][1]
    keep <- keep & if (is.na(value)) is.na(dat[[field]]) else !is.na(dat[[field]]) & dat[[field]] == value
  }
  dat[keep, , drop = FALSE]
}

.wmep_styles <- function(first) {
  if (first) list(labels = c("Known first stage", "Fitted naive", "Fitted adjusted"),
    short = c("Known", "Naive", "Adjusted"), colors = c("#0072B2", "#777777", "#D55E00"),
    pch = c(16, 17, 15), lty = c(1, 2, 1)) else
      list(labels = c("Original oracle", "Original feasible", "Stabilized oracle", "Stabilized feasible"),
        short = c("Orig/oracle", "Orig/fit", "Stab/oracle", "Stab/fit"),
        colors = c("#0072B2", "#0072B2", "#D55E00", "#D55E00"),
        pch = c(16, 17, 16, 17), lty = c(1, 2, 1, 2))
}

.wmep_panel <- function(rows, styles, variable, se = NULL, title, ylab,
                        reference = NULL, coverage = FALSE, benchmark = NULL,
                        benchmark_se = "benchmark_geometry_mcse") {
  sizes <- sort(unique(rows$n)); series <- sort(unique(rows$series))
  offset <- if (length(series) == 1L) 0 else seq(-.19, .19, length.out = length(series))
  names(offset) <- series
  y <- rows[[variable]]
  error <- if (is.null(se)) rep(NA_real_, length(y)) else rows[[se]]
  lo <- if (coverage) rows$coverage_exact95_lower else y - 1.96 * error
  hi <- if (coverage) rows$coverage_exact95_upper else y + 1.96 * error
  by <- if (is.null(benchmark)) numeric() else rows[[benchmark]]
  bse <- if (is.null(benchmark)) numeric() else rows[[benchmark_se]]
  values <- c(y, lo, hi, by, by - 1.96*bse, by + 1.96*bse, reference)
  values <- values[is.finite(values)]
  ylim <- if (length(values)) range(values) else c(0, 1)
  if (coverage) ylim <- c(0, 1.28) else {
    span <- diff(ylim)
    if (span == 0) span <- max(.05, abs(ylim[1])*.2)
    ylim <- ylim + c(-.12, .34)*span
  }
  graphics::plot(NA_real_, NA_real_, xlim = c(.55, length(sizes)+.45), ylim = ylim,
    xaxt = "n", yaxt = if (coverage) "n" else "s", xlab = "Total evaluation sample size n",
    ylab = ylab, main = title, cex.main = .96, cex.lab = .81, cex.axis = .77)
  graphics::axis(1, at = seq_along(sizes), labels = format(sizes, scientific = FALSE, trim = TRUE), cex.axis = .77)
  if (coverage) graphics::axis(2, at = seq(0, 1, .25), las = 1, cex.axis = .77)
  graphics::abline(h = if (coverage) seq(0, 1, .25) else graphics::axTicks(2), col = "#EEEEEE", lwd = .6)
  if (length(reference)) graphics::abline(h = reference, col = "#777777", lty = 3, lwd = 1)
  for (s in series) {
    idx <- which(rows$series == s); idx <- idx[order(rows$n[idx])]
    x <- match(rows$n[idx], sizes) + offset[as.character(s)]
    if (!is.null(benchmark)) x <- x - .025
    good <- is.finite(y[idx]); bars <- good & is.finite(lo[idx]) & is.finite(hi[idx])
    graphics::segments(x[bars], lo[idx][bars], x[bars], hi[idx][bars],
      col = grDevices::adjustcolor(styles$colors[s], .65), lwd = 1.25)
    # NA values break lines; failed sample-size groups are never bridged.
    graphics::lines(x, y[idx], type = "b", col = styles$colors[s],
      pch = styles$pch[s], lty = styles$lty[s], cex = .75, lwd = 1.15)
    if (!is.null(benchmark)) {
      bx <- x + .05; v <- rows[[benchmark]][idx]; e <- rows[[benchmark_se]][idx]
      resolved <- rows$geometry_resolved[idx]
      bbar <- is.finite(v) & is.finite(e)
      graphics::segments(bx[bbar], v[bbar]-1.96*e[bbar], bx[bbar], v[bbar]+1.96*e[bbar],
        col = styles$colors[s], lwd = 2.8)
      graphics::lines(bx, v, col = styles$colors[s], lty = 3, lwd = .9)
      graphics::points(bx, v, pch = ifelse(resolved, 5, 4), col = styles$colors[s], cex = .8)
    }
  }
  graphics::legend("topright", legend = styles$labels[series], col = styles$colors[series],
    pch = styles$pch[series], lty = styles$lty[series], ncol = if (length(series)>2L) 2L else 1L,
    cex = .66, bty = "n", inset = .005)
  if (!any(is.finite(y))) graphics::mtext(if (all(rows$completed == 0))
    "All fits failed; no estimator summaries" else "No finite statistic for this panel",
    side = 3, line = .1, adj = 0, col = "#8B0000", cex = .68)
}

.wmep_pair_panel <- function(rows) {
  r <- rows[rows$comparison == "fitted_adjusted", , drop = FALSE]
  # Geometry cancels in this difference. Do not attach the baseline's geometry SE.
  r$series <- 1L; r$pair_limit <- r$adjusted_minus_naive_reported_limit
  r$zero_geometry_se <- 0; r$geometry_resolved <- TRUE
  .wmep_panel(r, list(labels = "Paired adjusted - naive", colors = "#009E73", pch = 16, lty = 1),
    "adjusted_minus_naive_root_n_variance", "adjusted_minus_naive_root_n_variance_mcse",
    "Change in reported variance on paired data", "Root-N variance change", reference = 0,
    benchmark = "pair_limit", benchmark_se = "zero_geometry_se")
  if (all(r$paired_variance_completed == 0)) graphics::mtext("No jointly successful pairs; limit only",
    side = 1, line = 2.7, cex = .65, col = "#8B0000")
}

.wmep_counts <- function(rows, styles, first) {
  graphics::plot.new(); graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  graphics::title(main = "Replications, analysis sizes and diagnostics", cex.main = .96)
  lines <- c("Counts are per comparison; paired records are not independent datasets.")
  for (n in sort(unique(rows$n))) {
    r <- rows[rows$n == n, , drop = FALSE]; r <- r[order(r$series), , drop = FALSE]
    sizes <- if (first && r$branch[1] == "independent_training")
      sprintf("evaluation n=%d; independent training m=%d", n, r$m[1]) else if (first)
      sprintf("n=%d; first stage fitted on these same observations", n) else
      sprintf("total n=%d; analysis N=%d", n, r$analysis_n[1])
    lines <- c(lines, sizes,
      paste(sprintf("%s: %d/%d successful", styles$short[r$series], r$completed, r$requested), collapse = "; "))
    if (first) {
      p <- r[r$comparison == "fitted_adjusted", ]
      lines <- c(lines, sprintf("Paired successes %d/%d; warnings across comparisons %d.",
        p$paired_variance_completed, p$paired_variance_requested, sum(r$warning_records)))
    } else lines <- c(lines, paste0("Clipped records: ",
      paste(sprintf("%s %d", styles$short[r$series], r$clipped_records), collapse = "; "),
      sprintf(". Warnings %d.", sum(r$warning_records))))
  }
  status <- sort(unique(rows$geometry_status))
  lines <- c(lines, paste("Geometry status:", paste(status, collapse = ", ")),
    sprintf("Unresolved/not-at-goal benchmark rows: %d/%d.", sum(!rows$geometry_resolved), nrow(rows)),
    "Summaries condition on success; failed groups stay visible.",
    if (first) "Naive fitted intervals omit first-stage uncertainty by design." else
      "Clipping counts refer to fit records, not distinct individuals.")
  if (!first && rows$branch[1] == "split") lines <- c(lines,
    sprintf("Nominal block-one fraction = %.8g; score dimensions d0=1, d1=2.", rows$fraction1[1]),
    "Potential-mean pages use their own N; the contrast uses total n.")
  wrapped <- unlist(lapply(lines, strwrap, width = 77), use.names = FALSE)
  size <- if (length(wrapped) > 18L) .61 else .69
  y <- seq(.97, .02, length.out = max(2L, length(wrapped)))
  graphics::text(0, y[seq_along(wrapped)], wrapped, adj = c(0, 1), cex = size)
}

.wmep_draw_page <- function(rows, facet, first) {
  styles <- .wmep_styles(first)
  graphics::par(mfrow = c(3, 2), mar = c(4, 4.3, 2.8, 1), oma = c(4.4, .2, 3.4, .2),
    family = "sans", las = 1, mgp = c(2.6, .65, 0), tcl = -.25)
  .wmep_panel(rows, styles, "bias", "bias_mcse", "Bias relative to the target", "Bias", reference = 0)
  .wmep_panel(rows, styles, "coverage", title = "Coverage of nominal 95% intervals",
    ylab = "Coverage", reference = .95, coverage = TRUE)
  .wmep_panel(rows, styles, "mean_root_n_variance", "mean_root_n_variance_mcse",
    "Mean reported variance and its own limit", "Root-N variance", benchmark = "reported_limit")
  .wmep_panel(rows, styles, "empirical_root_n_variance", "empirical_se",
    "Empirical sampling variance and benchmark", "Root-N variance", benchmark = "benchmark_root_n_variance")
  if (first) .wmep_pair_panel(rows) else
    .wmep_panel(rows, styles, "variance_ratio", "variance_ratio_mcse",
      "Mean reported / empirical sampling variance", "Variance ratio", reference = 1)
  .wmep_counts(rows, styles, first)
  branch <- c(same_sample = "Same-sample fitted scores", independent_training = "Independent score training",
    estimated_weights = "Fitted individual weights", gaussian_same_sample = "Qualified Gaussian scalar fitting",
    spline = "Feasible nonlinear splines", strata = "Finite random strata", split = "Independent potential-mean blocks")
  dimension_label <- if (!first && facet$branch == "split") "d0=1, d1=2" else sprintf("d=%d", facet$d)
  title <- sprintf("%s | %s | %s | M=%d", branch[facet$branch], facet$estimand, dimension_label, facet$M)
  if (!first && facet$branch == "split") title <- paste0(title, sprintf(" | block-one fraction %.6g", facet$fraction1))
  if (first && facet$branch == "independent_training") title <- paste0(title, sprintf(" | m/n=%.6g", facet$training_ratio))
  graphics::mtext(title, outer = TRUE, side = 3, line = 1.3, cex = 1.07, font = 2)
  graphics::mtext("Joined synthetic summaries only; asymptotic limits are not exact finite-sample predictions.",
    outer = TRUE, side = 3, line = .15, cex = .74)
  captions <- c(
    "Bias/variance/paired bars: +/- 1.96 Monte Carlo SE; empirical-variance and ratio SEs are asymptotic/delta approximations.",
    "Coverage bars: supplied exact pointwise 95% binomial intervals. No simultaneous or significance claim is made.",
    "Variance panels: filled marks = simulation summaries; offset diamonds = limits, with separate +/- 1.96 geometry MCSE bars; x = unresolved precision.",
    "N is analysis_n (n for first-stage and effect pages); paired-change limits have no geometry error because alpha cancels.")
  if (!first) captions[4] <- "N is analysis_n; split potential-mean pages use their block size. Missing MCSE means no bar, not zero uncertainty."
  for (i in seq_along(captions)) graphics::mtext(captions[i], outer = TRUE, side = 1,
    line = .15 + .9*(i-1), cex = .66)
}

wm_plot_extension_summary <- function(input, output, png = TRUE, source_file = NULL) {
  if (!is.logical(png) || length(png) != 1L || is.na(png)) stop("png must be TRUE or FALSE")
  if (is.null(source_file) || length(source_file) != 1L || !file.exists(source_file))
    stop("Supply this plotting script as source_file for provenance")
  input <- normalizePath(input, mustWork = TRUE)
  source_file <- normalizePath(source_file, mustWork = TRUE)
  if (file.exists(output) || dir.exists(output)) stop("Use a new output directory; existing outputs are never overwritten")
  companion <- paste0(input, ".metadata.rds")
  inputs <- c(input, source_file, if (file.exists(companion)) companion)
  before <- tools::md5sum(inputs)
  dat <- utils::read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
  prepared <- .wmep_prepare(dat); rows <- prepared$data; panels <- prepared$panels
  if (nrow(panels) > 200L) stop("More than 200 pages exceeds the plotting work guard")
  if (!dir.create(output, recursive = TRUE, showWarnings = FALSE)) stop("Cannot create output directory")
  output <- normalizePath(output, mustWork = TRUE)
  metadata_file <- file.path(output, "figure_metadata.rds")
  metadata <- list(status = "started", started = Sys.time(), input = input,
    source_file = source_file, hashes_before = before, hash_algorithm = "MD5 via tools::md5sum",
    input_companion_metadata = if (file.exists(companion)) companion else NULL,
    session = utils::sessionInfo(), schema = prepared$schema, panels = panels,
    summary_rows = dat, panel_rows = vector("list", nrow(panels)),
    panel_diagnostics = vector("list", nrow(panels)),
    pdf_pages_completed = 0L, png_pages_completed = 0L,
    scope = paste("Plots already-joined summaries only; no estimator or geometry computation.",
      "Statistics condition on successful comparisons; failures and clipping are retained.",
      "Approximate MC error and numerical geometry error are separate; exact coverage intervals are pointwise.",
      "Figures do not certify asymptotic theorems or claim significance."))
  saveRDS(metadata, metadata_file)
  tryCatch({
    pdf_file <- file.path(output, "extension_summary.pdf")
    grDevices::pdf(pdf_file, width = 11.5, height = 10.5, onefile = TRUE, useDingbats = FALSE)
    tryCatch(for (i in seq_len(nrow(panels))) {
      r <- .wmep_rows(rows, panels[i, , drop = FALSE])
      .wmep_draw_page(r, panels[i, , drop = FALSE], prepared$first_stage)
      metadata$panel_rows[[i]] <- r
      metadata$panel_diagnostics[[i]] <- list(all_failed = all(r$completed == 0),
        comparison_rows = nrow(r), requested_comparisons = sum(r$requested),
        failed_comparisons = sum(r$failed),
        all_failed_comparison_groups = which(r$completed == 0),
        analysis_sizes = unique(r[c("n", "analysis_n")]),
        clipped_records = if (prepared$first_stage) NULL else sum(r$clipped_records),
        counts_scope = "Comparison counts can share underlying simulated datasets")
      metadata$pdf_pages_completed <- i
    }, finally = grDevices::dev.off())
    if (png) for (i in seq_len(nrow(panels))) {
      facet <- panels[i, , drop = FALSE]
      name <- sprintf("page%03d_%s_%s_d%d_M%d.png", i, facet$branch, facet$estimand, facet$d, facet$M)
      # Sequential page IDs keep different split fractions distinct in filenames.
      grDevices::png(file.path(output, name), width = 11.5, height = 10.5, units = "in", res = 240)
      tryCatch(.wmep_draw_page(metadata$panel_rows[[i]], facet, prepared$first_stage),
        finally = grDevices::dev.off())
      metadata$png_pages_completed <- i
    }
    if (!identical(before, tools::md5sum(inputs))) stop("Input, companion provenance or plotting code changed while rendering")
    figure_files <- list.files(output, "\\.(pdf|png)$", full.names = TRUE)
    metadata$figure_md5 <- tools::md5sum(figure_files)
    metadata$status <- "completed"; metadata$completed <- Sys.time()
    saveRDS(metadata, metadata_file)
    cat("Created", nrow(panels), "PDF pages and", metadata$png_pages_completed,
        "PNG pages from", nrow(dat), "joined rows; no estimator runs.\n")
    invisible(metadata)
  }, error = function(e) {
    metadata$status <- "failed"; metadata$error <- conditionMessage(e); metadata$stopped <- Sys.time()
    files <- list.files(output, "\\.(pdf|png)$", full.names = TRUE)
    metadata$partial_figure_md5 <- tools::md5sum(files)
    saveRDS(metadata, metadata_file)
    stop(conditionMessage(e), call. = FALSE)
  })
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) %in% 2:3 || (length(args) == 3L && args[3] != "--no-png"))
    stop("usage: joined_summary.csv new_output_directory [--no-png]")
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(script) != 1L) stop("Run the CLI with Rscript")
  wm_plot_extension_summary(args[1], args[2], png = length(args) == 2L, source_file = script)
}
