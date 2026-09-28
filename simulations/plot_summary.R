#!/usr/bin/env Rscript
# Rscript simulations/plot_summary.R joined_summary.csv new_output_directory
# Base-R PDF plus 300-dpi PNG figures from aggregated synthetic results only.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: joined_summary.csv new_output_directory")
input <- normalizePath(args[1L], mustWork = TRUE)
output <- args[2L]
if (file.exists(output) || dir.exists(output)) stop("Use a new output directory")
dat <- utils::read.csv(input, stringsAsFactors = FALSE)
required <- c("design", "estimand", "d", "M", "n", "method", "correction",
  "inference_status", "requested", "completed", "failed", "variance_ratio",
  "variance_ratio_mcse_delta", "coverage", "coverage_mcse", "bias", "bias_mcse",
  "benchmark_target_bias", "benchmark_scope")
if (!nrow(dat) || !all(required %in% names(dat))) stop("A joined simulation summary is required")
if (any(!dat$method %in% c("self_normalized", "stabilized")) ||
    any(!dat$correction %in% c("oracle", "polynomial", "raw")) ||
    any(!dat$design %in% c("strong", "weak")) ||
    any(!is.finite(dat$n) | dat$n <= 0)) stop("Unsupported plotting configuration")
panels <- unique(dat[c("design", "estimand", "d", "M")])
panels <- panels[order(panels$design, panels$estimand, panels$d, panels$M), , drop = FALSE]
colors <- c(self_normalized = "#0072B2", stabilized = "#D55E00")
symbols <- c(oracle = 16, polynomial = 17, raw = 1)
line_types <- c(oracle = 1, polynomial = 2, raw = 3)
method_labels <- c(self_normalized = "Original", stabilized = "Stabilized")
correction_labels <- c(oracle = "oracle", polynomial = "polynomial", raw = "raw point only")

draw_panel <- function(rows, variable, se_variable, title, ylab, reference = NULL,
                       inference_only = FALSE, show_xlab = FALSE) {
  if (inference_only) rows <- rows[rows$correction != "raw", , drop = FALSE]
  valid <- is.finite(rows[[variable]])
  x <- sort(unique(rows$n))
  if (!length(x)) x <- 1
  xlim <- if (length(x) > 1L) range(x) * c(.93, 1.07) else x * c(.75, 1.33)
  y <- rows[[variable]]
  se <- rows[[se_variable]]
  radius <- ifelse(is.finite(se), 1.96 * se, 0)
  limits <- c(y[valid] - radius[valid], y[valid] + radius[valid], reference)
  if (identical(variable, "bias")) limits <- c(limits, rows$benchmark_target_bias)
  limits <- limits[is.finite(limits)]
  ylim <- if (length(limits)) range(limits) else c(0, 1)
  if (diff(ylim) <= 0) ylim <- ylim + c(-1, 1) * max(.05, abs(ylim[1L]) * .1)
  else ylim <- ylim + c(-1, 1) * .09 * diff(ylim)
  if (variable == "coverage") ylim <- pmin(1, pmax(0, ylim))
  plot(NA_real_, NA_real_, xlim = xlim, ylim = ylim, log = "x", xaxt = "n",
       xlab = if (show_xlab) "Evaluation sample size" else "", ylab = ylab, main = title,
       cex.main = .95, cex.lab = .9, cex.axis = .8)
  axis(1, at = x, labels = format(x, scientific = FALSE, trim = TRUE), cex.axis = .8)
  abline(h = axTicks(2), col = "#EEEEEE", lwd = .6)
  if (!is.null(reference)) abline(h = reference, col = "#777777", lty = 3, lwd = .9)
  if (identical(variable, "bias") && any(rows$design == "weak")) {
    limit <- unique(rows$benchmark_target_bias[rows$method == "self_normalized"])
    if (length(limit) == 1L && is.finite(limit)) {
      abline(h = limit, col = colors["self_normalized"], lty = 4, lwd = 1.2)
      legend("bottomleft", legend = "Original-rule probability-limit bias",
             col = colors["self_normalized"], lty = 4, lwd = 1.2, cex = .60, bty = "n")
    }
  }
  series <- unique(rows[c("method", "correction")])
  labels <- character(); legend_colors <- character(); legend_pch <- integer(); legend_lty <- integer()
  for (j in seq_len(nrow(series))) {
    method <- series$method[j]; correction <- series$correction[j]
    keep <- rows$method == method & rows$correction == correction & valid
    r <- rows[keep, , drop = FALSE]
    if (!nrow(r)) next
    r <- r[order(r$n), , drop = FALSE]
    value <- r[[variable]]; error <- r[[se_variable]]
    lo <- value - 1.96 * error; hi <- value + 1.96 * error
    if (variable == "coverage") {lo <- pmax(0, lo); hi <- pmin(1, hi)}
    has_se <- is.finite(lo) & is.finite(hi)
    segments(r$n[has_se], lo[has_se], r$n[has_se], hi[has_se],
             col = grDevices::adjustcolor(colors[method], .65), lwd = 1.1)
    lines(r$n, value, type = "b", col = colors[method], pch = symbols[correction],
          lty = line_types[correction], lwd = 1.1, cex = .75)
    labels <- c(labels, paste(method_labels[method], correction_labels[correction]))
    legend_colors <- c(legend_colors, colors[method])
    legend_pch <- c(legend_pch, symbols[correction]); legend_lty <- c(legend_lty, line_types[correction])
  }
  if (length(labels)) legend("topright", legend = labels, col = legend_colors,
    pch = legend_pch, lty = legend_lty, lwd = 1, cex = .58, bty = "n", ncol = 2)
  else text(sqrt(prod(xlim)), mean(ylim), "No successful results for this panel", cex = .8)
}

draw_page <- function(facet) {
  take <- rep(TRUE, nrow(dat))
  for (field in names(facet)) take <- take & dat[[field]] == facet[[field]]
  rows <- dat[take, , drop = FALSE]
  par(mfrow = c(3, 1), mar = c(3.8, 4.3, 2.8, 1), oma = c(2.7, .2, 2.5, .2),
      family = "sans", las = 1, mgp = c(2.6, .65, 0), tcl = -.25)
  draw_panel(rows, "variance_ratio", "variance_ratio_mcse_delta",
    "Reported sampling variance relative to empirical variance", "Variance ratio", 1, TRUE)
  draw_panel(rows, "coverage", "coverage_mcse", "Coverage of nominal 95% intervals",
    "Coverage", .95, TRUE)
  draw_panel(rows, "bias", "bias_mcse", "Bias relative to the true causal target",
    "Bias", 0, FALSE, TRUE)
  mtext(sprintf("%s centering | %s | d = %d | M = %d",
    if (facet$design == "strong") "Strong" else "Weighted", facet$estimand,
    facet$d, facet$M), outer = TRUE, side = 3, line = .8, cex = 1.05, font = 2)
  total <- sum(rows$requested); failed <- sum(rows$failed)
  scope <- if (facet$design == "weak")
    "Original-rule intervals/variance ratios are assumption-failure diagnostics. " else ""
  mtext(paste0(scope, "Raw rules have no interval or variance claim."),
        outer = TRUE, side = 1, line = .35, cex = .61)
  mtext(sprintf("Bars: estimate +/- 1.96 Monte Carlo SE (approximate); summaries condition on success. Failed records: %d/%d.",
    failed, total), outer = TRUE, side = 1, line = 1.35, cex = .6)
}

if (!dir.create(output, recursive = TRUE, showWarnings = FALSE)) stop("Cannot create output directory")
pdf_file <- file.path(output, "simulation_summary.pdf")
grDevices::pdf(pdf_file, width = 7.2, height = 9.6, onefile = TRUE, useDingbats = FALSE)
tryCatch(for (i in seq_len(nrow(panels))) draw_page(panels[i, , drop = FALSE]),
         finally = grDevices::dev.off())
for (i in seq_len(nrow(panels))) {
  facet <- panels[i, , drop = FALSE]
  file <- file.path(output, sprintf("%s_%s_d%d_M%d.png", facet$design,
    facet$estimand, facet$d, facet$M))
  grDevices::png(file, width = 7.2, height = 9.6, units = "in", res = 300)
  tryCatch(draw_page(facet), finally = grDevices::dev.off())
}
saveRDS(list(input = input, input_md5 = tools::md5sum(input), created = Sys.time(),
  session = utils::sessionInfo(), panels = panels,
  interpretation = "Approximate Monte Carlo error bars on aggregated synthetic results; finite-sample scope labels retained"),
  file.path(output, "figure_metadata.rds"))
cat("Created one multipage PDF and", nrow(panels), "300-dpi PNG figures.\n")
