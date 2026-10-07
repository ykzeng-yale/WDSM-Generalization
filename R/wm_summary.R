# Presentation only: retain the estimator's existing uncertainty calculation.
.wm_result_summary <- function(object, wrapper = FALSE) {
  fit <- if (wrapper) object$fit else object
  fields <- c("estimate", "raw_estimate", "correction", "se", "variance",
              "root_n_variance")
  statistics <- vapply(fields, function(field) {
    value <- object[[field]]
    if (is.null(value)) NA_real_ else value
  }, numeric(1))
  graph <- fit$graph
  dimensions <- if (!is.null(graph$scores)) {
    stats::setNames(ncol(graph$scores), paste0("arm", graph$donor_arm))
  } else if (!is.null(graph$scores0)) {
    c(arm0 = ncol(graph$scores0),
      if (!is.null(graph$scores1)) c(arm1 = ncol(graph$scores1)))
  } else NULL
  donor_usage <- NULL
  if (!is.null(graph$edges)) {
    edges <- graph$edges
    arms <- sort(unique(edges$arm))
    donor_usage <- do.call(rbind, lapply(arms, function(arm) {
      take <- edges$arm == arm
      reuse <- tabulate(edges$donor[take], nbins = fit$n)
      data.frame(donor_arm = arm, queries = length(unique(edges$query[take])),
        selected_pairs = sum(take), donors_used = sum(reuse > 0L),
        maximum_reuse = max(reuse))
    }))
    rownames(donor_usage) <- NULL
  }
  call <- object$call
  if (is.null(call)) call <- object$fit_call
  status <- if (wrapper) object$inference$status else fit$info$inference_status
  contract <- if (wrapper) object$inference$contract else fit$info$contract
  structure(list(call = call, n = object$n, M = object$M,
    estimand = object$estimand, method = fit$method, dimensions = dimensions,
    statistics = statistics, conf.int = if (wrapper) object$conf.int else NULL,
    conf.level = if (wrapper) object$inference$conf.level else NULL,
    inference_status = status, contract = contract,
    donor_usage = donor_usage,
    graph_transform = graph$transform, tie_rule = graph$tie_rule),
    class = "summary_wm")
}

#' Summarize weighted matching results
#'
#' Compact presentation of an existing weighted matching fit. These methods
#' neither fit models nor calculate a new variance or confidence interval.
#'
#' @param object A result from \code{wm_match}, \code{wm_fit},
#'   \code{wm_wdsm_fit} or \code{wm_model_fit}. The \code{wm_match} methods
#'   also accept its potential-mean and independent-block subclasses.
#' @param x A fit for a print method, or a \code{summary_wm} result.
#' @param digits Integer number of significant digits to print, from 1 to 22.
#' @param ... Additional arguments (ignored).
#' @return Summary methods return a \code{summary_wm} list with the retained
#'   estimate, raw estimate, correction, SE and variances in \code{statistics};
#'   sample size, M, target, donor rule and available matching dimensions;
#'   retained inference status/contract; and graph donor-use counts when there
#'   is a single retained graph. Fitted wrappers additionally retain their
#'   original \code{conf.int} and \code{conf.level}. Print methods invisibly
#'   return their input.
#' @details
#' \code{summary(fit)$donor_usage} counts selected edges by donor arm; it is
#' not a covariate balance assessment or a weighted target sample size.
#' Independent-block fits have separate graphs in \code{fit$components}; their
#' summary does not merge local row identifiers or invent pooled graph counts.
#' An unavailable SE or interval is kept unavailable. Supplied-map fits do not
#' gain an interval merely by being summarized. For fitted wrappers the interval
#' is the one already produced by their selected inference route, at its recorded
#' confidence level. Summary methods do not verify identification, population
#' assumptions, numerical-root transfer or sampling coverage.
#' @seealso \code{\link{wm_match}}, \code{\link{wm_fit}},
#'   \code{\link{wm_wdsm_fit}}, \code{\link{wm_model_fit}},
#'   \code{\link{wm_bootstrap}}
#' @examples
#' x <- seq(0, 1, length.out = 16)
#' z <- rep(0:1, 8)
#' fit <- wm_fit(1 + 3 * x + 2 * z, z, 1 + x, x, M = 1)
#' print(fit)
#' summary(fit)$statistics
#' summary(fit)$donor_usage
#' @name wm_summary
NULL

#' @rdname wm_summary
#' @export
summary.wm_match <- function(object, ...) .wm_result_summary(object)

#' @rdname wm_summary
#' @export
summary.wm_wdsm_fit <- function(object, ...) .wm_result_summary(object, TRUE)

#' @rdname wm_summary
#' @export
summary.wm_model_fit <- function(object, ...) .wm_result_summary(object, TRUE)

#' @rdname wm_summary
#' @export
print.wm_match <- function(x, digits = 4L, ...) {
  print(summary.wm_match(x), digits = digits)
  invisible(x)
}

#' @rdname wm_summary
#' @export
print.wm_wdsm_fit <- function(x, digits = 4L, ...) {
  print(summary.wm_wdsm_fit(x), digits = digits)
  invisible(x)
}

#' @rdname wm_summary
#' @export
print.wm_model_fit <- function(x, digits = 4L, ...) {
  print(summary.wm_model_fit(x), digits = digits)
  invisible(x)
}

#' @rdname wm_summary
#' @export
print.summary_wm <- function(x, digits = 4L, ...) {
  if (!is.numeric(digits) || is.complex(digits) || !is.null(dim(digits)) ||
      length(digits) != 1L || !is.finite(digits) || digits != floor(digits) ||
      digits < 1L || digits > 22L) {
    stop("digits must be an integer from 1 to 22.", call. = FALSE)
  }
  cat("Weighted Matching\n")
  cat("Target:", x$estimand, " | donor rule:", x$method, "\n")
  cat("N:", x$n, " | M:", x$M, "\n")
  if (length(x$dimensions)) {
    cat("Matching dimensions:", paste(names(x$dimensions), x$dimensions,
      sep = "=", collapse = ", "), "\n")
  }
  print(x$statistics[c("estimate", "raw_estimate", "correction", "se")],
        digits = digits)
  if (!is.null(x$conf.int) && all(is.finite(x$conf.int))) {
    if (!is.null(x$conf.level)) cat(100 * x$conf.level, "% interval:\n", sep = "")
    else cat("Retained interval:\n")
    print(x$conf.int, digits = digits)
  }
  if (!is.null(x$inference_status)) cat("Inference:", x$inference_status, "\n")
  if (is.finite(x$statistics[["se"]])) {
    cat("Uncertainty retains the fitted procedure's stated assumptions.\n")
  } else cat("Sampling SE/interval unavailable for this fit.\n")
  if (!is.null(x$donor_usage)) {
    cat("Donor use (edge counts):\n")
    print(x$donor_usage, row.names = FALSE)
  }
  invisible(x)
}
