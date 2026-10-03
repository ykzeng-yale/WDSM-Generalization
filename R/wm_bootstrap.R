#' Bootstrap Replication for Weighted Matching
#'
#' Use contribution multipliers or fixed-reuse multinomial replication through
#' one interface, plus qualified Gaussian laws. The original point and weights
#' are retained; nested quadratic replication uses covariance-only local graphs.
#'
#' @param object For contribution replication, a \code{wm_match} result with
#'   variance estimation enabled or successful \code{wm_fitted_inference} with
#'   available complete rows in scalar_psm or a qualified zero-reciprocal scope,
#'   or complete corrected Gaussian variance in regular_joint_d_gt2, or prepared
#'   nested quadratic conditional-law inputs with replication_available = TRUE.
#'   For fixed_reuse, an ordinary self-normalized \code{wm_match} result.
#' @param B Positive integer number of draws, at least two for fixed_reuse.
#'   With supplied counts,
#'   infer their number of columns when B is omitted; explicit B must agree.
#' @param seed Optional nonnegative integer seed. An explicit seed restores
#'   the caller's \code{.Random.seed}, including its prior absence. The
#'   Box--Muller normal generator does not support this option because its
#'   hidden cached state cannot be restored; use \code{seed = NULL} with that
#'   generator or select another normal generator before this call.
#' @param conf.level Confidence level strictly between zero and one.
#' @param interval For contribution replication, \code{"normal"} uses the
#'   analytic conditional contribution variance and \code{"basic"} uses empirical
#'   root-n draw quantiles. For fixed_reuse, \code{"normal"} uses the original
#'   B-divisor replicate variance; basic is unsupported. \code{"none"} omits
#'   the interval under either method. Nested quadratic mixtures default to basic
#'   and reject an explicitly requested normal interval.
#' @param chunk_size Positive integer maximum number of multiplier values
#'   generated together for contribution replication. Ignored for fixed_reuse.
#'   For regular_joint_d_gt2 this bounds scalar Gaussian draw generation.
#'   No contributions-by-B matrix is allocated.
#' @param counts Required for fixed_reuse; optional for fitted contribution
#'   replication in zero-reciprocal/scalar_psm scopes; regular_joint_d_gt2 rejects
#'   them. Counts are finite nonnegative integers
#'   in original row order, with n rows and each column summing to n.
#'   Do not supply seed with frozen counts. Their draw law remains unverified.
#'   Fixed_reuse requires positive multiplicity in both arms in every column.
#' @param method \code{"contribution"} (default) or \code{"fixed_reuse"}.
#' @param refit Optional function of a full-n multiplicity vector returning
#'   mean0 and, for PATE, mean1 predictions. Only used with fixed_reuse.
#'   NULL retains the original predictions; supplied callbacks perform any
#'   intended nuisance refitting without changing graph, weights or reuse.
#' @return A list with \code{root_n_draws}, multiplier or raw mixture draws on the
#'   root-n scale; \code{draws}, estimates shifted by those draws divided by
#'   \code{sqrt(n)}; exact \code{conditional_root_n_variance} and
#'   \code{conditional_variance} in the applicable Gaussian/count branches
#'   (both NA for nested quadratic mixtures); Monte Carlo sample variances
#'   \code{monte_carlo_root_n_variance} and \code{monte_carlo_variance}; and
#'   \code{conf.int}, \code{interval}, \code{conf.level}, \code{B}, \code{n},
#'   \code{seed}, \code{estimate}, and \code{method}. Monte Carlo variances
#'   are \code{NA} when \code{B = 1}. Fixed-reuse results retain the
#'   \code{wm_bootstrap_refit} fields, including \code{variance}, \code{se},
#'   \code{variance_divisor = "B"} and prediction-refitting status.
#' @details
#' For normalized row contributions \eqn{r_i}, unordered-edge
#' contributions \eqn{e_j}, and optional independent-training contributions
#' \eqn{t_a}, a root-n draw is
#' \deqn{T^* = (\sum_i r_i \xi_i + \sum_j e_j \zeta_j +
#'                 \sum_a t_a \omega_a)/\sqrt{n},}
#' with mutually independent standard-normal multipliers. Its conditional
#' variance is \eqn{(\sum_i r_i^2 + \sum_j e_j^2 + \sum_a t_a^2)/n}.
#' Here \eqn{n} remains the evaluation sample size. Training coefficients
#' must already include any sample-size scaling; for a training sample of
#' size \eqn{m}, the independent-training adjustment supplies \eqn{n/m}
#' times its centered scalar influence. A missing or \code{NULL}
#' \code{contributions$training} is an empty vector. Components are not
#' centered separately. The fitted contributions already contain any
#' supplied nuisance influence correction. This procedure inherits the
#' inference assumptions of the fitted estimator and does not refit or
#' rematch observations.
#'
#' The normal interval uses no Monte Carlo quantiles, so increasing B does
#' not change it. With \code{seed = NULL}, draws advance the caller's RNG
#' stream. Within each replicate, multipliers are drawn in row, edge, then
#' training order. An absent training component leaves the RNG order unchanged.
#'
#' Successful fitted inference uses centered augmented_rows and the full observed n
#' for both PATE and PATT. Supplied multinomial columns give fixed-contribution
#' roots sum_i (N_i-1) augmented_row_i / sqrt(n), without refits or rematching.
#' The conditional variance under the stated Gaussian or multinomial law is
#' the centered-row square mean, not the empirical variance of supplied draws.
#' For d>=2 this requires the declared distinct_rarity, full_x, common_field
#' or patt zero-reciprocal branch, complete joint influence and full sensitivity.
#' All-scalar corrected graphs require PATE full_x or PATT patt with declared
#' zero/current-centering transport, correct full-X predictions and at least one
#' stated design-fitted, exact outcome-fitted-root or complete Gaussian sampling-law
#' alternative, retaining its weight bounds or integrated moments and current-row laws.
#' Population and numerical-root premises remain unverified; mixed dimensions
#' and scalar scopes outside these contracts are unsupported.
#' Source point, graph, dimensions and covariance identities are checked;
#' population assumptions and empirical row Lindeberg conditions remain required.
#'
#' PATE regular_joint_d_gt2 fitted objects instead use direct conditional Gaussian
#' roots with their complete reciprocal-corrected variance through this same
#' contribution dispatch. Supplied counts are rejected; augmented-row squares alone
#' omit the reciprocal term. The qualified bounded known-weight equal-d>2 baseline
#' regular-sheet sampling and feasible plug-in premises remain required. Finite-B
#' quantiles and Monte Carlo variances retain simulation error; the normal interval
#' uses analytic variance and does not depend on B. This path does not reproduce H
#' jointly, rematch, refit or silently floor the base/contrast covariance.
#'
#' With method = "fixed_reuse", supplied counts are passed to the unchanged
#' wm_bootstrap_refit calculation for supplied matching maps of any dimension.
#' The original weighted reuse and donor weights remain fixed, while an optional
#' callback replaces the mean predictions. PATE divides by the replicate
#' all-row weight total; PATT by the replicate treated-weight total. Variance
#' uses divisor B. Sampling validity requires the applicable refit expansion
#' and sampling/replication variance agreement. A fitted-inference object is
#' not automatically extracted or converted to prediction-refit replication.
#' In particular, fixed-reuse replication of raw estimated PSM with refit = NULL
#' does not add the fitted-score drift or its covariance. A unified interface
#' alone does not establish fitted-score variance validity.
#' Compatible same-count raw weighted-least-squares prediction callbacks retain the
#' separate complete-system and actual-row conditions. For deterministic polynomially
#' growing B, empirical nonlinear variance consistency additionally requires the
#' stated prediction derivative/Hessian moment above order two and simultaneous
#' selected-root, numerical-error and success conditions; see wm_bootstrap_refit.
#' Prepared nested quadratic inputs use the literal untruncated Gaussian mixture:
#' full contrast QR generates h, its evaluator returns m(h), s2(h), and the root
#' is m(h)+sqrt(s2(h))*eta. All B slots and innovations are retained. Any failed
#' slot invalidates aggregate inference; no cutoff, zero completion or redraw is
#' supported. The original point is unchanged. Counts, refit and fixed_reuse are
#' rejected. B=1 supplies draws only. B-1 sample variance over all raw roots/n
#' includes variation in m(h); basic intervals use uncentered root quantiles.
#' Exact finite-data conditional variance fields remain NA. Predeclared increasing
#' B and continuous unique limiting quantiles are required; empirical variance
#' additionally requires at-most-polynomial B. These are not runtime cutoffs.
#' Numerical RMS/point-error and population premises remain unverified; see the
#' manual for explicit scope, allocation caps and numerical qualifications.
#' @seealso \code{\link{wm_match}}, \code{\link{wm_bootstrap_refit}}
#' @export
wm_bootstrap <- function(object, B = 999L, seed = NULL, conf.level = 0.95,
                         interval = c("normal", "basic", "none"),
                         chunk_size = 65536L, counts = NULL,
                         method = c("contribution", "fixed_reuse"), refit = NULL) {
  method <- match.arg(method)
  if (inherits(object, "wm_fitted_inference") &&
      identical(object$status, "model_quadratic_contrast_inputs_ready")) {
    if (method != "contribution" || !is.null(refit)) {
      stop("Nested quadratic law replication does not support fixed_reuse or refit; the original point is retained.", call. = FALSE)
    }
    return(.wm_qc_bootstrap(object, B, seed, conf.level, interval, chunk_size, counts))
  }
  if (method == "fixed_reuse") {
    return(.wm_fixed_reuse_bootstrap(object, B, missing(B), seed, conf.level,
                                    interval, counts, refit))
  }
  if (!is.null(refit)) {
    stop("refit is supported only with method = 'fixed_reuse'.", call. = FALSE)
  }
  if (inherits(object, "wm_fitted_inference")) {
    if (identical(object$covariance_scope, "critical_planar")) {
      return(.wm_cp_bootstrap(object, B, seed, conf.level, interval, chunk_size, counts))
    }
    if (identical(object$covariance_scope, "scalar_psm")) {
      return(.wm_scalar_fitted_bootstrap(object, B, missing(B), seed, conf.level,
                                         interval, chunk_size, counts))
    }
    return(.wm_general_fitted_bootstrap(object, B, missing(B), seed, conf.level,
                                        interval, chunk_size, counts))
  }
  if (!is.null(counts)) {
    stop("Supplied contribution counts require successful wm_fitted_inference.",
         call. = FALSE)
  }
  positive_integer <- function(x) {
    is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x) &&
      x >= 1 && x <= .Machine$integer.max && x == floor(x)
  }
  finite_scalar <- function(x) {
    is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x)
  }
  finite_vector <- function(x) {
    is.numeric(x) && !is.complex(x) && is.null(dim(x)) && all(is.finite(x))
  }
  if (!inherits(object, "wm_match") || !is.list(object)) {
    stop("object must be a wm_match result", call. = FALSE)
  }
  if (!positive_integer(B)) {
    stop("B must be a positive integer", call. = FALSE)
  }
  if (!positive_integer(chunk_size)) {
    stop("chunk_size must be a positive integer", call. = FALSE)
  }
  if (!finite_scalar(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must lie strictly between zero and one", call. = FALSE)
  }
  interval <- match.arg(interval)
  if (!is.null(seed) &&
      (!finite_scalar(seed) || seed < 0 || seed > .Machine$integer.max ||
       seed != floor(seed))) {
    stop("seed must be NULL or a nonnegative integer", call. = FALSE)
  }
  if (!positive_integer(object$n) || !finite_scalar(object$estimate)) {
    stop("object must contain a positive integer n and finite estimate",
         call. = FALSE)
  }
  if (!finite_scalar(object$root_n_variance) || object$root_n_variance < 0) {
    stop("object must have finite nonnegative root_n_variance; enable variance estimation",
         call. = FALSE)
  }
  if (!is.list(object$contributions) ||
      !finite_vector(object$contributions$row) ||
      length(object$contributions$row) != object$n ||
      !finite_vector(object$contributions$edge)) {
    stop("object must contain finite row and edge contribution vectors, with n rows",
         call. = FALSE)
  }
  training <- object$contributions$training
  if (is.null(training)) training <- numeric()
  if (!finite_vector(training)) {
    stop("training contributions must be a finite numeric vector when present",
         call. = FALSE)
  }

  n <- object$n
  B <- as.integer(B)
  chunk_size <- as.integer(chunk_size)
  # Keep the exact decomposition and upstream training scaling: separate
  # centering changes the stabilized PATE variance, and actual rows omit its
  # edge correction. Independent training coefficients come after both.
  coefficients <- c(object$contributions$row, object$contributions$edge, training)
  scaled <- coefficients / sqrt(n)
  contribution_variance <- sum(scaled * scaled)
  tolerance <- 100 * .Machine$double.eps *
    max(contribution_variance, object$root_n_variance, .Machine$double.xmin)
  if (!is.finite(contribution_variance) ||
      abs(contribution_variance - object$root_n_variance) > tolerance) {
    stop("root_n_variance does not agree with the row, edge and training contributions",
         call. = FALSE)
  }

  if (!is.null(seed)) {
    # set.seed() clears a Box-Muller cache not represented by .Random.seed.
    # Reject before changing any state instead of silently losing that cache.
    if (identical(RNGkind()[2L], "Box-Muller")) {
      stop("an explicit seed cannot preserve the Box-Muller cache; use seed = NULL or another normal RNG",
           call. = FALSE)
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    previous_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv,
                                      inherits = FALSE) else NULL
    on.exit({
      if (had_seed) {
        assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  root_n_draws <- numeric(B)
  count <- length(scaled)
  starts <- seq.int(1, count, by = chunk_size)
  for (b in seq_len(B)) {
    value <- 0
    for (start in starts) {
      index <- seq.int(start, min(count, start + as.double(chunk_size) - 1))
      value <- value + sum(scaled[index] * stats::rnorm(length(index)))
    }
    root_n_draws[b] <- value
  }
  draws <- object$estimate + root_n_draws / sqrt(n)
  exact_root_n_variance <- object$root_n_variance
  exact_variance <- exact_root_n_variance / n
  mc_root_n_variance <- if (B > 1L) stats::var(root_n_draws) else NA_real_
  alpha <- 1 - conf.level
  conf.int <- switch(interval,
    normal = object$estimate + c(-1, 1) *
      stats::qnorm(alpha / 2, lower.tail = FALSE) * sqrt(exact_variance),
    basic = object$estimate - as.numeric(stats::quantile(
      root_n_draws, probs = c(1 - alpha / 2, alpha / 2), names = FALSE)) / sqrt(n),
    none = NULL
  )
  if (!is.null(conf.int)) names(conf.int) <- c("lower", "upper")
  list(
    root_n_draws = root_n_draws, draws = draws,
    conditional_root_n_variance = exact_root_n_variance,
    conditional_variance = exact_variance,
    monte_carlo_root_n_variance = mc_root_n_variance,
    monte_carlo_variance = mc_root_n_variance / n,
    conf.int = conf.int, interval = interval, conf.level = conf.level,
    B = B, n = n, seed = seed, estimate = object$estimate,
    method = "Gaussian contribution multipliers"
  )
}

# Explicit common-interface dispatch only; the refit algebra stays in its
# independently reviewed implementation and no count columns are generated.
.wm_fixed_reuse_bootstrap <- function(object, B, B_missing, seed, conf.level,
                                      interval, counts, refit) {
  if (inherits(object, "wm_fitted_inference") ||
      !inherits(object, "wm_match") || !is.list(object)) {
    stop("fixed_reuse requires an ordinary self_normalized wm_match result.",
         call. = FALSE)
  }
  interval <- match.arg(interval, c("normal", "basic", "none"))
  if (interval == "basic") {
    stop("fixed_reuse supports interval = 'normal' or 'none'; basic is unsupported.",
         call. = FALSE)
  }
  if (!is.null(seed)) {
    stop("Do not supply seed with frozen counts; this path uses no RNG.", call. = FALSE)
  }
  if (!is.matrix(counts)) {
    stop("fixed_reuse requires a supplied n-by-B count matrix with B >= 2.",
         call. = FALSE)
  }
  if (B_missing) B <- ncol(counts)
  if (!is.numeric(B) || is.complex(B) || length(B) != 1L || !is.finite(B) ||
      B < 2 || B > .Machine$integer.max || B != floor(B) || B != ncol(counts)) {
    stop("B must be an integer >= 2 equal to the supplied count columns.",
         call. = FALSE)
  }
  if (!is.null(rownames(counts)) &&
      !identical(rownames(counts), as.character(seq_len(nrow(counts))))) {
    stop("Count row labels must be original observation indices 1,...,n in order.",
         call. = FALSE)
  }
  result <- wm_bootstrap_refit(object, counts, refit = refit, conf.level = conf.level)
  if (interval == "none") result$conf.int <- NULL
  result
}

# The reviewed scalar complete rows supply a one-step contribution law. This
# adapter reuses the unchanged wm_match Gaussian generator; it never fits or
# invokes matching, and supplied counts never touch the RNG.
.wm_scalar_fitted_bootstrap <- function(object, B, B_missing, seed, conf.level,
                                        interval, chunk_size, counts) {
  positive_integer <- function(x) {
    is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x) &&
      x >= 1 && x <= .Machine$integer.max && x == floor(x)
  }
  finite_scalar <- function(x) {
    is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x)
  }
  finite_vector <- function(x) {
    is.numeric(x) && !is.complex(x) && is.null(dim(x)) && all(is.finite(x))
  }
  if (!is.list(object) || !identical(object$covariance_scope, "scalar_psm") ||
      !identical(object$status, "conditional_scalar_psm") ||
      !isTRUE(object$available) || !isTRUE(object$conditional_inference_available)) {
    stop("Require successful scalar_psm wm_fitted_inference with available complete rows.",
         call. = FALSE)
  }
  n <- object$n
  general_scalar <- !is.null(object$scalar_model)
  source_bound <- if (general_scalar) {
    inherits(object$source_object, "wm_match") && identical(object$source_object, object$fit)
  } else {
    inherits(object$source_object, "wm_scalar_logistic_match") &&
      identical(object$source_object$fit, object$fit) &&
      identical(object$source_object$estimate, object$estimate)
  }
  if (!positive_integer(n) || !finite_scalar(object$estimate) ||
      !is.character(object$estimand) || length(object$estimand) != 1L ||
      !object$estimand %in% c("PATE", "PATT") ||
      !inherits(object$fit, "wm_match") ||
      !identical(object$fit$n, n) || !identical(object$fit$estimate, object$estimate) ||
      !identical(object$fit$estimand, object$estimand) || !source_bound) {
    stop("Scalar contribution inference must retain the same raw point and full observed n.",
         call. = FALSE)
  }
  if (general_scalar) .wm_scalar_psm_general_rows(object)
  row <- object$augmented_rows
  if (!finite_vector(row) || length(row) != n ||
      !finite_scalar(object$root_n_variance) || object$root_n_variance <= 0) {
    stop("Scalar inference must contain n finite augmented rows and positive root-n variance.",
         call. = FALSE)
  }
  centered <- row - mean(row)
  V <- sum((centered / sqrt(n))^2)
  tolerance <- 100 * .Machine$double.eps *
    max(V, object$root_n_variance, .Machine$double.xmin)
  if (!is.finite(V) || V <= 0 ||
      abs(V - object$root_n_variance) > tolerance ||
      abs(mean(row)) > 100 * .Machine$double.eps *
        max(sqrt(V), .Machine$double.xmin)) {
    stop("Scalar root-n variance does not agree with its centered augmented rows.",
         call. = FALSE)
  }
  if (is.null(counts)) {
    proxy <- structure(list(n = n, estimate = object$estimate,
      root_n_variance = V, contributions = list(row = centered, edge = numeric())),
      class = c("wm_match", "list"))
    result <- wm_bootstrap(proxy, B = B, seed = seed, conf.level = conf.level,
      interval = interval, chunk_size = chunk_size)
    result$method <- "Gaussian scalar fitted-score contribution multipliers"
    result$draw_input <- "generated_gaussian_multipliers"
    result$supplied_draw_law_verified <- NA
  } else {
    if (!is.null(seed)) {
      stop("Do not supply seed with frozen counts; this path uses no RNG.", call. = FALSE)
    }
    if (!is.matrix(counts) || !is.numeric(counts) || is.complex(counts) ||
        nrow(counts) != n || !positive_integer(ncol(counts)) ||
        anyNA(counts) || any(!is.finite(counts)) ||
        any(counts < 0 | counts != floor(counts)) || any(colSums(counts) != n)) {
      stop("counts must be a finite nonnegative integer n-by-B matrix with each column summing to full n.",
           call. = FALSE)
    }
    if (!is.null(rownames(counts)) &&
        !identical(rownames(counts), as.character(seq_len(n)))) {
      stop("Count row labels must be original observation indices 1,...,n in order.",
           call. = FALSE)
    }
    if (B_missing) B <- ncol(counts)
    if (!positive_integer(B) || B != ncol(counts)) {
      stop("B must be a positive integer equal to the supplied count columns.", call. = FALSE)
    }
    if (!positive_integer(chunk_size)) stop("chunk_size must be a positive integer", call. = FALSE)
    if (!finite_scalar(conf.level) || conf.level <= 0 || conf.level >= 1) {
      stop("conf.level must lie strictly between zero and one", call. = FALSE)
    }
    interval <- match.arg(interval, c("normal", "basic", "none"))
    B <- as.integer(B)
    chunk_size <- as.integer(chunk_size)
    starts <- seq.int(1L, n, by = chunk_size)
    scaled <- centered / sqrt(n)
    roots <- numeric(B)
    for (b in seq_len(B)) {
      for (start in starts) {
        index <- seq.int(start, min(n, start + as.double(chunk_size) - 1))
        roots[b] <- roots[b] + sum(scaled[index] * (counts[index, b] - 1))
      }
    }
    draws <- object$estimate + roots / sqrt(n)
    exact_variance <- V / n
    mc_variance <- if (B > 1L) stats::var(roots) else NA_real_
    alpha <- 1 - conf.level
    confidence_interval <- switch(interval,
      normal = object$estimate + c(-1, 1) *
        stats::qnorm(alpha / 2, lower.tail = FALSE) * sqrt(exact_variance),
      basic = object$estimate - as.numeric(stats::quantile(roots,
        probs = c(1 - alpha / 2, alpha / 2), names = FALSE)) / sqrt(n),
      none = NULL)
    if (!is.null(confidence_interval)) names(confidence_interval) <- c("lower", "upper")
    result <- list(root_n_draws = roots, draws = draws,
      conditional_root_n_variance = V, conditional_variance = exact_variance,
      monte_carlo_root_n_variance = mc_variance, monte_carlo_variance = mc_variance / n,
      conf.int = confidence_interval, interval = interval, conf.level = conf.level,
      B = B, n = n, seed = NULL, estimate = object$estimate,
      method = "Supplied full-n scalar contribution counts",
      draw_input = "supplied_full_n_multinomial_counts",
      supplied_counts = counts, supplied_draw_law_verified = FALSE)
  }
  result$estimand <- object$estimand
  result$covariance_scope <- "scalar_psm"
  result$source_inference <- object
  result$population_assumptions_verified <- FALSE
  result$original_refit_bootstrap <- FALSE
  result$replication_contract <- paste(
    "Fixed centered augmented contributions and full observed n, under the scalar",
    "inference object's population and numerical-root premises. Gaussian multipliers",
    "are independent standard normal; supplied counts require independent columns from",
    "multinomial(n;1/n,...,1/n). Supplied matrices do not verify that draw law.",
    "Conditional variance is the centered-row square mean, not empirical draw variance.",
    "No nuisance refit, matching rerun or replicate treated-weight denominator.")
  result
}

# Public qualified fitted inference retains its corrected source fit and complete
# joint covariance blocks. Reconstruct their finite-data identities before
# treating the augmented rows as bootstrap coefficients; a row-only label is
# not sufficient provenance. This does not check population assumptions.
.wm_general_fitted_rows <- function(object, covariance_state = FALSE) {
  if (!is.logical(covariance_state) || length(covariance_state) != 1L ||
      is.na(covariance_state)) stop("Invalid covariance-state request.", call. = FALSE)
  scope <- object$covariance_scope
  reciprocal_scope <- identical(scope, "regular_joint_d_gt2")
  if (!is.list(object) || anyDuplicated(names(object)) ||
      !is.character(scope) || length(scope) != 1L || is.na(scope) ||
      !scope %in% c("distinct_rarity", "full_x", "common_field", "patt",
                    "regular_joint_d_gt2") ||
      !identical(object$status, paste0("conditional_", scope)) ||
      !isTRUE(object$available) || !isTRUE(object$numerically_available) ||
      !isTRUE(object$conditional_inference_available) ||
      !is.numeric(object$reciprocal_subtraction) ||
      (!reciprocal_scope && !identical(as.vector(object$reciprocal_subtraction), 0)) ||
      !is.list(object$contract) ||
      !all(c("target", "fitted_law", "covariance", "transport", "derivatives",
             "smooth_weights", "inference") %in% names(object$contract))) {
    stop("Require successful qualified public wm_fitted_inference.",
         call. = FALSE)
  }
  if (reciprocal_scope && !covariance_state) {
    stop("regular_joint_d_gt2 requires corrected-covariance Gaussian replication, not row/count multipliers.",
         call. = FALSE)
  }
  state <- .wm_fv_state(object$fit)
  n <- state$n
  if (!identical(object$n, object$fit$n) || object$n != n ||
      !identical(object$M, object$fit$M) ||
      !identical(object$estimate, object$fit$estimate) ||
      !identical(object$estimand, object$fit$estimand)) {
    stop("Fitted contribution inference must retain the same source point, M and full observed n.",
         call. = FALSE)
  }
  dimensions <- c(potential0 = ncol(object$fit$graph$scores0),
    if (state$pATE) c(potential1 = ncol(object$fit$graph$scores1)))
  scalar <- all(dimensions == 1L)
  if (!identical(object$dimensions, dimensions) ||
      (!all(dimensions >= 2L) && !(scalar && scope %in% c("full_x", "patt"))) ||
      (state$pATE && scope == "patt") || (!state$pATE && scope != "patt") ||
      (scope == "common_field" && dimensions[1L] != dimensions[2L]) ||
      (reciprocal_scope && (!state$pATE || length(dimensions) != 2L ||
        any(dimensions <= 2L) || dimensions[1L] != dimensions[2L]))) {
    stop("Fitted contribution scope and used matching dimensions disagree with the source fit.",
         call. = FALSE)
  }
  parameters <- object$parameter_names
  if (!is.character(parameters) || !length(parameters) || anyNA(parameters) ||
      any(!nzchar(parameters)) || anyDuplicated(parameters)) {
    stop("Fitted contributions require the complete named parameter order.", call. = FALSE)
  }
  p <- length(parameters)
  influence <- .wm_fi_matrix(object$nuisance_influence, n, parameters,
                             "stored nuisance_influence")
  named <- function(x, label) {
    value <- .wm_fv_number(x, p, label)
    if (!identical(names(value), parameters)) {
      stop(label, " must retain the exact parameter names/order.", call. = FALSE)
    }
    value
  }
  smooth <- named(object$smooth_sensitivity, "smooth_sensitivity")
  graph <- named(object$graph_sensitivity, "graph_sensitivity")
  total <- named(object$total_sensitivity, "total_sensitivity")
  derivative <- .wm_reciprocal_list(object$smooth_derivative, "smooth_derivative")
  if (!identical(derivative$parameter_names, parameters) ||
      !identical(derivative$estimate, object$estimate)) {
    stop("Smooth derivative must remain bound to the source point and parameter order.",
         call. = FALSE)
  }
  .wm_reciprocal_agree(named(derivative$sensitivity, "smooth derivative"), smooth,
                       "Stored smooth derivative")
  .wm_reciprocal_agree(named(derivative$weight_sensitivity, "weight sensitivity") +
    named(derivative$prediction_sensitivity, "prediction sensitivity"), smooth,
    "Complete smooth derivative")
  if (scalar && any(named(derivative$weight_sensitivity, "weight sensitivity") != 0)) {
    stop("Scalar full-X fitted contributions require supplied known weights and zero weight sensitivity.",
         call. = FALSE)
  }
  if (reciprocal_scope &&
      any(named(derivative$weight_sensitivity, "weight sensitivity") != 0)) {
    stop("regular_joint_d_gt2 requires supplied known weights and zero weight sensitivity.",
         call. = FALSE)
  }
  .wm_reciprocal_agree(total, smooth + graph, "Complete fitted slope")
  gamma <- mean(if (state$pATE) object$fit$weights else
    object$fit$data$Z * object$fit$weights)
  .wm_reciprocal_agree(object$raw_target_weight_mean, gamma,
                       "Raw target-weight normalization")
  transports <- .wm_reciprocal_list(object$graph_transport, "graph_transport")
  if (!identical(names(transports), names(dimensions))) {
    stop("Graph transport directions must agree with the source matching maps.", call. = FALSE)
  }
  j <- error <- matrix(0, length(dimensions), p)
  for (k in seq_along(dimensions)) {
    value <- .wm_reciprocal_list(transports[[k]], "graph transport direction")
    if (!isTRUE(value$numerically_available) ||
        !is.character(value$mode) || length(value$mode) != 1L ||
        is.na(value$mode) || !value$mode %in% c("zero", "estimate")) {
      stop("All stored graph transport directions must be available.", call. = FALSE)
    }
    if (scalar && value$mode != "zero") {
      stop("Scalar full-X fitted contributions require declared zero graph transport.",
           call. = FALSE)
    }
    j[k, ] <- named(value$graph_drift, "graph drift")
    error[k, ] <- named(value$quadrature_error_bound, "graph quadrature bound")
    if (any(error[k, ] < 0)) stop("Graph quadrature bounds must be nonnegative.", call. = FALSE)
    if (value$mode == "zero") {
      allowed <- if (scope == "full_x") "full_x" else if (scalar) "current_centering" else
        c("fixed_map", "current_centering")
      if (!is.character(value$basis) || length(value$basis) != 1L ||
          is.na(value$basis) || !value$basis %in% allowed ||
          any(c(j[k, ], error[k, ]) != 0)) {
        stop("Declared zero graph transport is inconsistent with its scope.", call. = FALSE)
      }
    } else if (scope == "full_x" || !identical(value$n, n) ||
               value$d != dimensions[k] || value$M != object$M ||
               value$donor_arm != k - 1L ||
               !identical(value$parameter_names, parameters) ||
               !is.list(value$map_binding)) {
      stop("Estimated graph transport must retain the source direction and parameter binding.",
           call. = FALSE)
    }
  }
  signs <- if (state$pATE) c(-1, 1) else -1
  .wm_reciprocal_agree(graph, as.vector(crossprod(signs, j)) / gamma,
                       "Normalized graph sensitivity")
  .wm_reciprocal_agree(named(object$graph_sensitivity_quadrature_error_bound,
    "graph sensitivity quadrature bound"), colSums(error) / gamma,
    "Normalized graph quadrature bound")
  # No coefficient fit, graph construction, transport estimation or RNG call.
  expected <- .wm_wdsm_fitted_variance(object$fit, influence, smooth, graph,
                                      scope, object$conf.level)
  for (field in c("base_rows", "augmented_rows", "V0", "C", "cross_term",
                  "nuisance_variance", "root_n_variance", "variance", "se")) {
    value <- .wm_fv_number(object[[field]], length(expected[[field]]), field)
    if (field %in% c("base_rows", "augmented_rows")) {
      .wm_fi_row_names(names(value), n, paste(field, "indices"))
    }
    .wm_reciprocal_agree(value, expected[[field]], paste("Complete fitted", field))
  }
  # This diagnostic is a cancellation residual in outcome-squared units.
  # Re-centering the retained influence can change its last bits without
  # changing any complete covariance block. Compare using covariance-block
  # magnitudes, without a currency-unit floor, rather than the near-zero residual.
  formula_error <- .wm_fv_number(object$covariance_formula_error, 1L,
                                  "covariance_formula_error")
  formula_scale <- max(abs(c(expected$V0, expected$cross_term,
    expected$nuisance_variance, if (reciprocal_scope)
      c(expected$diagonal_root_n_variance, expected$reciprocal_subtraction))))
  if (formula_error < 0 || !is.finite(formula_scale) || formula_scale <= 0) {
    stop("Complete fitted covariance_formula_error requires a nonnegative residual and a finite positive covariance-block scale.",
         call. = FALSE)
  }
  .wm_reciprocal_agree(formula_error / formula_scale,
    expected$covariance_formula_error / formula_scale,
    "Complete fitted covariance_formula_error")
  if (!identical(names(object$C), parameters) ||
      !is.matrix(object$Sigma) || !is.numeric(object$Sigma) ||
      is.complex(object$Sigma) || !identical(dim(object$Sigma), c(p, p)) ||
      !identical(dimnames(object$Sigma), list(parameters, parameters)) ||
      any(!is.finite(object$Sigma))) {
    stop("Complete covariance blocks must retain the exact parameter axes.", call. = FALSE)
  }
  .wm_reciprocal_agree(object$Sigma, expected$Sigma, "Complete joint covariance")
  .wm_reciprocal_agree(influence, expected$nuisance_influence, "Centered joint influence")
  if (reciprocal_scope) {
    for (field in c("reciprocal_subtraction", "reciprocal_numerator",
                    "diagonal_root_n_variance", "augmented_row_root_n_variance",
                    "unregularized_root_n_variance", "unregularized_variance",
                    "analysis_target_weight_mean", "weight_scale")) {
      value <- .wm_fv_number(object[[field]], 1L, field)
      .wm_reciprocal_agree(value, expected[[field]], paste("Corrected fitted", field))
    }
    pairs <- object$reciprocal_pairs
    if (!is.data.frame(pairs) || !identical(names(pairs), names(expected$reciprocal_pairs)) ||
        nrow(pairs) != nrow(expected$reciprocal_pairs)) {
      stop("Stored reciprocal pair structure is inconsistent with the source graph.", call. = FALSE)
    }
    for (field in names(pairs)) {
      .wm_reciprocal_agree(pairs[[field]], expected$reciprocal_pairs[[field]],
                           paste("Reciprocal pair", field))
    }
    if (!identical(object$variance_floor, NULL) || !identical(object$floor_active, FALSE) ||
        !identical(object$reciprocal_units, expected$reciprocal_units)) {
      stop("Corrected fitted inference must retain unfloored, declared reciprocal units.", call. = FALSE)
    }
  }
  if (!isTRUE(expected$available)) stop("Complete fitted variance is unavailable.", call. = FALSE)
  row <- object$augmented_rows
  mean_scale <- if (reciprocal_scope) sqrt(expected$augmented_row_root_n_variance) else
    sqrt(expected$root_n_variance)
  if (abs(mean(row)) > 100 * .Machine$double.eps *
      max(mean_scale, .Machine$double.xmin)) {
    stop("Augmented rows must already be centered.", call. = FALSE)
  }
  if (covariance_state) return(list(rows = row - mean(row), inference = expected))
  row - mean(row)
}

.wm_general_fitted_bootstrap <- function(object, B, B_missing, seed, conf.level,
                                         interval, chunk_size, counts) {
  if (identical(object$covariance_scope, "regular_joint_d_gt2")) {
    if (!is.null(counts)) {
      stop("regular_joint_d_gt2 rejects supplied counts: row/count multipliers omit reciprocal covariance.",
           call. = FALSE)
    }
    state <- .wm_general_fitted_rows(object, covariance_state = TRUE)
    return(.wm_corrected_fitted_gaussian(object, state$inference, B, seed,
                                        conf.level, interval, chunk_size))
  }
  row <- .wm_general_fitted_rows(object)
  n <- object$n
  V <- sum((row / sqrt(n))^2)
  if (is.null(counts)) {
    proxy <- structure(list(n = n, estimate = object$estimate,
      root_n_variance = V, contributions = list(row = row, edge = numeric())),
      class = c("wm_match", "list"))
    result <- wm_bootstrap(proxy, B = B, seed = seed, conf.level = conf.level,
      interval = interval, chunk_size = chunk_size)
    result$method <- "Gaussian complete fitted contribution multipliers"
    result$draw_input <- "generated_gaussian_multipliers"
    result$supplied_draw_law_verified <- NA
  } else {
    positive_integer <- function(x) {
      is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x) &&
        x >= 1 && x <= .Machine$integer.max && x == floor(x)
    }
    if (!is.null(seed)) {
      stop("Do not supply seed with frozen counts; this path uses no RNG.", call. = FALSE)
    }
    if (!is.matrix(counts) || !is.numeric(counts) || is.complex(counts) ||
        nrow(counts) != n || !positive_integer(ncol(counts)) ||
        anyNA(counts) || any(!is.finite(counts)) ||
        any(counts < 0 | counts != floor(counts)) || any(colSums(counts) != n)) {
      stop("counts must be a finite nonnegative integer n-by-B matrix with each column summing to full n.",
           call. = FALSE)
    }
    if (!is.null(rownames(counts)) &&
        !identical(rownames(counts), as.character(seq_len(n)))) {
      stop("Count row labels must be original observation indices 1,...,n in order.", call. = FALSE)
    }
    if (B_missing) B <- ncol(counts)
    if (!positive_integer(B) || B != ncol(counts)) {
      stop("B must be a positive integer equal to the supplied count columns.", call. = FALSE)
    }
    if (!positive_integer(chunk_size)) stop("chunk_size must be a positive integer", call. = FALSE)
    if (!is.numeric(conf.level) || is.complex(conf.level) || length(conf.level) != 1L ||
        !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
      stop("conf.level must lie strictly between zero and one", call. = FALSE)
    }
    interval <- match.arg(interval, c("normal", "basic", "none"))
    B <- as.integer(B)
    starts <- seq.int(1L, n, by = as.integer(chunk_size))
    scaled <- row / sqrt(n)
    roots <- numeric(B)
    for (b in seq_len(B)) for (start in starts) {
      index <- seq.int(start, min(n, start + as.double(chunk_size) - 1))
      roots[b] <- roots[b] + sum(scaled[index] * (counts[index, b] - 1))
    }
    draws <- object$estimate + roots / sqrt(n)
    mc <- if (B > 1L) stats::var(roots) else NA_real_
    alpha <- 1 - conf.level
    ci <- switch(interval,
      normal = object$estimate + c(-1, 1) *
        stats::qnorm(alpha / 2, lower.tail = FALSE) * sqrt(V / n),
      basic = object$estimate - as.numeric(stats::quantile(roots,
        probs = c(1 - alpha / 2, alpha / 2), names = FALSE)) / sqrt(n),
      none = NULL)
    if (!is.null(ci)) names(ci) <- c("lower", "upper")
    result <- list(root_n_draws = roots, draws = draws,
      conditional_root_n_variance = V, conditional_variance = V / n,
      monte_carlo_root_n_variance = mc, monte_carlo_variance = mc / n,
      conf.int = ci, interval = interval, conf.level = conf.level,
      B = B, n = n, seed = NULL, estimate = object$estimate,
      method = "Supplied full-n complete fitted contribution counts",
      draw_input = "supplied_full_n_multinomial_counts",
      supplied_counts = counts, supplied_draw_law_verified = FALSE)
  }
  result$estimand <- object$estimand
  result$covariance_scope <- object$covariance_scope
  result$source_inference <- object
  result$population_assumptions_verified <- FALSE
  result$original_refit_bootstrap <- FALSE
  result$replication_contract <- paste(
    "Fixed complete centered augmented contributions and full observed n, under the",
    "fitted inference object's zero-reciprocal population and joint-influence premises.",
    "Conditional Gaussian multipliers are independent standard normal; supplied counts",
    "require independent multinomial(n;1/n,...,1/n) columns. Supplied arrays do not",
    "verify that law. Exact conditional variance is the centered-row square mean,",
    "not empirical draw variance. Sampling transfer additionally needs the empirical",
    "row variance and Lindeberg conditions. No nuisance refit, rematching or replicate",
    "treated-weight denominator; this does not validate original-refit replication.")
  result
}

# Direct marginal Gaussian calibration inside the common fitted-object dispatch.
# Its variance is the full corrected contrast, not an augmented-row square mean.
.wm_corrected_fitted_gaussian <- function(object, checked, B, seed,
                                         conf.level, interval, chunk_size) {
  scalar <- function(x) is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
    length(x) == 1L && is.finite(x)
  positive_integer <- function(x) scalar(x) && x >= 1 &&
    x <= .Machine$integer.max && x == floor(x)
  if (!positive_integer(B)) stop("B must be a positive integer.", call. = FALSE)
  if (!positive_integer(chunk_size)) stop("chunk_size must be a positive integer.", call. = FALSE)
  if (!scalar(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must lie strictly between zero and one.", call. = FALSE)
  }
  interval <- match.arg(interval, c("normal", "basic", "none"))
  if (!is.null(seed) && (!scalar(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != floor(seed))) {
    stop("seed must be NULL or a nonnegative integer.", call. = FALSE)
  }
  V <- checked$root_n_variance
  if (!isTRUE(checked$available) || !scalar(V) || V <= 0) {
    stop("Complete reciprocal-corrected fitted variance is unavailable.", call. = FALSE)
  }
  if (!is.null(seed)) {
    if (identical(RNGkind()[2L], "Box-Muller")) {
      stop("An explicit seed cannot preserve the Box-Muller cache; use seed = NULL or another normal RNG.",
           call. = FALSE)
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    previous_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv,
                                      inherits = FALSE) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
  }
  B <- as.integer(B)
  roots <- numeric(B)
  for (start in seq.int(1L, B, by = as.integer(chunk_size))) {
    index <- seq.int(start, min(B, start + as.double(chunk_size) - 1))
    roots[index] <- sqrt(V) * stats::rnorm(length(index))
  }
  n <- object$n
  draws <- object$estimate + roots / sqrt(n)
  mc <- if (B > 1L) stats::var(roots) else NA_real_
  alpha <- 1 - conf.level
  ci <- switch(interval,
    normal = object$estimate + c(-1, 1) *
      stats::qnorm(alpha / 2, lower.tail = FALSE) * sqrt(V / n),
    basic = object$estimate - as.numeric(stats::quantile(roots,
      c(1 - alpha / 2, alpha / 2), names = FALSE)) / sqrt(n),
    none = NULL)
  if (any(!is.finite(c(roots, draws, ci))) || (B > 1L && !is.finite(mc))) {
    stop("Corrected fitted Gaussian arithmetic exceeded numerical range.", call. = FALSE)
  }
  if (!is.null(ci)) names(ci) <- c("lower", "upper")
  list(root_n_draws = roots, draws = draws,
    conditional_root_n_variance = V, conditional_variance = V / n,
    monte_carlo_root_n_variance = mc, monte_carlo_variance = mc / n,
    conf.int = ci, interval = interval, conf.level = conf.level,
    B = B, n = n, seed = seed, estimate = object$estimate,
    method = "Gaussian replication with reciprocal-corrected fitted-map variance",
    draw_input = "generated_gaussian_corrected_variance",
    supplied_draw_law_verified = NA, estimand = object$estimand,
    covariance_scope = object$covariance_scope, source_inference = object,
    population_assumptions_verified = FALSE, original_refit_bootstrap = FALSE,
    replication_contract = paste("Conditional marginal N(0,Vhat) with the complete",
      "reciprocal-corrected fitted contrast variance and full observed n. Qualified",
      "known-W equal-d>2 regular baseline-sheet sampling and feasible plug-in premises",
      "remain unverified. This is variance-calibrated Gaussian replication, not row/count",
      "multipliers, rematching or nuisance-refit validity. No joint nuisance draws or",
      "actual sampling-variance UI claim. Normal intervals use analytic variance; finite-B",
      "quantiles/sample variances retain Monte Carlo error (sample divisor B-1)."))
}
