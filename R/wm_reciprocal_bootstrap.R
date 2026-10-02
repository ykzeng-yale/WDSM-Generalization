#' Gaussian Replication with Reciprocal-Corrected Matching Variance
#'
#' Generate centered Gaussian roots using the fixed-map variance returned by
#' \code{wm_reciprocal_inference()}. This is a scalar Gaussian approximation,
#' distinct from row multipliers and original WDSM nuisance-refit replication.
#'
#' @param object An available \code{wm_reciprocal_inference} result.
#' @param B Positive integer number of draws.
#' @param seed Optional nonnegative integer seed, preserving the caller's RNG
#'   state, including its prior absence. Explicit seeds with the Box--Muller
#'   normal generator are rejected because its hidden cache cannot be restored.
#' @param conf.level Confidence level strictly between zero and one.
#' @param interval \code{"normal"} for a normal interval using the exact
#'   conditional variance, \code{"basic"} for simulated quantiles, or
#'   \code{"none"} for no interval.
#' @return A list of centered \code{root_n_draws}, shifted \code{draws}, exact
#'   conditional root-n and estimator variances, Monte Carlo variances using
#'   divisor B-1 (NA for B=1), confidence interval and procedure metadata.
#' @details
#' The root is \eqn{T^*=\sqrt{\widehat V}G}, where G is standard normal and
#' \eqn{\widehat V} is the corrected root-n variance, including any explicitly
#'   requested variance floor. The estimator draw is estimate plus
#'   \eqn{T^*/\sqrt n}. The function neither rematches nor refits.
#'
#' Sampling validity requires the separate fixed-map sampling CLT, consistent
#' reciprocal variance, positive variance lower bound and a negligible root-n
#' feasible-point remainder. The exact conditional Gaussian variance alone is
#' not a proof of sampling validity. Actual feasible sampling-variance agreement
#' additionally needs mean-square transfer or sufficient uniform integrability.
#' Generic estimated PS/PG maps, estimated weights and nuisance-refitting
#' replication are not justified by this procedure. The normal interval does
#' not depend on B. Any variance floor is inherited and reported, never added
#' silently here. A nonpositive unfloored estimate has no available interval.
#' Asymptotic validity with a floor requires a deterministic positive sequence
#' tending to zero, or a separate proof that the floor is asymptotically inactive;
#' an arbitrary fixed positive floor does not satisfy this theorem automatically.
#' @seealso \code{\link{wm_reciprocal_inference}}, \code{\link{wm_bootstrap}},
#'   \code{\link{wm_bootstrap_refit}}
#' @export
wm_reciprocal_bootstrap <- function(object, B = 999L, seed = NULL,
                                    conf.level = 0.95,
                                    interval = c("normal", "basic", "none")) {
  scalar <- function(x) is.numeric(x) && !is.complex(x) &&
    is.null(dim(x)) && length(x) == 1L && is.finite(x)
  positive_integer <- function(x) scalar(x) && x >= 1 &&
    x <= .Machine$integer.max && x == floor(x)
  if (!inherits(object, "wm_reciprocal_inference") || !is.list(object))
    stop("object must be a wm_reciprocal_inference result.", call. = FALSE)
  if (!positive_integer(B)) stop("B must be a positive integer.", call. = FALSE)
  if (!scalar(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("conf.level must lie strictly between zero and one.", call. = FALSE)
  interval <- match.arg(interval)
  if (!is.null(seed) && (!scalar(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != floor(seed)))
    stop("seed must be NULL or a nonnegative integer.", call. = FALSE)
  if (!isTRUE(object$available))
    stop("Reciprocal variance is unavailable; inspect the signed variance or explicitly choose a floor upstream.",
         call. = FALSE)
  if (!positive_integer(object$n) || !scalar(object$estimate) ||
      !scalar(object$gamma) || object$gamma <= 0 ||
      !scalar(object$used_numerator_variance) || object$used_numerator_variance <= 0 ||
      !scalar(object$root_n_variance) || object$root_n_variance <= 0 ||
      !scalar(object$variance) || object$variance <= 0)
    stop("object must contain finite positive variance and denominator fields.", call. = FALSE)
  reconstructed <- object$used_numerator_variance / object$gamma^2
  tolerance <- 100 * .Machine$double.eps *
    max(reconstructed, object$root_n_variance, .Machine$double.xmin)
  if (!is.finite(reconstructed) ||
      abs(reconstructed - object$root_n_variance) > tolerance ||
      abs(object$variance - object$root_n_variance / object$n) > tolerance / object$n)
    stop("Stored reciprocal variances disagree across scales.", call. = FALSE)
  if (!is.null(seed)) {
    if (identical(RNGkind()[2L], "Box-Muller"))
      stop("An explicit seed cannot preserve the Box-Muller cache; use seed = NULL or another normal RNG.",
           call. = FALSE)
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
  roots <- sqrt(object$root_n_variance) * stats::rnorm(B)
  draws <- object$estimate + roots / sqrt(object$n)
  mc_root <- if (B > 1L) stats::var(roots) else NA_real_
  alpha <- 1 - conf.level
  ci <- switch(interval,
    normal = object$estimate + c(-1, 1) * stats::qnorm(alpha / 2,
      lower.tail = FALSE) * sqrt(object$variance),
    basic = object$estimate - as.numeric(stats::quantile(roots,
      c(1 - alpha / 2, alpha / 2), names = FALSE)) / sqrt(object$n),
    none = NULL)
  if (any(!is.finite(c(roots, draws, ci))) ||
      (B > 1L && !is.finite(mc_root)))
    stop("Gaussian replication arithmetic exceeded numerical range.", call. = FALSE)
  if (!is.null(ci)) names(ci) <- c("lower", "upper")
  list(root_n_draws = roots, draws = draws,
       conditional_root_n_variance = object$root_n_variance,
       conditional_variance = object$variance,
       monte_carlo_root_n_variance = mc_root,
       monte_carlo_variance = mc_root / object$n,
       conf.int = ci, interval = interval, conf.level = conf.level,
       estimate = object$estimate, n = object$n, B = B, seed = seed,
       variance_floor = object$variance_floor, floor_active = object$floor_active,
       method = "Gaussian replication with reciprocal-corrected fixed-map variance",
       inference_contract = object$inference_contract)
}
