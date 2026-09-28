#' Gaussian Contribution Bootstrap for Weighted Matching
#'
#' Draw Gaussian multipliers on the normalized row, unordered-edge, and
#' optional independent-training contributions of a fitted \code{wm_match}
#' object. Matching graphs,
#' predictions, weights, and nuisance influence corrections remain fixed.
#'
#' @param object A \code{wm_match} result with variance estimation enabled.
#' @param B Positive integer number of multiplier draws.
#' @param seed Optional nonnegative integer seed. An explicit seed restores
#'   the caller's \code{.Random.seed}, including its prior absence. The
#'   Box--Muller normal generator does not support this option because its
#'   hidden cached state cannot be restored; use \code{seed = NULL} with that
#'   generator or select another normal generator before this call.
#' @param conf.level Confidence level strictly between zero and one.
#' @param interval \code{"normal"} uses the exact conditional Gaussian
#'   variance; \code{"basic"} uses simulated root-n quantiles; \code{"none"}
#'   omits the interval.
#' @param chunk_size Positive integer maximum number of multiplier values
#'   generated together. No contributions-by-B matrix is allocated.
#' @return A list with \code{root_n_draws}, centered multiplier draws on the
#'   root-n scale; \code{draws}, estimates shifted by those draws divided by
#'   \code{sqrt(n)}; exact \code{conditional_root_n_variance} and
#'   \code{conditional_variance}; Monte Carlo sample variances
#'   \code{monte_carlo_root_n_variance} and \code{monte_carlo_variance}; and
#'   \code{conf.int}, \code{interval}, \code{conf.level}, \code{B}, \code{n},
#'   \code{seed}, \code{estimate}, and \code{method}. Monte Carlo variances
#'   are \code{NA} when \code{B = 1}.
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
#' @seealso \code{\link{wm_match}}
#' @export
wm_bootstrap <- function(object, B = 999L, seed = NULL, conf.level = 0.95,
                         interval = c("normal", "basic", "none"),
                         chunk_size = 65536L) {
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
