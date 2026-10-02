#' Fixed-reuse multinomial replication with optional prediction refitting
#'
#' Reproduce the original WDSM replication calculation for self-normalized
#' matching in any matching dimension, retaining distinct arm-specific graphs.
#' The supplied count matrix makes paired comparisons reproducible.
#' A \code{source_random} graph retains its explicitly unsupported
#' sampling-inference scope; its replicate intervals are empirical calculations.
#' @export
wm_bootstrap_refit <- function(object, counts, refit = NULL, conf.level = 0.95) {
  if (!inherits(object, "wm_match") || !is.list(object) ||
      !identical(object$method, "self_normalized") ||
      !object$estimand %in% c("PATE", "PATT")) {
    stop("object must be a self_normalized wm_match result.", call. = FALSE)
  }
  n <- object$n
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 2 ||
      n != floor(n) || n > .Machine$integer.max) {
    stop("object must contain a valid sample size.", call. = FALSE)
  }
  if (!is.matrix(counts) || !is.numeric(counts) || is.complex(counts) ||
      nrow(counts) != n || ncol(counts) < 2L || anyNA(counts) ||
      any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts)) || any(colSums(counts) != n)) {
    stop("counts must be an n-by-B multinomial count matrix with B >= 2 and column sums n.",
         call. = FALSE)
  }
  if (!is.null(refit) && !is.function(refit)) {
    stop("refit must be NULL or a function of one multiplicity vector.", call. = FALSE)
  }
  if (!is.numeric(conf.level) || is.complex(conf.level) ||
      length(conf.level) != 1L || !is.finite(conf.level) ||
      conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must lie strictly between zero and one.", call. = FALSE)
  }
  if (!is.null(refit) && !isTRUE(object$info$corrected)) {
    stop("Raw matching uses fixed zero predictions; omit refit.", call. = FALSE)
  }
  Y <- .wm_numeric_vector(object$data$Y, n, "object outcomes")
  Z <- .wm_numeric_vector(as.numeric(object$data$Z), n, "object treatment")
  if (!all(Z %in% c(0, 1)) || length(unique(Z)) != 2L) {
    stop("object must contain both treatment arms.", call. = FALSE)
  }
  w <- .wm_numeric_vector(object$analysis_weights, n, "object analysis weights", TRUE)
  incoming <- object$loads$incoming
  if (!is.matrix(incoming) || !is.numeric(incoming) || is.complex(incoming) ||
      !identical(dim(incoming), c(as.integer(n), 2L)) ||
      anyNA(incoming) || any(!is.finite(incoming)) || any(incoming < 0)) {
    stop("object must contain finite nonnegative n-by-2 incoming loads.", call. = FALSE)
  }
  K <- incoming[cbind(seq_len(n), Z + 1L)]
  pate <- identical(object$estimand, "PATE")
  validate_predictions <- function(predictions) {
    if (!is.list(predictions)) stop("refit must return a list of mean predictions.")
    q0 <- .wm_numeric_vector(predictions$mean0, n, "mean0")
    q1 <- if (pate) .wm_numeric_vector(predictions$mean1, n, "mean1") else numeric(n)
    list(mean0 = q0, mean1 = q1)
  }
  original <- validate_predictions(object$predictions)
  # All quantities below use the same common weight normalization as wm_match.
  # Keep K fixed: multiplying query weights and rebuilding K is a different rule.
  calculate <- function(multiplicity, prediction) {
    q0 <- prediction$mean0
    if (pate) {
      q1 <- prediction$mean1
      numerator <- sum(multiplicity * w * (q1 - q0)) +
        sum(multiplicity * (2 * Z - 1) * (w + K) *
              (Y - ifelse(Z == 1, q1, q0)))
      denominator <- sum(multiplicity * w)
    } else {
      numerator <- sum(multiplicity * Z * w * (Y - q0)) -
        sum(multiplicity * (1 - Z) * K * (Y - q0))
      denominator <- sum(multiplicity * Z * w)
    }
    result <- numerator / denominator
    if (!is.finite(denominator) || denominator <= 0 || !is.finite(result)) {
      stop("Nonfinite replicate or invalid normalization total.")
    }
    result
  }
  estimate <- .wm_numeric_vector(object$estimate, 1L, "object estimate")
  reconstructed <- calculate(rep(1, n), original)
  if (abs(reconstructed - estimate) > 1e-10 * (1 + abs(estimate))) {
    stop("Original graph loads and predictions do not reconstruct the point estimate.",
         call. = FALSE)
  }
  if (any(colSums(counts[Z == 0, , drop = FALSE]) == 0) ||
      any(colSums(counts[Z == 1, , drop = FALSE]) == 0)) {
    stop("Both arms must have positive multiplicity in every replicate.", call. = FALSE)
  }
  B <- ncol(counts)
  draws <- numeric(B)
  for (b in seq_len(B)) {
    draws[b] <- tryCatch({
      multiplicity <- counts[, b]
      prediction <- if (is.null(refit)) original else
        validate_predictions(refit(multiplicity))
      calculate(multiplicity, prediction)
    }, error = function(e) {
      stop(sprintf("%s fixed-reuse replicate %d/%d failed: %s",
                   object$estimand, b, B, conditionMessage(e)), call. = FALSE)
    })
  }
  variance <- mean((draws - mean(draws))^2)
  se <- sqrt(variance)
  conf.int <- estimate + c(-1, 1) * stats::qnorm((1 + conf.level) / 2) * se
  if (any(!is.finite(c(variance, se, conf.int)))) {
    stop("Nonfinite replication uncertainty.", call. = FALSE)
  }
  names(conf.int) <- c("lower", "upper")
  result <- list(estimate = estimate, draws = draws,
       root_n_draws = sqrt(n) * (draws - estimate),
       variance = variance, se = se, conf.int = conf.int,
       variance_divisor = "B", B = B, n = n, conf.level = conf.level,
       method = "Fixed-reuse multinomial replication",
       predictions = if (is.null(refit)) "fixed" else "refitted",
       inference_contract = paste(
         "Replicate algebra only. Sampling validity requires the applicable",
         "nuisance-refit expansion and sampling/replication variance agreement.",
         "Graph, original weights and weighted reuse remain fixed; supplied",
         "influence corrections are not added to refitted replicates.",
         "The proved compatible prediction-only route uses the same-count raw weighted-LS",
         "embedding under complete-system/actual-row conditions. Fixed B retains Monte Carlo",
         "variation. Deterministic polynomial-growing-B empirical variance validity needs",
         "the stated prediction derivative/Hessian moment above order two and simultaneous",
         "nearby-root selection, numerical-error and success control; it does not assert",
         "nonlinear conditional-second-moment convergence."))
  if (identical(object$graph$tie_rule, "Source-compatible random squared-distance boundary")) {
    result$sampling_inference_available <- FALSE
    result$inference_contract <- paste("Empirical fixed-reuse replicate algebra for",
      "source_random only. Its tie policy, positive boundary tolerance and",
      "unrestricted discrete ties have no supported sampling-inference theorem.",
      "Graph, supplied weights and original weighted reuse remain fixed;",
      "replicate variance and intervals are descriptive calculations.")
  }
  result
}
