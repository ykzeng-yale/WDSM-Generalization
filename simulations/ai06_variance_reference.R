# Source-equation comparator: Abadie and Imbens (2006), published Section 4.
# Unit weights and one common supplied matching map only. No fitted-score
# adjustment, survey-weight extension, or package-inference claim is implicit.
ai06_variance_reference <- function(Y, Z, scores, M = 1L, J = 1L,
                                    estimand = c("PATE", "PATT")) {
  estimand <- match.arg(estimand)
  vector_ok <- function(x, n) is.numeric(x) && !is.complex(x) &&
    is.null(dim(x)) && length(x) == n && all(is.finite(x))
  if (!vector_ok(Y, length(Y)) || length(Y) < 4L) {
    stop("Y must be a finite numeric vector with at least four rows.")
  }
  n <- length(Y)
  if ((!is.numeric(Z) && !is.logical(Z)) || is.complex(Z) ||
      !is.null(dim(Z)) || length(Z) != n || anyNA(Z) ||
      any(!Z %in% c(0, 1)) || length(unique(Z)) != 2L) {
    stop("Z must contain both binary treatment arms and have length n.")
  }
  Z <- as.integer(Z)
  if (!is.matrix(scores) || !is.numeric(scores) || is.complex(scores) ||
      nrow(scores) != n || ncol(scores) < 1L || any(!is.finite(scores))) {
    stop("scores must be a finite numeric matrix with n rows.")
  }
  for (count in list(M = M, J = J)) {
    if (!vector_ok(count, 1L) || count < 1 || count != floor(count) || count > n) {
      stop("M and J must be positive integer counts no larger than n.")
    }
  }
  M <- as.integer(M)
  J <- as.integer(J)
  arms <- list(which(Z == 0L), which(Z == 1L))
  if (any(lengths(arms) <= J)) stop("Each arm needs at least J + 1 observations.")
  donor_arms <- if (estimand == "PATE") 1:2 else 1L
  if (any(lengths(arms)[donor_arms] < M)) stop("Each used donor arm needs M observations.")

  # Independent direct ordering; no wm_match graph or distance helper is reused.
  nearest <- function(i, candidates, count) {
    delta <- sweep(scores[candidates, , drop = FALSE], 2L, scores[i, ], "-")
    if (any(!is.finite(delta))) stop("Coordinate difference overflow; rescale scores.")
    # Per-candidate scaled norms avoid squaring tiny or huge coordinates.
    scale <- apply(abs(delta), 1L, max)
    distance <- scale * sqrt(rowSums((delta / ifelse(scale > 0, scale, 1))^2))
    if (any(!is.finite(distance))) stop("Distance arithmetic overflow; rescale scores.")
    candidates[order(distance, candidates, method = "radix")[seq_len(count)]]
  }
  auxiliary <- matrix(NA_integer_, nrow = n, ncol = J)
  sigma2 <- numeric(n)
  for (i in seq_len(n)) {
    candidates <- arms[[Z[i] + 1L]]
    auxiliary[i, ] <- nearest(i, candidates[candidates != i], J)
    sigma2[i] <- (J / (J + 1)) * (Y[i] - mean(Y[auxiliary[i, ]]))^2
  }
  queries <- if (estimand == "PATE") seq_len(n) else arms[[2L]]
  donors <- matrix(NA_integer_, nrow = length(queries), ncol = M)
  contrast <- numeric(length(queries))
  for (q in seq_along(queries)) {
    i <- queries[q]
    donors[q, ] <- nearest(i, arms[[2L - Z[i]]], M)
    contrast[q] <- (2L * Z[i] - 1L) * (Y[i] - mean(Y[donors[q, ]]))
  }
  reuse <- tabulate(as.vector(donors), nbins = n)
  estimate <- mean(contrast)
  target_count <- length(queries)
  if (estimand == "PATE") {
    conditional_sum <- sum((1 + reuse / M)^2 * sigma2)
    marginal_correction <- sum(((reuse / M)^2 +
      ((2 * M - 1) / M) * (reuse / M)) * sigma2)
  } else {
    conditional_sum <- sum((Z - (1 - Z) * reuse / M)^2 * sigma2)
    marginal_correction <- sum((1 - Z) * reuse * (reuse - 1) / M^2 * sigma2)
  }
  # Contrast squares already contain outcome noise; the correction above
  # replaces that noise with the complete residual component, not its sum.
  contrast_sum <- sum((contrast - estimate)^2)
  conditional_variance <- conditional_sum / target_count^2
  marginal_variance <- (contrast_sum + marginal_correction) / target_count^2
  quantities <- c(sigma2, estimate, conditional_variance, marginal_variance,
                  n * conditional_variance, n * marginal_variance)
  if (any(!is.finite(quantities))) stop("Outcome/variance arithmetic overflow; rescale outcomes.")
  list(estimate = estimate, estimand = estimand, n = n,
       target_count = target_count, M = M, J = J, dimension = ncol(scores),
       query = queries, donors = donors, reuse_count = reuse,
       auxiliary_donors = auxiliary, conditional_outcome_variance = sigma2,
       raw_contrast = contrast,
       conditional_variance = conditional_variance,
       marginal_variance = marginal_variance,
       conditional_root_n_variance = n * conditional_variance,
       marginal_root_n_variance = n * marginal_variance,
       conditional_root_target_count_variance = target_count * conditional_variance,
       marginal_root_target_count_variance = target_count * marginal_variance,
       conditional_se = sqrt(conditional_variance), marginal_se = sqrt(marginal_variance),
       contrast_component = contrast_sum / target_count^2,
       marginal_noise_correction = marginal_correction / target_count^2,
       contract = paste("AI06 common fixed map, unit weights, matching with replacement;",
         "M opposite-arm matches, J other same-arm variance neighbors;",
         "Euclidean supplied coordinates, exact-distance then original-row ties.",
         "Conditional residual and population-total variances have distinct targets.",
         "Validity requires the source assumptions and appropriate matching-bias control;",
         "no fitted-score adjustment or weighted extension is supplied."))
}
