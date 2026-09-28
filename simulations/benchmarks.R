# Deterministic mark enumeration and corrected-estimator variance benchmarks.
# Source this file from the simulation runner. No random numbers are generated.

.wm_bm_scalar <- function(x, name, lower = -Inf, upper = Inf,
                           strict_lower = FALSE, strict_upper = FALSE) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
      (if (strict_lower) x <= lower else x < lower) ||
      (if (strict_upper) x >= upper else x > upper)) {
    stop(name, " is outside its permitted finite scalar range", call. = FALSE)
  }
  as.numeric(x)
}

.wm_bm_integer <- function(M) {
  M <- .wm_bm_scalar(M, "M", 1, .Machine$integer.max)
  if (M != floor(M)) stop("M must be a positive integer", call. = FALSE)
  as.integer(M)
}

.wm_bm_pair <- function(x, name, positive = FALSE, nonnegative = FALSE) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || length(x) != 2L ||
      any(!is.finite(x)) || (positive && any(x <= 0)) ||
      (nonnegative && any(x < 0))) {
    stop(name, " must contain two finite values in the permitted range", call. = FALSE)
  }
  as.numeric(x)
}

# A and B condition on the root weight; C[,ell+1] conditions on ell shared
# nonroot donors. Private donor sets are independent conditional on that sum.
wm_two_point_mark_moments <- function(M, weights, prob_high,
                                      max_terms = 2e6) {
  M <- .wm_bm_integer(M)
  weights <- .wm_bm_pair(weights, "weights", positive = TRUE)
  p <- .wm_bm_scalar(prob_high, "prob_high", 0, 1)
  max_terms <- .wm_bm_scalar(max_terms, "max_terms", 1)
  terms <- as.double(M) * (M + 1) * (M + 2) / 3
  if (terms > max_terms) {
    stop("Discrete mark enumeration exceeds max_terms; increase the explicit guard",
         call. = FALSE)
  }
  # Fractions are invariant to this donor-law-only normalization.
  w <- weights / max(weights)
  if (any(w <= 0)) stop("Donor weight normalization underflowed", call. = FALSE)
  counts <- 0:(M - 1L)
  probabilities <- stats::dbinom(counts, M - 1L, p)
  other_sum <- (M - 1L - counts) * w[1L] + counts * w[2L]
  A <- B <- numeric(2L)
  C <- matrix(0, nrow = 2L, ncol = M)
  for (root in 1:2) {
    share <- w[root] / (w[root] + other_sum)
    A[root] <- sum(probabilities * share)
    B[root] <- sum(probabilities * share^2)
    for (ell in 0:(M - 1L)) {
      shared_counts <- 0:ell
      shared_probability <- stats::dbinom(shared_counts, ell, p)
      shared_sum <- (ell - shared_counts) * w[1L] + shared_counts * w[2L]
      private_n <- M - 1L - ell
      private_counts <- 0:private_n
      private_probability <- stats::dbinom(private_counts, private_n, p)
      private_sum <- (private_n - private_counts) * w[1L] + private_counts * w[2L]
      conditional <- vapply(shared_sum, function(s) {
        sum(private_probability * w[root] / (w[root] + s + private_sum))
      }, numeric(1L))
      C[root, ell + 1L] <- sum(shared_probability * conditional^2)
    }
  }
  colnames(C) <- paste0("beta_", 0:(M - 1L))
  list(weights = weights, probabilities = c(1 - p, p),
       A = A, B = B, C = C, overlap = 0:(M - 1L),
       precision = "Finite binomial enumeration, floating-point arithmetic; no Monte Carlo")
}

wm_beta_1d <- function(M) {
  M <- .wm_bm_integer(M)
  if (M > 1e6) stop("Requested beta vector exceeds the allocation guard", call. = FALSE)
  beta <- 2 * seq_len(M)
  beta[M] <- 3 * M / 2
  names(beta) <- paste0("beta_", 0:(M - 1L))
  beta
}

.wm_bm_geometry <- function(M, beta, beta_cov) {
  if (!is.numeric(beta) || is.complex(beta) || !is.null(dim(beta)) ||
      length(beta) != M || any(!is.finite(beta)) || any(beta < 0)) {
    stop("beta must be a nonnegative finite vector of length M", call. = FALSE)
  }
  beta <- as.numeric(beta)
  names(beta) <- paste0("beta_", 0:(M - 1L))
  if (!is.null(beta_cov)) {
    if (!is.matrix(beta_cov) || !is.numeric(beta_cov) || is.complex(beta_cov) ||
        !identical(dim(beta_cov), c(M, M)) || any(!is.finite(beta_cov))) {
      stop("beta_cov must be a finite M-by-M covariance matrix", call. = FALSE)
    }
    scale <- max(abs(beta_cov), .Machine$double.xmin)
    tolerance <- 1e-10 * scale
    if (max(abs(beta_cov - t(beta_cov))) > tolerance) {
      stop("beta_cov must be symmetric", call. = FALSE)
    }
    beta_cov <- beta_cov / 2 + t(beta_cov) / 2
    eigenvalues <- eigen(beta_cov, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(eigenvalues)) || min(eigenvalues) < -tolerance) {
      stop("beta_cov must be positive semidefinite", call. = FALSE)
    }
  }
  if (!is.finite(sum(beta))) stop("The sum of beta exceeds numerical range", call. = FALSE)
  list(beta = beta, beta_cov = beta_cov, alpha = sum(beta),
       alpha_jensen_bound_satisfied = sum(beta) >= as.double(M)^2)
}

.wm_bm_reciprocal <- function(M, q) {
  # Nonnegative log-scale summands avoid factorial overflow.
  value <- 0
  for (a in 0:(M - 1L)) for (b in 0:(M - 1L)) {
    value <- value + exp(lchoose(a + b, a) + a * log1p(-q) + b * log(q))
  }
  q * (1 - q) * value
}

.wm_bm_finish <- function(M, geometry, table, coefficients, marks,
                          assumptions, diagnostics = list()) {
  colnames(coefficients) <- names(geometry$beta)
  rownames(coefficients) <- paste(table$method, table$estimand, sep = ":")
  established <- table$variance_scope != "not_established"
  if (any(!is.finite(as.matrix(table[c("target", "probability_limit", "bias")])))) {
    stop("Target or probability-limit arithmetic exceeded numerical range", call. = FALSE)
  }
  if (any(!is.finite(table$variance_intercept[established])) ||
      any(!is.finite(coefficients[established, , drop = FALSE]))) {
    stop("Variance benchmark arithmetic exceeded numerical range", call. = FALSE)
  }
  table$root_n_variance <- NA_real_
  table$geometry_mcse <- NA_real_
  for (i in seq_len(nrow(table))) {
    coefficient <- coefficients[i, ]
    if (anyNA(coefficient) || is.na(table$variance_intercept[i])) next
    table$root_n_variance[i] <- table$variance_intercept[i] +
      sum(coefficient * geometry$beta)
    if (!is.null(geometry$beta_cov)) {
      uncertainty <- as.numeric(crossprod(coefficient, geometry$beta_cov %*% coefficient))
      if (!is.finite(uncertainty)) {
        stop("Geometric uncertainty arithmetic exceeded numerical range", call. = FALSE)
      }
      table$geometry_mcse[i] <- sqrt(max(0, uncertainty))
    } else if (all(coefficient == 0)) {
      table$geometry_mcse[i] <- 0
    }
  }
  if (any(!is.finite(table$root_n_variance[established])) ||
      any(table$root_n_variance[established] < 0)) {
    stop("The supplied marks and geometry produced an invalid variance benchmark",
         call. = FALSE)
  }
  list(table = table, beta_coefficients = coefficients, beta = geometry$beta,
       beta_cov = geometry$beta_cov, alpha = geometry$alpha,
       alpha_jensen_bound_satisfied = geometry$alpha_jensen_bound_satisfied,
       mark_moments = marks, M = M, assumptions = assumptions,
       diagnostics = diagnostics,
       precision = if (is.null(geometry$beta_cov))
         "Mark enumeration is deterministic; geometric uncertainty was not supplied"
       else "Geometry MCSE propagates the supplied covariance of beta estimates")
}

# S is independent of treatment and weight/error marks. The effect depends
# only on S, has the supplied mean/variance, and each error is centered given
# its arm and weight. These conditions eliminate heterogeneity/residual cross
# terms; they are not a formula for general covariate-dependent mark laws.
wm_benchmark_strong <- function(M, beta, beta_cov = NULL, q = 0.4,
                                 weight0 = c(0.5, 1.5), weight1 = c(1, 3),
                                 prob0 = 0.5, prob1 = 0.5,
                                 variance0 = c(0.25, 1), variance1 = c(0.5, 1.5),
                                 effect_variance = 1 / 12, effect_mean = 1.5) {
  M <- .wm_bm_integer(M)
  q <- .wm_bm_scalar(q, "q", 0, 1, TRUE, TRUE)
  effect_variance <- .wm_bm_scalar(effect_variance, "effect_variance", 0)
  effect_mean <- .wm_bm_scalar(effect_mean, "effect_mean")
  geometry <- .wm_bm_geometry(M, beta, beta_cov)
  marks <- list(wm_two_point_mark_moments(M, weight0, prob0),
                wm_two_point_mark_moments(M, weight1, prob1))
  variances <- list(.wm_bm_pair(variance0, "variance0", nonnegative = TRUE),
                    .wm_bm_pair(variance1, "variance1", nonnegative = TRUE))
  arm_probability <- c(1 - q, q)
  rho <- nu <- direct <- cross <- diagonal <- numeric(2L)
  pair <- matrix(0, 2L, M)
  for (z in 1:2) {
    a <- marks[[z]]
    rho[z] <- sum(a$probabilities * a$weights)
    nu[z] <- sum(a$probabilities * a$weights^2)
    weighted_variance <- a$probabilities * variances[[z]]
    direct[z] <- sum(weighted_variance * a$weights^2)
    cross[z] <- sum(weighted_variance * a$weights * a$A)
    diagonal[z] <- sum(weighted_variance * a$B)
    pair[z, ] <- colSums(a$C * weighted_variance)
  }
  gamma <- sum(arm_probability * rho)
  gamma_t <- q * rho[2L]
  heterogeneity <- sum(arm_probability * nu) * effect_variance
  t <- direct / rho^2
  sn_intercept <- heterogeneity
  st_intercept <- heterogeneity
  sn_coefficient <- st_coefficient <- numeric(M)
  for (z in 1:2) {
    other <- 3L - z
    r <- arm_probability[other] / arm_probability[z]
    sn_intercept <- sn_intercept + arm_probability[z] *
      (direct[z] + 2 * r * rho[other] * M * cross[z] + r * nu[other] * M * diagonal[z])
    sn_coefficient <- sn_coefficient + arm_probability[z] * r^2 * rho[other]^2 * pair[z, ]
    st_intercept <- st_intercept + arm_probability[z] *
      (direct[z] + 2 * r * rho[other] * rho[z] * t[z] + t[z] * r * nu[other] / M)
    st_coefficient <- st_coefficient + arm_probability[z] * t[z] * r^2 * rho[other]^2 / M^2
  }
  r0 <- q / (1 - q)
  treated_direct <- q * (nu[2L] * effect_variance + direct[2L])
  sn_t_intercept <- treated_direct + (1 - q) * r0 * nu[2L] * M * diagonal[1L]
  sn_t_coefficient <- (1 - q) * r0^2 * rho[2L]^2 * pair[1L, ]
  st_t_intercept <- treated_direct + (1 - q) * t[1L] * r0 * nu[2L] / M
  st_t_coefficient <- rep((1 - q) * t[1L] * r0^2 * rho[2L]^2 / M^2, M)
  table <- data.frame(method = rep(c("self_normalized", "stabilized"), each = 2L),
                      estimand = rep(c("PATE", "PATT"), 2L),
                      target = effect_mean, probability_limit = effect_mean, bias = 0,
                      variance_scope = "bias_corrected",
                      variance_intercept = c(sn_intercept / gamma^2,
                        sn_t_intercept / gamma_t^2, st_intercept / gamma^2,
                        st_t_intercept / gamma_t^2), stringsAsFactors = FALSE)
  coefficients <- rbind(sn_coefficient / gamma^2, sn_t_coefficient / gamma_t^2,
                         st_coefficient / gamma^2, st_t_coefficient / gamma_t^2)
  .wm_bm_finish(M, geometry, table, coefficients, marks,
    assumptions = c("S independent of arm/weight/error marks; effect depends only on S",
      "Errors conditionally centered given arm and individual weight",
      "Variance benchmarks concern oracle or valid feasible corrected estimators, not raw stabilized imputation"),
    diagnostics = list(q = q, rho = rho, nu = nu, gamma = gamma,
                        conditional_eta_second_moment = t,
                        heterogeneity_numerator = heterogeneity))
}

.wm_bm_donor_mean <- function(M, weights, prob_high, outcome) {
  k <- 0:M
  w <- weights / max(weights)
  high_share <- k * w[2L] / ((M - k) * w[1L] + k * w[2L])
  sum(stats::dbinom(k, M, prob_high) *
        ((1 - high_share) * outcome[1L] + high_share * outcome[2L]))
}

# Constant weighted arm means, independent uniform matching coordinates.
# For nondefault parameters the target is the observed weighted-mean contrast;
# causal identification must be checked separately. Delta is identically zero
# after subtracting that constant contrast, so no delta*kappa term is dropped.
wm_benchmark_weak <- function(M, beta, beta_cov = NULL, q = 8 / 17,
                               weight0 = c(1, 2), weight1 = c(1, 3),
                               prob0 = 1 / 3, prob1 = 1 / 4,
                               outcome0 = c(0, 1), outcome1 = c(0, 1)) {
  M <- .wm_bm_integer(M)
  q <- .wm_bm_scalar(q, "q", 0, 1, TRUE, TRUE)
  geometry <- .wm_bm_geometry(M, beta, beta_cov)
  marks <- list(wm_two_point_mark_moments(M, weight0, prob0),
                wm_two_point_mark_moments(M, weight1, prob1))
  outcomes <- list(.wm_bm_pair(outcome0, "outcome0"),
                   .wm_bm_pair(outcome1, "outcome1"))
  probabilities <- c(prob0, prob1)
  arm_probability <- c(1 - q, q)
  rho <- nu <- m <- t <- kappa <- donor_mean <- numeric(2L)
  for (z in 1:2) {
    a <- marks[[z]]
    rho[z] <- sum(a$probabilities * a$weights)
    nu[z] <- sum(a$probabilities * a$weights^2)
    m[z] <- sum(a$probabilities * a$weights * outcomes[[z]]) / rho[z]
    eta <- a$weights * (outcomes[[z]] - m[z]) / rho[z]
    t[z] <- sum(a$probabilities * eta^2)
    kappa[z] <- sum(a$probabilities * a$weights * eta)
    donor_mean[z] <- .wm_bm_donor_mean(M, a$weights, probabilities[z], outcomes[[z]])
  }
  gamma <- sum(arm_probability * rho)
  gamma_t <- q * rho[2L]
  target <- m[2L] - m[1L]
  pate_limit <- (q * rho[2L] * (m[2L] - donor_mean[1L]) +
                  (1 - q) * rho[1L] * (donor_mean[2L] - m[1L])) / gamma
  patt_limit <- m[2L] - donor_mean[1L]
  st_intercept <- sum(arm_probability * rho^2 * t)
  st_coefficient <- numeric(M)
  for (z in 1:2) {
    other <- 3L - z
    r <- arm_probability[other] / arm_probability[z]
    st_intercept <- st_intercept + arm_probability[z] *
      (2 * r * rho[other] * rho[z] * t[z] + t[z] * r * nu[other] / M)
    st_coefficient <- st_coefficient + arm_probability[z] * t[z] * r^2 * rho[other]^2 / M^2
  }
  reciprocal <- -2 * .wm_bm_reciprocal(M, q) * kappa[1L] * kappa[2L] / M^2
  diagonal_intercept <- st_intercept / gamma^2
  st_intercept <- st_intercept + reciprocal
  r0 <- q / (1 - q)
  st_t_intercept <- q * rho[2L]^2 * t[2L] + (1 - q) * t[1L] * r0 * nu[2L] / M
  st_t_coefficient <- rep((1 - q) * t[1L] * r0^2 * rho[2L]^2 / M^2, M)
  table <- data.frame(method = rep(c("self_normalized", "stabilized"), each = 2L),
                      estimand = rep(c("PATE", "PATT"), 2L), target = target,
                      probability_limit = c(pate_limit, patt_limit, target, target),
                      bias = c(pate_limit - target, patt_limit - target, 0, 0),
                      variance_scope = c(rep("not_established", 2L), rep("bias_corrected", 2L)),
                      variance_intercept = c(NA_real_, NA_real_, st_intercept / gamma^2,
                                             st_t_intercept / gamma_t^2),
                      stringsAsFactors = FALSE)
  coefficients <- rbind(rep(NA_real_, M), rep(NA_real_, M),
                         st_coefficient / gamma^2, st_t_coefficient / gamma_t^2)
  .wm_bm_finish(M, geometry, table, coefficients, marks,
    assumptions = c("Matching coordinates independent of arm and two-point outcome/weight marks",
      "Constant weighted arm means; target is their contrast (default design identifies causal zero)",
      "No original-rule variance assertion under failed strong centering"),
    diagnostics = list(q = q, rho = rho, nu = nu, mean = m, gamma = gamma,
      conditional_eta_second_moment = t, kappa = kappa,
      donor_weighted_means = donor_mean,
      stabilized_pate_actual_row_second_moment = diagonal_intercept +
        sum(st_coefficient * geometry$beta) / gamma^2,
      stabilized_pate_reciprocal_correction = reciprocal / gamma^2))
}
