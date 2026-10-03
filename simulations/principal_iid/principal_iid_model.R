# Source-derived iid mixture model, exact selection sampler and target calculation.
# Sourcing this file defines functions only. See README.md for the scientific recipe.

.principal_iid_with_seed <- function(seed, fun) {
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed))
    stop("seed must be one nonnegative R integer.")
  old_kind <- RNGkind()
  if (identical(old_kind[[2L]], "Box-Muller"))
    stop("Box-Muller caller RNG has hidden cached normals; use a separate Inversion session.")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))
  fun()
}

principal_iid_model <- function() {
  prefix <- .principal_iid_with_seed(20220930L, function() {
    # Do not shorten these to six columns: all 70 draws precede all 140 draws.
    strata <- matrix(rnorm(10L * 7L, mean = 0, sd = 0.35),
                     nrow = 10L, ncol = 7L)
    components <- matrix(rnorm(20L * 7L, mean = 0, sd = 0.15),
                         nrow = 20L, ncol = 7L)
    list(strata = strata, components = components,
         rng_end = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  })
  labels <- data.frame(H = seq_len(200L),
    stratum = rep(seq_len(10L), each = 20L),
    component = rep(seq_len(20L), times = 10L))
  means <- prefix$strata[labels$stratum, seq_len(6L), drop = FALSE] +
    prefix$components[labels$component, seq_len(6L), drop = FALSE]
  colnames(means) <- paste0("X", seq_len(6L))
  rownames(means) <- as.character(labels$H)
  values <- list(
    schema = "principal_iid_model_v1", labels = labels, means = means,
    component_probability = rep(1 / 200, 200L),
    strata_offsets_7 = prefix$strata,
    component_offsets_7 = prefix$components,
    prefix_seed = 20220930L, prefix_rng_end = prefix$rng_end,
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    reference_R = "4.4.2",
    treatment = data.frame(overlap = c("GoodOverlap", "PoorOverlap"),
      intercept = c(log(35 / 80), log(20 / 80)), kappa = c(0.6, 2)),
    treatment_coefficient = log(c(1.1, 1.25, 1.5, 1.75, 2, 2.5)),
    treatment_interaction = log(1.1),
    selection_coefficient = log(c(1.05, 1.1, 1.15, 1.1, 1.05, 1.1)),
    selection_intercept = log(0.005 / 0.995),
    selection_treatment = log(0.9),
    outcome_coefficient = c(2.5, -2, 1.75, -1.25, 1.5, 1.1),
    source_sha256 = c(
      GoodOverlap = "8ee8accb56911964c535d49a409607c14132fc685857887f1a07304fc7613955",
      PoorOverlap = "beccd37262a80670a7ba5dbe3bde91f13f8bfbfcd294b366ce2a2cdaef7dee81",
      sampler = "a955b73849b430d62e51bd67886f18e27c79486a52908aaa01821ccf2b61a906"))
  out <- list2env(values, parent = emptyenv())
  class(out) <- "principal_iid_model_v1"
  lockEnvironment(out, bindings = TRUE)
  out
}

.principal_iid_check_model <- function(model) {
  if (!inherits(model, "principal_iid_model_v1") || !is.environment(model) ||
      !environmentIsLocked(model) ||
      !identical(model$schema, "principal_iid_model_v1") ||
      !identical(dim(model$means), c(200L, 6L)) ||
      any(!is.finite(model$means)))
    stop("Supply the locked principal_iid_model() object or its saved RDS.")
  invisible(model)
}

.principal_iid_X <- function(X) {
  if (!is.matrix(X) || !is.numeric(X) || ncol(X) != 6L ||
      any(!is.finite(X))) stop("X must be a finite numeric six-column matrix.")
  X
}

.principal_iid_overlap <- function(model, overlap) {
  .principal_iid_check_model(model)
  if (length(overlap) != 1L || !is.character(overlap) || is.na(overlap) ||
      !overlap %in% model$treatment$overlap)
    stop("overlap must be GoodOverlap or PoorOverlap.")
  model$treatment[match(overlap, model$treatment$overlap), , drop = FALSE]
}

principal_iid_known_mean <- function(model, X) {
  .principal_iid_check_model(model)
  X <- .principal_iid_X(X)
  linear <- drop(X %*% model$outcome_coefficient)
  interaction <- X[, 1L] * X[, 2L]
  # Preserve the original outcome construction's algebra and shared E0.
  mu0 <- 0.3 * (linear + 2.5 * interaction)
  delta <- 1 + 0.2 * (linear + 1.5 * interaction)
  list(mu0 = mu0, mu1 = mu0 + delta, delta = delta)
}

principal_iid_delta <- function(model, X) {
  principal_iid_known_mean(model, X)$delta
}

principal_iid_mu <- function(model, X, z) {
  means <- principal_iid_known_mean(model, X)
  if (!is.numeric(z) || anyNA(z) || !all(z %in% c(0, 1)) ||
      !length(z) %in% c(1L, length(means$mu0)))
    stop("z must be 0/1, of length one or nrow(X).")
  if (length(z) == 1L) return(if (z == 1) means$mu1 else means$mu0)
  ifelse(z == 1, means$mu1, means$mu0)
}

principal_iid_propensity <- function(model, X, overlap) {
  setting <- .principal_iid_overlap(model, overlap)
  X <- .principal_iid_X(X)
  lp <- setting$intercept + setting$kappa *
    (drop(X %*% model$treatment_coefficient) +
     model$treatment_interaction * X[, 1L] * X[, 2L])
  plogis(lp)
}

principal_iid_selection <- function(model, X, Z) {
  .principal_iid_check_model(model)
  X <- .principal_iid_X(X)
  if (!is.numeric(Z) || length(Z) != nrow(X) || anyNA(Z) ||
      !all(Z %in% c(0, 1))) stop("Z must contain one 0/1 per row.")
  lp <- model$selection_intercept +
    drop(X %*% model$selection_coefficient) + model$selection_treatment * Z
  pi <- plogis(lp)
  W <- 1 / pi
  if (any(!is.finite(W)) || any(pi <= 0) || any(pi > 1))
    stop("Nonfinite supplied weight: report the sample failure; do not clip/redraw.")
  list(lp = lp, pi = pi, W = W)
}

principal_iid_sample <- function(model, overlap, n, seed, batch_size = 4096L) {
  .principal_iid_overlap(model, overlap)
  for (value in list(n, batch_size)) {
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
        value < 1 || value > .Machine$integer.max || value != floor(value))
      stop("n and batch_size must be positive R integers.")
  }
  n <- as.integer(n)
  batch_size <- as.integer(batch_size)
  .principal_iid_with_seed(seed, function() {
    rng_start <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    a <- model$selection_coefficient
    log_h <- drop(model$means %*% a)
    h_probability <- exp(log_h - max(log_h))
    h_probability <- h_probability / sum(h_probability)
    X <- matrix(NA_real_, nrow = n, ncol = 6L,
                dimnames = list(NULL, paste0("X", seq_len(6L))))
    H <- Z <- integer(n)
    filled <- 0L
    proposals <- all_accepted <- batches <- 0
    while (filled < n) {
      proposed_h <- sample.int(200L, size = batch_size, replace = TRUE,
                               prob = h_probability)
      proposed_x <- matrix(rnorm(batch_size * 6L), nrow = batch_size,
                           ncol = 6L) +
        model$means[proposed_h, , drop = FALSE] +
        matrix(a, nrow = batch_size, ncol = 6L, byrow = TRUE)
      proposed_z <- rbinom(batch_size, size = 1L,
        prob = principal_iid_propensity(model, proposed_x, overlap))
      lp_pi <- model$selection_intercept + drop(proposed_x %*% a) +
        model$selection_treatment * proposed_z
      acceptance <- exp(model$selection_treatment * proposed_z) * plogis(-lp_pi)
      if (any(!is.finite(acceptance)) || any(acceptance < 0 | acceptance > 1))
        stop("Invalid tilted acceptance probability; no sample substitution.")
      accepted <- which(runif(batch_size) < acceptance)
      take <- min(length(accepted), n - filled)
      if (take > 0L) {
        source_rows <- accepted[seq_len(take)]
        target_rows <- filled + seq_len(take)
        X[target_rows, ] <- proposed_x[source_rows, , drop = FALSE]
        H[target_rows] <- proposed_h[source_rows]
        Z[target_rows] <- proposed_z[source_rows]
        filled <- filled + take
      }
      proposals <- proposals + batch_size
      all_accepted <- all_accepted + length(accepted)
      batches <- batches + 1
    }
    means <- principal_iid_known_mean(model, X)
    E0 <- rnorm(n)
    E1 <- rnorm(n)
    Y0 <- means$mu0 + E0
    Y1 <- Y0 + means$delta + E1
    selection <- principal_iid_selection(model, X, Z)
    data <- data.frame(X, Z = as.integer(Z), W = selection$W,
      Y = ifelse(Z == 1L, Y1, Y0), Y0 = Y0, Y1 = Y1,
      mu0 = means$mu0, mu1 = means$mu1, delta = means$delta,
      e = principal_iid_propensity(model, X, overlap), pi = selection$pi,
      H = H, E0 = E0, E1 = E1, check.names = FALSE)
    if (any(!is.finite(as.matrix(data))))
      stop("Nonfinite generated row: retain failure, without clipping or redraw.")
    list(schema = "principal_iid_sample_v1", overlap = overlap, n = n,
      seed = as.integer(seed), rng_kind = model$rng_kind,
      prefix_seed = model$prefix_seed, source_sha256 = model$source_sha256,
      data = data, rng_start = rng_start,
      rng_end = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
      diagnostics = list(sampler = "exact_normal_tilt_rejection",
        batch_size = batch_size, proposals = proposals, batches = batches,
        accepted_in_generated_batches = all_accepted, accepted_used = n,
        unused_final_batch_successes = all_accepted - n,
        supplied_weight = "1/pi_Z(X)", first_n_accepted = TRUE))
  })
}

principal_iid_pate <- function(model) {
  .principal_iid_check_model(model)
  1 + 0.2 * sum(model$outcome_coefficient * colMeans(model$means)) +
    0.3 * mean(model$means[, 1L] * model$means[, 2L])
}

.principal_iid_gh_normal <- function(q) {
  # Probabilists' orthonormal Hermite Jacobi matrix: diagonal 0,
  # off diagonal sqrt(j). Eigenvalues/first eigenvector-row squared
  # integrate the standard-normal probability law directly.
  jacobi <- matrix(0, q, q)
  off <- sqrt(seq_len(q - 1L))
  jacobi[cbind(seq_len(q - 1L), seq.int(2L, q))] <- off
  jacobi[cbind(seq.int(2L, q), seq_len(q - 1L))] <- off
  eig <- eigen(jacobi, symmetric = TRUE)
  index <- order(eig$values)
  nodes <- eig$values[index]
  weights <- eig$vectors[1L, index]^2
  weights <- weights / sum(weights)
  if (any(!is.finite(nodes)) || any(!is.finite(weights)) || any(weights < 0))
    stop("Nonfinite Gauss-Hermite rule.")
  list(nodes = nodes, weights = weights)
}

.principal_iid_patt_order <- function(model, overlap, rule, order) {
  setting <- .principal_iid_overlap(model, overlap)
  t <- model$treatment_coefficient
  b <- model$outcome_coefficient
  kappa <- setting$kappa
  sd_T <- kappa * sqrt(sum(t[3:6]^2))
  stein <- 0.2 * kappa * sum(b[3:6] * t[3:6])
  weights_12 <- rep(rule$weights, each = order) *
    rep(rule$weights, times = order)
  component <- data.frame(overlap = overlap, order = order, H = seq_len(200L),
                           numerator = NA_real_, denominator = NA_real_)
  for (h in seq_len(200L)) {
    m <- model$means[h, ]
    x1 <- rep(m[1L] + rule$nodes, each = order)
    x2 <- rep(m[2L] + rule$nodes, times = order)
    lp_12 <- setting$intercept + kappa *
      (t[1L] * x1 + t[2L] * x2 + model$treatment_interaction * x1 * x2)
    T_nodes <- kappa * sum(t[3:6] * m[3:6]) + sd_T * rule$nodes
    probability <- plogis(outer(lp_12, T_nodes, "+"))
    conditional_e <- drop(probability %*% rule$weights)
    conditional_derivative <- drop((probability * (1 - probability)) %*%
                                     rule$weights)
    delta_h <- 1 + 0.2 * (b[1L] * x1 + b[2L] * x2 +
                          sum(b[3:6] * m[3:6])) + 0.3 * x1 * x2
    component$numerator[h] <- sum(weights_12 *
      (delta_h * conditional_e + stein * conditional_derivative))
    component$denominator[h] <- sum(weights_12 * conditional_e)
  }
  numerator <- mean(component$numerator)
  denominator <- mean(component$denominator)
  value <- numerator / denominator
  list(component = component,
       row = data.frame(overlap = overlap, order = order,
         numerator = numerator, denominator = denominator, PATT = value,
         finite = all(is.finite(component$numerator)) &&
           all(is.finite(component$denominator)) &&
           all(component$denominator > 0) && is.finite(value)))
}

principal_iid_truth <- function(model, on_order = NULL) {
  .principal_iid_check_model(model)
  if (!is.null(on_order) && !is.function(on_order))
    stop("on_order must be NULL or a progress-output function.")
  orders <- c(32L, 64L, 96L, 128L)
  rows <- components <- list()
  for (q in orders) {
    rule <- .principal_iid_gh_normal(q)
    for (overlap in model$treatment$overlap) {
      result <- .principal_iid_patt_order(model, overlap, rule, q)
      rows[[length(rows) + 1L]] <- result$row
      components[[length(components) + 1L]] <- result$component
      if (!is.null(on_order)) on_order(do.call(rbind, rows),
                                       do.call(rbind, components))
    }
  }
  order_table <- do.call(rbind, rows)
  order_table$successive_difference <- NA_real_
  targets <- list()
  for (overlap in model$treatment$overlap) {
    indices <- which(order_table$overlap == overlap)
    values <- order_table$PATT[indices]
    differences <- abs(diff(values))
    order_table$successive_difference[indices[-1L]] <- differences
    last_two_max <- max(differences[2:3])
    admitted <- all(order_table$finite[indices]) && is.finite(last_two_max) &&
      last_two_max <= 1e-7
    target <- if (admitted) values[4L] else NA_real_
    tau_A <- principal_iid_pate(model)
    targets[[length(targets) + 1L]] <- data.frame(overlap = overlap,
      estimand = c("PATE", "PATT"), target = c(tau_A, target),
      admitted = c(is.finite(tau_A), admitted),
      truth_method = c("analytic_fixed_mixture", "normal_GH_3D_Stein_order128"),
      numerical_tolerance = c(NA_real_, 1e-7),
      last_two_difference_max = c(NA_real_, last_two_max),
      sensitivity_lower = c(tau_A, target) - 1e-6,
      sensitivity_upper = c(tau_A, target) + 1e-6,
      sensitivity_is_rigorous_enclosure = FALSE)
  }
  targets <- do.call(rbind, targets)
  list(schema = "principal_iid_truth_v1", orders = order_table,
       components = do.call(rbind, components), targets = targets,
       admitted = all(targets$admitted),
       accuracy_scope = paste("Operational last-two-order convergence <= 1e-7;",
         "not a rigorous quadrature error bound. Preserve +/-1e-6 sensitivity."))
}
