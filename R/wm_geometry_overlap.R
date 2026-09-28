# Experimental positive boundary-rank integral for complete neighbor-set overlap.
# This is separate from wm_geometry and from the alpha-only Poisson-count sum.

.wm_overlap_log_add <- function(x, y) {
  answer <- pmax(x, y)
  finite <- is.finite(answer)
  answer[finite] <- answer[finite] +
    log(exp(x[finite] - answer[finite]) + exp(y[finite] - answer[finite]))
  answer
}

.wm_overlap_surface_fraction <- function(radius, other_radius, separation, d) {
  n <- length(radius)
  answer <- numeric(n)
  positive <- radius > 0
  if (d == 1L) {
    answer[positive] <- (as.numeric(abs(radius[positive] - separation[positive]) <
                                      other_radius[positive]) +
                          as.numeric(radius[positive] + separation[positive] <
                                       other_radius[positive])) / 2
    return(answer)
  }
  concentric <- positive & separation == 0
  answer[concentric] <- as.numeric(radius[concentric] < other_radius[concentric])
  take <- positive & separation > 0
  if (any(take)) {
    r <- radius[take]; s <- other_radius[take]; delta <- separation[take]
    cutoff <- ((r - s) * (r + s) + delta^2) / (2 * r * delta)
    value <- as.numeric(cutoff <= -1)
    partial <- cutoff > -1 & cutoff < 1
    small <- stats::pbeta(cutoff[partial]^2, 1 / 2, (d - 1) / 2,
                          lower.tail = FALSE) / 2
    value[partial] <- ifelse(cutoff[partial] >= 0, small, 1 - small)
    answer[take] <- value
  }
  answer
}

.wm_overlap_ball_parts <- function(radius1, radius2, separation, d) {
  n <- length(radius1)
  if (length(radius2) != n || length(separation) != n ||
      any(!is.finite(c(radius1, radius2, separation))) ||
      any(c(radius1, radius2, separation) < 0) ||
      any(pmax(radius1, radius2) <= 0)) {
    stop("Boundary-rank balls need finite nonnegative radii and positive maximum radius.",
         call. = FALSE)
  }
  scale <- pmax(radius1, radius2)
  r <- radius1 / scale; s <- radius2 / scale; delta <- separation / scale
  vr <- r^d; vs <- s^d
  if (any((r > 0 & vr == 0) | (s > 0 & vs == 0))) {
    stop("Boundary-rank ball volume underflowed.", call. = FALSE)
  }
  intersection <- pmin(vr, vs)
  disjoint <- delta >= r + s
  intersection[disjoint] <- 0
  partial <- !disjoint & delta > abs(r - s)
  if (any(partial)) {
    rr <- r[partial]; ss <- s[partial]; dd <- delta[partial]
    x <- (dd + (rr - ss) * (rr + ss) / dd) / 2
    y <- dd - x
    volume_fraction <- function(t) {
      t <- pmax(-1, pmin(1, t))
      small <- stats::pbeta(t^2, 1 / 2, (d + 1) / 2,
                            lower.tail = FALSE) / 2
      ifelse(t >= 0, small, 1 - small)
    }
    intersection[partial] <- vr[partial] * volume_fraction(x / rr) +
      vs[partial] * volume_fraction(y / ss)
  }
  tolerance <- 1000 * .Machine$double.eps
  if (any(intersection < -tolerance | intersection > pmin(vr, vs) + tolerance)) {
    stop("Boundary-rank ball intersection violated its geometric bounds.", call. = FALSE)
  }
  # Only a geometric round-off correction. No estimated constant is projected.
  intersection <- pmax(0, pmin(pmin(vr, vs), intersection))
  first <- vr - intersection
  second <- vs - intersection
  union <- intersection + first + second
  log_ball <- d * log(pi) / 2 - lgamma(d / 2 + 1)
  volume_scale <- log_ball + d * log(scale)
  frac_r <- .wm_overlap_surface_fraction(r, s, delta, d)
  frac_s <- .wm_overlap_surface_fraction(s, r, delta, d)
  surface_r <- surface_s <- rep(-Inf, n)
  positive_r <- r > 0
  positive_s <- s > 0
  surface_r[positive_r] <- log(d) + log_ball +
    (d - 1) * (log(scale[positive_r]) + log(r[positive_r]))
  surface_s[positive_s] <- log(d) + log_ball +
    (d - 1) * (log(scale[positive_s]) + log(s[positive_s]))
  list(I = volume_scale + log(intersection),
       A = volume_scale + log(first), B = volume_scale + log(second),
       U = volume_scale + log(union),
       s = surface_r + log(frac_r), t = surface_s + log(frac_s),
       p = surface_r + log1p(-frac_r), q = surface_s + log1p(-frac_s))
}

.wm_overlap_log_monomial <- function(parts, coefficient, powers) {
  n <- length(parts$U)
  if (coefficient == 0) return(rep(-Inf, n))
  if (coefficient < 0 || any(powers < 0)) {
    stop("Invalid boundary-rank polynomial.", call. = FALSE)
  }
  answer <- rep(log(coefficient), n)
  for (name in names(powers)) {
    if (powers[[name]] > 0) answer <- answer + powers[[name]] * parts[[name]]
  }
  answer
}

.wm_overlap_log_polynomial <- function(parts, M, ell, family) {
  a <- M - 1L - ell
  monomial <- function(coefficient, powers) {
    .wm_overlap_log_monomial(parts, coefficient, powers)
  }
  if (family == 1L) return(monomial(1, c(I = ell, A = a, B = a)))
  result <- rep(-Inf, length(parts$U))
  add <- function(value) result <<- .wm_overlap_log_add(result, value)
  if (family == 2L) {
    if (ell > 0L) add(monomial(ell, c(I = ell - 1L, s = 1L, A = a, B = a)))
    if (a > 0L) add(monomial(a, c(I = ell, A = a - 1L, B = a, p = 1L)))
  } else if (family == 3L) {
    if (ell > 1L) add(monomial(ell * (ell - 1),
      c(I = ell - 2L, s = 1L, t = 1L, A = a, B = a)))
    if (ell > 0L && a > 0L) {
      add(monomial(ell * a, c(I = ell - 1L, s = 1L, q = 1L, A = a, B = a - 1L)))
      add(monomial(ell * a, c(I = ell - 1L, t = 1L, p = 1L, A = a - 1L, B = a)))
    }
    if (a > 0L) add(monomial(a^2, c(I = ell, A = a - 1L, B = a - 1L, p = 1L, q = 1L)))
  } else if (family == 4L && ell > 0L) {
    add(monomial(ell, c(I = ell - 1L, A = a, B = a)))
  }
  result
}

.wm_overlap_active <- function(M, d) {
  ell <- 0:(M - 1L)
  a <- M - 1L - ell
  active <- cbind(root_root = rep(TRUE, M), one_nonroot = rep(M > 1L, M),
                 distinct_nonroot = ell > 1L | a > 0L, shared_nonroot = ell > 0L)
  # In 1d root-balls are either nested or disjoint; intermediate shared/private
  # root-root combinations are structurally zero, not unresolved MC events.
  if (d == 1L) active[ell > 0L & a > 0L, "root_root"] <- FALSE
  active
}

.wm_overlap_angular_terms <- function(normals, M, d) {
  dimensions <- c(2 * d, 2 * d + 1, 2 * d + 2, 3 * d)
  active <- .wm_overlap_active(M, d)
  values <- matrix(0, nrow(normals), 4 * M)
  norm_at <- function(x) sqrt(rowSums(x * x))
  for (family in seq_len(4L)) {
    if (!any(active[, family])) next
    dimension <- dimensions[family]
    coordinate <- normals[, seq_len(dimension), drop = FALSE]
    length_at <- norm_at(coordinate)
    if (any(!is.finite(length_at)) || any(length_at <= 0)) {
      stop("Boundary-rank sphere generation encountered an invalid norm.", call. = FALSE)
    }
    coordinate <- coordinate / length_at
    x <- coordinate[, seq_len(d), drop = FALSE]
    y <- coordinate[, d + seq_len(d), drop = FALSE]
    rx <- norm_at(x); ry <- norm_at(y)
    if (family == 1L) {
      radius1 <- rx; radius2 <- ry; include <- rep(TRUE, nrow(coordinate))
    } else if (family == 2L) {
      radius1 <- coordinate[, 2 * d + 1L]; radius2 <- ry
      include <- radius1 > rx
    } else if (family == 3L) {
      radius1 <- coordinate[, 2 * d + 1L]; radius2 <- coordinate[, 2 * d + 2L]
      include <- radius1 > rx & radius2 > ry
    } else {
      u <- coordinate[, 2 * d + seq_len(d), drop = FALSE]
      radius1 <- norm_at(u - x); radius2 <- norm_at(u - y)
      include <- radius1 > rx & radius2 > ry
    }
    if (!any(include)) next
    parts <- .wm_overlap_ball_parts(radius1[include], radius2[include],
                                    norm_at(x - y)[include], d)
    log_area <- log(2) + dimension * log(pi) / 2 - lgamma(dimension / 2)
    for (ell in 0:(M - 1L)) {
      if (!active[ell + 1L, family]) next
      a <- M - 1L - ell
      k <- 2 * M - ell
      log_value <- .wm_overlap_log_polynomial(parts, M, ell, family) +
        lgamma(k) + log_area - log(d) - lgamma(ell + 1) - 2 * lgamma(a + 1) -
        k * parts$U + if (family == 2L) log(2) else 0
      if (any(is.na(log_value)) || any(log_value == Inf) ||
          any(log_value > log(.Machine$double.xmax))) {
        stop("Boundary-rank integrand exceeded numerical range.", call. = FALSE)
      }
      value <- exp(log_value)
      if (any(is.finite(log_value) & value == 0)) {
        stop("A positive boundary-rank integrand underflowed.", call. = FALSE)
      }
      values[include, (family - 1L) * M + ell + 1L] <- value
    }
  }
  values
}

.wm_overlap_mc <- function(M, d, draws, chunk_size) {
  families <- colnames(.wm_overlap_active(M, d))
  active <- as.vector(.wm_overlap_active(M, d))
  dimensions <- c(2 * d, 2 * d + 1, 2 * d + 2, 3 * d)
  max_dimension <- max(dimensions[colSums(matrix(active, M, 4L)) > 0])
  count <- 0L
  average <- numeric(4 * M)
  deviation <- matrix(0, 4 * M, 4 * M)
  nonzero <- maximum <- numeric(4 * M)
  beta_nonzero <- beta_maximum <- numeric(M)
  alpha_mean <- alpha_deviation <- 0
  value_limit <- sqrt(.Machine$double.xmax) / (16 * M * sqrt(draws))
  while (count < draws) {
    batch <- min(chunk_size, draws - count)
    normals <- matrix(stats::rnorm(batch * max_dimension), batch, byrow = TRUE)
    value <- .wm_overlap_angular_terms(normals, M, d)
    if (any(!is.finite(value)) || any(value > value_limit)) {
      stop("Boundary-rank sample moments exceeded numerical range.", call. = FALSE)
    }
    nonzero <- nonzero + colSums(value > 0)
    maximum <- pmax(maximum, apply(value, 2L, max))
    beta_value <- value[, seq_len(M), drop = FALSE]
    for (family in 2:4) beta_value <- beta_value + value[, (family - 1L) * M + seq_len(M), drop = FALSE]
    beta_nonzero <- beta_nonzero + colSums(beta_value > 0)
    beta_maximum <- pmax(beta_maximum, apply(beta_value, 2L, max))
    batch_average <- colMeans(value)
    centered <- sweep(value, 2L, batch_average, "-")
    change <- batch_average - average
    updated <- count + batch
    deviation <- deviation + crossprod(centered) +
      tcrossprod(change) * (count * as.double(batch) / updated)
    average <- average + change * (batch / updated)
    alpha_value <- rowSums(value)
    alpha_batch <- mean(alpha_value)
    alpha_change <- alpha_batch - alpha_mean
    alpha_deviation <- alpha_deviation + sum((alpha_value - alpha_batch)^2) +
      alpha_change^2 * (count * as.double(batch) / updated)
    alpha_mean <- alpha_mean + alpha_change * (batch / updated)
    count <- updated
  }
  overlap <- as.character(0:(M - 1L))
  term_names <- unlist(lapply(families, function(f) paste(f, overlap, sep = ":")), use.names = FALSE)
  term_covariance <- deviation / (draws - 1) / draws
  dimnames(term_covariance) <- list(term_names, term_names)
  term_estimates <- matrix(average, M, 4L, dimnames = list(overlap, families))
  collapse <- matrix(0, 4 * M, M)
  collapse[cbind(seq_len(4 * M), rep(seq_len(M), 4L))] <- 1
  covariance <- crossprod(collapse, term_covariance %*% collapse)
  dimnames(covariance) <- list(overlap, overlap)
  list(beta = rowSums(term_estimates), covariance = covariance,
       beta_mcse = stats::setNames(sqrt(pmax(0, diag(covariance))), overlap),
       alpha_mcse = sqrt(max(0, alpha_deviation / (draws - 1) / draws)),
       term_estimates = term_estimates, term_covariance = term_covariance,
       term_mcse = matrix(sqrt(pmax(0, diag(term_covariance))), M, 4L,
                          dimnames = list(overlap, families)),
       nonzero = stats::setNames(beta_nonzero, overlap),
       maximum = stats::setNames(beta_maximum, overlap),
       term_nonzero = matrix(nonzero, M, 4L, dimnames = list(overlap, families)),
       term_maximum = matrix(maximum, M, 4L, dimnames = list(overlap, families)))
}

.wm_geometry_overlap_dispatch <- function(M, d, draws, seed, chunk_size, exact_1d) {
  positive_integer <- function(x) is.numeric(x) && !is.complex(x) &&
    length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x) &&
    x <= .Machine$integer.max
  for (name in c("M", "d", "draws", "chunk_size")) {
    if (!positive_integer(get(name))) stop(name, " must be one positive integer.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || is.complex(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max || seed != floor(seed))) {
    stop("seed must be NULL or one nonnegative integer.", call. = FALSE)
  }
  M <- as.integer(M); d <- as.integer(d); draws <- as.integer(draws)
  requested_chunk <- as.integer(chunk_size)
  if (d == 1L && exact_1d) return(wm_geometry(M, d, draws, seed, chunk_size))
  if (draws < 2L) stop("Boundary-rank Monte Carlo requires draws >= 2.", call. = FALSE)
  if (M > 64L || d > 128L) stop("Boundary-rank work cap requires M <= 64 and d <= 128.", call. = FALSE)
  active <- .wm_overlap_active(M, d)
  dimensions <- c(2 * d, 2 * d + 1, 2 * d + 2, 3 * d)
  geometric_work <- as.double(draws) * sum(dimensions[colSums(active) > 0])
  covariance_work <- as.double(draws) * (4 * as.double(M))^2
  if (geometric_work > 50000000 || covariance_work > 100000000) {
    stop("Boundary-rank integration exceeds the bounded work cap (50 million sphere coordinates or 100 million covariance products).",
         call. = FALSE)
  }
  max_dimension <- max(dimensions[colSums(active) > 0])
  chunk_size <- as.integer(min(requested_chunk, draws, floor(1000000 / max(max_dimension, 4 * M))))
  if (!is.null(seed)) {
    if (identical(RNGkind()[2L], "Box-Muller")) {
      stop("An explicit seed cannot preserve the Box-Muller cache; use seed = NULL or another normal RNG.", call. = FALSE)
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    })
    set.seed(as.integer(seed))
  }
  result <- .wm_overlap_mc(M, d, draws, chunk_size)
  term_precision <- matrix("monte_carlo_estimate", M, 4L, dimnames = dimnames(result$term_estimates))
  term_precision[!active] <- "structurally_zero"
  term_precision[active & result$term_nonzero == 0] <- "unresolved_zero_hits"
  term_precision[active & result$term_nonzero > 0 & result$term_mcse == 0] <- "unresolved_zero_empirical_variance"
  unresolved_term <- apply(matrix(grepl("^unresolved", term_precision), M, 4L), 1L, any)
  component_precision <- ifelse(result$nonzero == 0, "unresolved_zero_hits",
    ifelse(result$beta_mcse == 0, "unresolved_zero_empirical_variance",
      ifelse(unresolved_term, "unresolved_boundary_term", "monte_carlo_estimate")))
  structure(list(M = M, d = d, overlap = 0:(M - 1L), beta = result$beta,
    alpha = sum(result$beta), beta_mcse = result$beta_mcse, alpha_mcse = result$alpha_mcse,
    covariance = result$covariance, method = "boundary_rank_monte_carlo",
    precision_status = if (any(grepl("^unresolved", component_precision))) "unresolved_components" else "monte_carlo_estimate",
    component_precision = unname(component_precision), draws = draws, requested_draws = draws,
    seed = seed, chunk_size = chunk_size, requested_chunk_size = requested_chunk,
    definition = paste("Complete order-M overlap beta: positive root/root, twice root/nonroot,",
      "distinct-boundary and shared-boundary Mecke integrals; Gamma(2M-l)/(d l! ((M-1-l)!)^2)",
      "times the corresponding sphere-area-weighted angular expectations. alpha=sum(beta)."),
    term_estimates = result$term_estimates, term_mcse = result$term_mcse,
    term_covariance = result$term_covariance, term_precision = term_precision,
    normalization = data.frame(family = colnames(active), sphere_dimension = dimensions,
      multiplicity = c(1, 2, 1, 1),
      log_sphere_area = log(2) + dimensions * log(pi) / 2 - lgamma(dimensions / 2)),
    diagnostics = list(nonzero_draws = result$nonzero, maximum_draw_value = result$maximum,
      term_nonzero_draws = result$term_nonzero, term_maximum_draw_value = result$term_maximum,
      active_terms = active, rng_used = TRUE, coupled_components = TRUE, coupled_boundary_terms = TRUE,
      jensen_lower_bound = as.double(M)^2, below_jensen = sum(result$beta) < as.double(M)^2,
      geometric_coordinate_work = geometric_work, covariance_product_work = covariance_work,
      gaussian_allocation_cap = 1000000, sphere_coordinate_work_cap = 50000000,
      covariance_product_work_cap = 100000000,
      validation_status = "experimental_method_requires_independent_numerical_validation",
      note = paste("All boundary terms and overlap components share normalized Gaussian prefixes.",
        "term_covariance and covariance are covariance matrices OF the Monte Carlo estimates.",
        "Empirical MCSE is not an error guarantee; zero-hit and zero-variance active terms remain unresolved.",
        "No constant is projected onto Jensen's bound."))),
    class = c("wm_geometry", "list"))
}

#' Experimental boundary-rank integration of matching overlap constants
#' @export
wm_geometry_overlap <- function(M, d, draws = 4096L, seed = NULL, chunk_size = 256L) {
  .wm_geometry_overlap_dispatch(M, d, draws, seed, chunk_size, exact_1d = TRUE)
}
