# Source-only deterministic angle reference.
# The caller supplies env$.wm_alpha_integrand from wm_geometry_alpha.R.
# Adaptive-quadrature errors remain diagnostic, not certified bounds.
alpha_angles <- function(M, d, relative_tolerance = 1e-7) {
  evaluations <- 0L
  largest_inner_error <- 0
  F <- function(r, s, separation) {
    evaluations <<- evaluations + length(r)
    if (evaluations > 1000000L) stop("angle-reference evaluation work cap exceeded")
    env$.wm_alpha_integrand(r, s, separation, M, d)
  }
  if (d == 1L) {
    integrand <- function(theta) {
      r <- cos(theta)
      s <- sin(theta)
      (F(r, s, abs(r - s)) + F(r, s, r + s)) / pi
    }
  } else {
    integrand <- function(theta) vapply(theta, function(t) {
      if (t == 0 || t == pi / 2) return(0)
      r <- cos(t)
      s <- sin(t)
      inner <- stats::integrate(function(phi) {
        separation <- sqrt(pmax(0, 1 - sin(2 * t) * cos(phi)))
        density <- if (d == 2L) rep(1 / pi, length(phi)) else
          exp((d - 2) * log(sin(phi)) - lbeta((d - 1) / 2, 1 / 2))
        F(rep(r, length(phi)), rep(s, length(phi)), separation) * density
      }, 0, pi, rel.tol = relative_tolerance / 10, abs.tol = 1e-10,
      subdivisions = 200L, stop.on.error = TRUE)
      largest_inner_error <<- max(largest_inner_error, inner$abs.error)
      theta_density <- 2 * exp((d - 1) * (log(sin(t)) + log(cos(t))) -
                                  lbeta(d / 2, d / 2))
      inner$value * theta_density
    }, numeric(1L))
  }
  # Symmetry removes the radius-order kink at pi/4 from the interval interior.
  outer <- stats::integrate(integrand, 0, pi / 4,
    rel.tol = relative_tolerance, abs.tol = 1e-9,
    subdivisions = 200L, stop.on.error = TRUE)
  list(alpha = 2 * outer$value, outer_absolute_error = 2 * outer$abs.error,
       maximum_observed_inner_error = largest_inner_error, evaluations = evaluations,
       error_scope = "Observed nested adaptive-quadrature diagnostics, not a certified global bound")
}
