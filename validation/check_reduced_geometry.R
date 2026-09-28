# Deterministic angle quadrature independently checks the sphere sampler.
# Nested quadrature error diagnostics are not rigorous certified error bounds.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript")
script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
env <- new.env(parent = globalenv())
sys.source(file.path(root, "R", "wm_geometry.R"), envir = env)
sys.source(file.path(root, "R", "wm_geometry_alpha.R"), envir = env)

source(file.path(dirname(script), "geometry_angle_reference.R"))

results <- list()
for (M in c(1L, 2L, 3L, 5L)) {
  reference <- alpha_angles(M, 1L, relative_tolerance = 1e-9)
  exact <- M^2 + M / 2
  stopifnot(abs(reference$alpha - exact) < 1e-8 * exact)
  results[[length(results) + 1L]] <- data.frame(M = M, d = 1L,
    exact = exact, angle_alpha = reference$alpha, mc_alpha = NA_real_,
    mcse = 0, discrepancy_in_mcse = NA_real_,
    outer_error = reference$outer_absolute_error,
    largest_observed_inner_error = reference$maximum_observed_inner_error,
    evaluations = reference$evaluations)
}
for (d in c(2L, 3L, 5L)) for (M in c(1L, 3L)) {
  reference <- alpha_angles(M, d)
  result <- env$wm_geometry_alpha(M, d, draws = 16384L,
                                 seed = 930000L + 100L * d + M)
  discrepancy <- (result$alpha - reference$alpha) / result$alpha_mcse
  # This fixed-seed discrepancy check is a numerical diagnostic, not proof
  # that a Monte Carlo interval always covers.
  stopifnot(is.finite(discrepancy), abs(discrepancy) < 7,
            result$diagnostics$maximum_draw_value <= result$diagnostics$integrand_upper_bound)
  results[[length(results) + 1L]] <- data.frame(M = M, d = d,
    exact = NA_real_, angle_alpha = reference$alpha, mc_alpha = result$alpha,
    mcse = result$alpha_mcse, discrepancy_in_mcse = discrepancy,
    outer_error = reference$outer_absolute_error,
    largest_observed_inner_error = reference$maximum_observed_inner_error,
    evaluations = reference$evaluations)
  if (M == 1L) {
    original <- env$wm_geometry(M, d, draws = 256L, seed = 832L)
    reduced <- env$wm_geometry_alpha(M, d, draws = 256L, seed = 832L)
    stopifnot(isTRUE(all.equal(original$alpha, reduced$alpha, tolerance = 1e-12)),
              isTRUE(all.equal(original$alpha_mcse, reduced$alpha_mcse, tolerance = 1e-12)))
  }
}
print(do.call(rbind, results), row.names = FALSE, digits = 10)
cat("PASS: exact 1D quadrature, identical M1 sphere draws, and reduced-alpha MC versus two-angle quadrature.\n")
print(sessionInfo())
