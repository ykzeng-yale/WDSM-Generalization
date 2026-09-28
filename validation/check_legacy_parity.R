# The legacy public API and generalized core agree when scores, score scaling,
# fitted bias correction, donor rule, and donor ordering are aligned.
# Run: Rscript code-release/validation/check_legacy_parity.R
# This validates point estimates, NOT equivalence of the distinct bootstraps.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this validation with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
env <- new.env(parent = globalenv())
for (f in sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE))) sys.source(f, env)
equal <- function(a, b) {
  comparison <- all.equal(a, b, tolerance = 1e-11, check.attributes = TRUE)
  if (!isTRUE(comparison)) stop(paste(comparison, collapse = "; "))
}
set.seed(928521)
n <- 160L
X <- data.frame(x1 = runif(n, -1, 1), x2 = rnorm(n), x3 = runif(n))
Z <- rep(0:1, n / 2)
weights <- exp(0.2 * X$x1 + 0.15 * Z)
ps <- plogis(0.6 * X$x1 - 0.3 * X$x2)
pg <- cbind(0.4 + X$x2 + 0.7 * X$x3, 1 + 0.2 * X$x1 - 0.4 * X$x3)
Y <- pg[cbind(seq_len(n), Z + 1L)] + 0.3 * X$x2^2 + rnorm(n, sd = 0.2)
checked <- 0L
for (estimand in c("PATE", "PATT")) for (M in c(1L, 3L)) {
  public <- if (estimand == "PATE") env$wdsmatchATE else env$wdsmatchATT
  # Keep the public API's default random/tolerance tie policy. This continuous
  # fixture must prove that policy is inactive before claiming parity.
  old <- public(Y, X, Z, weights, M = M, ps = ps, pg = pg,
                use.bias.correction = TRUE, varest = FALSE)
  stopifnot(all(vapply(old$diagnostics$ties$arms, function(a)
    a$boundary_tie_recipients == 0L && a$randomized_recipients == 0L, logical(1))))
  scores <- env$estimate_scores(Y, X, Z, weights, ps, pg,
                                estimand = estimand, use.bias.correction = TRUE)
  args <- list(Y = Y, Z = Z, weights = weights, scores0 = scores$D0,
               M = M, estimand = estimand, method = "self_normalized",
               mean0 = scores$q0, variance = FALSE)
  if (estimand == "PATE") {
    args$scores1 <- scores$D1
    args$mean1 <- scores$q1
  }
  new <- do.call(env$wm_match, args)
  allocation <- env$wdsm_make_matches(Z, scores$D0, scores$D1, M, estimand)
  neighbors <- allocation$matches_0
  if (estimand == "PATE") neighbors[Z == 0] <- allocation$matches_1[Z == 0]
  stopifnot(identical(new$graph$neighbors, neighbors))
  equal(new$estimate, old$estimate)
  K <- env$wdsm_reuse(Z, weights, allocation$matches_0, allocation$matches_1, estimand)
  equal(unname(rowSums(new$loads$incoming) * new$weight_scale), as.numeric(K))
  checked <- checked + 1L
}

# Exact ties demonstrate the policy boundary. Align the legacy internal graph
# explicitly; do not silently replace the public API's random tie convention.
Z <- rep(0:1, each = 4L)
D <- matrix(0, 8L, 2L)
fixed <- env$wdsm_make_matches(Z, D, D, M = 2, tie_rule = "row_order", tie_tolerance = 0)
core <- env$wm_match(1:8, Z, rep(1, 8), D, M = 2, variance = FALSE)
neighbors <- fixed$matches_0
neighbors[Z == 0] <- fixed$matches_1[Z == 0]
stopifnot(identical(core$graph$neighbors, neighbors))
random <- env$wdsm_make_matches(Z, D, D, M = 2)
stopifnot(random$tie_diagnostics$rule == "random", all(vapply(random$tie_diagnostics$arms,
  function(a) a$randomized_recipients > 0L, logical(1))))

# The default legacy tolerance can also randomize unequal but very close
# squared distances. The core has no tolerance and deterministically ranks
# the exact Euclidean distances. This is intentionally not a parity claim.
D <- cbind(c(1, 1 + 1e-15, 3, 4, 0, 10, 20, 30), 0)
near <- env$wdsm_make_matches(Z, D, D, M = 1)
core <- env$wm_match(1:8, Z, rep(1, 8), D, M = 1, variance = FALSE)
stopifnot(core$graph$neighbors[[5L]] == 1L,
          near$tie_diagnostics$arms[["0"]]$boundary_tie_recipients > 0L,
          near$tie_diagnostics$arms[["0"]]$exact_boundary_tie_recipients == 0L)
cat("PASS:", checked, "legacy public point-estimate/graph/load parity cases;",
    "exact and tolerance tie policies checked separately.\n")
print(sessionInfo())
