# Small deterministic calculation checks; not a Monte Carlo experiment.
source("R/utils.R")
source("R/wm_core.R")
source("R/wm_scalar_logistic_match.R")
checks <- character()
check <- function(value, label) {
  if (!isTRUE(value)) stop(label)
  checks <<- c(checks, label)
}
near <- function(x, y, label, tolerance = 2e-10) {
  check(identical(dim(x), dim(y)) && length(x) == length(y) &&
    all(is.finite(c(x, y))) &&
    max(abs(x - y)) <= tolerance * max(1, abs(x), abs(y)), label)
}
fails <- function(expr, label) check(inherits(try(expr, silent = TRUE), "try-error"), label)
T <- seq(-1.1, 1.1, length.out = 24)
V <- rep(c(-0.8, 0.2, 0.6, -0.3), 6)
L <- cbind(intercept = 1, index = T, auxiliary = V)
Z <- rep(c(0, 1, 1, 0, 1, 0), 4)
U <- plogis(T)
Y <- 0.5 + 0.8 * Z + T^2 - 0.3 * V + rep(c(-0.2, 0.1, 0.3), 8)
model <- list(w0 = function(u) 2 + 0.2 * u, w1 = function(u) 1.5 + 0.1 * u,
              dw0 = function(u) rep(0.2, length(u)),
              dw1 = function(u) rep(0.1, length(u)), domain = c(0.001, 0.999))
W <- ifelse(Z == 1, model$w1(U), model$w0(U))
start_time <- proc.time()[["elapsed"]]
had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
if (had_rng) rng_before <- get(".Random.seed", envir = .GlobalEnv)
fits <- lapply(c("PATE", "PATT"), function(target)
  wm_scalar_logistic_match(Y, Z, W, L, M = 3, estimand = target, weight_model = model))
names(fits) <- c("PATE", "PATT")
for (target in names(fits)) {
  x <- fits[[target]]
  check(identical(x$status, "point_computed"), paste(target, "point success"))
  check(identical(x$inference$status, "conditional_formula"), paste(target, "conditional formula"))
  check(!x$fit$info$corrected && is.na(x$fit$variance), paste(target, "raw fit preserved"))
  check(identical(x$numerical_root_bridge, "not_certified") &&
          !x$inference$population_assumptions_verified &&
          !x$inference$bootstrap_validity_claimed, paste(target, "scope retained"))
  e <- x$propensity$probability
  direct <- wm_match(Y, Z, W, matrix(e, ncol = 1), M = 3,
                     estimand = target, variance = FALSE)
  check(identical(x$fit, direct), paste(target, "exact existing engine identity"))
  check(x$propensity$diagnostics$rank == ncol(L), paste(target, "all coefficients fitted"))
  # Independent direct edge reconstruction, retaining finite-M donor ratios.
  contrast <- numeric(length(Y))
  query <- if (target == "PATE") seq_along(Y) else which(Z == 1)
  for (i in query) {
    pool <- which(Z != Z[i])
    selected <- pool[order(abs(e[pool] - e[i]), pool)[1:3]]
    imputed <- sum(W[selected] * Y[selected]) / sum(W[selected])
    contrast[i] <- (2 * Z[i] - 1) * (Y[i] - imputed)
  }
  point <- sum(W[query] * contrast[query]) / sum(W[query])
  near(point, x$estimate, paste(target, "independent weighted edge point"))
  # Independent per-query NW means versus the bounded-memory batch routine.
  predictions <- sapply(0:1, function(z) vapply(e, function(at) {
    d <- which(Z == z)
    kernel <- vapply(e[d], function(score) max(0, 1 - abs(score - at) /
                                               x$inference$predictions$bandwidth), numeric(1))
    sum(kernel * Y[d]) / sum(kernel)
  }, numeric(1)))
  near(unname(predictions), unname(x$inference$predictions$prediction),
       paste(target, "independent own-arm kernel means"))
  # Raw-unit formula evaluated with explicit row sums, outside the helper.
  n <- length(Y); sw <- sum(W); a <- n / sw; rho <- sum(W * Z) / sw
  residual <- Y - predictions[cbind(seq_len(n), Z + 1)]
  d <- predictions[, 2] - predictions[, 1] - point
  J <- Reduce("+", lapply(seq_len(n), function(i)
    W[i] * e[i] * (1 - e[i]) * tcrossprod(L[i, ]))) / n
  B <- Reduce("+", lapply(seq_len(n), function(i)
    W[i]^2 * (Z[i] - e[i])^2 * tcrossprod(L[i, ]))) / n
  Mmoment <- function(z, F) {
    i <- which(Z == z); p <- if (z == 1) e[i] else 1 - e[i]
    sum(W[i] * F[i] * residual[i]^2 / p) / sw
  }
  Qmoment <- function(z, F) {
    i <- which(Z == z); p <- if (z == 1) e[i] else 1 - e[i]
    Reduce("+", lapply(seq_along(i), function(j)
      W[i[j]] * F[i[j]] * residual[i[j]] * L[i[j], ] / p[j])) / sw
  }
  Ev <- function(F) Reduce("+", lapply(seq_len(n), function(i)
    W[i] * F[i] * L[i, ])) / sw
  w0 <- model$w0(e); w1 <- model$w1(e); f <- e * (1 - e); M <- 3
  F0 <- e + e^2 * (1 - e) * model$dw0(e) / (M * w0)
  F1 <- 1 - e - e * (1 - e)^2 * model$dw1(e) / (M * w1)
  if (target == "PATE") {
    v0 <- a * (Mmoment(1, w1/e + w1*(1-e)^2/(2*M*e) + w0*(1-e)/M) +
      Mmoment(0, w0/(1-e) + w0*e^2/(2*M*(1-e)) + w1*e/M) +
      sum(W * (e*w1+(1-e)*w0) * d^2)/sw)
    b <- -Qmoment(0, F0) - Qmoment(1, F1)
    C <- solve(J, Qmoment(1, w1*(1-e)) + Qmoment(0, w0*e) + Ev(f*(w1-w0)*d))
  } else {
    v0 <- a/rho^2 * (Mmoment(1, w1*e) +
      Mmoment(0, w1*e/M + (1+1/(2*M))*w0*e^2/(1-e)) + sum(W*e*w1*d^2)/sw)
    b <- -Qmoment(0, F0)/rho
    C <- solve(J, Qmoment(1, f*w1) + Qmoment(0, e^2*w0) + Ev(f*w1*d))/rho
  }
  Sigma <- solve(J) %*% B %*% solve(J)
  expected <- v0 + 2*sum(b*C) + drop(crossprod(b, Sigma %*% b))
  near(v0, x$inference$blocks$V0, paste(target, "raw-unit V0"))
  near(unname(b), unname(x$inference$blocks$b), paste(target, "raw-unit drift"))
  near(unname(C), unname(x$inference$blocks$C), paste(target, "raw-unit covariance"))
  near(unname(Sigma), unname(x$inference$blocks$Sigma), paste(target, "raw-unit sandwich"))
  near(expected, x$inference$root_n_variance, paste(target, "complete variance formula"))
  near(expected/n, x$inference$variance, paste(target, "total-n scaling"))
  # Rescale the complete weight model; no propensity refit or fresh graph.
  scaled_fit <- x$fit; scaled_fit$weights <- 7 * W
  scaled_model <- model
  for (name in c("w0","w1","dw0","dw1")) {
    fun <- model[[name]]
    scaled_model[[name]] <- local({ original <- fun; function(u) 7 * original(u) })
  }
  scaled <- .wm_scalar_moments(scaled_fit, L, e, scaled_model,
                              x$inference$predictions$bandwidth, 0.95)
  near(scaled$root_n_variance, expected, paste(target, "common weight scale"))
  # Same exact predictors under a nonsingular affine design transform.
  A <- rbind(c(1,0,0), c(0.3,1.2,0.2), c(-0.2,-0.1,0.9))
  transformed <- .wm_scalar_moments(x$fit, L %*% t(A), e, model,
                                    x$inference$predictions$bandwidth, 0.95)
  near(transformed$root_n_variance, expected, paste(target, "affine scalar invariance"))
  near(unname(transformed$blocks$b), unname(drop(A %*% b)), paste(target, "affine full b"))
  near(unname(transformed$blocks$C), unname(drop(solve(t(A), C))),
       paste(target, "affine full C"))
}
fails(wm_scalar_logistic_match(Y,Z,W,L,weight_model=NULL), "missing weight model")
bad <- model; bad$domain <- c(0,1)
fails(.wm_scalar_weight_model(bad), "invalid callback domain")
bad <- model; bad$w0 <- 1
fails(.wm_scalar_weight_model(bad), "nonfunction callback")
fails(wm_scalar_logistic_match(Y,Z,W,L,M=0,weight_model=model), "invalid M")
fails(wm_scalar_logistic_match(Y,Z,W,L,weight_model=model,maxit=1i), "complex maxit")
fails(wm_scalar_logistic_match(Y,Z,W,L,weight_model=model,score_tolerance=1i),
      "complex tolerance")
fails(wm_scalar_logistic_match(Y,Z,W,cbind(L,L[,2]),weight_model=model),
      "singular design")
fails(wm_scalar_logistic_match(Y,Z,W,L[,2:3],weight_model=model), "missing intercept")
# One additional fit checks point preservation after a declared inference failure.
outside <- model; outside$domain <- c(0.999, 0.9999)
failed_inference <- wm_scalar_logistic_match(Y,Z,W,L,M=3,estimand="PATE",
                                           weight_model=outside)
check(identical(failed_inference$fit, fits$PATE$fit), "inference failure preserves full point")
check(failed_inference$inference$status == "unavailable" &&
      all(is.na(failed_inference$inference$conf.int)), "unavailable interval retained")
fails(.wm_scalar_nw(seq(0.1,0.9,length.out=24), Y, Z, 1e-20),
      "zero kernel denominator retained")
check(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE) == had_rng &&
      (!had_rng || identical(rng_before, get(".Random.seed", envir=.GlobalEnv))),
      "no RNG used")

# Final reporting guards and declared inference scope.
model$label <- "deterministic_known_weight_fixture"
ordinary <- .wm_scalar_interval(1,2,24,0.95)
check(ordinary$status=="conditional_formula","ordinary interval")
tail <- .wm_scalar_interval(1,2,24,1-.Machine$double.eps/2)
check(tail$status=="conditional_formula" && all(is.finite(tail$conf.int)),
      "near-one level has finite stable tail")
collapsed <- .wm_scalar_interval(1e8,1e-18,100,0.95)
check(collapsed$status=="unrepresentable_interval" && all(is.na(collapsed$conf.int)),
      "rounded collapsed interval withheld")
check(collapsed$root_n_variance==1e-18 && collapsed$computed_sampling_variance>0,
      "collapsed interval retains variance evidence")
underflow <- .wm_scalar_interval(0,.Machine$double.xmin*.Machine$double.eps,24,0.95)
check(underflow$status=="unrepresentable_sampling_variance" &&
      all(is.na(underflow$conf.int)), "positive variance underflow withheld")
negative <- .wm_scalar_interval(0,-1,24,0.95)
check(negative$status=="nonpositive_root_variance" && negative$root_n_variance == -1,
      "negative raw variance retained")
tiny_level <- .wm_scalar_interval(1,2,24,1e-300)
check(tiny_level$status=="unrepresentable_interval_scale",
      "unrepresentable central probability withheld")
# One new two-column fit verifies the added inference-dimension failure path.
old <- fits$PATE
low_dimension <- wm_scalar_logistic_match(old$fit$data$Y,old$fit$data$Z,old$fit$weights,
  L[,1:2],M=3,weight_model=model)
check(low_dimension$status=="point_computed" && is.finite(low_dimension$estimate),
      "lower-dimensional point retained")
check(low_dimension$inference$status=="unavailable" &&
      grepl("three design columns",low_dimension$inference$reason),
      "unproved lower-dimensional inference withheld")
check(is.na(low_dimension$fit$variance) && !low_dimension$fit$info$corrected,
      "lower-dimensional raw engine unchanged")
check(low_dimension$bandwidth_rule=="n^(-1/5)" &&
      identical(low_dimension$weight_model$label,model$label),
      "unavailable inference retains provenance")
# Custom rate is a computation input, not an automatically verified sequence.
custom <- .wm_scalar_moments(old$fit,L,old$propensity$probability,model,0.6,0.95)
check(custom$predictions$bandwidth_rule=="user_supplied" &&
      !custom$predictions$user_bandwidth_rate_verified,"custom bandwidth not certified")

check(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE) == had_rng &&
      (!had_rng || identical(rng_before,get(".Random.seed",envir=.GlobalEnv))),
      "no RNG through final reporting guards")
cat(sprintf("%d deterministic checks PASS; %.3f seconds; one24-row fixture,4 fits,0 random draws.\n",
            length(checks),proc.time()[["elapsed"]]-start_time))
