# Synthetic data only. Callers set a reproducible RNG stream before generation.
wm_simulate <- function(n, d, design = c("strong", "weak")) {
  design <- match.arg(design)
  valid <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= 1 && x == floor(x) && x <= .Machine$integer.max
  if (!valid(n) || n < 2 || !valid(d) || n * d > 1e7)
    stop("n and d must be positive integers within the allocation guard")
  n <- as.integer(n); d <- as.integer(d)
  S <- matrix(stats::runif(n * d), n, d)
  if (design == "strong") {
    Z <- stats::rbinom(n, 1, .4)
    U <- stats::rbinom(n, 1, .5)
    W <- ifelse(Z == 0, .5 + U, 1 + 2 * U)
    variance <- ifelse(Z == 0, .25 + .75 * U, .5 + U)
    mean0 <- 1 + S[, 1] + rowSums(S^2) / (2 * d)
    mean1 <- mean0 + 1 + S[, 1]
    Y <- ifelse(Z == 0, mean0, mean1) + sqrt(variance) *
      sample(c(-1, 1), n, replace = TRUE)
    rho0 <- rep(1, n); rho1 <- rep(2, n)
    target <- c(PATE = 1.5, PATT = 1.5)
    degree <- 2L
  } else {
    Z <- stats::rbinom(n, 1, 8 / 17)
    U <- stats::rbinom(n, 1, ifelse(Z == 0, 1 / 3, 1 / 4))
    W <- ifelse(Z == 0, 1 + U, 1 + 2 * U)
    Y <- as.numeric(U)
    mean0 <- mean1 <- rep(.5, n)
    rho0 <- rep(4 / 3, n); rho1 <- rep(3 / 2, n)
    target <- c(PATE = 0, PATT = 0)
    degree <- 0L
  }
  list(Y = Y, Z = Z, weights = W, scores = S, U = U,
       mean0 = mean0, mean1 = mean1, rho0 = rho0, rho1 = rho1,
       target = target, degree = degree, design = design)
}

# An indexed L'Ecuyer-CMRG stream makes batches invariant to worker count/order.
wm_sim_stream <- function(base_seed, replication) {
  if (!is.numeric(base_seed) || length(base_seed) != 1L || !is.finite(base_seed) ||
      base_seed < 0 || base_seed > .Machine$integer.max || base_seed != floor(base_seed) ||
      !is.numeric(replication) || length(replication) != 1L || !is.finite(replication) ||
      replication < 1 || replication != floor(replication)) stop("invalid seed or replication")
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  previous <- if (had_seed) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  previous_kind <- RNGkind()
  if (identical(previous_kind[2L], "Box-Muller")) stop("use a noncached normal RNG before indexing streams")
  on.exit({
    do.call(RNGkind, as.list(previous_kind))
    if (had_seed) assign(".Random.seed", previous, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(as.integer(base_seed), kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  stream <- .Random.seed
  if (replication > 1) for (i in seq_len(replication - 1)) stream <- parallel::nextRNGStream(stream)
  stream
}
