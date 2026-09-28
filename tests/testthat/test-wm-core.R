wm_fixture <- function() {
  list(Y = c(1, 3, 8, 4, 7, 13), Z = c(0, 0, 0, 1, 1, 1),
       weights = c(1, 2, 4, 3, 5, 7),
       scores0 = matrix(c(0, 1, 2, 0.5, 1.5, 3), ncol = 1L),
       M = 2L, mean0 = c(0, 1, 3, 0.4, 2, 5),
       mean1 = c(2, 4, 6, 3, 5, 9))
}

wm_fixture_neighbors <- function() {
  list(c(4L, 5L), c(4L, 5L), c(5L, 6L),
       c(1L, 2L), c(2L, 3L), c(3L, 2L))
}

# Independent exhaustive oracle: scan every row for cell/arm eligibility,
# compute each robust Euclidean norm separately, and sort all donor distances.
wm_brute_neighbors <- function(args) {
  n <- length(args$Y)
  ans <- vector("list", n)
  fold <- if (is.null(args$fold_id)) rep(1, n) else args$fold_id
  stratum <- if (is.null(args$strata)) rep(1, n) else args$strata
  pate <- is.null(args$estimand) || args$estimand == "PATE"
  queries <- if (pate) seq_len(n) else which(args$Z == 1)
  for (i in queries) {
    donors <- which(args$Z != args$Z[i] & fold == fold[i] & stratum == stratum[i])
    scores <- if (args$Z[i] == 1 || is.null(args$scores1)) args$scores0 else args$scores1
    distance <- vapply(donors, function(j) {
      delta <- scores[j, ] - scores[i, ]
      scale <- max(abs(delta))
      if (scale == 0) 0 else scale * sqrt(sum((delta / scale)^2))
    }, numeric(1))
    ans[[i]] <- donors[order(distance, donors, method = "radix")[seq_len(args$M)]]
  }
  ans
}

test_that("partial selection agrees with exhaustive graphs across dimensions and cells", {
  set.seed(92831)
  n <- 96L
  for (dimension in c(1L, 2L, 3L, 6L)) for (M in c(1L, 3L)) {
    score <- matrix(runif(n * dimension, -2, 2), n, dimension)
    # Rounding creates repeated rows and exact distance boundaries. The large
    # case uses dyadic values so the exhaustive norm has unambiguous ties.
    for (coordinates in list(score, round(score), round(score) * 2^660)) {
      args <- list(Y = rnorm(n), Z = rep(0:1, n / 2), weights = 1 + runif(n),
                   scores0 = coordinates, M = M, mean0 = runif(n), mean1 = runif(n),
                   fold_id = rep(c("a.b", "a"), each = n / 2),
                   strata = rep(rep(c("c", "b.c"), each = n / 4), 2))
      for (estimand in c("PATE", "PATT")) for (method in c("self_normalized", "stabilized")) {
        args$estimand <- estimand
        args$method <- method
        args$mean1 <- if (estimand == "PATE") rep(0.8, n) else NULL
        args$rho0 <- if (method == "stabilized") rep(1.4, n) else NULL
        args$rho1 <- if (method == "stabilized" && estimand == "PATE") rep(1.6, n) else NULL
        fit <- do.call(wm_match, args)
        expected <- wm_brute_neighbors(args)
        expect_identical(fit$graph$neighbors, expected,
                         info = paste(dimension, M, estimand, method))
        expect_identical(fit$graph$edges$query,
                         rep(if (estimand == "PATE") seq_len(n) else which(args$Z == 1), each = M))
        expect_identical(fit$graph$edges$donor, unlist(expected, use.names = FALSE))
      }
    }
  }
})

test_that("full donor selection, unequal directional maps and boundary ties are exact", {
  args <- wm_fixture()
  args$M <- 3L
  args$scores1 <- cbind(args$scores0, c(3, 1, 4, 2, 6, 5))
  expect_identical(do.call(wm_match, args)$graph$neighbors, wm_brute_neighbors(args))
  args$M <- 2L
  args$scores0 <- matrix(c(0, 0, 0, 1, 1, 1))
  args$scores1 <- NULL
  expect_identical(do.call(wm_match, args)$graph$neighbors,
                   c(rep(list(c(4L, 5L)), 3), rep(list(c(1L, 2L)), 3)))
  args$M <- .Machine$integer.max
  expect_error(do.call(wm_match, args), "Insufficient donors for query row 1")
})

test_that("self-normalized estimates and loads agree with a hand graph", {
  d <- wm_fixture()
  fit <- do.call(wm_match, d)
  nn <- wm_fixture_neighbors()
  expect_identical(fit$graph$neighbors, nn)
  y0 <- y1 <- d$Y
  K <- numeric(6)
  for (i in 1:6) {
    j <- nn[[i]]
    lambda <- d$weights[j] / sum(d$weights[j])
    mu <- if (d$Z[i] == 1) d$mean0 else d$mean1
    prediction <- mu[i] + sum(lambda * (d$Y[j] - mu[j]))
    if (d$Z[i] == 1) y0[i] <- prediction else y1[i] <- prediction
    for (k in seq_along(j)) K[j[k]] <- K[j[k]] + d$weights[i] * lambda[k]
  }
  tau <- sum(d$weights * (y1 - y0)) / sum(d$weights)
  eps <- d$Y - ifelse(d$Z == 1, d$mean1, d$mean0)
  rows <- (d$weights * (d$mean1 - d$mean0 - tau) +
             (2 * d$Z - 1) * (d$weights + K) * eps) / mean(d$weights)
  expect_equal(fit$estimate, tau)
  expect_equal(fit$imputed, cbind(Y0 = y0, Y1 = y1))
  expect_equal(rowSums(fit$loads$incoming) * fit$weight_scale, K)
  expect_equal(sum(K[d$Z == 0]), sum(d$weights[d$Z == 1]))
  expect_equal(sum(K[d$Z == 1]), sum(d$weights[d$Z == 0]))
  expect_equal(fit$contributions$row, rows)
  expect_equal(fit$root_n_variance, mean(rows^2))
  expect_equal(fit$variance, mean(rows^2) / 6)
  expect_equal(fit$numerator_identity$difference, 0, tolerance = 1e-12)
  expect_equal(sum(fit$contributions$actual), 0, tolerance = 1e-12)
})

test_that("PATT needs no treated mean and has the correct two-component rows", {
  d <- wm_fixture()
  d$mean1 <- NULL
  d$estimand <- "PATT"
  fit <- do.call(wm_match, d)
  nn <- wm_fixture_neighbors()
  K <- numeric(6)
  imputed0 <- d$Y
  for (i in 4:6) {
    j <- nn[[i]]
    lambda <- d$weights[j] / sum(d$weights[j])
    imputed0[i] <- d$mean0[i] + sum(lambda * (d$Y[j] - d$mean0[j]))
    for (k in seq_along(j)) K[j[k]] <- K[j[k]] + d$weights[i] * lambda[k]
  }
  tau <- sum(d$weights[4:6] * (d$Y[4:6] - imputed0[4:6])) / sum(d$weights[4:6])
  rows <- ifelse(d$Z == 1, d$weights * (d$Y - d$mean0 - tau),
                 -K * (d$Y - d$mean0)) / mean(d$Z * d$weights)
  expect_equal(fit$estimate, tau)
  expect_equal(fit$contributions$row, rows)
  expect_equal(fit$root_n_variance, mean(rows^2))
  expect_true(all(vapply(fit$graph$neighbors[1:3], is.null, logical(1))))
  expect_equal(fit$numerator_identity$difference, 0, tolerance = 1e-12)
})

test_that("stabilized PATE uses the positive row and unordered-edge formula", {
  d <- wm_fixture()
  d$method <- "stabilized"
  d$rho0 <- c(2, 2.5, 3, 2.2, 2.8, 3.2)
  d$rho1 <- c(4, 4.5, 5, 4.2, 4.8, 5.2)
  fit <- do.call(wm_match, d)
  nn <- wm_fixture_neighbors()
  rho <- ifelse(d$Z == 1, d$rho1, d$rho0)
  mu <- ifelse(d$Z == 1, d$mean1, d$mean0)
  eta <- d$weights * (d$Y - mu) / rho
  C <- R <- numeric(6)
  y0 <- y1 <- d$Y
  for (i in 1:6) {
    j <- nn[[i]]
    prediction <- if (d$Z[i] == 1) d$mean0[i] else d$mean1[i]
    prediction <- prediction + mean(eta[j])
    if (d$Z[i] == 1) y0[i] <- prediction else y1[i] <- prediction
    for (donor in j) {
      C[donor] <- C[donor] + d$weights[i] / 2
      R[donor] <- R[donor] + rho[i] / 2
    }
  }
  tau <- sum(d$weights * (y1 - y0)) / sum(d$weights)
  phi <- d$weights * (d$mean1 - d$mean0 - tau) +
    (2 * d$Z - 1) * (rho + R) * eta
  H <- numeric()
  for (i in 1:3) for (j in 4:6) {
    i_to_j <- j %in% nn[[i]]
    j_to_i <- i %in% nn[[j]]
    if (i_to_j || j_to_i) {
      H[paste(i, j, sep = ":")] <-
        ((-eta[i] * (d$weights[j] - rho[j]) * j_to_i) +
           eta[j] * (d$weights[i] - rho[i]) * i_to_j) / 2
    }
  }
  actual <- d$weights * (d$mean1 - d$mean0 - tau) +
    (2 * d$Z - 1) * (rho + C) * eta
  key <- paste(fit$graph$unordered$i, fit$graph$unordered$j, sep = ":")
  expect_equal(fit$estimate, tau)
  expect_equal(fit$contributions$row, phi / mean(d$weights))
  expect_equal(fit$contributions$edge, unname(H[key]) / mean(d$weights))
  expect_equal(fit$contributions$actual, actual / mean(d$weights))
  expect_equal(sum(phi) + sum(H), sum(actual), tolerance = 1e-12)
  expected <- (sum(phi^2) + sum(H^2)) / (6 * mean(d$weights)^2)
  expect_equal(fit$root_n_variance, expected)
  expect_gt(abs(fit$root_n_variance - mean(fit$contributions$actual^2)), 1e-4)
  expect_equal(fit$numerator_identity$row_edge_sum, 0, tolerance = 1e-12)
  expect_equal(fit$numerator_identity$difference, 0, tolerance = 1e-12)
  expect_equal(nrow(fit$graph$unordered), length(H))
  expect_identical(fit$graph$unordered$i_to_j,
                   vapply(seq_len(nrow(fit$graph$unordered)), function(k)
                     fit$graph$unordered$j[k] %in% nn[[fit$graph$unordered$i[k]]], logical(1)))
  expect_identical(fit$graph$unordered$j_to_i,
                   vapply(seq_len(nrow(fit$graph$unordered)), function(k)
                     fit$graph$unordered$i[k] %in% nn[[fit$graph$unordered$j[k]]], logical(1)))
})

test_that("stabilized PATT squares its treated and control components", {
  d <- wm_fixture()
  d$mean1 <- NULL
  d$estimand <- "PATT"
  d$method <- "stabilized"
  d$rho0 <- c(2, 2.5, 3, 2.2, 2.8, 3.2)
  fit <- do.call(wm_match, d)
  nn <- wm_fixture_neighbors()
  eta <- d$weights * (d$Y - d$mean0) / d$rho0
  missing0 <- vapply(4:6, function(i) d$mean0[i] + mean(eta[nn[[i]]]), numeric(1))
  tau <- sum(d$weights[4:6] * (d$Y[4:6] - missing0)) / sum(d$weights[4:6])
  C <- numeric(6)
  for (i in 4:6) for (j in nn[[i]]) C[j] <- C[j] + d$weights[i] / 2
  rows <- ifelse(d$Z == 1, d$weights * (d$Y - d$mean0 - tau), -C * eta) /
    mean(d$Z * d$weights)
  expect_equal(fit$estimate, tau)
  expect_equal(fit$contributions$row, rows)
  expect_length(fit$contributions$edge, 0)
  expect_equal(fit$root_n_variance, mean(rows^2))
})

test_that("common rescaling preserves estimates and inference in both rules", {
  for (method in c("self_normalized", "stabilized")) {
    d <- wm_fixture()
    d$method <- method
    if (method == "stabilized") {
      d$rho0 <- rep(3, 6)
      d$rho1 <- rep(4, 6)
    }
    fit <- do.call(wm_match, d)
    d$weights <- d$weights * 1e100
    if (method == "stabilized") {
      d$rho0 <- d$rho0 * 1e100
      d$rho1 <- d$rho1 * 1e100
    }
    scaled <- do.call(wm_match, d)
    expect_equal(scaled$estimate, fit$estimate)
    expect_equal(scaled$root_n_variance, fit$root_n_variance)
    expect_equal(scaled$contributions, fit$contributions)
    expect_equal(scaled$graph$neighbors, fit$graph$neighbors)
  }
})

test_that("ties, arbitrary score dimension, and cell restrictions are explicit", {
  d <- wm_fixture()
  d$M <- 1L
  fit <- do.call(wm_match, d)
  expect_identical(fit$graph$neighbors[[4]], 1L)
  expect_identical(fit$graph$neighbors[[5]], 2L)
  expect_equal(fit$graph$edges$share, rep(1, nrow(fit$graph$edges)))
  d$M <- 2L
  d$scores1 <- cbind(d$scores0, 2 * d$scores0, -d$scores0)
  fit <- do.call(wm_match, d)
  expect_identical(fit$graph$neighbors, wm_fixture_neighbors())
  expect_equal(ncol(fit$graph$scores1), 3)
  restricted <- wm_match(1:8, rep(0:1, 4), rep(1, 8), matrix(1:8), M = 1,
                         variance = FALSE, fold_id = rep(1:2, each = 4),
                         strata = rep(c("a", "a", "b", "b"), 2))
  expect_identical(restricted$graph$neighbors,
                   lapply(c(2L, 1L, 4L, 3L, 6L, 5L, 8L, 7L), identity))
  expect_error(wm_match(1:8, rep(0:1, 4), rep(1, 8), matrix(1:8), M = 2,
                         variance = FALSE, fold_id = rep(1:2, each = 4),
                         strata = rep(c("a", "a", "b", "b"), 2)), "Insufficient donors")
})

test_that("raw point estimates and normalized nuisance corrections are distinct", {
  d <- wm_fixture()
  raw <- d
  raw$mean0 <- raw$mean1 <- NULL
  raw$variance <- FALSE
  fit_raw <- do.call(wm_match, raw)
  fit <- do.call(wm_match, d)
  expect_equal(fit$raw_estimate, fit_raw$estimate)
  expect_true(is.na(fit_raw$root_n_variance))
  expect_identical(fit_raw$info$inference_status, "point_estimate_only")
  U <- cbind(c(1, -2, 3, 0, 1, -1), c(2, 1, 0, -2, 4, 1))
  g <- c(0.3, -0.2)
  d$nuisance_influence <- U
  d$sensitivity <- g
  adjusted <- do.call(wm_match, d)
  rows <- fit$contributions$row + as.vector(U %*% g)
  rows <- rows - mean(rows)
  expect_equal(adjusted$estimate, fit$estimate)
  expect_equal(adjusted$contributions$row, rows)
  expect_equal(adjusted$root_n_variance, mean(rows^2))
  expect_equal(adjusted$contributions$actual, fit$contributions$actual)
})

test_that("restriction labels preserve exact pair identity without string collisions", {
  args <- list(Y = 1:4, Z = c(0, 1, 0, 1), weights = rep(1, 4),
               scores0 = matrix(c(0, 100, 100, 0)), M = 1, variance = FALSE)
  args$fold_id <- c("a.b", "a.b", "a", "a")
  args$strata <- c("c", "c", "b.c", "b.c")
  fit <- do.call(wm_match, args)
  expected <- lapply(c(2L, 1L, 4L, 3L), identity)
  expect_identical(fit$graph$neighbors, expected)
  expect_equal(length(unique(fit$graph$cell)), 2)
  args$fold_id <- c(1, 1, 1 + .Machine$double.eps, 1 + .Machine$double.eps)
  args$strata <- NULL
  fit <- do.call(wm_match, args)
  expect_identical(fit$graph$neighbors, expected)
  expect_equal(length(unique(fit$graph$cell)), 2)
  args$strata <- args$fold_id
  args$fold_id <- NULL
  expect_identical(do.call(wm_match, args)$graph$neighbors, expected)
})

test_that("invalid inputs fail before silently changing the estimator", {
  d <- wm_fixture()
  bad <- d; bad$weights[1] <- 0
  expect_error(do.call(wm_match, bad), "strictly positive")
  bad <- d; bad$Y[1] <- Inf
  expect_error(do.call(wm_match, bad), "finite numeric")
  bad <- d; bad$M <- 1.5
  expect_error(do.call(wm_match, bad), "positive integer")
  bad <- d; bad$M <- 1 + 1i
  expect_error(do.call(wm_match, bad), "positive integer")
  bad <- d; bad$scores0 <- as.vector(d$scores0)
  expect_error(do.call(wm_match, bad), "numeric matrix")
  bad <- d; bad$mean1 <- NULL
  expect_error(do.call(wm_match, bad), "both mean0 and mean1")
  bad <- d; bad$mean0 <- bad$mean1 <- NULL
  expect_error(do.call(wm_match, bad), "requires the relevant")
  bad <- d; bad$estimand <- "PATT"
  expect_error(do.call(wm_match, bad), "omit arm-1")
  bad <- d; bad$method <- "stabilized"
  expect_error(do.call(wm_match, bad), "requires rho")
  bad$rho0 <- bad$rho1 <- rep(3, 6)
  bad$scores1 <- cbind(d$scores0, d$scores0)
  expect_error(do.call(wm_match, bad), "same supplied score")
  bad <- d; bad$weights <- c(1e-300, 1e300, rep(1, 4))
  expect_error(do.call(wm_match, bad), "underflowed")
  bad <- d; bad$method <- "stabilized"; bad$weights <- rep(1e-300, 6)
  bad$rho0 <- bad$rho1 <- rep(1e300, 6)
  expect_error(do.call(wm_match, bad), "Normalized rho")
  bad <- d; bad$nuisance_influence <- matrix(1, 6, 2)
  bad$sensitivity <- 1
  expect_error(do.call(wm_match, bad), "length 2")
  bad <- d; bad$fold_id <- c(NA, rep(1, 5))
  expect_error(do.call(wm_match, bad), "nonmissing cell")
})

test_that("distance calculation avoids avoidable squared-distance overflow", {
  d <- wm_fixture()
  d$scores0 <- d$scores0 * 1e200
  fit <- do.call(wm_match, d)
  expect_true(all(is.finite(fit$graph$edges$distance)))
  expect_identical(lapply(fit$graph$neighbors, sort), lapply(wm_fixture_neighbors(), sort))
})
