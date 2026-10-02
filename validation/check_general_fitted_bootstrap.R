# Args: package root, optional receipt RDS, optional isolated installed library,
# optional explicitly supplied original-pilot results directory.
# Default: curated public synthetic fixture in validation/fixtures; no private
# workspace required. No new fit, graph, data or RNG: reparameterize saved
# full-X coefficients and assemble influence on genuine retained d6 graphs.
# An additive saved scalar/full-X section uses the same common public API,
# literal stored predictions/complete PS+prediction+scale influences/counts,
# and actual d1 graphs. No subject/graph/model is regenerated.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1L] else ".", mustWork = TRUE)
installed <- length(args) >= 3L && nzchar(args[3L])
if (installed) {
  lib <- normalizePath(args[3L], mustWork = TRUE)
  ns <- loadNamespace("wdsmatch", lib.loc = lib)
  stopifnot(identical(normalizePath(getNamespaceInfo(ns, "path")),
    normalizePath(file.path(lib, "wdsmatch"))))
  wm_bootstrap <- getExportedValue("wdsmatch", "wm_bootstrap")
  wm_fitted_inference <- getExportedValue("wdsmatch", "wm_fitted_inference")
  source_api <- new.env(parent = .GlobalEnv)
  for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f, local = source_api)
  stopifnot(identical(formals(wm_bootstrap), formals(source_api$wm_bootstrap)),
    identical(body(wm_bootstrap), body(source_api$wm_bootstrap)),
    identical(formals(wm_fitted_inference), formals(source_api$wm_fitted_inference)),
    identical(body(wm_fitted_inference), body(source_api$wm_fitted_inference)))
} else for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
original_pilot <- length(args) >= 4L && nzchar(args[4L])
if (original_pilot) {
  pilot <- normalizePath(args[4L], mustWork = TRUE)
  inputs <- lapply(c("GoodOverlap_retrospective_CorCor_PATE",
    "GoodOverlap_prospective_MisCor_PATT"), function(directory) {
      paths <- file.path(pilot, directory, c("data.rds", "weighted.rds"))
      list(label = directory, data = readRDS(paths[1L]),
        saved = readRDS(paths[2L]), reference = NULL, source_paths = paths)
    })
} else {
  fixture_path <- file.path(root, "validation", "fixtures", "general_fitted_full_x_saved.rds")
  fixture <- readRDS(fixture_path)
  stopifnot(identical(fixture$schema, "wm_saved_full_x_v1"),
    identical(fixture$synthetic, TRUE), length(fixture$cases) == 2L)
  inputs <- fixture$cases
}
had_rng <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
if (had_rng) rng_before <- get(".Random.seed", .GlobalEnv)
labels <- character(); objects <- list(); source_pins <- character()
check <- function(ok, label) {
  if (!isTRUE(ok)) stop(label, call. = FALSE)
  labels <<- c(labels, label)
}
near <- function(x, y, label, tolerance = 2e-10) {
  check(identical(dim(x), dim(y)) && length(x) == length(y) &&
    all(is.finite(c(x, y))) && all(abs(x - y) <= tolerance * pmax(1, abs(x), abs(y))), label)
}
fails <- function(expression, label) check(inherits(try(expression, silent = TRUE), "try-error"), label)
forbidden <- function(...) stop("Forbidden fitting, graph building or random draw.")
wm_match <- wm_fit <- wm_wdsm_fit <- wdsm_fit_ps <- wm_scalar_logistic_match <- forbidden
lm <- glm <- runif <- rnorm <- rmultinom <- sample <- forbidden
start <- proc.time()[["elapsed"]]
quadratic <- function(X) {
  B <- cbind(intercept = 1, X)
  for (j in 1:6) for (k in j:6) {
    B <- cbind(B, X[, j] * X[, k])
    colnames(B)[ncol(B)] <- paste0(colnames(X)[j], ":", colnames(X)[k])
  }
  B
}
raw_coefficients <- function(theta, center, scale) {
  out <- theta; out[] <- 0
  out[1] <- theta[1] - sum(theta[2:7] * center / scale)
  out[2:7] <- theta[2:7] / scale
  index <- 7L
  for (j in 1:6) for (k in j:6) {
    index <- index + 1L
    q <- theta[index] / (scale[j] * scale[k])
    out[index] <- q
    out[1] <- out[1] + q * center[j] * center[k]
    out[j + 1L] <- out[j + 1L] - q * center[k]
    out[k + 1L] <- out[k + 1L] - q * center[j]
  }
  out
}
if (!original_pilot) source_pins["validation/fixtures/general_fitted_full_x_saved.rds"] <-
  digest::digest(fixture_path, file = TRUE, algo = "sha256", serialize = FALSE)
for (input in inputs) {
  directory <- input$label
  if (original_pilot) source_pins[input$source_paths] <- vapply(input$source_paths,
    function(p) digest::digest(p, file = TRUE, algo = "sha256", serialize = FALSE), character(1))
  data <- input$data; saved <- input$saved
  check(saved$bases$X6_means$ok, paste(directory, "saved correction available"))
  prediction <- saved$bases$X6_means$value
  X <- as.matrix(data[paste0("X", 1:6)])
  center <- colMeans(X); scale <- sqrt(colMeans(sweep(X, 2, center, "-")^2))
  Braw <- quadratic(X)
  Bstd <- quadratic(sweep(sweep(X, 2, center, "-"), 2, scale, "/"))
  n <- nrow(X); arms <- if (saved$estimand == "PATE") 0:1 else 0L
  parameters <- unlist(lapply(arms, function(z) paste0("mean", z, ":", colnames(Braw))))
  influence <- D0 <- matrix(0, n, length(parameters), dimnames = list(NULL, parameters))
  D1 <- if (length(arms) == 2L) D0 else NULL
  for (arm in arms) {
    mu <- prediction[[paste0("mean", arm)]]
    theta <- prediction$coefficients[[paste0("arm", arm)]]
    beta <- raw_coefficients(theta, center, scale)
    near(as.vector(Bstd %*% theta), mu, paste(directory, arm, "saved standardized predictions"))
    near(as.vector(Braw %*% beta), mu, paste(directory, arm, "exact raw-X reparameterization"))
    rows <- data$A == arm
    w <- data$survey_weight / max(data$survey_weight)
    residual <- data$Y - mu
    score <- Braw * (rows * w * residual)
    H <- crossprod(Braw, Braw * (rows * w)) / n
    check(max(abs(colMeans(score))) < 1e-8 * max(1, max(abs(H))),
      paste(directory, arm, "stored weighted normal equations"))
    index <- match(paste0("mean", arm, ":", colnames(Braw)), parameters)
    influence[, index] <- score %*% solve(H)
    if (arm == 0L) D0[, index] <- Braw else D1[, index] <- Braw
  }
  if (!original_pilot) check(identical(influence, input$reference$influence) &&
    identical(D0, input$reference$D0) && identical(D1, input$reference$D1),
    paste(directory, "literal saved influence and derivative arrays"))
  counts <- matrix(1, n, 4L)
  i <- which(data$A == 0L)[1]; j <- which(data$A == 1L)[1]
  counts[i, 2] <- 2; counts[j, 2] <- 0
  counts[i, 3] <- 0; counts[j, 3] <- 2
  counts[, 4] <- 0; counts[1, 4] <- n
  for (M in c(1L, 3L, 5L)) {
    key <- paste0("WM_X6_M", M); check(saved$points[[key]]$ok, paste(key, "saved point available"))
    fit <- saved$points[[key]]$value; before <- serialize(fit, NULL)
    label <- paste(directory, key)
    check(fit$M == M && ncol(fit$graph$scores0) == 6L, paste(label, "actual saved graph dimensions"))
    near(fit$predictions$mean0, prediction$mean0, paste(label, "same original control predictions"))
    if (saved$estimand == "PATE") near(fit$predictions$mean1, prediction$mean1,
      paste(label, "same original treated predictions"))
    infer <- wm_fitted_inference(fit, influence, D0, D1,
      covariance_scope = if (saved$estimand == "PATE") "full_x" else "patt",
      transport0 = if (saved$estimand == "PATT") list(mode = "zero", basis = "current_centering") else NULL)
    reference_ok <- original_pilot || all(vapply(names(input$reference$metrics[[key]]), function(name) {
      x <- infer[[name]]; y <- input$reference$metrics[[key]][[name]]
      identical(dim(x), dim(y)) && length(x) == length(y) &&
        all(is.finite(c(x, y))) && all(abs(x - y) <= 2e-10 * pmax(1, abs(x), abs(y)))
    }, logical(1)))
    check(infer$available && identical(infer$fit, fit) && reference_ok,
      paste(label, "genuine public constructor retains fit and saved metrics"))
    # Independent incoming loads and complete corrected same-row formula.
    edges <- fit$graph$edges; q <- edges$query; d <- edges$donor
    W <- fit$weights; Z <- fit$data$Z; Y <- fit$data$Y
    donor_total <- numeric(n)
    tab <- rowsum(W[d], q, reorder = FALSE)
    donor_total[as.integer(rownames(tab))] <- as.vector(tab)
    K <- numeric(n)
    tab <- rowsum(W[q] * W[d] / donor_total[q], d, reorder = FALSE)
    K[as.integer(rownames(tab))] <- as.vector(tab)
    mu0 <- fit$predictions$mean0; mu1 <- fit$predictions$mean1
    gamma <- mean(if (saved$estimand == "PATE") W else Z * W)
    if (saved$estimand == "PATE") {
      L <- (W * (mu1 - mu0 - fit$estimate) +
        (2 * Z - 1) * (W + K) * (Y - ifelse(Z == 1, mu1, mu0))) / gamma
      slope <- colSums(D0 * ((1 - Z) * K - Z * W) +
        D1 * ((1 - Z) * W - Z * K)) / sum(W)
    } else {
      L <- (Z * W * (Y - mu0 - fit$estimate) - (1 - Z) * K * (Y - mu0)) / gamma
      slope <- colSums(D0 * ((1 - Z) * K - Z * W)) / sum(Z * W)
    }
    L <- L - mean(L); U <- sweep(influence, 2, colMeans(influence), "-")
    complete <- L + as.vector(U %*% slope); complete <- complete - mean(complete)
    near(infer$base_rows, L, paste(label, "independent complete base rows"))
    near(infer$total_sensitivity, slope, paste(label, "independent full-X prediction slope"))
    near(infer$augmented_rows, complete, paste(label, "independent augmented row formula"))
    out <- wm_bootstrap(infer, counts = counts, conf.level = .9)
    roots <- as.vector(crossprod(complete, counts - 1)) / sqrt(n)
    near(out$root_n_draws, roots, paste(label, "literal full-n count roots"))
    near(out$draws, fit$estimate + roots / sqrt(n), paste(label, "unchanged point and total-n estimates"))
    near(out$conditional_root_n_variance, mean(complete^2), paste(label, "analytic row variance"))
    near(out$conditional_variance, mean(complete^2) / n, paste(label, "full observed n variance"))
    near(out$monte_carlo_root_n_variance, stats::var(roots), paste(label, "separate B-minus-one diagnostic"))
    near(out$conf.int, fit$estimate + c(-1, 1) * stats::qnorm(.95) * sqrt(mean(complete^2) / n),
      paste(label, "analytic normal interval"))
    check(!isTRUE(out$population_assumptions_verified) && !out$original_refit_bootstrap &&
      !out$supplied_draw_law_verified && identical(out$source_inference, infer),
      paste(label, "scientific premises unverified and original source retained"))
    labeled <- counts; rownames(labeled) <- as.character(seq_len(n))
    near(wm_bootstrap(infer, counts = labeled, B = 4L, chunk_size = 37L)$root_n_draws,
      roots, paste(label, "row-index binding and bounded chunks"))
    singleton <- wm_bootstrap(infer, counts = counts[, 4, drop = FALSE], interval = "none")
    check(singleton$B == 1L && is.na(singleton$monte_carlo_variance) && is.null(singleton$conf.int),
      paste(label, "collapsed arm and B1 contribution remains defined"))
    check(identical(serialize(fit, NULL), before), paste(label, "saved graph and point unchanged"))
    objects[[label]] <- infer
  }
}
object <- objects[[1]]; n <- object$n
counts <- matrix(1, n, 2L)
guard <- function(change, label) {bad <- object; bad <- change(bad); fails(wm_bootstrap(bad, counts = counts), label)}
for (flag in c("available", "numerically_available", "conditional_inference_available"))
  guard(function(x) {x[[flag]] <- FALSE; x}, paste("reject unavailable", flag))
guard(function(x) {x$status <- "transport_unavailable"; x}, "reject unavailable status")
guard(function(x) {x$reciprocal_subtraction <- 1; x}, "reject reciprocal correction")
guard(function(x) {x$estimate <- x$estimate + 1; x}, "reject changed point")
guard(function(x) {x$n <- x$n - 1L; x}, "reject changed total n")
guard(function(x) {x$M <- 99L; x}, "reject changed M")
guard(function(x) {x$dimensions[] <- 2L; x}, "reject dimension label only")
guard(function(x) {x$fit$data$Y[1] <- x$fit$data$Y[1] + 1; x}, "reject changed source outcomes")
guard(function(x) {x$parameter_names <- rev(x$parameter_names); x}, "reject parameter order")
guard(function(x) {colnames(x$nuisance_influence) <- rev(x$parameter_names); x}, "reject influence axis")
guard(function(x) {rownames(x$nuisance_influence) <- rev(as.character(seq_len(n))); x}, "reject influence row reorder")
guard(function(x) {x$total_sensitivity[1] <- x$total_sensitivity[1] + 1; x}, "reject missing full slope")
guard(function(x) {x$base_rows[1] <- x$base_rows[1] + 1; x}, "reject stale base rows")
guard(function(x) {x$augmented_rows[1] <- x$augmented_rows[1] + 1; x}, "reject stale complete rows")
guard(function(x) {x$Sigma[1, 1] <- x$Sigma[1, 1] + 1; x}, "reject changed covariance block")
guard(function(x) {x$cross_term <- x$cross_term + 1; x}, "reject omitted covariance term")
guard(function(x) {x$root_n_variance <- x$root_n_variance + 1; x}, "reject changed variance")
fake <- structure(object[c("n", "estimate", "augmented_rows", "root_n_variance", "covariance_scope")],
  class = c("wm_fitted_inference", "list"))
fails(wm_bootstrap(fake, counts = counts), "reject fabricated row-only fitted label")
counts <- matrix(1, n, 2)
fails(wm_bootstrap(object, counts = counts, seed = 1L), "reject seed with frozen counts")
fails(wm_bootstrap(object, counts = counts, B = 3L), "reject B mismatch")
badcounts <- counts; badcounts[1, 1] <- .5
fails(wm_bootstrap(object, counts = badcounts), "reject fractional counts")
badcounts <- counts; rownames(badcounts) <- rev(as.character(seq_len(n)))
fails(wm_bootstrap(object, counts = badcounts), "reject count row reordering")
fails(wm_bootstrap(object, counts = counts, refit = forbidden), "reject contribution refit")
fails(wm_bootstrap(object, method = "fixed_reuse", counts = counts), "reject inferred object in fixed-reuse")
fails(wm_bootstrap(object, B = 0L), "Gaussian proxy rejects B before random draw")
check(identical(had_rng, exists(".Random.seed", .GlobalEnv, inherits = FALSE)) &&
  (!had_rng || identical(rng_before, get(".Random.seed", .GlobalEnv))), "RNG state unchanged")
# Additive checks on genuine saved scalar propensity graphs with full-X correction.
# Existing d6 cases/checks above are retained literally. No new fit or graph.
if (!original_pilot) {
  check(identical(fixture$scalar_design_extension, "wm_saved_scalar_design_full_x_v1") &&
    length(fixture$scalar_design_cases) == 2L, "saved scalar full-X extension schema")
  if (installed) {
    for (name in c(".wm_wdsm_fitted_variance", ".wm_general_fitted_rows")) {
      check(identical(formals(get(name, envir = ns)), formals(source_api[[name]])) &&
        identical(body(get(name, envir = ns)), body(source_api[[name]])),
        paste(name, "installed generic helper matches source"))
    }
  }
  scalar_objects <- list(); scalar_arguments <- list()
  for (input in fixture$scalar_design_cases) {
    label <- input$label; data <- input$data; saved <- input$saved
    n <- nrow(data); pate <- saved$estimand == "PATE"
    X <- as.matrix(data[paste0("X", 1:6)])
    Bpg <- cbind("(Intercept)" = 1, X, "X1:X2" = data$X1 * data$X2)
    Bps <- input$model$ps_design
    expected_design <- if (ncol(Bps) == 8L) Bpg else cbind("(Intercept)" = 1, X)
    check(identical(unname(Bps), unname(expected_design)) &&
      identical(colnames(Bps), colnames(expected_design)), paste(label, "literal known logistic design"))
    W <- data$survey_weight; Z <- data$A; Y <- data$Y
    Wps <- if (input$model$ps_weight_rule == "supplied_weights") W / max(W) else rep(1, n)
    check(input$model$ps_weight_rule %in% c("supplied_weights", "unit") &&
      identical(Wps, input$model$ps_fit_weights), paste(label, "actual original PS weight rule"))
    probability <- input$model$probability
    check(identical(probability, saved$bases$PS$value$probability), paste(label, "literal stored PS probability"))
    equation <- Bps * (Wps * (Z - probability))
    Hps <- crossprod(Bps, Bps * (Wps * probability * (1 - probability))) / n
    Ips <- equation %*% solve(Hps)
    De <- Bps * (probability * (1 - probability))
    center <- mean(probability); scale <- sqrt(mean((probability - center)^2))
    Ic <- probability - center + as.vector(Ips %*% colMeans(De))
    Is <- ((probability - center)^2 - scale^2) / (2 * scale) +
      as.vector(Ips %*% (colMeans(De * (probability - center)) / scale))
    check(max(abs(colMeans(equation))) < 1e-12, paste(label, "stored PS score equation"))
    near(as.vector(saved$bases$PS$value$scores), (probability - center) / scale,
      paste(label, "same original scalar standardization"))
    arms <- if (pate) 0:1 else 0L
    parameters <- c(paste0("ps:", colnames(Bps)), unlist(lapply(arms,
      function(z) paste0("mean", z, ":", colnames(Bpg)))), "ps:center", "ps:scale")
    U <- D0 <- matrix(0, n, length(parameters), dimnames = list(NULL, parameters))
    D1 <- if (pate) D0 else NULL
    U[, seq_len(ncol(Bps))] <- Ips
    U[, length(parameters) - 1L] <- Ic; U[, length(parameters)] <- Is
    prediction <- saved$bases$PS_means$value
    for (arm in arms) {
      mu <- prediction[[paste0("mean", arm)]]
      beta <- prediction$coefficients[[paste0("arm", arm)]]
      check(identical(as.vector(Bpg %*% beta), mu), paste(label, arm, "literal saved full-X predictions"))
      w <- W / max(W)
      score <- Bpg * ((Z == arm) * w * (Y - mu))
      H <- crossprod(Bpg, Bpg * ((Z == arm) * w)) / n
      check(max(abs(colMeans(score))) < 1e-12, paste(label, arm, "stored full-X normal equations"))
      index <- match(paste0("mean", arm, ":", colnames(Bpg)), parameters)
      U[, index] <- score %*% solve(H)
      if (arm == 0L) D0[, index] <- Bpg else D1[, index] <- Bpg
    }
    check(identical(U, input$reference$influence) && identical(D0, input$reference$D0) &&
      identical(D1, input$reference$D1), paste(label, "literal complete PS/prediction/scale IF and derivatives"))
    counts <- input$counts
    check(is.matrix(counts) && nrow(counts) == n && ncol(counts) == 200L &&
      all(counts >= 0 & counts == floor(counts)) && all(colSums(counts) == n),
      paste(label, "original 200 full-n count arrays"))
    for (M in c(1L, 3L, 5L)) {
      key <- paste0("WM_PS_M", M); fit <- saved$points[[key]]$value
      before <- serialize(fit, NULL); item <- paste(label, key)
      check(saved$points[[key]]$ok && fit$M == M && ncol(fit$graph$scores0) == 1L &&
        (!pate || ncol(fit$graph$scores1) == 1L), paste(item, "actual saved scalar graph/M"))
      check(identical(fit$data$Y, Y) && identical(fit$data$Z, Z) &&
        identical(fit$weights, W) &&
        identical(unname(fit$graph$scores0), unname(saved$bases$PS$value$scores)) &&
        (!pate || identical(unname(fit$graph$scores1), unname(saved$bases$PS$value$scores))),
        paste(item, "literal original data/weight/scalar map binding"))
      arguments <- list(fit = fit, nuisance_influence = U, mean_derivative0 = D0,
        mean_derivative1 = D1, covariance_scope = if (pate) "full_x" else "patt",
        transport0 = if (pate) NULL else list(mode = "zero", basis = "current_centering"))
      infer <- do.call(wm_fitted_inference, arguments)
      check(infer$available && identical(infer$fit, fit) &&
        identical(infer$estimate, fit$estimate) && infer$n == n,
        paste(item, "common fitted API preserves scalar point/graph/full n"))
      check(identical(infer$nuisance_influence, input$reference$centered_influence) &&
        all(vapply(names(input$reference$metrics[[key]]), function(k)
          identical(infer[[k]], input$reference$metrics[[key]][[k]]), logical(1))),
        paste(item, "all nine saved scalar variance/slope components exact"))
      edges <- fit$graph$edges; q <- edges$query; d <- edges$donor
      donor_total <- numeric(n); tab <- rowsum(W[d], q, reorder = FALSE)
      donor_total[as.integer(rownames(tab))] <- as.vector(tab)
      K <- numeric(n); tab <- rowsum(W[q] * W[d] / donor_total[q], d, reorder = FALSE)
      K[as.integer(rownames(tab))] <- as.vector(tab)
      near(edges$share, W[d] / donor_total[q], paste(item, "independent finite-M donor fractions"))
      near(fit$loads$incoming[cbind(seq_len(n), Z + 1L)] * fit$weight_scale, K,
        paste(item, "independent literal incoming loads"))
      mu0 <- prediction$mean0; mu1 <- prediction$mean1
      gamma <- mean(if (pate) W else Z * W)
      if (pate) {
        L <- (W * (mu1 - mu0 - fit$estimate) + (2 * Z - 1) * (W + K) *
          (Y - ifelse(Z == 1, mu1, mu0))) / gamma
        b <- colSums(D0 * ((1 - Z) * K - Z * W) +
          D1 * ((1 - Z) * W - Z * K)) / sum(W)
      } else {
        L <- (Z * W * (Y - mu0 - fit$estimate) - (1 - Z) * K * (Y - mu0)) / gamma
        b <- colSums(D0 * ((1 - Z) * K - Z * W)) / sum(Z * W)
      }
      L <- L - mean(L); Ui <- sweep(U, 2, colMeans(U), "-")
      complete <- L + as.vector(Ui %*% b); complete <- complete - mean(complete)
      C <- as.vector(crossprod(Ui, L)) / n; Sigma <- crossprod(Ui) / n
      V <- mean(complete^2)
      near(infer$base_rows, L, paste(item, "independent scalar complete base rows"))
      near(infer$total_sensitivity, b, paste(item, "independent scalar full prediction slope"))
      near(infer$augmented_rows, complete, paste(item, "independent scalar complete augmented rows"))
      near(infer$C, C, paste(item, "complete scalar same-row covariance"))
      near(infer$Sigma, Sigma, paste(item, "complete scalar joint covariance"))
      near(infer$root_n_variance, mean(L^2) + 2 * sum(b * C) +
        as.numeric(crossprod(b, Sigma %*% b)), paste(item, "scalar covariance identity"))
      check(all(infer$graph_sensitivity == 0) && all(infer$smooth_derivative$weight_sensitivity == 0) &&
        infer$reciprocal_subtraction == 0 && grepl("design", infer$contract$fitted_law, fixed = TRUE) &&
        !infer$assumptions_verified && !infer$application_verified,
        paste(item, "restricted known-W/design/full-X premises unverified"))
      out <- wm_bootstrap(infer, counts = counts, conf.level = .9)
      roots <- as.vector(crossprod(complete, counts - 1)) / sqrt(n)
      near(roots, input$reference$count_roots[[key]], paste(item, "prior saved scalar root arrays"))
      near(out$root_n_draws, roots, paste(item, "same 200 full-n count roots through common API"))
      near(out$draws, fit$estimate + roots / sqrt(n), paste(item, "unchanged point plus full-n scalar roots"))
      near(out$conditional_root_n_variance, V, paste(item, "analytic scalar row variance"))
      near(out$conditional_variance, V / n, paste(item, "PATE/PATT full observed n variance"))
      near(out$monte_carlo_root_n_variance, stats::var(roots), paste(item, "separate scalar B-minus-one diagnostic"))
      near(out$conf.int, fit$estimate + c(-1, 1) * stats::qnorm(.95) * sqrt(V / n),
        paste(item, "scalar analytic normal interval"))
      check(identical(out$source_inference, infer) && !out$population_assumptions_verified &&
        !out$original_refit_bootstrap && !out$supplied_draw_law_verified,
        paste(item, "no automatic scalar callback/coverage/draw-law validity"))
      zero <- arguments; zero$weight_derivative <- matrix(0, n, ncol(U), dimnames = list(NULL, parameters))
      check(identical(do.call(wm_fitted_inference, zero), infer), paste(item, "exact zero DW compatibility"))
      check(identical(serialize(fit, NULL), before), paste(item, "original scalar fit serialization preserved"))
      scalar_objects[[item]] <- infer; scalar_arguments[[item]] <- arguments
      objects[[item]] <- infer
    }
  }
  argument <- scalar_arguments[[1L]]; scalar <- scalar_objects[[1L]]; n <- scalar$n
  for (scope in c("distinct_rarity", "common_field")) {
    bad <- argument; bad$covariance_scope <- scope
    fails(do.call(wm_fitted_inference, bad), paste("reject scalar", scope))
  }
  bad <- argument; bad$fit$graph$scores1 <- cbind(bad$fit$graph$scores1, bad$fit$graph$scores1)
  fails(do.call(wm_fitted_inference, bad), "reject mixed scalar/higher matching dimensions")
  bad <- argument; bad$weight_derivative <- matrix(0, n, ncol(argument$nuisance_influence),
    dimnames = list(NULL, colnames(argument$nuisance_influence))); bad$weight_derivative[1, 1] <- 1e-100
  fails(do.call(wm_fitted_inference, bad), "reject tiny nonzero scalar supplied-W derivative")
  bad <- argument; bad$transport0 <- list(mode = "zero", basis = "current_centering")
  fails(do.call(wm_fitted_inference, bad), "reject supplied scalar PATE transport")
  patt <- scalar_objects[[grep("patt", names(scalar_objects), fixed = TRUE)[1L]]]
  patt_argument <- scalar_arguments[[grep("patt", names(scalar_arguments), fixed = TRUE)[1L]]]
  for (spec in list(NULL, list(mode = "zero", basis = "fixed_map"), list(mode = "estimate"),
    list(mode = "zero", basis = "current_centering", extra = 1))) {
    bad <- patt_argument; bad$transport0 <- spec
    fails(do.call(wm_fitted_inference, bad), "reject unsupported scalar PATT transport")
  }
  counts <- matrix(1, n, 2L)
  scalar_guard <- function(change, label) {
    bad <- change(scalar); fails(wm_bootstrap(bad, counts = counts), label)
  }
  scalar_guard(function(x) {x$covariance_scope <- "common_field"; x$status <- "conditional_common_field"; x},
    "scalar binder rejects common_field label")
  scalar_guard(function(x) {x$fit$graph$scores1 <- cbind(x$fit$graph$scores1, x$fit$graph$scores1); x$dimensions[2] <- 2L; x},
    "scalar binder rejects mixed dimensions")
  scalar_guard(function(x) {x$smooth_derivative$weight_sensitivity[1] <- 1e-100; x},
    "scalar binder rejects tiny nonzero weight sensitivity")
  scalar_guard(function(x) {x$graph_transport$potential0$mode <- "estimate"; x},
    "scalar binder rejects estimated transport")
  scalar_guard(function(x) {x$augmented_rows[1] <- x$augmented_rows[1] + 1; x}, "scalar binder rejects stale rows")
  scalar_guard(function(x) {x$estimate <- x$estimate + 1; x}, "scalar binder rejects changed point")
  bad <- patt; bad$graph_transport$potential0$basis <- "fixed_map"
  fails(wm_bootstrap(bad, counts = matrix(1, bad$n, 2L)), "scalar binder rejects PATT fixed_map label")
  fails(wm_bootstrap(scalar, counts = counts, seed = 1L), "scalar counts reject seed")
  fails(wm_bootstrap(scalar, counts = counts, B = 3L), "scalar counts reject B mismatch")
  check(identical(had_rng, exists(".Random.seed", .GlobalEnv, inherits = FALSE)) &&
    (!had_rng || identical(rng_before, get(".Random.seed", .GlobalEnv))), "scalar extension RNG state unchanged")
}

receipt <- list(status = "PASS", check_count = length(labels), labels = labels,
  mode = if (installed) "installed namespace" else "source", objects = objects,
  installed_api_formals_and_bodies_match_source = if (installed) TRUE else NA,
  input_mode = if (original_pilot) "explicit original-pilot audit" else "public saved synthetic fixture",
  source_pins = source_pins, elapsed_seconds = proc.time()[["elapsed"]] - start,
  scientific_scope = if (original_pilot) "Saved actual full-X d6 corrected graphs, M1/3/5, PATE full_x and PATT patt; explicit original-pilot audit." else "Saved actual d6 full-X and d1 design-fitted PS graphs with full-X correction, M1/3/5, PATE full_x/PATT patt. Complete same-stack IF, full-n supplied count arithmetic and original points only. No coefficients fitted, graphs built, new data, draws or calibration; population premises, coverage and original-refit validity unverified. Other covariance/transport scopes unexercised.")
if (length(args) >= 2L && nzchar(args[2L])) saveRDS(receipt, args[2L], version = 2)
cat("PASS", length(labels), "saved-data algebra/binding checks; zero fits, graphs, RNG or new data.\n")
