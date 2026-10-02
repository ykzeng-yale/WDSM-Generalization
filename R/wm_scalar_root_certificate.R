# Candidate: exact-data root, donor-set and point-arithmetic certificates.
# GMP rational endpoints; precision controls enclosure width, not validity.

.wm_cert_arithmetic <- function(bits) {
  Q <- gmp::as.bigq
  pow2 <- function(exponent) {
    positive <- exponent >= 0
    out <- Q(rep(1, length(exponent)))
    if (any(positive)) out[positive] <- Q(gmp::as.bigz(2)^exponent[positive])
    if (any(!positive)) out[!positive] <- Q(1, gmp::as.bigz(2)^(-exponent[!positive]))
    out
  }
  # Never use .bigq2mpfr's rnd.mode to claim rational enclosures.
  roundq <- function(x, upper) {
    out <- x
    nz <- which(x != 0)
    if (!length(nz)) return(out)
    a <- abs(x[nz])
    exponent <- nchar(as.character(gmp::numerator(a), b = 2)) -
      nchar(as.character(gmp::denominator(a), b = 2))
    exponent <- exponent - as.integer(a < pow2(exponent))
    shift <- bits - 1L - exponent
    scaled <- a * pow2(shift)
    num <- gmp::numerator(scaled)
    den <- gmp::denominator(scaled)
    integer <- num %/% den
    round_up <- if (upper) x[nz] > 0 else x[nz] < 0
    integer <- integer + as.integer(round_up & (num %% den != 0))
    value <- Q(integer) * pow2(-shift)
    value[x[nz] < 0] <- -value[x[nz] < 0]
    out[nz] <- value
    out
  }
  pack <- function(lo, hi = lo) list(lo = roundq(lo, FALSE), hi = roundq(hi, TRUE))
  exact <- function(x) { x <- Q(x); list(lo = x, hi = x) }
  neg <- function(a) list(lo = -a$hi, hi = -a$lo)
  add <- function(a, b) pack(a$lo + b$lo, a$hi + b$hi)
  sub <- function(a, b) add(a, neg(b))
  low <- function(a, b) { k <- a > b; a[k] <- b[k]; a }
  high <- function(a, b) { k <- a < b; a[k] <- b[k]; a }
  mul <- function(a, b) {
    values <- list(a$lo * b$lo, a$lo * b$hi, a$hi * b$lo, a$hi * b$hi)
    lo <- hi <- values[[1L]]
    for (j in 2:4) { lo <- low(lo, values[[j]]); hi <- high(hi, values[[j]]) }
    pack(lo, hi)
  }
  div <- function(a, b) {
    if (any(b$lo <= 0 & b$hi >= 0)) stop("Certificate division interval includes zero.")
    mul(a, pack(1 / b$hi, 1 / b$lo))
  }
  total <- function(a) pack(sum(a$lo), sum(a$hi))
  absmax <- function(a) high(abs(a$lo), abs(a$hi))
  at <- function(a, i) list(lo = a$lo[i], hi = a$hi[i])
  mp <- function(x) {
    ans <- Rmpfr::.bigq2mpfr(x, precB = bits)
    if (any(!is.finite(ans)) || any(Rmpfr::getPrec(ans) != bits)) {
      stop("Endpoint is not representable at certificate precision.")
    }
    nz <- which(x != 0)
    if (length(nz) && any(ans[nz] == 0)) stop("A nonzero endpoint underflowed during MPFR conversion.")
    if (length(nz)) {
      erange <- Rmpfr::.mpfr_erange()
      exponent <- Rmpfr::frexpMpfr(ans[nz])$e
      # .mpfr2bigq constructs both a p-bit integer and a power-of-two scale.
      shift <- pmax(0, -exponent) + bits
      if (erange[[1L]] > min(exponent - bits, 1L) ||
          erange[[2L]] < max(shift + 1L, exponent + shift)) {
        stop("MPFR exponent range is too narrow for exact reverse conversion.")
      }
    }
    if (any(ans[x == 0] != 0) ||
        (length(nz) && any(Rmpfr::.mpfr2bigq(ans[nz]) != x[nz]))) {
      stop("Endpoint did not roundtrip exactly through MPFR.")
    }
    ans
  }
  exponential <- function(a) {
    erange <- Rmpfr::.mpfr_erange()
    if (erange[[1L]] > -(bits + 2500L) || erange[[2L]] < bits + 2500L ||
        any(a$lo < -1000 | a$hi > 1000)) {
      stop("Exp enclosure requires arguments in [-1000,1000] and a safe MPFR exponent range.")
    }
    evaluate <- function(x, upper) {
      value <- exp(mp(x))
      if (any(!is.finite(value)) || any(value <= 0) ||
          any(Rmpfr::getPrec(value) != bits)) stop("Unexpected MPFR exp result.")
      q <- Rmpfr::.mpfr2bigq(value)
      padding <- abs(q) * pow2(1L - bits)
      if (upper) q + padding else q - padding
    }
    pack(evaluate(a$lo, FALSE), evaluate(a$hi, TRUE))
  }
  logistic <- function(a) div(exact(1), add(exact(1), exponential(neg(a))))
  display <- function(x, upper) {
    # The exact rational string is authoritative; double conversion is outward.
    x <- roundq(x, upper)
    Rmpfr::toNum(mp(x), rnd.mode = if (upper) "U" else "D")
  }
  record <- function(a) {
    lower <- display(a$lo, FALSE); upper <- display(a$hi, TRUE)
    list(lower = lower, upper = upper, numeric_display_finite = all(is.finite(c(lower, upper))),
      lower_exact = as.character(a$lo), upper_exact = as.character(a$hi))
  }
  list(Q = Q, pow2 = pow2, roundq = roundq, pack = pack, exact = exact,
       neg = neg, add = add, sub = sub, mul = mul, div = div, total = total,
       absmax = absmax, at = at, exp = exponential, logistic = logistic,
       display = display, record = record, low = low, high = high)
}

wm_scalar_root_certificate <- function(object, design, center = NULL,
                                       radius = 1e-4, precision = 128L) {
  if (!requireNamespace("gmp", quietly = TRUE) ||
      !requireNamespace("Rmpfr", quietly = TRUE)) {
    stop("Numerical certification requires the optional gmp and Rmpfr packages.", call. = FALSE)
  }
  if (!is.numeric(precision) || length(precision) != 1L || !is.finite(precision) ||
      precision != floor(precision) || precision < 64 || precision > 4096) {
    stop("precision must be an integer from 64 to 4096.", call. = FALSE)
  }
  if (!is.numeric(radius) || length(radius) != 1L || !is.finite(radius) || radius <= 0) {
    stop("radius must be finite and positive.", call. = FALSE)
  }
  if (!inherits(object, "wm_scalar_logistic_match") ||
      !identical(object$status, "point_computed") || !inherits(object$fit, "wm_match")) {
    stop("Supply a computed raw wm_scalar_logistic_match object.", call. = FALSE)
  }
  fit <- object$fit
  n <- fit$n
  Y <- .wm_numeric_vector(fit$data$Y, n, "stored Y")
  Z <- fit$data$Z
  W <- .wm_numeric_vector(fit$weights, n, "stored weights", positive = TRUE)
  probability <- .wm_numeric_vector(object$propensity$probability, n, "stored probability")
  if (any(probability <= 0 | probability >= 1) || !all(Z %in% c(0, 1)) ||
      length(Z) != n || length(unique(Z)) != 2L ||
      !fit$estimand %in% c("PATE", "PATT") ||
      !identical(fit$method, "self_normalized") || !identical(fit$info$corrected, FALSE) ||
      !identical(fit$estimate, object$estimate) || !is.finite(fit$estimate) ||
      length(unique(fit$graph$cell)) != 1L ||
      !identical(as.numeric(fit$graph$scores0), probability) ||
      (fit$estimand == "PATE" && !identical(as.numeric(fit$graph$scores1), probability))) {
    stop("Stored object is not an unrestricted raw same-probability scalar match.", call. = FALSE)
  }
  design <- .wm_score_matrix(design, n, "design")
  p <- ncol(design)
  if (p >= n || any(design[, 1L] != 1)) stop("Supply the fitting design with leading intercept.")
  center_origin <- if (is.null(center)) "QR of stored probability logits; arbitrary proposal" else
    "User supplied binary coefficient vector; arbitrary proposal"
  if (is.null(center)) center <- tryCatch(drop(qr.solve(design, stats::qlogis(probability))),
      error = function(e) stop("Could not propose a center: ", conditionMessage(e), call. = FALSE))
  center <- .wm_numeric_vector(center, p, "center")
  A <- .wm_cert_arithmetic(as.integer(precision))
  Q <- A$Q
  query <- if (fit$estimand == "PATE") seq_len(n) else which(Z == 1L)
  neighbors <- fit$graph$neighbors
  if (length(neighbors) != n || !is.numeric(fit$M) || length(fit$M) != 1L ||
      !is.finite(fit$M) || fit$M < 1 || fit$M != floor(fit$M)) stop("Invalid stored donor graph.")
  for (i in query) {
    j <- neighbors[[i]]
    if (!is.numeric(j) || length(j) != fit$M || anyNA(j) ||
        any(j < 1 | j > n | j != floor(j)) || anyDuplicated(j) || any(Z[j] == Z[i])) {
      stop("Invalid selected donor set at row ", i, ".", call. = FALSE)
    }
  }
  edges <- fit$graph$edges
  if (!is.data.frame(edges) || !all(c("query", "donor", "arm") %in% names(edges)) ||
      !all(vapply(edges[c("query", "donor", "arm")], function(x)
        is.numeric(x) && !is.complex(x) && all(is.finite(x)) && all(x == floor(x)), logical(1))) ||
      !identical(as.numeric(edges$query), as.numeric(rep(query, each = fit$M))) ||
      !identical(as.numeric(edges$donor), as.numeric(unlist(neighbors[query], use.names = FALSE))) ||
      !identical(as.numeric(edges$arm), as.numeric(1L - Z[rep(query, each = fit$M)]))) {
    stop("Stored edge and neighbor representations disagree.", call. = FALSE)
  }
  # Each input value is interpreted exactly as its stored binary64 value.
  binding <- list(design = design, Z = Z, weights = W, probability = probability,
    neighbors = neighbors, edges = edges[c("query", "donor", "arm")],
    Y = Y, estimate = fit$estimate, M = fit$M, estimand = fit$estimand)
  base <- list(status = "not_certified", root_contained = FALSE,
    donor_sets_certified = FALSE, point_root_error_bound = NULL,
    inputs = binding, center = center, center_origin = center_origin,
    radius = radius, precision_bits = as.integer(precision),
    arithmetic = "GMP exact rational operations and exact dyadic outward rounding; padded MPFR exp",
    runtime = list(Rmpfr = as.character(utils::packageVersion("Rmpfr")),
      gmp = as.character(utils::packageVersion("gmp")), MPFR = as.character(Rmpfr::mpfrVersion()),
      MPFR_exponent_range = Rmpfr::.mpfr_erange()),
    population_assumptions_verified = FALSE, asymptotic_success_rate_claimed = FALSE,
    variance_arithmetic_certified = FALSE, bootstrap_validity_claimed = FALSE)
  calculate <- function() {
    # Approximate linear algebra merely proposes a coordinate system.
    eta_double <- drop(design %*% center)
    e_double <- stats::plogis(eta_double)
    H <- crossprod(design, design * ((W / max(W)) * e_double * (1 - e_double))) / n
    C <- backsolve(chol(H), diag(p))
    if (any(!is.finite(C))) stop("Nonfinite preconditioner proposal.")
    L <- lapply(seq_len(p), function(k) Q(design[, k]))
    eta <- Q(rep(0, n))
    R <- vector("list", p)
    for (k in seq_len(p)) {
      eta <- eta + L[[k]] * Q(center[k])
      R[[k]] <- Q(rep(0, n))
      for (j in seq_len(p)) R[[k]] <- R[[k]] + L[[j]] * Q(C[j, k])
    }
    row_norm <- Reduce(`+`, lapply(R, abs)) # exact L1 upper bound on Euclidean row norm
    weight <- Q(W) / Q(max(W))
    eta_interval <- A$pack(eta)
    e0 <- A$logistic(eta_interval)
    gradient <- lapply(R, function(column) A$div(A$total(A$mul(
      A$exact(weight * column), A$sub(e0, A$exact(Z)))), A$exact(n)))
    gradient_bound <- A$roundq(sum(Reduce(c, lapply(gradient, A$absmax))), TRUE)
    maximum_eta <- A$roundq(abs(eta) + Q(radius) * row_norm, TRUE)
    exp_negative <- A$exp(A$neg(A$exact(maximum_eta)))
    denominator <- A$add(A$exact(1), exp_negative)
    derivative <- A$div(exp_negative, A$mul(denominator, denominator))
    G <- vector("list", p * p)
    for (j in seq_len(p)) for (k in seq_len(p)) {
      G[[(j - 1L) * p + k]] <- A$div(A$total(A$mul(
        A$exact(weight * derivative$lo * R[[j]] * R[[k]]), A$exact(1))), A$exact(n))
    }
    row_bounds <- Q(rep(0, p))
    for (j in seq_len(p)) {
      off <- Q(0)
      for (k in setdiff(seq_len(p), j)) off <- off + A$absmax(G[[(j - 1L) * p + k]])
      row_bounds[j] <- G[[(j - 1L) * p + j]]$lo - off
    }
    curvature <- A$roundq(min(row_bounds), FALSE)
    base$preconditioner <- C
    base$gradient_norm_upper <- A$record(A$exact(gradient_bound))
    base$curvature_lower <- A$record(A$exact(curvature))
    base$boundary_margin <- A$record(A$pack(curvature * Q(radius) - gradient_bound))
    if (curvature <= 0 || gradient_bound >= curvature * Q(radius)) {
      base$reason <- "Verified curvature/strict outward-gradient condition did not pass."
      return(base)
    }
    delta <- A$roundq(gradient_bound / curvature, TRUE)
    base$root_contained <- TRUE
    base$status <- "root_contained_graph_unresolved"
    base$transformed_root_error_upper <- A$record(A$exact(delta))
    base$coefficient_error_upper <- A$record(A$pack(delta * sum(abs(Q(C)))))
    lip <- A$roundq(row_norm * delta / 4, TRUE)
    root_probability <- A$pack(A$low(A$high(e0$lo - lip, Q(rep(0, n))), Q(rep(1, n))),
      A$high(A$low(e0$hi + lip, Q(rep(1, n))), Q(rep(0, n))))
    pq <- Q(probability)
    epsilon <- max(A$absmax(A$sub(A$exact(pq), root_probability)))
    base$root_probability <- A$record(root_probability)
    base$stored_probability_error_upper <- A$record(A$exact(epsilon))
    orders <- lapply(0:1, function(z) {
      rows <- which(Z == z); rows[order(probability[rows], rows, method = "radix")]
    })
    minimum_gap <- NULL
    unresolved <- integer()
    for (i in query) {
      selected <- neighbors[[i]]
      donors <- orders[[2L - Z[i]]]
      if (length(selected) == length(donors)) next
      split <- findInterval(probability[i], probability[donors])
      left <- split
      right <- split + 1L
      while (left > 0L && donors[left] %in% selected) left <- left - 1L
      while (right <= length(donors) && donors[right] %in% selected) right <- right + 1L
      excluded <- c(if (left > 0L) donors[left], if (right <= length(donors)) donors[right])
      gap <- min(abs(pq[excluded] - pq[i])) - max(abs(pq[selected] - pq[i]))
      if (is.null(minimum_gap) || gap < minimum_gap) minimum_gap <- gap
      if (gap <= 4 * epsilon) unresolved <- c(unresolved, i)
    }
    base$minimum_stored_gap <- if (is.null(minimum_gap)) NULL else A$record(A$pack(minimum_gap))
    base$required_gap_margin <- A$record(A$pack(4 * epsilon))
    base$unresolved_queries <- unresolved
    base$donor_sets_certified <- !length(unresolved)
    # Evaluate the raw matching contrast with intervals, preserving donor ratios.
    contribution <- list(lo = Q(rep(0, length(query))), hi = Q(rep(0, length(query))))
    Yq <- Q(Y)
    for (k in seq_along(query)) {
      i <- query[k]; j <- neighbors[[i]]
      imputed <- A$div(A$exact(sum(weight[j] * Yq[j])), A$exact(sum(weight[j])))
      term <- A$mul(A$exact(weight[i] * (2 * Z[i] - 1)), A$sub(A$exact(Yq[i]), imputed))
      contribution$lo[k] <- term$lo; contribution$hi[k] <- term$hi
    }
    exact_point <- A$div(A$total(contribution), A$exact(sum(weight[query])))
    point_error <- max(A$absmax(A$sub(A$exact(fit$estimate), exact_point)))
    base$exact_current_graph_point <- A$record(exact_point)
    base$point_arithmetic_error_upper <- A$record(A$exact(point_error))
    if (base$donor_sets_certified) {
      base$status <- "root_and_donor_sets_certified"
      base$point_root_error_bound <- base$point_arithmetic_error_upper
    }
    base
  }
  out <- tryCatch(calculate(), error = function(e) {
    base$reason <- conditionMessage(e)
    base
  })
  structure(out, class = c("wm_scalar_root_certificate", "list"))
}
