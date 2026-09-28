# Nuisance fits for known, fixed matching maps. These routines deliberately
# retain every requested basis term and fail when that basis is not estimable.

.wm_fit_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x != floor(x) || x > .Machine$integer.max) {
    stop(name, " must be a single integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

.wm_fit_arm_argument <- function(x, name, arm, estimand) {
  allowed <- if (estimand == "PATT") 1L else c(1L, 2L)
  if (!is.numeric(x) || is.complex(x) ||
      !(length(x) %in% allowed) || anyNA(x)) {
    stop(name, " must be numeric with length one",
         if (estimand == "PATE") " or two (control, treated)" else "",
         ".", call. = FALSE)
  }
  x[if (length(x) == 1L) 1L else arm + 1L]
}

.wm_fit_scores <- function(x, name, n) {
  if (is.numeric(x) && is.null(dim(x))) x <- matrix(x, ncol = 1L)
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || nrow(x) != n ||
      ncol(x) < 1L || any(!is.finite(x))) {
    stop(name, " must be a finite numeric n-by-d matrix (or a scalar-score ",
         "vector) with d >= 1.", call. = FALSE)
  }
  x
}

.wm_fit_labels <- function(x, name, n) {
  if (is.factor(x)) x <- as.character(x)
  if ((!is.numeric(x) && !is.character(x) && !is.logical(x)) ||
      !is.null(dim(x)) || length(x) != n || anyNA(x) ||
      (is.numeric(x) && any(!is.finite(x)))) {
    stop(name, " must contain one nonmissing label per row.", call. = FALSE)
  }
  # Preserve adjacent numeric values: formatting them as character can merge
  # distinct cells. Training masks below use equality on these original values.
  x
}

.wm_fit_support <- function(support, scores, name) {
  d <- ncol(scores)
  if (d == 1L && is.numeric(support) && is.null(dim(support)) &&
      length(support) == 2L) support <- matrix(support, nrow = 1L)
  if (!is.matrix(support) || !is.numeric(support) || is.complex(support) ||
      !identical(dim(support), c(d, 2L)) ||
      any(!is.finite(support)) || any(support[, 1L] >= support[, 2L]) ||
      any(!is.finite(support[, 2L] - support[, 1L]))) {
    stop(name, " must be a supplied d-by-2 finite matrix of known lower and ",
         "upper support endpoints, with lower < upper.", call. = FALSE)
  }
  for (j in seq_len(d)) {
    if (any(scores[, j] < support[j, 1L] |
            scores[, j] > support[j, 2L])) {
      stop(name, ": every query and donor score must lie inside the declared ",
           "support; column ", j, " contains an out-of-support value.",
           call. = FALSE)
    }
  }
  colnames(support) <- c("lower", "upper")
  support
}

.wm_fit_basis_size <- function(log_size, n_training, n_evaluation, d,
                               max_basis, max_basis_elements, context) {
  if (!is.finite(log_size) || log_size > log(max_basis) + 1e-12) {
    stop(context, ": requested basis exceeds max_basis = ", max_basis,
         "; choose an explicitly smaller basis or increase the guard.",
         call. = FALSE)
  }
  k <- round(exp(log_size))
  if (k > n_training) {
    stop(context, ": insufficient training donors: ", n_training,
         " rows for ", k, " requested basis columns.", call. = FALSE)
  }
  # Bounds the two design matrices and the stored tensor/exponent index.
  if (k * (n_training + n_evaluation + d) > max_basis_elements) {
    stop(context, ": requested basis allocation exceeds max_basis_elements = ",
         max_basis_elements, ".", call. = FALSE)
  }
  as.integer(k)
}

.wm_fit_polynomial_exponents <- function(d, degree, k) {
  ans <- matrix(0L, nrow = k, ncol = d)
  row <- 0L
  current <- integer(d)
  visit <- function(position, remaining) {
    if (position == d) {
      current[position] <<- remaining
      row <<- row + 1L
      ans[row, ] <<- current
    } else {
      for (value in 0:remaining) {
        current[position] <<- value
        visit(position + 1L, remaining - value)
      }
    }
    invisible(NULL)
  }
  for (total in 0:degree) visit(1L, total)
  colnames(ans) <- paste0("score", seq_len(d))
  ans
}

.wm_fit_basis_spec <- function(scores, regression, degree, spline_order,
                               mesh_exponent, moment_order, support,
                               donor_rule, n_training, n_evaluation,
                               max_basis, max_basis_elements, context) {
  d <- ncol(scores)
  if (regression == "polynomial") {
    degree <- .wm_fit_integer(degree, "degree", 0L)
    k <- .wm_fit_basis_size(lchoose(as.double(d) + degree, degree), n_training,
                            n_evaluation, d, max_basis, max_basis_elements,
                            context)
    # Recursion depth is dimension. High-dimensional constant fits need none.
    if (degree == 0L) {
      exponents <- matrix(0L, nrow = 1L, ncol = d)
      colnames(exponents) <- paste0("score", seq_len(d))
    } else {
      if (d > 100L) {
        stop(context, ": polynomial dimension above 100 is unsupported by ",
             "the exact reference basis constructor.", call. = FALSE)
      }
      exponents <- .wm_fit_polynomial_exponents(d, degree, k)
    }
    return(list(type = regression, dimension = d, degree = degree,
                exponents = exponents, n_columns = k, support = NULL))
  }
  r <- .wm_fit_integer(spline_order, "spline_order", 2L)
  if (!is.numeric(mesh_exponent) || length(mesh_exponent) != 1L ||
      !is.finite(mesh_exponent) || mesh_exponent <= 0) {
    stop("mesh_exponent must be an explicitly supplied positive finite number.",
         call. = FALSE)
  }
  if (!is.numeric(moment_order) || length(moment_order) != 1L ||
      is.na(moment_order) || moment_order <= 2) {
    stop("moment_order must explicitly declare p > 2 or Inf.", call. = FALSE)
  }
  t_d <- max(0, 0.5 - 1 / d)
  lower <- t_d / (r - 1)
  moment_upper <- if (is.infinite(moment_order)) 1 / d else
    (1 - 2 / moment_order) / d
  upper <- min((1 - 2 * t_d) / (d + 2), moment_upper)
  if (donor_rule == "stabilized") {
    lower <- max(lower, 1 / (4 * r))
    upper <- min(upper, 1 / (2 * d))
  }
  if (!(mesh_exponent > lower && mesh_exponent < upper)) {
    stop(context, ": mesh_exponent must satisfy the strict sufficient rate ",
         "interval (", signif(lower, 6), ", ", signif(upper, 6),
         ") for this dimension, order, moment declaration, and donor rule.",
         call. = FALSE)
  }
  kappa <- ceiling(n_training^mesh_exponent)
  k <- .wm_fit_basis_size(d * log(kappa + r - 1), n_training,
                          n_evaluation, d, max_basis, max_basis_elements,
                          context)
  if (!is.finite(kappa) || kappa > .Machine$integer.max) {
    stop(context, ": unsupported spline mesh size.", call. = FALSE)
  }
  kappa <- as.integer(kappa)
  interiors <- if (kappa > 1L) seq_len(kappa - 1L) / kappa else numeric()
  knots <- c(rep(0, r), interiors, rep(1, r))
  q <- kappa + r - 1L
  # First coordinate varies fastest, matching the tensor multiplication below.
  tensor_index <- matrix(0L, nrow = k, ncol = d)
  for (j in seq_len(d)) {
    tensor_index[, j] <- rep(rep(seq_len(q), each = q^(j - 1L)),
                             length.out = k)
  }
  list(type = regression, dimension = d, order = r, kappa = kappa,
       mesh_exponent = mesh_exponent, moment_order = moment_order,
       rate_interval = c(lower = lower, upper = upper),
       support = support, knots = rep(list(knots), d),
       tensor_index = tensor_index, n_columns = k,
       coordinate_transform = "fixed known support to [0,1]")
}

.wm_fit_basis <- function(scores, spec) {
  n <- nrow(scores)
  if (spec$type == "polynomial") {
    ans <- matrix(1, nrow = n, ncol = spec$n_columns)
    for (k in seq_len(spec$n_columns)) {
      for (j in which(spec$exponents[k, ] != 0L)) {
        ans[, k] <- ans[, k] * scores[, j]^spec$exponents[k, j]
      }
    }
  } else {
    lower <- spec$support[, 1L]
    width <- spec$support[, 2L] - lower
    scaled <- sweep(sweep(scores, 2L, lower, "-"), 2L, width, "/")
    ans <- matrix(1, nrow = n, ncol = 1L)
    for (j in seq_len(spec$dimension)) {
      one <- splines::splineDesign(spec$knots[[j]], scaled[, j],
                                   ord = spec$order, outer.ok = FALSE)
      old_columns <- ncol(ans)
      next_basis <- matrix(0, nrow = n,
                            ncol = old_columns * ncol(one))
      for (column in seq_len(ncol(one))) {
        block <- (column - 1L) * old_columns + seq_len(old_columns)
        next_basis[, block] <- ans * one[, column]
      }
      ans <- next_basis
    }
  }
  if (any(!is.finite(ans))) {
    stop("The requested basis contains nonfinite values; respecify the fixed ",
         "coordinates or basis without changing them implicitly.",
         call. = FALSE)
  }
  colnames(ans) <- paste0("b", seq_len(ncol(ans)))
  ans
}

.wm_fit_qr <- function(design, response, weights, qr_tol, context) {
  zero <- which(colSums(abs(design)) == 0)
  if (length(zero)) {
    stop(context, ": zero training basis column(s): ",
         paste(zero, collapse = ", "),
         ". No terms were dropped.", call. = FALSE)
  }
  weight_scale <- max(weights)
  root_weight <- sqrt(weights / weight_scale)
  weighted_design <- design * root_weight
  weighted_response <- response * root_weight
  decomposition <- qr(weighted_design, tol = qr_tol, LAPACK = FALSE)
  if (decomposition$rank < ncol(design)) {
    stop(context, ": rank-deficient training basis (rank ",
         decomposition$rank, " of ", ncol(design),
         "). No terms were dropped and no ridge penalty was added.",
         call. = FALSE)
  }
  coefficients <- as.vector(qr.coef(decomposition, weighted_response))
  if (any(!is.finite(coefficients))) {
    stop(context, ": nonfinite least-squares coefficients.", call. = FALSE)
  }
  names(coefficients) <- colnames(design)
  list(coefficients = coefficients, rank = decomposition$rank,
       qr_pivot = decomposition$pivot, qr_tol = qr_tol,
       weight_scale = weight_scale)
}

#' Fit nuisance regressions and match using fixed known score maps
#'
#' @details
#' Common-score stabilized PATE requires bounded outcomes. Stabilized PATT
#' permits finite p greater than two under its separate residual and
#' treated-outcome moment, sampling-equivalence and uniform nuisance-consistency
#' conditions. The concrete spline procedure additionally requires declared
#' smoothness, the strict mesh window, positive normalizer bounds and bounded
#' conditional residual variance. Its finite-moment squared-row variance and
#' Gaussian multiplier result do not extend common-score PATE to unbounded
#' outcomes. Score maps, weights and target identification obey the fixed-map
#' contract in \code{\link{wm_match}}.
#'
#' @export
wm_fit <- function(Y, Z, weights, scores0, scores1 = NULL, M = 3L,
                   estimand = c("PATE", "PATT"),
                   method = c("self_normalized", "stabilized"),
                   regression = c("polynomial", "spline"), degree = 1L,
                   spline_order = NULL, mesh_exponent = NULL,
                   moment_order = NULL, support0 = NULL, support1 = NULL,
                   rho_bounds = NULL, fold_id = NULL, strata = NULL,
                   variance = TRUE, max_basis = 10000L,
                   max_basis_elements = 1e7, qr_tol = 1e-10) {
  fit_call <- match.call()
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  regression <- match.arg(regression)
  if (!is.numeric(Y) || is.complex(Y) || !is.null(dim(Y)) ||
      length(Y) < 2L) {
    stop("Y must be a finite numeric vector with at least two rows.",
         call. = FALSE)
  }
  n <- length(Y)
  Y <- .wm_numeric_vector(Y, n, "Y")
  weights <- .wm_numeric_vector(weights, n, "weights", positive = TRUE)
  if (any(weights / max(weights) <= 0)) {
    stop("Common weight normalization underflowed; the weight range is too wide.",
         call. = FALSE)
  }
  if ((!is.numeric(Z) && !is.logical(Z)) || is.complex(Z) ||
      !is.null(dim(Z)) || length(Z) != n || anyNA(Z) ||
      any(!Z %in% c(0, 1)) || length(unique(Z)) != 2L) {
    stop("Z must have length n, contain only 0 and 1, and include both arms.",
         call. = FALSE)
  }
  Z <- as.integer(Z)
  M <- .wm_fit_integer(M, "M")
  max_basis <- .wm_fit_integer(max_basis, "max_basis")
  if (!is.numeric(max_basis_elements) || is.complex(max_basis_elements) ||
      length(max_basis_elements) != 1L || !is.finite(max_basis_elements) ||
      max_basis_elements < 1) {
    stop("max_basis_elements must be a positive finite allocation limit.",
         call. = FALSE)
  }
  if (!is.numeric(qr_tol) || is.complex(qr_tol) || length(qr_tol) != 1L ||
      !is.finite(qr_tol) || qr_tol <= 0 || qr_tol >= 1) {
    stop("qr_tol must be finite and strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.logical(variance) || length(variance) != 1L || is.na(variance)) {
    stop("variance must be TRUE or FALSE.", call. = FALSE)
  }
  scores0 <- .wm_fit_scores(scores0, "scores0", n)
  if (estimand == "PATT" && (!is.null(scores1) || !is.null(support1))) {
    stop("PATT fits only the control regression; omit scores1 and support1.",
         call. = FALSE)
  }
  score1_supplied <- !is.null(scores1)
  if (estimand == "PATE") {
    scores1 <- if (is.null(scores1)) scores0 else
      .wm_fit_scores(scores1, "scores1", n)
  }
  same_maps <- estimand == "PATE" &&
    identical(dim(scores0), dim(scores1)) &&
    isTRUE(all(scores0 == scores1))
  if (method == "stabilized" && estimand == "PATE" && !same_maps) {
    stop("Stabilized PATE requires the same supplied score matrix in both arms.",
         call. = FALSE)
  }
  folds <- if (is.null(fold_id)) rep("all", n) else
    .wm_fit_labels(fold_id, "fold_id", n)
  restrictions <- if (is.null(strata)) rep("all", n) else
    .wm_fit_labels(strata, "strata", n)
  if (!is.null(fold_id) && length(unique(folds)) < 2L) {
    stop("fold_id must define at least two nonempty fixed, data-independent folds.",
         call. = FALSE)
  }
  if (regression == "spline" && method == "stabilized" && is.null(fold_id)) {
    stop("Stabilized spline fitting requires a supplied fixed fold_id with ",
         "at least two folds; matching is wholly within evaluation folds.",
         call. = FALSE)
  }
  arms <- if (estimand == "PATE") 0:1 else 0L
  supports <- list(NULL, NULL)
  degree_arm <- order_arm <- exponent_arm <- moment_arm <- vector("list", 2L)
  if (regression == "spline") {
    if (is.null(spline_order) || is.null(mesh_exponent) ||
        is.null(moment_order)) {
      stop("Spline fitting requires explicit spline_order, mesh_exponent, ",
           "and moment_order.", call. = FALSE)
    }
    supports[[1L]] <- .wm_fit_support(support0, scores0, "support0")
    if (estimand == "PATE") {
      if (is.null(support1) && same_maps) support1 <- support0
      supports[[2L]] <- .wm_fit_support(support1, scores1, "support1")
    }
    for (arm in arms) {
      order_arm[[arm + 1L]] <- .wm_fit_arm_argument(
        spline_order, "spline_order", arm, estimand)
      exponent_arm[[arm + 1L]] <- .wm_fit_arm_argument(
        mesh_exponent, "mesh_exponent", arm, estimand)
      moment_arm[[arm + 1L]] <- .wm_fit_arm_argument(
        moment_order, "moment_order", arm, estimand)
    }
  } else {
    if (!is.null(spline_order) || !is.null(mesh_exponent) ||
        !is.null(moment_order) || !is.null(support0) || !is.null(support1)) {
      stop("spline_order, mesh_exponent, moment_order and support arguments ",
           "apply only to regression = 'spline'.", call. = FALSE)
    }
    for (arm in arms) {
      degree_arm[[arm + 1L]] <- .wm_fit_arm_argument(
        degree, "degree", arm, estimand)
    }
  }
  bounds <- NULL
  if (method == "stabilized") {
    if (is.numeric(rho_bounds) && !is.complex(rho_bounds) &&
        is.null(dim(rho_bounds)) && length(rho_bounds) == 2L) {
      bounds <- matrix(rep(rho_bounds, length(arms)), ncol = 2L, byrow = TRUE)
    } else if (is.matrix(rho_bounds) && is.numeric(rho_bounds) &&
               !is.complex(rho_bounds) &&
               identical(dim(rho_bounds), c(length(arms), 2L))) {
      bounds <- rho_bounds
    }
    if (is.null(bounds) || any(!is.finite(bounds)) ||
        any(bounds[, 1L] <= 0) || any(bounds[, 1L] >= bounds[, 2L])) {
      stop("rho_bounds must supply positive finite lower < upper bounds ",
           "(a pair, or one matrix row per fitted arm). The caller declares ",
           "these strictly contain the true conditional weight means.",
           call. = FALSE)
    }
    dimnames(bounds) <- list(paste0("arm", arms), c("lower", "upper"))
  } else if (!is.null(rho_bounds)) {
    stop("rho_bounds applies only to method = 'stabilized'.", call. = FALSE)
  }

  means <- list(rep(NA_real_, n), rep(NA_real_, n))
  rhos <- list(rep(NA_real_, n), rep(NA_real_, n))
  models <- list()
  model_id <- 0L
  stratum_levels <- unique(restrictions)
  for (stratum_id in seq_along(stratum_levels)) {
    stratum <- stratum_levels[stratum_id]
    stratum_rows <- which(restrictions == stratum)
    fold_levels <- unique(folds[stratum_rows])
    for (fold_index in seq_along(fold_levels)) {
      fold <- fold_levels[fold_index]
      evaluation <- stratum_rows[folds[stratum_rows] == fold]
      training_pool <- if (is.null(fold_id)) stratum_rows else
        stratum_rows[folds[stratum_rows] != fold]
      for (arm in arms) {
        context <- paste0("stratum '", stratum, "', fold '", fold,
                          "', arm ", arm)
        if (any(Z[evaluation] == 1L - arm) &&
            sum(Z[evaluation] == arm) < M) {
          stop(context, ": insufficient evaluation donors: ",
               sum(Z[evaluation] == arm), " but M = ", M, ".",
               call. = FALSE)
        }
        training <- training_pool[Z[training_pool] == arm]
        if (!length(training)) {
          stop(context, ": no training donors in the required arm and stratum.",
               call. = FALSE)
        }
        scores <- if (arm == 0L) scores0 else scores1
        spec <- .wm_fit_basis_spec(
          scores, regression, degree_arm[[arm + 1L]],
          order_arm[[arm + 1L]], exponent_arm[[arm + 1L]],
          moment_arm[[arm + 1L]], supports[[arm + 1L]], method,
          length(training), length(evaluation), max_basis,
          max_basis_elements, context)
        training_design <- .wm_fit_basis(scores[training, , drop = FALSE], spec)
        evaluation_design <- .wm_fit_basis(scores[evaluation, , drop = FALSE], spec)
        mean_fit <- .wm_fit_qr(training_design, Y[training], weights[training],
                                qr_tol, paste0(context, " mean fit"))
        predicted_mean <- as.vector(evaluation_design %*% mean_fit$coefficients)
        if (any(!is.finite(predicted_mean))) {
          stop(context, ": nonfinite outcome predictions.", call. = FALSE)
        }
        means[[arm + 1L]][evaluation] <- predicted_mean
        mean_fit$weighted <- TRUE
        rho_fit <- NULL
        if (method == "stabilized") {
          rho_fit <- .wm_fit_qr(training_design, weights[training],
                                 rep(1, length(training)), qr_tol,
                                 paste0(context, " rho fit"))
          raw_rho <- as.vector(evaluation_design %*% rho_fit$coefficients)
          if (any(!is.finite(raw_rho))) {
            stop(context, ": nonfinite conditional-weight predictions.",
                 call. = FALSE)
          }
          arm_bounds <- bounds[arm + 1L, ]
          clipped <- raw_rho < arm_bounds[1L] | raw_rho > arm_bounds[2L]
          rhos[[arm + 1L]][evaluation] <- pmin(arm_bounds[2L],
                                             pmax(arm_bounds[1L], raw_rho))
          rho_fit$weighted <- FALSE
          rho_fit$bounds <- arm_bounds
          rho_fit$clipped_count <- sum(clipped)
          rho_fit$clipped_rows <- evaluation[clipped]
          rho_fit$prediction_range_before_clipping <- range(raw_rho)
        }
        model_id <- model_id + 1L
        models[[model_id]] <- list(
          arm = arm, stratum = stratum, fold = fold,
          training_rows = training, evaluation_rows = evaluation,
          n_training = length(training), basis = spec,
          mean = mean_fit, rho = rho_fit)
        names(models)[model_id] <- paste0("arm", arm, "_stratum",
                                           stratum_id, "_fold", fold_index)
      }
    }
  }
  result <- wm_match(
    Y = Y, Z = Z, weights = weights, scores0 = scores0,
    scores1 = if (estimand == "PATE" && score1_supplied) scores1 else NULL,
    M = M, estimand = estimand, method = method,
    mean0 = means[[1L]], mean1 = if (estimand == "PATE") means[[2L]] else NULL,
    rho0 = if (method == "stabilized") rhos[[1L]] else NULL,
    rho1 = if (method == "stabilized" && estimand == "PATE") rhos[[2L]] else NULL,
    fold_id = fold_id, strata = strata, variance = variance)
  result$fit_call <- fit_call
  result$nuisance <- list(
    regression = regression, models = models,
    predictions = list(mean0 = means[[1L]],
                       mean1 = if (estimand == "PATE") means[[2L]] else NULL,
                       rho0 = if (method == "stabilized") rhos[[1L]] else NULL,
                       rho1 = if (method == "stabilized" && estimand == "PATE")
                         rhos[[2L]] else NULL),
    fold_id = if (is.null(fold_id)) NULL else folds,
    strata = if (is.null(strata)) NULL else restrictions,
    known_support = if (regression == "spline") supports[arms + 1L] else NULL,
    declared_moment_order = if (regression == "spline") moment_order else NULL,
    rho_bounds = bounds,
    rho_clipped_count = sum(vapply(models, function(x) {
      if (is.null(x$rho)) 0L else x$rho$clipped_count
    }, numeric(1L))),
    max_basis = max_basis, max_basis_elements = max_basis_elements,
    qr_tol = qr_tol,
    inference_contract = list(
      score_maps = "fixed known maps; no estimated-map uncertainty is included",
      weights = "known positive individual weights",
      moment_declaration = if (regression == "spline")
        "Caller-supplied residual moment order; Inf means enough finite moments, not a data check"
        else "Finite moments and regular correct polynomial models are caller assumptions",
      conditional_variance = "Bounded conditional residual variance is required for spline rates",
      outcome_model = if (method == "self_normalized")
        "Strong graph-field residual centering and the stated outcome-regression assumptions"
        else "Correct weighted outcome and ordinary conditional-weight mean functions",
      normalizer_bounds = if (method == "stabilized")
        "Caller declares rho_bounds strictly contain the true conditional-weight means" else NULL,
      bounded_outcomes = if (method == "stabilized")
        "Common-score stabilized PATE requires bounded outcomes; PATT permits finite p>2 moments under its separate nuisance and variance-transfer conditions, with bounded conditional residual variance for spline rates"
        else NULL,
      folds = if (!is.null(fold_id))
        "Caller declares a fixed data-independent partition; fits exclude each evaluation fold and matches stay inside it"
        else "Same-sample nuisance fits",
      strata = if (!is.null(strata))
        paste("Fits and matching stay in each fixed stratum with positive limiting proportions;",
              "spline arms use one common declared rectangle across strata,",
              "with positive donor density throughout that rectangle in each stratum")
        else NULL))
  class(result) <- c("wm_fit", class(result))
  result
}
