# Focused deterministic source validation, not a simulation/performance study.
# Args: package source root, preserved original integration-pilot results root.
# Uses two original saved datasets x PS/full-X at M3, with their existing counts.
# No fixture, data, fit, graph, bootstrap count or random number is generated.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
package_root <- normalizePath(args[1L], mustWork = TRUE)
pilot_root <- normalizePath(args[2L], mustWork = TRUE)
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
wm <- new.env(parent = .GlobalEnv)
for (path in sort(list.files(file.path(package_root, "R"), pattern = "\\.R$", full.names = TRUE))) {
  sys.source(path, envir = wm)
}
sys.source(file.path(package_root, "simulations", "wm_saved_full_x_stack.R"), envir = wm)
forbidden <- function(...) stop("Forbidden fit, rematch or RNG call.", call. = FALSE)
for (name in c("wm_match", "wm_fit", "wm_wdsm_fit", ".wm_ns_ols",
               ".wm_wdsm_nuisance_stack", ".wm_wdsm_prediction_stack",
               "wm_scalar_logistic_match", "wm_bootstrap", "wm_bootstrap_refit")) {
  assign(name, forbidden, envir = wm)
}
checks <- 0L
near <- function(a, b, label, tolerance = 2e-10) {
  if (length(a) != length(b) || any(!is.finite(c(a, b))) ||
      any(abs(a - b) > tolerance * pmax(1, abs(a), abs(b)))) {
    stop(label, call. = FALSE)
  }
  checks <<- checks + 1L
  invisible(max(abs(a - b) / pmax(1, abs(a), abs(b))))
}
quadratic <- function(x) {
  out <- cbind(intercept = 1, x)
  for (j in seq_len(ncol(x))) for (k in j:ncol(x)) {
    out <- cbind(out, x[, j] * x[, k])
    colnames(out)[ncol(out)] <- paste0(colnames(x)[j], ":", colnames(x)[k])
  }
  out
}
count_point <- function(fit, m, q0, q1) {
  z <- fit$data$Z; y <- fit$data$Y; w <- fit$analysis_weights
  K <- fit$loads$incoming[cbind(seq_along(z), z + 1L)]
  if (fit$estimand == "PATE") {
    numerator <- sum(m * w * (q1 - q0)) +
      sum(m * (2 * z - 1) * (w + K) * (y - ifelse(z == 1, q1, q0)))
    denominator <- sum(m * w)
  } else {
    numerator <- sum(m * z * w * (y - q0)) - sum(m * (1 - z) * K * (y - q0))
    denominator <- sum(m * z * w)
  }
  numerator / denominator
}
paths <- character(); input_md5 <- character(); summaries <- list()
case_labels <- c("GoodOverlap_retrospective_CorCor_PATE", "GoodOverlap_prospective_MisCor_PATT")
for (case_label in case_labels) {
  case_root <- file.path(pilot_root, case_label)
  case_paths <- file.path(case_root, c("data.rds", "weighted.rds", "counts.rds"))
  paths <- c(paths, case_paths); input_md5 <- c(input_md5, unname(tools::md5sum(case_paths)))
  data <- readRDS(case_paths[1L]); weighted <- readRDS(case_paths[2L])
  counts <- readRDS(case_paths[3L])$counts
  stopifnot(identical(ncol(counts), weighted$B), all(colSums(counts) == nrow(data)))
  original_data <- list(Y = data$Y, Z = data$A, weights = data$survey_weight)
  for (family in c("PS", "X6")) {
    key <- paste0("WM_", family, "_M3")
    stopifnot(isTRUE(weighted$points[[key]]$ok))
    fit <- weighted$points[[key]]$value
    fit_before <- serialize(fit, NULL, version = 2L)
    means <- weighted$bases[[paste0(family, "_means")]]$value
    required_means <- if (fit$estimand == "PATE") c("mean0", "mean1") else "mean0"
    if (family == "X6") {
      raw_x <- as.matrix(data[paste0("X", 1:6)])
      # The saved graph intentionally drops column names. Restore only the
      # recorded coordinate metadata on a local copy, never the source fit.
      source_standardized_x <- fit$graph$scores0
      stopifnot(ncol(source_standardized_x) == ncol(raw_x))
      colnames(source_standardized_x) <- colnames(raw_x)
    }
    models <- list()
    for (label in required_means) {
      arm <- if (label == "mean0") "arm0" else "arm1"
      D <- if (family == "PS") {
        weighted$bases$DSM$value$stack$inputs[[if (label == "mean0") "pg0_design" else "pg1_design"]]
      } else quadratic(source_standardized_x)
      if (family == "X6") {
        stopifnot(identical(colnames(D), names(means$coefficients[[arm]])))
      }
      models[[label]] <- list(design = D, coefficients = means$coefficients[[arm]])
    }
    if (family == "PS") {
      ps_input <- list(design = weighted$bases$DSM$value$stack$inputs$ps_design,
        probability = weighted$bases$PS$value$probability,
        weighting = if (weighted$sampling_design == "retrospective") "probability" else "unit")
      stack <- wm$wm_saved_full_x_stack(original_data, fit, models, family = "PS", ps = ps_input, wm = wm)
      stopifnot(all(!stack$parameter_values_available[stack$blocks$ps]))
    } else {
      stack <- wm$wm_saved_full_x_stack(original_data, fit, models, family = "full_X", raw_X = raw_x, wm = wm)
    }
    # An independent local equation evaluator. For missing logistic coefficients,
    # the recorded e fixes its local predictor; no new root or fit is selected.
    equation_at <- function(change, scores = FALSE) {
      block <- stack$blocks; n <- nrow(data)
      center <- stack$coordinate_center + change[block$center]
      variance <- stack$coordinate_variance + change[block$variance]
      if (family == "PS") {
        e <- stats::plogis(stats::qlogis(ps_input$probability) +
          as.vector(ps_input$design %*% change[block$ps]))
        coordinate <- matrix(e, n, 1L)
      } else coordinate <- raw_x
      centered <- sweep(coordinate, 2L, center, "-")
      if (scores) return(sweep(centered, 2L, sqrt(variance), "/"))
      out <- matrix(0, n, length(change))
      if (family == "PS") {
        psw <- if (ps_input$weighting == "probability") data$survey_weight / max(data$survey_weight) else rep(1, n)
        out[, block$ps] <- ps_input$design * (psw * (data$A - e))
      }
      out[, block$center] <- centered
      out[, block$variance] <- sweep(centered^2, 2L, variance, "-")
      for (label in required_means) {
        model <- stack$mean_models_raw[[label]]
        arm <- if (label == "mean0") 0L else 1L
        q <- as.vector(model$design %*% (model$coefficients + change[block[[label]]]))
        iw <- data$survey_weight / max(data$survey_weight) * (data$A == arm)
        out[, block[[label]]] <- model$design * (iw * (data$Y - q))
      }
      out
    }
    p <- length(stack$parameter_names); zero <- numeric(p)
    near(equation_at(zero), stack$estimating_equations, "Complete equations at saved parameter")
    jacobian_fd <- matrix(0, p, p); max_score_error <- 0
    for (j in seq_len(p)) {
      value <- if (is.finite(stack$parameter[j])) stack$parameter[j] else 0
      # Relative steps are needed for the inverse-square-root variance map;
      # an absolute step creates O((h/v)^2) central-quotient truncation.
      h <- 2e-5 * if (j %in% stack$blocks$variance) value else (1 + abs(value))
      change <- zero; change[j] <- h
      jacobian_fd[, j] <- (colMeans(equation_at(change)) - colMeans(equation_at(-change))) / (2 * h)
      score_fd <- (equation_at(change, TRUE) - equation_at(-change, TRUE)) / (2 * h)
      max_score_error <- max(max_score_error, near(score_fd, stack$score_derivative[, , j],
        paste("All coordinate derivative slots:", stack$parameter_names[j]), tolerance = 2e-7))
    }
    jacobian_error <- near(jacobian_fd, stack$jacobian, "Complete analytic Jacobian", tolerance = 2e-7)
    near(stack$estimating_equations + stack$nuisance_influence %*% t(stack$jacobian),
      matrix(0, nrow(data), p), "Influence orientation Psi+IF J'=0", tolerance = 2e-8)
    near(stack$nuisance_covariance,
      crossprod(sweep(stack$nuisance_influence, 2L, colMeans(stack$nuisance_influence), "-")) / nrow(data),
      "All same-row covariance slots")
    stopifnot(identical(colnames(stack$nuisance_influence), stack$parameter_names),
      length(stack$blocks$center) == ncol(fit$graph$scores0),
      length(stack$blocks$variance) == ncol(fit$graph$scores0))

    # Existing public common inference, independent finite-row/slope identity.
    inference <- do.call(wm$wm_fitted_inference, stack$inference_arguments)
    stopifnot(isTRUE(inference$available), identical(inference$estimate, fit$estimate),
      identical(serialize(fit, NULL, version = 2L), fit_before),
      identical(serialize(inference$fit, NULL, version = 2L), fit_before),
      identical(stack$assumptions_verified, FALSE), identical(inference$assumptions_verified, FALSE))
    z <- fit$data$Z; y <- fit$data$Y; w <- fit$analysis_weights
    K <- fit$loads$incoming[cbind(seq_along(z), z + 1L)]
    q0 <- fit$predictions$mean0; q1 <- fit$predictions$mean1
    slope <- colSums(stack$mean0_derivative * ((1 - z) * K - z * w))
    if (fit$estimand == "PATE") {
      slope <- slope + colSums(stack$mean1_derivative * ((1 - z) * w - z * K))
      rows <- (w * (q1 - q0 - fit$estimate) +
        (2 * z - 1) * (w + K) * (y - ifelse(z == 1, q1, q0))) / mean(w)
    } else rows <- (z * w * (y - q0 - fit$estimate) - (1 - z) * K * (y - q0)) / mean(z * w)
    slope <- slope / fit$denominator
    near(slope, inference$total_sensitivity, "Finite-M complete slope")
    centered_if <- sweep(stack$nuisance_influence, 2L, colMeans(stack$nuisance_influence), "-")
    direct_rows <- rows - mean(rows) + as.vector(centered_if %*% slope)
    direct_rows <- direct_rows - mean(direct_rows)
    near(direct_rows, inference$augmented_rows, "Existing common-API complete rows")
    near(mean(direct_rows^2) / nrow(data), inference$variance, "Direct V/n finite-row formula")

    # Existing B200 predictions/counts only: no WLS, new counts or rematching.
    count_comparisons <- 0L
    for (b in seq_len(ncol(counts))) {
      saved <- weighted$refit_diagnostics[[b]]$fits[[family]]
      if (!isTRUE(saved$ok)) stop("Preserved source count model is unavailable.", call. = FALSE)
      q_std <- q_raw <- list(mean0 = NULL, mean1 = NULL)
      for (label in required_means) {
        arm <- if (label == "mean0") "arm0" else "arm1"
        beta <- saved$coefficients[[arm]]
        q_std[[label]] <- as.vector(models[[label]]$design %*% beta)
        beta_raw <- if (family == "X6") as.vector(stack$raw_polynomial_transform %*% beta) else beta
        q_raw[[label]] <- as.vector(stack$mean_models_raw[[label]]$design %*% beta_raw)
        near(q_raw[[label]], q_std[[label]], "Saved count raw/source prediction equivalence")
      }
      stopifnot(identical(digest::digest(q_std, algo = "sha256"), saved$prediction_sha256))
      near(count_point(fit, counts[, b], q_std$mean0, q_std$mean1),
           weighted$summaries[[key]]$draws[b], "Preserved original source count point")
      near(count_point(fit, counts[, b], q_raw$mean0, q_raw$mean1),
           weighted$summaries[[key]]$draws[b], "Raw-chart count prediction point equivalence")
      count_comparisons <- count_comparisons + 1L
    }
    summaries[[paste(case_label, family, sep = "/")]] <- data.frame(
      case = case_label, family = family, n = nrow(data), p = p,
      jacobian_scaled_error = jacobian_error, score_scaled_error = max_score_error,
      existing_count_columns = count_comparisons, checks = checks,
      scientific_premises = "UNKNOWN", stringsAsFactors = FALSE)
  }
}
stopifnot(identical(unname(tools::md5sum(paths)), input_md5),
  !exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
print(do.call(rbind, summaries), row.names = FALSE)
cat("PASS: four saved PS/full-X PATE/PATT cases; source inputs unchanged; no fit/rematch/RNG.\n")
cat("Software arithmetic only. Source sampling/CLT/callback/coverage premises remain UNKNOWN.\n")
