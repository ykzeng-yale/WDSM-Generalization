validation_inputs <- function() {
  list(Y = c(1, 3, 2, 6, 4, 8, 5, 9), X = data.frame(x = 1:8),
       Z = rep(c(0, 1), 4), weights = rep(1, 8), M = 1,
       ps = seq(0.15, 0.85, length.out = 8),
       pg = cbind(psi0 = sin(1:8), psi1 = cos(1:8)),
       use.bias.correction = FALSE, varest = FALSE)
}

validation_fit <- function(..., estimand = "PATE") {
  args <- validation_inputs()
  changes <- list(...)
  args[names(changes)] <- changes
  do.call(if (estimand == "PATE") wdsmatchATE else wdsmatchATT, args)
}

validation_run <- function(counts, estimand = "PATE", sampling = "retrospective",
                           alpha = 0.05) {
  args <- validation_inputs()
  args$varest <- TRUE
  args$boots <- 2L
  args$alpha <- alpha
  args$sampling <- sampling
  args$estimand <- estimand
  args$call <- quote(quote(wdsmatchATE()))
  args[c("model.ps", "model.pg")] <- list(NULL, NULL)
  args$bootstrap_counts <- counts
  do.call(wdsmatch:::wdsm_run, args)
}

test_that("valid fixed-score inputs reach both public estimators", {
  expect_s3_class(validation_fit(), "wdsmatch")
  expect_s3_class(validation_fit(estimand = "PATT"), "wdsmatch")
  one_score <- matrix(sin(1:8), ncol = 1)
  expect_s3_class(validation_fit(pg = one_score, estimand = "PATT"), "wdsmatch")
  expect_equal(validation_fit(pg = one_score, estimand = "PATT")$estimate,
               validation_fit(estimand = "PATT")$estimate)
})

test_that("outcomes, treatment, and weights require aligned finite vectors", {
  args <- validation_inputs()
  args$weights <- NULL
  expect_error(do.call(wdsmatchATE, args), "weights.*required")
  expect_error(validation_fit(weights = NULL), "weights.*required")
  for (w in list(rep(1, 7), c(0, rep(1, 7)), c(-1, rep(1, 7)),
                 c(Inf, rep(1, 7)), c(NA_real_, rep(1, 7)),
                 matrix(1, 8, 1), rep("1", 8))) {
    expect_error(validation_fit(weights = w), "weights.*finite, positive")
  }
  for (y in list(c(NA_real_, 2:8), c(Inf, 2:8), matrix(1:8, 8, 1), letters[1:8])) {
    expect_error(validation_fit(Y = y), "Y.*finite numeric vector")
  }
  for (z in list(rep(0, 7), c(2, rep(0, 7)), c(NA_real_, rep(0, 7)),
                 rep("0", 8), matrix(rep(c(0, 1), 4), 8, 1))) {
    expect_error(validation_fit(Z = z), "Z.*binary")
  }
})

test_that("covariate rows, data, and names are checked before fitting", {
  expect_error(validation_fit(X = list(x = 1:8)), "X.*matrix or data frame")
  for (x in list(data.frame(x = 1:7), data.frame(x = c(NA_real_, 2:8)),
                 data.frame(x = c(Inf, 2:8)), data.frame(x = factor(1:8)))) {
    expect_error(validation_fit(X = x), "X.*finite numeric row")
  }
  for (bad_name in c("Y", "Z", "bad`name")) {
    x <- data.frame(x = 1:8)
    names(x) <- bad_name
    expect_error(validation_fit(X = x), "X column names")
  }
  duplicated <- data.frame(x = 1:8, y = 8:1)
  names(duplicated) <- c("x", "x")
  expect_error(validation_fit(X = duplicated), "X column names")
})

test_that("M and donor availability are validated for the target estimand", {
  for (m in list(0, -1, 1.5, Inf, NA_real_, c(1, 3), "3")) {
    expect_error(validation_fit(M = m), "M.*positive integer")
  }
  for (z in list(rep(0, 8), rep(1, 8))) {
    expect_error(validation_fit(Z = z), "Insufficient treatment-arm")
  }
  expect_error(validation_fit(M = 5), "opposite-arm donors")
  expect_error(validation_fit(M = 5, estimand = "PATT"), "opposite-arm donors")
  one_treated <- c(1, rep(0, 7))
  expect_error(validation_fit(Z = one_treated, M = 3), "opposite-arm donors")
  expect_s3_class(validation_fit(Z = one_treated, M = 3, estimand = "PATT"),
                  "wdsmatch")
})

test_that("inference controls reject invalid alpha and replication counts", {
  for (a in list(0, 1, -0.1, Inf, NA_real_, c(0.05, 0.1), "0.05")) {
    expect_error(validation_fit(alpha = a), "alpha.*strictly between")
  }
  for (b in list(0, 1, 2.5, Inf, NA_real_, c(2, 3), "2")) {
    expect_error(validation_fit(varest = TRUE, boots = b), "boots.*at least two")
  }
  expect_error(validation_fit(varest = NA), "varest.*TRUE or FALSE")
  expect_error(validation_fit(use.bias.correction = 1), "varest.*TRUE or FALSE")
  expect_error(validation_fit(sampling = "unknown"), "arg.*should be one of")
})

test_that("supplied score dimensions and fitted formula responses are checked", {
  for (ps in list(rep(0.5, 7), c(0, rep(0.5, 7)), c(1, rep(0.5, 7)),
                  c(NA_real_, rep(0.5, 7)), matrix(0.5, 8, 1))) {
    expect_error(validation_fit(ps = ps), "ps.*finite probability")
  }
  for (pg in list(sin(1:8), matrix(1, 7, 2), matrix(1, 8, 3),
                  matrix(NA_real_, 8, 2), matrix("score", 8, 2))) {
    expect_error(validation_fit(pg = pg), "pg.*finite numeric matrix")
  }
  expect_error(validation_fit(pg = matrix(sin(1:8), 8, 1)), "two columns for PATE")
  expect_error(validation_fit(ps = NULL, model.ps = Y ~ x), "response Z")
  expect_error(validation_fit(ps = NULL, model.ps = ~ x), "response Z")
  expect_error(validation_fit(ps = NULL, model.ps = "Z ~ x"), "response Z")
  expect_error(validation_fit(pg = NULL, model.pg = Z ~ x), "response Y")
  expect_error(validation_fit(ps = NULL, model.ps = Z ~ missing_predictor),
               "predictors must be columns of X")
})

test_that("malformed multinomial counts cannot produce partial inference", {
  counts <- matrix(1, 8, 2)
  expect_error(validation_run(counts[, 1, drop = FALSE]), "n-by-boots multinomial counts")
  wrong_total <- counts
  wrong_total[1, 1] <- 0
  expect_error(validation_run(wrong_total), "n-by-boots multinomial counts")
  fractional <- counts
  fractional[1:2, 1] <- c(0.5, 1.5)
  expect_error(validation_run(fractional), "n-by-boots multinomial counts")
  negative <- counts
  negative[1:2, 1] <- c(-1, 3)
  expect_error(validation_run(negative), "n-by-boots multinomial counts")
  nonfinite <- counts
  nonfinite[1, 1] <- NA_real_
  expect_error(validation_run(nonfinite), "n-by-boots multinomial counts")
})

test_that("an empty-arm draw stops with replicate context in all four paths", {
  counts <- matrix(1, 8, 2)
  # Replicate 1 is valid. Replicate 2 contains only controls; inference
  # must stop, not omit the failed draw or substitute the point estimate.
  counts[, 2] <- rep(c(2, 0), 4)
  for (estimand in c("PATE", "PATT")) {
    for (sampling in c("retrospective", "prospective")) {
      expect_error(validation_run(counts, estimand, sampling),
        paste0(sampling, " ", estimand,
               " bootstrap replicate 2/2 failed: both treatment arms must have positive bootstrap multiplicity"),
        fixed = TRUE)
    }
  }
})

test_that("printed interval levels follow alpha and retain legacy compatibility", {
  counts <- cbind(rep(1, 8), c(2, 0, 2, 0, 0, 2, 0, 2))
  fit <- validation_run(counts, alpha = 0.1)
  expect_output(print(fit), "90% Wald CI", fixed = TRUE)
  expect_output(summary(fit), "90% Wald CI", fixed = TRUE)
  invisible(capture.output(visible <- withVisible(print(fit))))
  expect_identical(visible$visible, FALSE)
  legacy <- fit
  legacy$alpha <- legacy$interval.type <- NULL
  expect_output(print(legacy), "95% CI", fixed = TRUE)
})
