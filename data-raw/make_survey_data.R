# Regenerate the synthetic example from the package root.
# Requires the optional development package usethis; it is not needed to use wdsmatch.
set.seed(2025)
n <- 800
X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
X4 <- rnorm(n); X5 <- rnorm(n); X6 <- rnorm(n)
ps_true <- plogis(0.3 + 0.6*X1 + 0.4*X2 - 0.3*X3 + 0.2*X1*X2)
Z <- rbinom(n, 1, ps_true)
mu0 <- 1 + X1 + 0.5*X2 - 0.3*X3 + 0.2*X4 + 0.3*X1*X2
mu1 <- mu0 + 0.8 + 0.2*X1
Y0 <- mu0 + rnorm(n, sd = 0.5)
Y1 <- mu1 + rnorm(n, sd = 0.5)
Y <- Z * Y1 + (1 - Z) * Y0
sw_prob <- plogis(-2 + 0.3*Z + 0.2*X1 + 0.15*X2)
S <- rbinom(n, 1, sw_prob)
survey_obs <- data.frame(
  Y = Y[S == 1], Z = Z[S == 1],
  X1 = X1[S == 1], X2 = X2[S == 1], X3 = X3[S == 1],
  X4 = X4[S == 1], X5 = X5[S == 1], X6 = X6[S == 1],
  survey_weight = 1 / sw_prob[S == 1]
)
rownames(survey_obs) <- NULL
usethis::use_data(survey_obs, overwrite = TRUE)
