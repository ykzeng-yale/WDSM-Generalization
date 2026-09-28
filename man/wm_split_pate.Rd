\name{wm_split_pate}
\alias{wm_split_pate}
\title{Stabilized PATE from Independent Evaluation Blocks}
\usage{
wm_split_pate(Y, Z, weights, scores0, scores1, block_id,
  mean0, mean1, rho0, rho1, M = 3L, strata = NULL, variance = TRUE)
}
\arguments{
\item{Y, Z, weights}{Finite outcome, binary treatment and positive individual-weight vectors.}
\item{scores0, scores1}{Fixed matching matrices for the two potential means. Their positive dimensions may differ.}
\item{block_id}{Fixed data-independent evaluation-block labels, zero or one. Each block must contain both treatment groups and sufficient donors.}
\item{mean0, mean1, rho0, rho1}{Supplied predictions at every original row. Each block uses only the predictions for its assigned potential mean.}
\item{M}{Positive fixed donor count.}
\item{strata}{Optional finite matching-stratum labels.}
\item{variance}{Whether to calculate contribution-based variance.}
}
\description{Estimates the treated potential mean in block one and the control potential mean in block zero, then subtracts them. The independent-block construction permits different continuous matching maps for the stabilized potential means.}
\details{
The two blocks contain independent iid observations with positive limiting size fractions. Each block has its own weight denominator. Supplied nuisance predictions must be known, external, or separately trained with a proved negligible remainder. Training them on the full evaluation data requires a separate argument and is not automatically justified by this interface.

If a block of size \eqn{N_a} has normalized potential-mean contribution \eqn{L_{a,i}}, its pooled row is \eqn{(2a-1)nL_{a,i}/N_a}. Thus the sampling variance is \eqn{\widehat V_1/N_1+\widehat V_0/N_0}. Gaussian multipliers on the pooled rows automatically have this variance. No row or donor is shared between blocks.

This is a distinct estimator from common-score full-sample stabilized PATE and may lose precision. It does not assert validity of a full-sample stabilized PATE estimator with different arm-specific maps. See \code{wm_potential} for identification, moment, geometry and nuisance requirements.
}
\value{A \code{wm_split_pate} object inheriting from \code{wm_match}, with the PATE estimate, variances, standard error, block labels, both complete component fits (including their original row indices) and pooled contribution rows.}
\seealso{\code{\link{wm_potential}}, \code{\link{wm_bootstrap}}}
