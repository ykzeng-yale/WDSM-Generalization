\name{wm_summary}
\alias{wm_summary}
\alias{print.wm_match}
\alias{summary.wm_match}
\alias{print.wm_wdsm_fit}
\alias{summary.wm_wdsm_fit}
\alias{print.wm_model_fit}
\alias{summary.wm_model_fit}
\alias{print.summary_wm}
\title{Summarize Weighted Matching Results}
\usage{
\method{summary}{wm_match}(object, ...)
\method{summary}{wm_wdsm_fit}(object, ...)
\method{summary}{wm_model_fit}(object, ...)
\method{print}{wm_match}(x, digits = 4L, ...)
\method{print}{wm_wdsm_fit}(x, digits = 4L, ...)
\method{print}{wm_model_fit}(x, digits = 4L, ...)
\method{print}{summary_wm}(x, digits = 4L, ...)
}
\arguments{
  \item{object}{A result from \code{wm_match}, \code{wm_fit},
  \code{wm_wdsm_fit} or \code{wm_model_fit}. The \code{wm_match} methods
  also accept its potential-mean and independent-block subclasses.}
  \item{x}{A fit for a print method, or a \code{summary_wm} result.}
  \item{digits}{Integer number of significant digits to print, from 1 to 22.}
  \item{...}{Additional arguments (ignored).}
}
\description{
Compact presentation of an existing weighted matching fit. These methods
neither fit models nor calculate a new variance or confidence interval.
}
\details{
\code{summary(fit)$donor_usage} counts selected edges by donor arm; it is
not a covariate balance assessment or a weighted target sample size.
Independent-block fits have separate graphs in \code{fit$components}; their
summary does not merge local row identifiers or invent pooled graph counts.
An unavailable SE or interval is kept unavailable. Supplied-map fits do not
gain an interval merely by being summarized. For fitted wrappers the interval
is the one already produced by their selected inference route, at its recorded
confidence level. Summary methods do not verify identification, population
assumptions, numerical-root transfer or sampling coverage.
}
\value{
Summary methods return a \code{summary_wm} list with the retained
estimate, raw estimate, correction, SE and variances in \code{statistics};
sample size, M, target, donor rule and available matching dimensions;
retained inference status/contract; and graph donor-use counts when there
is a single retained graph. Fitted wrappers additionally retain their
original \code{conf.int} and \code{conf.level}. Print methods invisibly
return their input.
}
\seealso{\code{\link{wm_match}}, \code{\link{wm_fit}},
\code{\link{wm_wdsm_fit}}, \code{\link{wm_model_fit}},
\code{\link{wm_bootstrap}}}
\examples{
x <- seq(0, 1, length.out = 16)
z <- rep(0:1, 8)
fit <- wm_fit(1 + 3 * x + 2 * z, z, 1 + x, x, M = 1)
print(fit)
summary(fit)$statistics
summary(fit)$donor_usage
}
