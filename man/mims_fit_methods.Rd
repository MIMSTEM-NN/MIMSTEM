\name{mims_fit_methods}
\alias{print.mims_fit}
\alias{summary.mims_fit}
\alias{print.summary.mims_fit}
\alias{coef.mims_fit}
\alias{confint.mims_fit}
\alias{as.data.frame.mims_fit}
\alias{plot.mims_fit}
\title{Methods for fitted MIMS objects}
\usage{
\method{print}{mims_fit}(x, ..., digits = 4)

\method{summary}{mims_fit}(object, tau = NULL,
  effect_type = c("risk", "RMTL"), estimand = NULL, ...)

\method{print}{summary.mims_fit}(x, ..., digits = 4)

\method{coef}{mims_fit}(object, tau = NULL,
  effect_type = c("risk", "RMTL"), estimand = NULL, ...)

\method{confint}{mims_fit}(object, parm = NULL, level = 0.95,
  tau = NULL, effect_type = c("risk", "RMTL"), ...)

\method{as.data.frame}{mims_fit}(x, row.names = NULL, optional = FALSE, ...)

\method{plot}{mims_fit}(x, estimand = NULL,
  effect_type = c("risk", "RMTL"), ci = TRUE, zero_line = TRUE,
  xlab = "Evaluation time", ylab = NULL, main = NULL, ...)
}
\arguments{
  \item{x, object}{A fitted object returned by \code{mims_tem()},
  \code{mims_s2tem_ge()}, or \code{mims_s2tem_pe()}.}
  \item{tau}{Optional fitted evaluation time. If omitted by \code{coef()} or
  \code{confint()}, the largest fitted evaluation time is used.}
  \item{effect_type}{Either \code{"risk"} or \code{"RMTL"}.}
  \item{estimand}{Optional character vector of estimands to retain.}
  \item{digits}{Number of printed significant digits.}
  \item{parm}{Optional estimand names for \code{confint()}.}
  \item{level}{Confidence level used to reconstruct Wald intervals from the
  stored bootstrap standard errors.}
  \item{row.names, optional}{Standard arguments for \code{as.data.frame()}.}
  \item{ci}{Logical; for a single plotted estimand, draw stored confidence
  limits when available.}
  \item{zero_line}{Logical; draw a horizontal reference line at zero.}
  \item{xlab, ylab, main}{Plot labels.}
  \item{...}{Additional arguments passed to the relevant method.}
}
\value{
\code{summary()} returns a \code{summary.mims_fit} object. \code{coef()}
returns a named numeric vector. \code{confint()} returns a two-column matrix.
\code{as.data.frame()} returns the tidy effect table. \code{plot()} returns the
plotted effect table invisibly.
}
\description{
S3 methods for inspecting, extracting, and plotting effect estimates from
MIMSTEM fits.
}
