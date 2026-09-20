\name{MIMSTEM-package}
\alias{MIMSTEM}
\alias{MIMSTEM-package}
\title{MIMSTEM: Marginal Interventional Multi-State Time-to-Event Mediation}
\description{
MIMSTEM implements three marginal interventional mediation estimators for
one or two sequential time-to-event mediators using multi-state models.
}
\details{
The public estimators are \code{mims_tem()}, \code{mims_s2tem_ge()}, and
\code{mims_s2tem_pe()}. The package supports individual-level nonparametric
bootstrap standard errors and Wald confidence intervals.

Synthetic data files illustrating the required column structure are installed
under \code{extdata}. They are demonstration data only and are not empirical
dissertation data. A package workflow vignette is available through
\code{vignette("MIMSTEM-workflow", package = "MIMSTEM")}.
}
\seealso{
\code{\link{mims_tem}}, \code{\link{mims_s2tem_ge}},
\code{\link{mims_s2tem_pe}}
}
