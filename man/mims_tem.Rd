\name{mims_tem}
\alias{mims_tem}
\title{MIMS-TEM for one time-to-event mediator}
\usage{
mims_tem(dat, tau = 1:10, basec = character(0), a = 1, astar = 0,
  bootstrap = TRUE, B = 200L, seed = 1L, conf_level = 0.95,
  stratified = FALSE, bh_method = "breslow", EMint = FALSE,
  scale = c("risk", "survival"), min_success_fraction = 0.8,
  keep_replicates = FALSE, keep_fit_objects = FALSE, verbose = TRUE)
}
\arguments{
  \item{dat}{Subject-level data frame with canonical MIMS-TEM variables.}
  \item{tau}{Positive evaluation times.}
  \item{basec}{Baseline covariate names.}
  \item{a, astar}{Exposure levels defining the contrast.}
  \item{bootstrap}{Whether to run nonparametric bootstrap inference.}
  \item{B}{Number of bootstrap replicates.}
  \item{seed}{Bootstrap seed.}
  \item{conf_level}{Confidence level for Wald intervals.}
  \item{stratified}{Use exposure-stratified resampling when TRUE.}
  \item{bh_method}{Cox tie method.}
  \item{EMint}{Include the supported exposure-by-mediator-time interaction.}
  \item{scale}{Effect display scale. Counterfactual endpoint columns remain risks.}
  \item{min_success_fraction}{Minimum successful bootstrap fraction.}
  \item{keep_replicates}{Keep bootstrap replicate array.}
  \item{keep_fit_objects}{Keep fitted point-estimator objects when bootstrap is FALSE.}
  \item{verbose}{Show progress.}
}
\value{An object inheriting from \code{mims_fit}.}
\description{Fits the validated single time-to-event mediator marginal interventional multi-state estimator.}
\examples{
\dontrun{
path <- system.file("extdata", "mims_example_tem.csv", package = "MIMSTEM")
dat <- utils::read.csv(path)
fit <- mims_tem(
  dat = dat,
  tau = c(2, 4, 6),
  basec = c("age", "sex"),
  bootstrap = FALSE
)
summary(fit, tau = 6)
}
}
