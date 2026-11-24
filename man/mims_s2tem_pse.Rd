\name{mims_s2tem_pse}
\alias{mims_s2tem_pse}
\title{MIMS-S2TEM-PSE for path-specific effects}
\usage{
mims_s2tem_pse(dat, tau = 1:10, basec = character(0), a = 1, astar = 0,
  bootstrap = TRUE, B = 200L, seed = 1L, conf_level = 0.95,
  stratified = FALSE, bh_method = "breslow", EMint = FALSE,
  scale = c("risk", "survival"), min_success_fraction = 0.8,
  keep_replicates = FALSE, verbose = TRUE, v_chunk_size = 4096L,
  max_grid_points = 200L,
  primary_rule = c("t2_overrides", "require_both"),
  t1_coupling = "rank_preserving", subject_chunk_size = 512L)
}
\arguments{
  \item{dat}{Subject-level data frame with canonical sequential-mediator variables.}
  \item{tau}{Positive evaluation times.}
  \item{basec}{Baseline covariate names.}
  \item{a, astar}{Exposure levels defining the contrast.}
  \item{bootstrap}{Whether to run nonparametric bootstrap inference.}
  \item{B}{Number of bootstrap replicates.}
  \item{seed}{Bootstrap seed.}
  \item{conf_level}{Confidence level for Wald intervals.}
  \item{stratified}{Use exposure-stratified resampling when TRUE.}
  \item{bh_method}{Cox tie method.}
  \item{EMint}{Include supported exposure-by-mediator-time interactions.}
  \item{scale}{Effect display scale.}
  \item{min_success_fraction}{Minimum successful bootstrap fraction.}
  \item{keep_replicates}{Keep bootstrap replicate array.}
  \item{verbose}{Show progress.}
  \item{v_chunk_size}{Internal numerical integration chunk size.}
  \item{max_grid_points}{Maximum adaptive grid size.}
  \item{primary_rule}{Mixed-history rule; \code{"t2_overrides"} is primary.}
  \item{t1_coupling}{T1 coupling; currently only \code{"rank_preserving"}.}
  \item{subject_chunk_size}{Subject block size for path-kernel calculations.}
}
\value{An object inheriting from \code{mims_fit}.}
\description{Fits the validated path-specific estimator. The main effect table reports TE, DE, overall indirect effect, three mediator-related path-specific intervention effects, their RMTL analogues, and mediation proportions.}
\details{
Mediation proportions include \code{PM_path_risk}, \code{PM_path_RMTL}, and
\code{prop_T1_risk}, \code{prop_T2_risk}, \code{prop_T1_T2_risk} and their RMTL analogues.
All use the total effect as denominator. The survival display uses the corresponding survival ratio names. Path-specific effects depend on the prespecified switch order.
Effect confidence intervals use bootstrap standard errors and Wald intervals;
ratio confidence intervals use first-order delta-method standard errors based
on paired bootstrap covariance of numerator and denominator. Ratios are NA
when the absolute total effect is at most 1e-6 and may lie outside [0, 1].
}
\examples{
\dontrun{
path <- system.file("extdata", "mims_example_s2tem.csv", package = "MIMSTEM")
dat <- utils::read.csv(path)
fit <- mims_s2tem_pse(
  dat = dat,
  tau = c(2, 4, 6),
  basec = c("age", "sex"),
  bootstrap = FALSE,
  primary_rule = "t2_overrides"
)
summary(fit, tau = 6)
}
}
