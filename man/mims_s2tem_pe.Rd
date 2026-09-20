\name{mims_s2tem_pe}
\alias{mims_s2tem_pe}
\title{MIMS-S2TEM-PE for path-specific effects}
\usage{
mims_s2tem_pe(dat, tau = 1:10, basec = character(0), a = 1, astar = 0,
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
\description{Fits the validated path-specific estimator. The main effect table reports TE, DE, and three mediator-related path-specific intervention effects.}
\examples{
\dontrun{
path <- system.file("extdata", "mims_example_s2tem.csv", package = "MIMSTEM")
dat <- utils::read.csv(path)
fit <- mims_s2tem_pe(
  dat = dat,
  tau = c(2, 4, 6),
  basec = c("age", "sex"),
  bootstrap = FALSE,
  primary_rule = "t2_overrides"
)
summary(fit, tau = 6)
}
}
