# MIMSTEM

`MIMSTEM` implements three marginal interventional multi-state mediation methods for time-to-event mediators:

- `mims_tem()`: one time-to-event mediator (MIMS-TEM).
- `mims_s2tem_ge()`: two sequential time-to-event mediators, overall effects (MIMS-S2TEM-GE).
- `mims_s2tem_pe()`: two sequential time-to-event mediators, path-specific decomposition (MIMS-S2TEM-PE).


## Frozen dissertation release

Version 1.0.0 is the frozen dissertation release. Its executable numerical core is unchanged from the regression-validated 0.1.x series. Any future edit to files under `R/` should be treated as a numerical-core change and should trigger the regression checks documented under `inst/validation/`.

The package maintainer contact in `DESCRIPTION` has been finalized as `wnn0924@126.com`. This release-metadata update does not alter the numerical core.

## Canonical data names

The current package deliberately keeps the variable names used by the validated dissertation code so that packaging does not alter the estimators.

MIMS-TEM expects `A`, `t_out`, `d_out` (or `yevent`), `t_med`, and `d_med`.

MIMS-S2TEM-GE and MIMS-S2TEM-PE expect `A`, `t_out`, `d_out` (or `yevent`), `t_med1`, `d_med1`, `t_med2`, and `d_med2`.

Baseline covariates are supplied by name through `basec`.


## Synthetic example data

Two deterministic synthetic CSV files are installed with the package. They demonstrate the required input structure and are not empirical dissertation data.

```r
tem_dat <- utils::read.csv(
  system.file("extdata", "mims_example_tem.csv", package = "MIMSTEM")
)

s2_dat <- utils::read.csv(
  system.file("extdata", "mims_example_s2tem.csv", package = "MIMSTEM")
)
```

A longer installed guide is available at `system.file("doc", "USER_GUIDE.md", package = "MIMSTEM")`.

A formal workflow vignette is also included:

```r
vignette("MIMSTEM-workflow", package = "MIMSTEM")
```

## Recommended workflow

Start with a point-estimate fit while checking data coding and model diagnostics:

```r
library(MIMSTEM)

fit_ge <- mims_s2tem_ge(
  dat = dat2,
  tau = c(5, 10),
  basec = c("age", "sex"),
  bootstrap = FALSE
)
summary(fit_ge, tau = 10)
```

After the point fit is satisfactory, request the dissertation-standard individual-level nonparametric bootstrap with Wald intervals:

```r
fit_ge_boot <- mims_s2tem_ge(
  dat = dat2,
  tau = c(5, 10),
  basec = c("age", "sex"),
  bootstrap = TRUE,
  B = 200,
  seed = 2026
)
summary(fit_ge_boot, tau = 10)
confint(fit_ge_boot, tau = 10)
```

For MIMS-S2TEM-PE, the primary mixed-history rule is `"t2_overrides"`:

```r
fit_pe <- mims_s2tem_pe(
  dat = dat2,
  tau = c(5, 10),
  basec = c("age", "sex"),
  primary_rule = "t2_overrides",
  bootstrap = FALSE
)

summary(fit_pe, tau = 10)
coef(fit_pe, tau = 10)
plot(fit_pe)
plot(fit_pe, estimand = "P_A_T1_T2_Y")
```

The main MIMS-S2TEM-PE effect table reports TE, DE, `P_A_T2_Y`, `P_A_T1_Y`, and `P_A_T1_T2_Y`. `IE_total` remains available in `$point` for the decomposition check

`IE_total = P_A_T2_Y + P_A_T1_Y + P_A_T1_T2_Y`.

The common GE/PE endpoints satisfy `psi0000 = psi00`, `psi1000 = psi10`, and `psi1111 = psi11` by construction of the shared backend.

## Result extraction

Every public estimator returns a `mims_fit` object. Common access patterns are:

```r
fit$point          # complete estimator output
fit$endpoints      # counterfactual endpoint curves
fit$effects        # tidy main-effect table
fit$diagnostics    # bootstrap/path diagnostics
fit$settings       # fitted settings

as.data.frame(fit)
coef(fit, tau = 10)
confint(fit, tau = 10)
summary(fit, tau = 10, estimand = c("TE", "DE"))
```

Bootstrap defaults to ordinary individual-level nonparametric resampling (`stratified = FALSE`) with bootstrap standard errors and Wald confidence intervals.

## Important preprocessing requirements

The package intentionally does not silently modify event ordering or missing
values. Same-time mediator/outcome rules, construction of first-event mediator
times, eligibility criteria, administrative censoring, and missing-data handling
must be completed before fitting and should be documented in the analysis.

A compact input data dictionary is installed at:

```r
system.file("doc", "DATA_DICTIONARY.md", package = "MIMSTEM")
```

## Release validation

The installed release checks can be located with:

```r
system.file("validation", "run_release_check.R", package = "MIMSTEM")
```

The validated numerical core is accompanied by a source-file SHA-256 manifest
under `inst/validation/CORE_SHA256.txt`. This is intended as an audit aid: if a
core R file is edited after validation, its hash changes and numerical
regression checks should be repeated before release.
