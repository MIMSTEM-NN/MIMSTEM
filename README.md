# MIMSTEM

`MIMSTEM` implements three marginal interventional multi-state mediation methods for time-to-event mediators:

- `mims_tem()`: one time-to-event mediator (MIMS-TEM).
- `mims_s2tem_ge()`: two sequential time-to-event mediators, overall effects (MIMS-S2TEM-GE).
- `mims_s2tem_pse()`: two sequential time-to-event mediators, path-specific decomposition (MIMS-S2TEM-PSE).


## Dissertation release

Version 1.0.0 includes the PSE naming update and mediation proportion inference. Record the Git commit SHA alongside the version for reproducibility. See `inst/validation/` for historical validation and the checks required for this revision.

## Canonical data names

The current package deliberately keeps the variable names used by the validated dissertation code so that packaging does not alter the estimators.

MIMS-TEM expects `A`, `t_out`, `d_out` (or `yevent`), `t_med`, and `d_med`.

MIMS-S2TEM-GE and MIMS-S2TEM-PSE expect `A`, `t_out`, `d_out` (or `yevent`), `t_med1`, `d_med1`, `t_med2`, and `d_med2`.

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

A formal workflow vignette is also included:

```r
vignette("MIMSTEM-workflow", package = "MIMSTEM")
```

## Recommended workflow

Start with a point-estimate fit while checking data coding and model diagnostics:

```r
library(MIMSTEM)

fit_ge <- mims_s2tem_ge(
  dat = s2_dat,
  tau = c(5, 10),
  basec = c("age", "sex"),
  bootstrap = FALSE
)
summary(fit_ge, tau = 10)
```

After the point fit is satisfactory, request the dissertation-standard individual-level nonparametric bootstrap with Wald intervals:

```r
fit_ge_boot <- mims_s2tem_ge(
  dat = s2_dat,
  tau = c(5, 10),
  basec = c("age", "sex"),
  bootstrap = TRUE,
  B = 200,
  seed = 2026
)
summary(fit_ge_boot, tau = 10)
confint(fit_ge_boot, tau = 10)
```

For MIMS-S2TEM-PSE, the primary mixed-history rule is `"t2_overrides"`:

```r
fit_pse <- mims_s2tem_pse(
  dat = s2_dat,
  tau = c(5, 10),
  basec = c("age", "sex"),
  primary_rule = "t2_overrides",
  bootstrap = FALSE
)

summary(fit_pse, tau = 10)
coef(fit_pse, tau = 10)
plot(fit_pse)
plot(fit_pse, estimand = "P_A_T1_T2_Y")
```

The main MIMS-S2TEM-PSE effect table reports TE, DE, `IE_total`, `P_A_T2_Y`, `P_A_T1_Y`, and `P_A_T1_T2_Y`, their RMTL analogues, and mediation proportions. `IE_total` is also available in `$point` for the decomposition check

`IE_total = P_A_T2_Y + P_A_T1_Y + P_A_T1_T2_Y`.

The common GE/PSE endpoints satisfy `psi0000 = psi00`, `psi1000 = psi10`, and `psi1111 = psi11` by construction of the shared backend.

## Mediation proportions

For TEM and GE, `PM_risk` and `PM_RMTL` are the indirect effect divided by the total effect on the corresponding scale. For PSE, `PM_path_risk` and `PM_path_RMTL` use the overall indirect effect as numerator. `prop_T1_risk`, `prop_T2_risk`, and `prop_T1_T2_risk` (and their `_RMTL` analogues) divide the respective path-specific effect by the total effect, not by the overall indirect effect. With `scale = "survival"`, the `_risk` ratio names become `_survival`; RMTL names are unchanged.

Ratios are available even with `bootstrap = FALSE`. With bootstrap inference, ordinary effects use bootstrap standard errors and Wald intervals; ratios use first-order delta-method standard errors from paired bootstrap covariance of the numerator and denominator and Wald intervals centered on the original-sample ratio. Ratios are not percentages and are not constrained to [0, 1]. Ratios are undefined (NA) when the absolute total effect is at most 1e-6. Inspect `fit$diagnostics$bootstrap$ratio_diagnostics` for denominator instability and paired replicate counts. Intervals are unavailable when the minimum bootstrap success requirement is not met.

The default `B = 200` is useful for demonstrations; increase it for formal inference (e.g., `B = 1000`). Path-specific effects depend on the prespecified mechanism-switching order and coupling assumptions.

## Result extraction

Every public estimator returns a `mims_fit` object. Common access patterns are:

```r
fit <- fit_ge_boot
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

## Release validation

The installed release checks can be located with:

```r
system.file("validation", "run_release_check.R", package = "MIMSTEM")
```

The current source snapshot is accompanied by a source-file SHA-256 manifest
under `inst/validation/CORE_SHA256.txt`. This is intended as an audit aid: if a
core R file is edited after validation, its hash changes and numerical
regression checks should be repeated before release.
