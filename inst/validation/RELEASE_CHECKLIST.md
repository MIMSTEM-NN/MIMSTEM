# MIMSTEM release checklist

This checklist separates numerical-method validation from package-engineering
checks. A release should not be tagged until all applicable items are complete.

## Frozen numerical core

- MIMS-TEM point estimator agrees with the frozen validated source implementation.
- MIMS-S2TEM-GE point estimator agrees with the frozen validated source implementation.
- MIMS-S2TEM-PE point estimator agrees with the frozen validated source implementation.
- GE/PE common endpoints agree exactly on the same fitted settings.
- PE telescoping identities are satisfied.
- Unified bootstrap sampling and bootstrap standard errors agree with the validated refactor.
- Any intended inference change (for example, percentile to Wald intervals) is documented separately from point-estimator validation.

## Package behavior

- Package installs from the source tarball.
- `library(MIMSTEM)` loads without warning or error.
- Public exports are exactly `mims_tem`, `mims_s2tem_ge`, and `mims_s2tem_pe`.
- Synthetic example files are installed and readable.
- `vignette("MIMSTEM-workflow", package = "MIMSTEM")` is available.
- Installed smoke test ends with `Smoke test passed.`
- `devtools::check()` / `R CMD check --as-cran` reports 0 errors and 0 warnings.
- Environment-only timestamp notes are recorded separately and not treated as package-code failures.

## Release metadata

- Maintainer email finalized in `DESCRIPTION`: `wnn0924@126.com`.
- Confirm author spelling and roles in `Authors@R`.
- Update `inst/CITATION` to the release version.
- Update `NEWS.md`.
- Record the R version/platform used for the final package check.
- Archive the source tarball used in the dissertation or supplementary materials.

## Dissertation reproducibility

- Report the exact MIMSTEM version in Methods or the software appendix.
- Archive analysis scripts separately from the package source.
- Preserve the preprocessing rules used to construct mediator and outcome times.
- Preserve the bootstrap seed and number of bootstrap replicates.
