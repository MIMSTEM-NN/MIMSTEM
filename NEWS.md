# MIMSTEM 1.0.0

- Dissertation release of MIMS-TEM, MIMS-S2TEM-GE, and MIMS-S2TEM-PSE.
- Standardized the public PSE function name to `mims_s2tem_pse()` and updated documentation and examples.
- Added mediation proportions and path-specific contributions relative to the total effect, with paired-bootstrap delta-method Wald inference.
- Retained the supplied point-estimation framework and the shared S2TEM backend.
- Updated validation scripts and source hashes; see validation records for checks still required.

# MIMSTEM 0.1.5

- Added a formal input-data dictionary, including required coding, sequential mediator ordering, same-time event preprocessing, and RMTL grid guidance.
- Added a release checklist separating numerical regression validation, package checks, release metadata, and dissertation reproducibility.
- Added an installed-package release-check script for API, example data, citation, vignette registration, and session information.
- Added SHA-256 hashes of the frozen executable R source files to support release auditing.
- No changes were made to the validated point estimators, path kernel, or bootstrap algorithms.

# MIMSTEM 0.1.4

- Added a formal package workflow vignette covering data structure, all three estimators, bootstrap inference, PE interpretation, diagnostics, and reproducibility.
- Added fast unit tests for PE telescoping identities, risk/survival display conventions, unified bootstrap reproducibility, Wald intervals, RNG restoration, and synthetic data structure.
- Added an installed-package check record documenting the clean 0-error/0-warning status achieved before this release.
- Added a manual installed-package smoke-test script using the synthetic example data and GE/PE identity checks.
- No changes were made to the validated point estimators or bootstrap algorithms.

# MIMSTEM 0.1.3

- Added deterministic synthetic example data files for MIMS-TEM and the two sequential-mediator methods.
- Added package-level documentation and runnable workflow examples (marked `dontrun` to keep routine package checks fast).
- Added an installed user guide covering data structure, point estimation, bootstrap inference, PE decomposition, diagnostics, and reproducibility.
- Added a package citation entry.
- No changes were made to the validated numerical estimators or bootstrap algorithms.

# MIMSTEM 0.1.2

- Removed non-ASCII characters from executable R source code.
- Qualified `head()` and `tail()` calls to satisfy R CMD check.
- No changes were made to the numerical estimators.

# MIMSTEM 0.1.1

- Added standard S3 result methods including `coef()`, `confint()`, `as.data.frame()`, and `plot()`.
- Added validation notes and improved result printing.
- No changes were made to the validated numerical estimators.

# MIMSTEM 0.1.0

- Initial package assembly of MIMS-TEM, MIMS-S2TEM-GE, and MIMS-S2TEM-PE.
