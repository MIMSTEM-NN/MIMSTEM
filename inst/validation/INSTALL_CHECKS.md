# Installed-package validation status

The installed package was checked on Windows with R 4.4.3 using
`devtools::check()` / `R CMD check --as-cran` conventions.

For MIMSTEM 0.1.2, the observed result was:

- 0 errors;
- 0 warnings;
- 1 environment-related note: `unable to verify current time` during the
  future-file-timestamp check.

The timestamp note is external to the package source and does not indicate a
numerical or namespace problem.

Installed-package numerical regression checks also confirmed exact equality
between the packaged GE/PE point estimators and the frozen verified refactors,
exact equality of the three shared GE/PE endpoints, and exact PE telescoping
identities.
