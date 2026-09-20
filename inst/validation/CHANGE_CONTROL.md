# MIMSTEM numerical-core change control

Version 1.0.0 is the frozen dissertation release.

## Files treated as numerical core

All executable files under `R/` are treated as part of the numerical core. The SHA-256 manifest in `CORE_SHA256.txt` records their frozen contents.

## Changes that require numerical regression validation

Repeat the frozen point-estimator and bootstrap regression checks if any executable `R/*.R` file changes, including apparently minor changes to numerical integration, grid construction, transition-model design matrices, path-history construction, bootstrap resampling, or result transformations.

For MIMS-S2TEM-GE and MIMS-S2TEM-PE, also re-check the common endpoint identities `psi0000 = psi00`, `psi1000 = psi10`, and `psi1111 = psi11`, together with the PE telescoping identities.

## Changes that do not by themselves require numerical regression validation

Pure edits to prose documentation, maintainer contact information, citation text, NEWS, or other non-executable metadata do not alter the numerical core. They should still be followed by package installation and `R CMD check`.

## Release rule

Do not replace the archived dissertation source tarball after a core-code change without repeating the regression validation and recording the new package version.
