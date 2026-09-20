# Numerical regression validation

The package core was refactored from the dissertation analysis scripts without
changing the estimands or numerical integration rules.

Before packaging, the following checks were completed on the development data:

- MIMS-TEM: point estimates, counterfactual risks/survivals, RMTL effects,
  Cox coefficients, time grid, and selected model terms agreed exactly with
  the validated source estimator (`max error = 0`).
- MIMS-S2TEM-GE: all point estimates and RMTL quantities agreed exactly with
  the validated source estimator (`max error = 0`).
- MIMS-S2TEM-PE: all five counterfactual risks, direct/indirect/total effects,
  all three mediator-related path-specific effects, RMTL/RMST quantities,
  controls, and diagnostics agreed exactly with the validated source estimator
  (`max error = 0`).
- GE and PE bootstrap point estimates and bootstrap standard errors agreed
  exactly with the source implementations in regression checks.
- Installed-package regression checks confirmed that GE and PE point results
  remained identical to the verified refactored source functions, and the PE
  common endpoints matched GE exactly:
  `psi0000 = psi00`, `psi1000 = psi10`, `psi1111 = psi11`.
- PE decomposition identities were exactly satisfied:
  `IE_total = P_A_T2_Y + P_A_T1_Y + P_A_T1_T2_Y` and
  `TE = DE + IE_total`.

These checks document numerical equivalence of the refactor. They are not a
substitute for method assumptions, model diagnostics, or application-specific
validation.
