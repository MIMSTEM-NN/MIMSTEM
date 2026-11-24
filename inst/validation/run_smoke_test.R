# MIMSTEM installed-package smoke test
# Run after installing the package. This script uses synthetic package data
# and does not run bootstrap inference.

library(MIMSTEM)

message("MIMSTEM version: ", as.character(packageVersion("MIMSTEM")))

tem_dat <- utils::read.csv(
    system.file("extdata", "mims_example_tem.csv", package = "MIMSTEM")
)
s2_dat <- utils::read.csv(
    system.file("extdata", "mims_example_s2tem.csv", package = "MIMSTEM")
)

fit_tem <- mims_tem(
    dat = tem_dat,
    tau = c(2, 4, 6),
    basec = c("age", "sex"),
    bootstrap = FALSE,
    verbose = FALSE
)

fit_ge <- mims_s2tem_ge(
    dat = s2_dat,
    tau = c(2, 4, 6),
    basec = c("age", "sex"),
    bootstrap = FALSE,
    verbose = FALSE
)

fit_pse <- mims_s2tem_pse(
    dat = s2_dat,
    tau = c(2, 4, 6),
    basec = c("age", "sex"),
    bootstrap = FALSE,
    primary_rule = "t2_overrides",
    verbose = FALSE
)

stopifnot(
    all(is.finite(fit_tem$point$TE)),
    all(is.finite(fit_ge$point$TE)),
    all(is.finite(fit_pse$point$TE))
)

endpoint_errors <- c(
    psi0000_vs_psi00 = max(abs(fit_pse$point$psi0000 - fit_ge$point$psi00)),
    psi1000_vs_psi10 = max(abs(fit_pse$point$psi1000 - fit_ge$point$psi10)),
    psi1111_vs_psi11 = max(abs(fit_pse$point$psi1111 - fit_ge$point$psi11))
)

decomposition_errors <- c(
    IE_total = max(abs(
        fit_pse$point$IE_total -
            fit_pse$point$P_A_T2_Y -
            fit_pse$point$P_A_T1_Y -
            fit_pse$point$P_A_T1_T2_Y
    )),
    TE = max(abs(
        fit_pse$point$TE - fit_pse$point$DE - fit_pse$point$IE_total
    ))
)

stopifnot(
    max(endpoint_errors) < 1e-10,
    max(decomposition_errors) < 1e-10
)

cat("\nMIMS-TEM final-time effects:\n")
print(summary(fit_tem, tau = 6)$table, row.names = FALSE)

cat("\nMIMS-S2TEM-GE final-time effects:\n")
print(summary(fit_ge, tau = 6)$table, row.names = FALSE)

cat("\nMIMS-S2TEM-PSE final-time effects:\n")
print(summary(fit_pse, tau = 6)$table, row.names = FALSE)

cat("\nGE/PSE endpoint errors:\n")
print(endpoint_errors)

cat("\nPSE decomposition errors:\n")
print(decomposition_errors)

for (suffix in c("risk", "RMTL")) {
    pp <- fit_pse$point
    pm <- pp[[paste0("PM_path_", suffix)]]
    parts <- pp[[paste0("prop_T1_", suffix)]] +
        pp[[paste0("prop_T2_", suffix)]] + pp[[paste0("prop_T1_T2_", suffix)]]
    stopifnot(isTRUE(all.equal(pm, parts, tolerance = 1e-10)))
    stopifnot(isTRUE(all.equal(pm, fit_ge$point[[paste0("PM_", suffix)]], tolerance = 1e-10)))
}

cat("\nSmoke test passed.\n")
