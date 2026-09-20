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

fit_pe <- mims_s2tem_pe(
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
    all(is.finite(fit_pe$point$TE))
)

endpoint_errors <- c(
    psi0000_vs_psi00 = max(abs(fit_pe$point$psi0000 - fit_ge$point$psi00)),
    psi1000_vs_psi10 = max(abs(fit_pe$point$psi1000 - fit_ge$point$psi10)),
    psi1111_vs_psi11 = max(abs(fit_pe$point$psi1111 - fit_ge$point$psi11))
)

decomposition_errors <- c(
    IE_total = max(abs(
        fit_pe$point$IE_total -
            fit_pe$point$P_A_T2_Y -
            fit_pe$point$P_A_T1_Y -
            fit_pe$point$P_A_T1_T2_Y
    )),
    TE = max(abs(
        fit_pe$point$TE - fit_pe$point$DE - fit_pe$point$IE_total
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

cat("\nMIMS-S2TEM-PE final-time effects:\n")
print(summary(fit_pe, tau = 6)$table, row.names = FALSE)

cat("\nGE/PE endpoint errors:\n")
print(endpoint_errors)

cat("\nPE decomposition errors:\n")
print(decomposition_errors)

cat("\nSmoke test passed.\n")
