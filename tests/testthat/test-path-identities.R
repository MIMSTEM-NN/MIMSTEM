test_that("PE endpoint differences telescope exactly", {
    q <- data.frame(
        tau = c(1, 2, 3),
        psi0000 = c(0.10, 0.20, 0.30),
        psi1000 = c(0.08, 0.17, 0.25),
        psi1100 = c(0.07, 0.15, 0.22),
        psi1110 = c(0.065, 0.14, 0.20),
        psi1111 = c(0.05, 0.12, 0.18)
    )

    z <- MIMSTEM:::.pe2_format_risks(q, scale = "risk")

    expect_equal(z$DE, z$psi1000 - z$psi0000, tolerance = 1e-14)
    expect_equal(z$P_A_T2_Y, z$psi1100 - z$psi1000, tolerance = 1e-14)
    expect_equal(z$P_A_T1_Y, z$psi1110 - z$psi1100, tolerance = 1e-14)
    expect_equal(z$P_A_T1_T2_Y, z$psi1111 - z$psi1110, tolerance = 1e-14)
    expect_equal(
        z$IE_total,
        z$P_A_T2_Y + z$P_A_T1_Y + z$P_A_T1_T2_Y,
        tolerance = 1e-14
    )
    expect_equal(z$TE, z$DE + z$IE_total, tolerance = 1e-14)
})

test_that("PE survival display reverses effect signs but not risk endpoints", {
    q <- data.frame(
        tau = c(1, 2),
        psi0000 = c(0.10, 0.20),
        psi1000 = c(0.08, 0.17),
        psi1100 = c(0.07, 0.15),
        psi1110 = c(0.06, 0.14),
        psi1111 = c(0.05, 0.12)
    )

    r <- MIMSTEM:::.pe2_format_risks(q, scale = "risk")
    s <- MIMSTEM:::.pe2_format_risks(q, scale = "survival")

    risk_cols <- c("psi0000", "psi1000", "psi1100", "psi1110", "psi1111")
    expect_equal(r[risk_cols], s[risk_cols], tolerance = 1e-14)
    expect_equal(s$DE, -r$DE, tolerance = 1e-14)
    expect_equal(s$IE_total, -r$IE_total, tolerance = 1e-14)
    expect_equal(s$TE, -r$TE, tolerance = 1e-14)
    expect_equal(s$DE_RMTL, r$DE_RMTL, tolerance = 1e-14)
})
