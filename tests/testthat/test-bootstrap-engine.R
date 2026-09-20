toy_point_estimator <- function(dat, tau, multiplier = 1) {
    mu <- mean(dat$y) * multiplier
    data.frame(
        tau = tau,
        DE = rep(mu, length(tau)),
        TE = rep(2 * mu, length(tau))
    )
}

test_that("unified bootstrap is reproducible and reports Wald intervals", {
    dat <- data.frame(
        id = seq_len(20),
        A = rep(c(0, 1), each = 10),
        y = seq_len(20) / 20
    )

    b1 <- MIMSTEM:::.bootstrap_mims_wald(
        dat = dat,
        estimator = toy_point_estimator,
        estimator_args = list(multiplier = 1),
        tau = c(1, 2),
        metrics = c("DE", "TE"),
        B = 20,
        seed = 123,
        conf_level = 0.95,
        stratified = FALSE,
        keep_replicates = TRUE,
        verbose = FALSE
    )

    b2 <- MIMSTEM:::.bootstrap_mims_wald(
        dat = dat,
        estimator = toy_point_estimator,
        estimator_args = list(multiplier = 1),
        tau = c(1, 2),
        metrics = c("DE", "TE"),
        B = 20,
        seed = 123,
        conf_level = 0.95,
        stratified = FALSE,
        keep_replicates = TRUE,
        verbose = FALSE
    )

    expect_equal(b1$replicates, b2$replicates, tolerance = 0)
    expect_equal(b1$intervals, b2$intervals, tolerance = 0)
    expect_identical(b1$diagnostics$ci_method, "wald_bootstrap_se")
    expect_false(b1$diagnostics$stratified)
    expect_equal(b1$diagnostics$B_ok, 20L)

    z <- stats::qnorm(0.975)
    expect_equal(
        b1$intervals$lower,
        b1$intervals$estimate - z * b1$intervals$se_boot,
        tolerance = 1e-14
    )
    expect_equal(
        b1$intervals$upper,
        b1$intervals$estimate + z * b1$intervals$se_boot,
        tolerance = 1e-14
    )
})

test_that("unified bootstrap restores the caller RNG state", {
    dat <- data.frame(
        A = rep(c(0, 1), each = 5),
        y = seq_len(10)
    )

    set.seed(999)
    seed_before <- .Random.seed

    invisible(MIMSTEM:::.bootstrap_mims_wald(
        dat = dat,
        estimator = toy_point_estimator,
        tau = 1,
        metrics = c("DE", "TE"),
        B = 5,
        seed = 17,
        verbose = FALSE
    ))

    expect_identical(.Random.seed, seed_before)
})
