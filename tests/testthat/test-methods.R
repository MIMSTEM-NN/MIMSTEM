make_mock_mims_fit <- function() {
    effects <- data.frame(
        tau = rep(c(1, 2), 3),
        estimand = rep(c("DE", "IE", "TE"), each = 2),
        estimate = c(-0.1, -0.2, -0.02, -0.04, -0.12, -0.24),
        se = rep(0.01, 6),
        lower = c(-0.12, -0.22, -0.04, -0.06, -0.14, -0.26),
        upper = c(-0.08, -0.18, 0, -0.02, -0.10, -0.22),
        effect_type = "risk",
        stringsAsFactors = FALSE
    )
    structure(
        list(
            method = "MIMS_S2TEM_GE",
            effects = effects,
            diagnostics = list(),
            settings = list(tau = c(1, 2), scale = "risk", bootstrap = TRUE)
        ),
        class = c("mims_s2tem_ge", "mims_fit")
    )
}

test_that("coef and confint use the final tau by default", {
    x <- make_mock_mims_fit()
    cc <- coef(x)
    expect_equal(unname(cc), c(-0.2, -0.04, -0.24))
    expect_equal(names(cc), c("DE", "IE", "TE"))

    ci <- confint(x)
    expect_equal(dim(ci), c(3, 2))
    expect_equal(rownames(ci), c("DE", "IE", "TE"))
})

test_that("summary can filter tau and estimand", {
    x <- make_mock_mims_fit()
    s <- summary(x, tau = 2, estimand = c("DE", "TE"))
    expect_equal(s$table$estimand, c("DE", "TE"))
    expect_true(all(s$table$tau == 2))
})

test_that("mims fits coerce to the tidy effect table", {
    x <- make_mock_mims_fit()
    d <- as.data.frame(x)
    expect_equal(d, x$effects)
})
