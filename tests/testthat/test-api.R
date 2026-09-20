test_that("public API functions are available", {
    expect_true(is.function(mims_tem))
    expect_true(is.function(mims_s2tem_ge))
    expect_true(is.function(mims_s2tem_pe))
})

test_that("missing requested baseline covariates fail early", {
    dd <- data.frame(A = c(0, 1), x = 1:2)
    expect_error(
        MIMSTEM:::.mims_validate_public_input(dd, 1, "not_here", 1, 0),
        "missing"
    )
})
