test_that("synthetic example data files are installed and structurally valid", {
    p1 <- system.file("extdata", "mims_example_tem.csv", package = "MIMSTEM")
    p2 <- system.file("extdata", "mims_example_s2tem.csv", package = "MIMSTEM")
    expect_true(nzchar(p1))
    expect_true(nzchar(p2))
    expect_true(file.exists(p1))
    expect_true(file.exists(p2))

    d1 <- utils::read.csv(p1)
    d2 <- utils::read.csv(p2)

    expect_true(all(c("id", "A", "t_out", "d_out", "t_med", "d_med") %in% names(d1)))
    expect_true(all(c(
        "id", "A", "t_out", "d_out", "t_med1", "d_med1", "t_med2", "d_med2"
    ) %in% names(d2)))

    expect_setequal(unique(d1$A), c(0, 1))
    expect_setequal(unique(d2$A), c(0, 1))
    expect_true(all(d1$d_out %in% c(0, 1)))
    expect_true(all(d1$d_med %in% c(0, 1)))
    expect_true(all(d2$d_out %in% c(0, 1)))
    expect_true(all(d2$d_med1 %in% c(0, 1)))
    expect_true(all(d2$d_med2 %in% c(0, 1)))
    expect_true(all(d2$d_med2 == 0 | d2$d_med1 == 1))
    expect_true(all(d2$d_med2 == 0 | d2$t_med2 > d2$t_med1))
})
