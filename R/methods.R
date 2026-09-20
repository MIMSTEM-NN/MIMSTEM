# S3 methods ---------------------------------------------------------------

.mims_method_label <- function(x) {
    switch(
        x$method,
        MIMS_TEM = "MIMS-TEM",
        MIMS_S2TEM_GE = "MIMS-S2TEM-GE",
        MIMS_S2TEM_PE = "MIMS-S2TEM-PE",
        x$method
    )
}

.mims_select_effects <- function(object, tau = NULL,
                                 effect_type = c("risk", "RMTL"),
                                 estimand = NULL) {
    effect_type <- match.arg(effect_type)
    tab <- object$effects[object$effects$effect_type == effect_type, , drop = FALSE]

    if (!is.null(tau)) {
        tau <- as.numeric(tau)
        if (length(tau) != 1L || !is.finite(tau)) {
            stop("`tau` must be one finite number.", call. = FALSE)
        }
        available <- sort(unique(tab$tau))
        if (!tau %in% available) {
            stop(
                "Requested tau is not in the fitted evaluation times: ",
                paste(available, collapse = ", "),
                call. = FALSE
            )
        }
        tab <- tab[tab$tau == tau, , drop = FALSE]
    }

    if (!is.null(estimand)) {
        estimand <- as.character(estimand)
        bad <- setdiff(estimand, unique(tab$estimand))
        if (length(bad)) {
            stop(
                "Requested estimand(s) not available: ",
                paste(bad, collapse = ", "),
                call. = FALSE
            )
        }
        tab <- tab[tab$estimand %in% estimand, , drop = FALSE]
    }

    tab
}

#' @export
print.mims_fit <- function(x, ..., digits = 4) {
    cat(.mims_method_label(x), " fit\n", sep = "")
    cat("Evaluation times: ", paste(x$settings$tau, collapse = ", "), "\n", sep = "")
    cat("Effect display scale: ", x$settings$scale, "\n", sep = "")

    if (isTRUE(x$settings$bootstrap)) {
        dg <- x$diagnostics$bootstrap
        cat(
            "Bootstrap: ", dg$B_ok, "/", dg$B_total,
            " successful; ", dg$ci_method,
            "; conf.level=", format(dg$conf_level), "\n",
            sep = ""
        )
    } else {
        cat("Bootstrap: not requested\n")
    }

    if (!is.null(x$diagnostics$path) && nrow(x$diagnostics$path) > 0L &&
        "fraction_clipped" %in% names(x$diagnostics$path)) {
        mx <- suppressWarnings(max(x$diagnostics$path$fraction_clipped, na.rm = TRUE))
        if (is.finite(mx) && mx > 0) {
            cat("Path-kernel diagnostic: maximum clipped fraction = ",
                format(mx, digits = digits), "\n", sep = "")
        }
    }

    tt <- max(x$effects$tau)
    tab <- .mims_select_effects(x, tau = tt, effect_type = "risk")
    cat("\n", x$settings$scale, "-scale effects at tau = ", tt, ":\n", sep = "")
    print(tab[c("estimand", "estimate", "se", "lower", "upper")],
          row.names = FALSE, digits = digits)
    invisible(x)
}

#' @export
summary.mims_fit <- function(object, tau = NULL,
                             effect_type = c("risk", "RMTL"),
                             estimand = NULL, ...) {
    tab <- .mims_select_effects(
        object,
        tau = tau,
        effect_type = effect_type,
        estimand = estimand
    )

    out <- list(
        method = object$method,
        effect_type = match.arg(effect_type),
        table = tab,
        bootstrap = if (!is.null(object$diagnostics$bootstrap)) object$diagnostics$bootstrap else NULL,
        settings = object$settings
    )
    class(out) <- "summary.mims_fit"
    out
}

#' @export
print.summary.mims_fit <- function(x, ..., digits = 4) {
    cat(
        switch(
            x$method,
            MIMS_TEM = "MIMS-TEM",
            MIMS_S2TEM_GE = "MIMS-S2TEM-GE",
            MIMS_S2TEM_PE = "MIMS-S2TEM-PE",
            x$method
        ),
        " summary (", x$effect_type, ")\n",
        sep = ""
    )
    print(x$table[c("tau", "estimand", "estimate", "se", "lower", "upper")],
          row.names = FALSE, digits = digits)
    invisible(x)
}

#' @export
coef.mims_fit <- function(object, tau = NULL,
                          effect_type = c("risk", "RMTL"),
                          estimand = NULL, ...) {
    if (is.null(tau)) tau <- max(object$effects$tau)
    tab <- .mims_select_effects(
        object,
        tau = tau,
        effect_type = effect_type,
        estimand = estimand
    )
    out <- tab$estimate
    names(out) <- tab$estimand
    out
}

#' @export
confint.mims_fit <- function(object, parm = NULL, level = 0.95,
                             tau = NULL,
                             effect_type = c("risk", "RMTL"), ...) {
    if (length(level) != 1L || !is.numeric(level) || !is.finite(level) ||
        level <= 0 || level >= 1) {
        stop("`level` must be one number strictly between 0 and 1.", call. = FALSE)
    }
    if (is.null(tau)) tau <- max(object$effects$tau)

    tab <- .mims_select_effects(
        object,
        tau = tau,
        effect_type = effect_type,
        estimand = parm
    )

    alpha <- 1 - level
    z <- stats::qnorm(1 - alpha / 2)
    lo <- tab$estimate - z * tab$se
    hi <- tab$estimate + z * tab$se

    out <- cbind(lo, hi)
    colnames(out) <- paste0(format(c(alpha / 2, 1 - alpha / 2) * 100, trim = TRUE), "%")
    rownames(out) <- tab$estimand
    out
}

#' @export
as.data.frame.mims_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
    out <- x$effects
    if (!is.null(row.names)) rownames(out) <- row.names
    out
}

#' @export
plot.mims_fit <- function(x,
                          estimand = NULL,
                          effect_type = c("risk", "RMTL"),
                          ci = TRUE,
                          zero_line = TRUE,
                          xlab = "Evaluation time",
                          ylab = NULL,
                          main = NULL,
                          ...) {
    effect_type <- match.arg(effect_type)
    tab <- .mims_select_effects(x, effect_type = effect_type, estimand = estimand)

    if (is.null(estimand)) {
        preferred <- if (x$method == "MIMS_S2TEM_PE") {
            c("TE", "DE", "P_A_T2_Y", "P_A_T1_Y", "P_A_T1_T2_Y")
        } else {
            c("TE", "DE", "IE")
        }
        if (effect_type == "RMTL") preferred <- paste0(preferred, "_RMTL")
        preferred <- preferred[preferred %in% unique(tab$estimand)]
        tab <- tab[tab$estimand %in% preferred, , drop = FALSE]
    }

    ests <- unique(tab$estimand)
    times <- sort(unique(tab$tau))
    mat <- sapply(ests, function(nm) {
        z <- tab[tab$estimand == nm, , drop = FALSE]
        z$estimate[match(times, z$tau)]
    })
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1L)
    colnames(mat) <- ests

    if (is.null(ylab)) {
        ylab <- if (effect_type == "RMTL") "RMTL effect" else paste0(x$settings$scale, "-scale effect")
    }
    if (is.null(main)) main <- paste0(.mims_method_label(x), ": ", effect_type, " effects")

    graphics::matplot(
        times, mat,
        type = "l", lty = seq_len(ncol(mat)),
        xlab = xlab, ylab = ylab, main = main,
        ...
    )

    if (isTRUE(zero_line)) graphics::abline(h = 0, lty = 3)

    # Confidence bands are added only when one estimand is shown; with several
    # curves, separate bands usually obscure the comparison.
    if (isTRUE(ci) && length(ests) == 1L) {
        z <- tab[tab$estimand == ests[1L], , drop = FALSE]
        z <- z[match(times, z$tau), , drop = FALSE]
        if (all(is.finite(z$lower)) && all(is.finite(z$upper))) {
            graphics::lines(times, z$lower, lty = 2)
            graphics::lines(times, z$upper, lty = 2)
        }
    }

    if (length(ests) > 1L) {
        graphics::legend("topright", legend = ests,
                         lty = seq_along(ests), bty = "n")
    }

    invisible(tab)
}
