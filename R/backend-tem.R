# ============================================================
# MIMS-TEM point estimator refactor
# Version: 2026-09-18
# Numerical target: interventional_cmaverse_marginal()
# ============================================================

.mims_tem_compute_rmtl <- function(tau, psi, method = c("trapezoid", "left")) {
    method <- match.arg(method)
    tau <- as.numeric(tau); psi <- as.numeric(psi)
    if (length(tau) != length(psi)) stop("tau and psi must have the same length.")
    if (any(!is.finite(tau))) stop("tau must be finite.")
    if (is.unsorted(tau, strictly = TRUE)) stop("tau must be strictly increasing.")
    n <- length(tau); out <- numeric(n)
    if (n == 0L) return(out)
    t_prev <- 0; p_prev <- 0
    for (k in seq_len(n)) {
        dt <- tau[k] - t_prev
        if (dt < 0) stop("tau must be increasing.")
        incr <- if (method == "left") p_prev * dt else 0.5 * (p_prev + psi[k]) * dt
        out[k] <- if (k == 1L) incr else out[k - 1L] + incr
        t_prev <- tau[k]; p_prev <- psi[k]
    }
    out
}

.mims_tem_add_rmtl <- function(out, integration = c("trapezoid", "left")) {
    integration <- match.arg(integration)
    need <- c("tau", "psi00", "psi10", "psi11")
    miss <- setdiff(need, names(out))
    if (length(miss)) stop("Missing columns for RMTL: ", paste(miss, collapse = ", "))
    tau <- out$tau
    out$rmtl00 <- .mims_tem_compute_rmtl(tau, out$psi00, method = integration)
    out$rmtl10 <- .mims_tem_compute_rmtl(tau, out$psi10, method = integration)
    out$rmtl11 <- .mims_tem_compute_rmtl(tau, out$psi11, method = integration)
    out$DE_RMTL <- out$rmtl10 - out$rmtl00
    out$IE_RMTL <- out$rmtl11 - out$rmtl10
    out$TE_RMTL <- out$rmtl11 - out$rmtl00
    out
}

.mims_tem_prepare <- function(dat) {
    dat <- as.data.frame(dat)
    if (!all(c("A", "t_out") %in% names(dat))) stop("dat must contain A and t_out (plus d_out or yevent).")
    if (!("yevent" %in% names(dat))) {
        if (!("d_out" %in% names(dat))) stop("Missing d_out or yevent.")
        dat$yevent <- as.integer(dat$d_out == 1L)
    }
    if (!all(c("t_med_eff", "mevent_eff") %in% names(dat))) {
        if (!all(c("t_med", "d_med") %in% names(dat))) stop("Missing t_med_eff/mevent_eff and raw t_med/d_med.")
        med_happens <- with(dat, d_med == 1L & is.finite(t_med) & is.finite(t_out) & (t_med < t_out))
        dat$mevent_eff <- as.integer(med_happens)
        dat$t_med_eff <- ifelse(dat$mevent_eff == 1L, dat$t_med, dat$t_out)
    }
    dat$A <- as.numeric(dat$A)
    dat
}

.mims_tem_H_at <- function(Hdf, t) {
    if (nrow(Hdf) == 0L) return(0)
    idx <- findInterval(t, Hdf$time)
    if (idx <= 0L) return(0)
    Hdf$Haz[idx]
}

.mims_tem_H_align <- function(Hdf, tgrid) {
    if (nrow(Hdf) == 0L) return(rep(0, length(tgrid)))
    idx <- findInterval(tgrid, Hdf$time)
    out <- numeric(length(tgrid)); ok <- idx > 0L
    out[ok] <- Hdf$Haz[idx[ok]]
    out
}

.mims_tem_safe_exp <- function(x) exp(pmin(pmax(x, -50), 50))
.mims_tem_escape_regex <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)

.mims_tem_pick_base_terms <- function(v, cols, all_basec = character(0)) {
    pat <- paste0("^", .mims_tem_escape_regex(v), ".*\\.[123]$")
    out <- grep(pat, cols, value = TRUE)
    longer_vars <- setdiff(all_basec[nchar(all_basec) > nchar(v)], v)
    if (length(longer_vars) > 0L && length(out) > 0L) {
        for (w in longer_vars) {
            pat_w <- paste0("^", .mims_tem_escape_regex(w), ".*\\.[123]$")
            out <- out[!grepl(pat_w, out)]
        }
    }
    out
}

.mims_tem_align_to_coef_cols <- function(df, coef_cols) {
    miss <- setdiff(coef_cols, colnames(df))
    if (length(miss)) for (nm in miss) df[[nm]] <- 0
    as.matrix(df[, coef_cols, drop = FALSE])
}

.fit_tem_backend <- function(dat, tau = 1:10, basec = character(0), a = 1, astar = 0,
                             EMint = FALSE, bh_method = "breslow", verbose = FALSE) {
    dat <- .mims_tem_prepare(dat)
    tau <- sort(unique(as.numeric(tau)))
    basec <- intersect(basec, names(dat))
    if (!length(tau) || any(!is.finite(tau))) stop("tau must contain finite evaluation times.")
    if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required.")
    if (!requireNamespace("mstate", quietly = TRUE)) stop("Package 'mstate' is required.")

    trans <- mstate::transMat(x = list(c(2,3), c(3), c()), names = c("A", "t_med_eff", "t_out"))
    covs_df <- c("A", "t_med_eff", basec)
    mstate_dat <- mstate::msprep(time = c(NA, "t_med_eff", "t_out"),
                                status = c(NA, "mevent_eff", "yevent"),
                                data = dat, trans = trans, keep = covs_df)
    mstate_dat <- mstate::expand.covs(mstate_dat, covs_df, append = TRUE, longnames = FALSE)

    base_terms <- character(0)
    if (length(basec) > 0L) {
        base_terms <- unique(unlist(lapply(basec, .mims_tem_pick_base_terms,
                                           cols = colnames(mstate_dat), all_basec = basec)))
    }
    terms <- c(paste0("A.", 1:3), "t_med_eff.3", base_terms)
    if (isTRUE(EMint)) terms <- c(terms, "A.3:t_med_eff.3")
    fml <- stats::as.formula(paste0("survival::Surv(Tstart, Tstop, status) ~ ",
                                    paste(terms, collapse = " + "), " + strata(trans)"))
    fit <- survival::coxph(fml, data = mstate_dat, method = bh_method, x = TRUE)
    b <- stats::coef(fit)
    if (length(b) == 0L) stop("coxph estimated no coefficients.")
    if (anyNA(b)) {
        if (isTRUE(verbose)) warning(sprintf("coxph produced %d NA coefficients; setting them to 0.", sum(is.na(b))))
        b[is.na(b)] <- 0
    }

    mf <- stats::model.frame(fit)
    need_cols <- setdiff(colnames(mf), c("(weights)", "(offset)"))
    need_cols <- setdiff(need_cols, "Surv(Tstart, Tstop, status)")
    need_cols <- intersect(need_cols, colnames(mstate_dat))
    newd0 <- mstate_dat[0, need_cols, drop = FALSE]
    newd0 <- newd0[rep(1, 3), , drop = FALSE]
    newd0[,] <- 0; newd0$trans <- 1:3; newd0$strata <- 1:3
    attr(newd0, "trans") <- trans; class(newd0) <- c("msdata", "data.frame")
    Haz0 <- mstate::msfit(fit, newdata = newd0, trans = trans)$Haz
    H01 <- subset(Haz0, trans == 1)[, c("time", "Haz")]
    H02 <- subset(Haz0, trans == 2)[, c("time", "Haz")]
    H03 <- subset(Haz0, trans == 3)[, c("time", "Haz")]
    H01 <- H01[!duplicated(H01$time), , drop = FALSE]
    H02 <- H02[!duplicated(H02$time), , drop = FALSE]
    H03 <- H03[!duplicated(H03$time), , drop = FALSE]
    if (nrow(H01) < 1L) stop("Transition-1 cumulative hazard grid is empty.")

    time_vec <- sort(unique(c(0, as.numeric(H01$time), as.numeric(H02$time), as.numeric(H03$time), tau)))
    time_vec <- time_vec[is.finite(time_vec)]; time_vec <- time_vec[time_vec >= 0]
    if (length(time_vec) < 2L) stop("Joint time grid is too short.")
    t_left <- time_vec[-length(time_vec)]; t_right <- time_vec[-1L]
    dt_base <- diff(time_vec)
    if (any(dt_base <= 0)) stop("Joint time grid is not strictly increasing.")
    H01_grid <- .mims_tem_H_align(H01, time_vec)
    H02_grid <- .mims_tem_H_align(H02, time_vec)
    H03_grid <- .mims_tem_H_align(H03, time_vec)
    h01_base <- diff(H01_grid) / dt_base
    H01_u <- H01_grid[-length(H01_grid)]
    H02_u <- H02_grid[-length(H02_grid)]
    H03_u <- H03_grid[-length(H03_grid)]

    coef_cols <- names(b); n <- nrow(dat); keep_vars <- unique(c("A", "t_med_eff", basec))
    proto_all <- dat[rep(seq_len(n), each = 3), keep_vars, drop = FALSE]
    proto_all$id <- rep(seq_len(n), each = 3)
    proto_all$trans <- factor(rep(1:3, times = n), levels = 1:3)
    proto_all$strata <- factor(rep(1:3, times = n), levels = 1:3)
    for (nm in basec) {
        if (nm %in% names(proto_all) && is.character(proto_all[[nm]])) {
            proto_all[[nm]] <- factor(proto_all[[nm]], levels = unique(dat[[nm]]))
        }
    }
    attr(proto_all, "trans") <- trans; class(proto_all) <- c("msdata", "data.frame")
    proto_exp <- mstate::expand.covs(proto_all, covs_df, append = TRUE, longnames = FALSE)
    Xall <- .mims_tem_align_to_coef_cols(proto_exp, coef_cols)
    X1 <- Xall[which(proto_all$trans == 1), , drop = FALSE]
    X2 <- Xall[which(proto_all$trans == 2), , drop = FALSE]
    X3 <- Xall[which(proto_all$trans == 3), , drop = FALSE]

    colA1 <- "A.1"; colA2 <- "A.2"; colA3 <- "A.3"; colM3 <- "t_med_eff.3"; colInt <- "A.3:t_med_eff.3"
    set_A <- function(X, which_col, Aval) {
        if (colA1 %in% colnames(X)) X[, colA1] <- 0
        if (colA2 %in% colnames(X)) X[, colA2] <- 0
        if (colA3 %in% colnames(X)) X[, colA3] <- 0
        if (which_col %in% colnames(X)) X[, which_col] <- Aval
        X
    }
    X3_const <- X3
    if (colM3 %in% colnames(X3_const)) X3_const[, colM3] <- 0
    if (colInt %in% colnames(X3_const)) X3_const[, colInt] <- 0
    betaM <- if (colM3 %in% names(b)) b[[colM3]] else 0
    betaInt <- if (colInt %in% names(b)) b[[colInt]] else 0

    eta1_astar <- as.numeric(set_A(X1, colA1, astar) %*% b)
    eta1_a <- as.numeric(set_A(X1, colA1, a) %*% b)
    eta2_astar <- as.numeric(set_A(X2, colA2, astar) %*% b)
    eta2_a <- as.numeric(set_A(X2, colA2, a) %*% b)
    eta3c_astar <- as.numeric(set_A(X3_const, colA3, astar) %*% b)
    eta3c_a <- as.numeric(set_A(X3_const, colA3, a) %*% b)

    list(version = "tem_backend_2026-09-18_v1", dat = dat, tau = tau, basec = basec,
         a = a, astar = astar, EMint = EMint, bh_method = bh_method, fit = fit,
         mstate_dat = mstate_dat, trans = trans, base_terms = base_terms,
         hazards = list(H01 = H01, H02 = H02, H03 = H03),
         grid = list(time_vec = time_vec, t_left = t_left, t_right = t_right,
                     h01_base = h01_base, H01_u = H01_u, H02_u = H02_u, H03_u = H03_u,
                     max_grid = max(t_right)),
         mult = list(e1_astar = .mims_tem_safe_exp(eta1_astar), e1_a = .mims_tem_safe_exp(eta1_a),
                     e2_astar = .mims_tem_safe_exp(eta2_astar), e2_a = .mims_tem_safe_exp(eta2_a),
                     e3c_astar = .mims_tem_safe_exp(eta3c_astar), e3c_a = .mims_tem_safe_exp(eta3c_a)),
         beta = list(betaM = betaM, betaInt = betaInt))
}

.predict_tem <- function(backend, scale = c("risk", "survival"), keep_fit_objects = TRUE, verbose = FALSE) {
    scale <- match.arg(scale)
    if (!is.list(backend) || !identical(backend$version, "tem_backend_2026-09-18_v1")) stop("Invalid TEM backend object.")
    tau <- backend$tau; a <- backend$a; astar <- backend$astar
    H01 <- backend$hazards$H01; H02 <- backend$hazards$H02; H03 <- backend$hazards$H03
    gr <- backend$grid; t_left <- gr$t_left; t_right <- gr$t_right; h01_base <- gr$h01_base
    H01_u <- gr$H01_u; H02_u <- gr$H02_u; H03_u <- gr$H03_u; max_grid <- gr$max_grid
    e1_astar <- backend$mult$e1_astar; e1_a <- backend$mult$e1_a
    e2_astar <- backend$mult$e2_astar; e2_a <- backend$mult$e2_a
    e3c_astar <- backend$mult$e3c_astar; e3c_a <- backend$mult$e3c_a
    betaM <- backend$beta$betaM; betaInt <- backend$beta$betaInt

    out <- data.frame(tau = tau, DE = NA_real_, IE = NA_real_, TE = NA_real_,
                      psi00 = NA_real_, psi10 = NA_real_, psi11 = NA_real_,
                      S00 = NA_real_, S10 = NA_real_, S11 = NA_real_)
    if (any(tau > max_grid) && isTRUE(verbose)) warning("tau exceeds the joint time grid.")

    for (kk in seq_along(tau)) {
        tt <- min(tau[kk], max_grid)
        H01_tt <- .mims_tem_H_at(H01, tt); H02_tt <- .mims_tem_H_at(H02, tt); H03_tt <- .mims_tem_H_at(H03, tt)
        idx <- which(t_left < tt)
        if (length(idx) == 0L) {
            psi00 <- psi10 <- psi11 <- 0
        } else {
            len_use <- pmin(tt, t_right[idx]) - t_left[idx]
            h01_use <- h01_base[idx]; tL_use <- t_left[idx]
            H01u_use <- H01_u[idx]; H02u_use <- H02_u[idx]; H03u_use <- H03_u[idx]
            dH03_use <- pmax(H03_tt - H03u_use, 0)

            P00_i <- exp(-(H01_tt * e1_astar + H02_tt * e2_astar))
            e3_mult_u <- .mims_tem_safe_exp((betaM + betaInt * astar) * tL_use)
            surv0_u <- exp(-(outer(H01u_use, e1_astar, "*") + outer(H02u_use, e2_astar, "*")))
            haz1_u <- outer(h01_use, e1_astar, "*")
            tail_surv <- exp(-outer(dH03_use * e3_mult_u, e3c_astar, "*"))
            P01 <- colSums(surv0_u * haz1_u * tail_surv * len_use)
            S00_i <- P00_i + P01

            P00g_i <- exp(-(H01_tt * e1_astar + H02_tt * e2_a))
            e3_mult_u_g <- .mims_tem_safe_exp((betaM + betaInt * a) * tL_use)
            surv0_u_g <- exp(-(outer(H01u_use, e1_astar, "*") + outer(H02u_use, e2_a, "*")))
            haz1_u_g <- outer(h01_use, e1_astar, "*")
            tail_surv_g <- exp(-outer(dH03_use * e3_mult_u_g, e3c_a, "*"))
            P01g <- colSums(surv0_u_g * haz1_u_g * tail_surv_g * len_use)
            S10_i <- P00g_i + P01g

            P10_i <- exp(-(H01_tt * e1_a + H02_tt * e2_a))
            surv0_u_a <- exp(-(outer(H01u_use, e1_a, "*") + outer(H02u_use, e2_a, "*")))
            haz1_u_a <- outer(h01_use, e1_a, "*")
            tail_surv_a <- exp(-outer(dH03_use * e3_mult_u_g, e3c_a, "*"))
            P01a <- colSums(surv0_u_a * haz1_u_a * tail_surv_a * len_use)
            S11_i <- P10_i + P01a

            psi00 <- mean(1 - S00_i); psi10 <- mean(1 - S10_i); psi11 <- mean(1 - S11_i)
        }
        IE <- psi11 - psi10; DE <- psi10 - psi00; TE <- psi11 - psi00
        if (scale == "survival") { IE <- -IE; DE <- -DE; TE <- -TE }
        out$IE[kk] <- IE; out$DE[kk] <- DE; out$TE[kk] <- TE
        out$psi00[kk] <- psi00; out$psi10[kk] <- psi10; out$psi11[kk] <- psi11
        out$S00[kk] <- 1 - psi00; out$S10[kk] <- 1 - psi10; out$S11[kk] <- 1 - psi11
    }

    if (isTRUE(keep_fit_objects)) {
        attr(out, "fit") <- backend$fit
        attr(out, "mstate_dat") <- backend$mstate_dat
    }
    attr(out, "time_grid") <- list(grid = backend$grid$time_vec, max_grid = backend$grid$max_grid,
                                    max_H01_time = if (nrow(H01) > 0L) max(H01$time, na.rm = TRUE) else NA_real_,
                                    max_H02_time = if (nrow(H02) > 0L) max(H02$time, na.rm = TRUE) else NA_real_,
                                    max_H03_time = if (nrow(H03) > 0L) max(H03$time, na.rm = TRUE) else NA_real_)
    attr(out, "base_terms") <- backend$base_terms
    .mims_tem_add_rmtl(out, integration = "trapezoid")
}

mims_tem_point <- function(dat, tau = 1:10, basec = character(0), a = 1, astar = 0,
                           EMint = FALSE, bh_method = "breslow", scale = c("risk", "survival"),
                           verbose = FALSE, keep_fit_objects = TRUE) {
    scale <- match.arg(scale)
    backend <- .fit_tem_backend(dat = dat, tau = tau, basec = basec, a = a, astar = astar,
                                EMint = EMint, bh_method = bh_method, verbose = verbose)
    .predict_tem(backend = backend, scale = scale, keep_fit_objects = keep_fit_objects, verbose = verbose)
}

verify_mims_tem_equivalence <- function(dat, legacy_fun, tau = 1:10, basec = character(0),
                                        a = 1, astar = 0, EMint = FALSE, bh_method = "breslow",
                                        scale = "risk", tolerance = 1e-12) {
    if (!is.function(legacy_fun)) stop("legacy_fun must be a function.")
    old <- legacy_fun(dat = dat, tau = tau, basec = basec, a = a, astar = astar,
                      EMint = EMint, bh_method = bh_method, scale = scale, verbose = FALSE)
    new <- mims_tem_point(dat = dat, tau = tau, basec = basec, a = a, astar = astar,
                          EMint = EMint, bh_method = bh_method, scale = scale,
                          verbose = FALSE, keep_fit_objects = TRUE)
    metrics <- c("tau", "psi00", "psi10", "psi11", "S00", "S10", "S11",
                 "DE", "IE", "TE", "rmtl00", "rmtl10", "rmtl11", "DE_RMTL", "IE_RMTL", "TE_RMTL")
    by_metric <- do.call(rbind, lapply(metrics, function(nm) {
        if (!nm %in% names(old) || !nm %in% names(new)) return(data.frame(metric = nm, max_abs_error = Inf))
        data.frame(metric = nm, max_abs_error = max(abs(old[[nm]] - new[[nm]]), na.rm = TRUE))
    }))
    old_grid <- attr(old, "time_grid", exact = TRUE); new_grid <- attr(new, "time_grid", exact = TRUE)
    grid_error <- max(abs(old_grid$grid - new_grid$grid),
                      abs(old_grid$max_grid - new_grid$max_grid),
                      abs(old_grid$max_H01_time - new_grid$max_H01_time),
                      abs(old_grid$max_H02_time - new_grid$max_H02_time),
                      abs(old_grid$max_H03_time - new_grid$max_H03_time), na.rm = TRUE)
    base_terms_equal <- identical(attr(old, "base_terms", exact = TRUE), attr(new, "base_terms", exact = TRUE))
    old_coef <- stats::coef(attr(old, "fit", exact = TRUE)); new_coef <- stats::coef(attr(new, "fit", exact = TRUE))
    coef_error <- if (identical(names(old_coef), names(new_coef))) max(abs(old_coef - new_coef), na.rm = TRUE) else Inf
    max_error <- max(by_metric$max_abs_error, grid_error, coef_error, na.rm = TRUE)
    list(passed = isTRUE(base_terms_equal) && is.finite(max_error) && max_error <= tolerance,
         max_error = max_error, by_metric = by_metric, grid_max_error = grid_error,
         coef_max_error = coef_error, base_terms_equal = base_terms_equal, old = old, new = new)
}
