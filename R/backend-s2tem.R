# imediation_estimator.R

# ============================================================
# Helpers for chapter 2 (two sequential time mediators)
# Required raw columns:
#   A, t_out, d_out(or yevent), t_med1, d_med1, t_med2, d_med2
#
# Construct:
#   mevent1_eff, t_med1_eff
#   mevent2_eff, t_med2_eff
#
# Definitions:
#   M1 effective: M1 occurs before death/censoring
#   M2 effective: M2 occurs after M1 and before death/censoring
# ============================================================
.as01 <- function(x, nm) {
    if (is.factor(x)) x <- as.character(x)
    x <- as.numeric(x)
    bad <- !is.na(x) & !(x %in% c(0, 1))
    if (any(bad)) stop(sprintf("%s must be coded as 0/1.", nm))
    as.integer(x)
}

.im2_prepare_eff_vars <- function(dat) {
    dat <- as.data.frame(dat)
    
    if (!all(c("A", "t_out") %in% names(dat))) {
        stop("dat must contain A and t_out, plus d_out or yevent.")
    }
    if (!("yevent" %in% names(dat))) {
        if (!("d_out" %in% names(dat))) stop("Missing d_out or yevent.")
        dat$yevent <- as.integer(dat$d_out == 1L)
    }
    
    need_raw <- c("t_med1", "d_med1", "t_med2", "d_med2")
    miss_raw <- setdiff(need_raw, names(dat))
    if (length(miss_raw) > 0) {
        stop("Missing raw sequential-mediator variables: ", paste(miss_raw, collapse = ", "))
    }
    
    dat$A <- .as01(dat$A, "A")
    dat$d_med1 <- .as01(dat$d_med1, "d_med1")
    dat$d_med2 <- .as01(dat$d_med2, "d_med2")
    
    if ("d_out" %in% names(dat)) {
        dat$d_out <- .as01(dat$d_out, "d_out")
    }
    if ("yevent" %in% names(dat)) {
        dat$yevent <- .as01(dat$yevent, "yevent")
    }
    
    # ---- mediator 1 effective ----
    med1_happens <- with(dat,
                         d_med1 == 1L &
                             is.finite(t_med1) & is.finite(t_out) &
                             (t_med1 < t_out)
    )
    dat$mevent1_eff <- as.integer(med1_happens)
    dat$t_med1_eff  <- ifelse(dat$mevent1_eff == 1L, dat$t_med1, dat$t_out)
    
    # ---- mediator 2 effective ----
    med2_happens <- with(dat,
                         d_med2 == 1L &
                             dat$mevent1_eff == 1L &
                             is.finite(t_med2) & is.finite(t_out) &
                             (t_med2 > dat$t_med1_eff) &
                             (t_med2 < t_out)
    )
    dat$mevent2_eff <- as.integer(med2_happens)
    dat$t_med2_eff  <- ifelse(dat$mevent2_eff == 1L, dat$t_med2, dat$t_out)
    
    # handy indicators for naive methods
    dat$death_before_m1 <- with(dat, as.integer(yevent == 1L & mevent1_eff == 0L))
    dat$death_between_m1_m2 <- with(dat, as.integer(yevent == 1L & mevent1_eff == 1L & mevent2_eff == 0L))
    dat$completed_chain <- dat$mevent2_eff
    
    dat
}
# ============================================================
# Utilities: trapezoid cumulative integral for RMTL
# psi(t) = cumulative death risk by tau
# RMTL(t_k) = \int_0^{t_k} psi(u) du
# ============================================================

.compute_rmtl_from_psi <- function(tau, psi, method = c("trapezoid", "left")) {
    method <- match.arg(method)
    
    tau <- as.numeric(tau)
    psi <- as.numeric(psi)
    
    if (length(tau) != length(psi)) stop("tau and psi must have the same length.")
    if (any(!is.finite(tau))) stop("tau must contain finite values.")
    if (is.unsorted(tau, strictly = TRUE)) stop("tau must be strictly increasing.")
    
    n <- length(tau)
    out <- numeric(n)
    if (n == 0) return(out)
    
    # Start integration at 0 with psi(0)=0.
    t_prev <- 0
    p_prev <- 0
    
    for (k in seq_len(n)) {
        dt <- tau[k] - t_prev
        if (dt < 0) stop("tau must be increasing.")
        
        if (method == "left") {
            incr <- p_prev * dt
        } else {
            incr <- 0.5 * (p_prev + psi[k]) * dt
        }
        
        out[k] <- if (k == 1) incr else out[k - 1] + incr
        t_prev <- tau[k]
        p_prev <- psi[k]
    }
    
    out
}

# ============================================================
# Add RMTL-scale quantities to an estimator output
# Requires columns: tau, psi00, psi10, psi11
# Does NOT output INT_RMTL
# ============================================================

add_rmtl_to_output <- function(out, integration = c("trapezoid", "left")) {
    integration <- match.arg(integration)
    
    need <- c("tau", "psi00", "psi10", "psi11")
    miss <- setdiff(need, names(out))
    if (length(miss) > 0) {
        stop("Missing columns required to compute RMTL: ", paste(miss, collapse = ", "))
    }
    
    tau <- out$tau
    
    out$rmtl00 <- .compute_rmtl_from_psi(tau, out$psi00, method = integration)
    out$rmtl10 <- .compute_rmtl_from_psi(tau, out$psi10, method = integration)
    out$rmtl11 <- .compute_rmtl_from_psi(tau, out$psi11, method = integration)
    
    if ("psi01" %in% names(out) && !all(is.na(out$psi01))) {
        out$rmtl01 <- .compute_rmtl_from_psi(tau, out$psi01, method = integration)
    }
    
    out$DE_RMTL <- out$rmtl10 - out$rmtl00
    out$IE_RMTL <- out$rmtl11 - out$rmtl10
    out$TE_RMTL <- out$rmtl11 - out$rmtl00
    
    out
}

#####01#####

# ============================================================
# Build design matrices for 5 transitions
# transitions:
#   1: start -> M1
#   2: start -> death
#   3: M1 -> M2
#   4: M1 -> death
#   5: M2 -> death
# ============================================================
.build_all_X_2med <- function(dat, trans, covs_df, coef_cols, basec) {
    .align_to_coef_cols <- function(df, coef_cols) {
        miss <- setdiff(coef_cols, colnames(df))
        if (length(miss) > 0) {
            for (nm in miss) df[[nm]] <- 0
        }
        as.matrix(df[, coef_cols, drop = FALSE])
    }
    
    dat <- as.data.frame(dat)
    
    n <- nrow(dat)
    keep_vars <- unique(c("A", "t_med1_eff", "t_med2_eff", basec))
    
    proto_all <- dat[rep(seq_len(n), each = 5), keep_vars, drop = FALSE]
    proto_all$id <- rep(seq_len(n), each = 5)
    proto_all$trans  <- factor(rep(1:5, times = n), levels = 1:5)
    proto_all$strata <- factor(rep(1:5, times = n), levels = 1:5)
    
    for (nm in basec) {
        if (nm %in% names(proto_all) && is.character(proto_all[[nm]])) {
            proto_all[[nm]] <- factor(proto_all[[nm]], levels = unique(dat[[nm]]))
        }
    }
    
    attr(proto_all, "trans") <- trans
    class(proto_all) <- c("msdata", "data.frame")
    
    proto_exp <- mstate::expand.covs(
        proto_all,
        covs = c("A", "t_med1_eff", "t_med2_eff", basec),
        append = TRUE,
        longnames = FALSE
    )
    
    Xall <- .align_to_coef_cols(proto_exp, coef_cols)
    
    idx1 <- which(proto_all$trans == 1)
    idx2 <- which(proto_all$trans == 2)
    idx3 <- which(proto_all$trans == 3)
    idx4 <- which(proto_all$trans == 4)
    idx5 <- which(proto_all$trans == 5)
    
    list(
        X1 = Xall[idx1, , drop = FALSE],  # start -> M1
        X2 = Xall[idx2, , drop = FALSE],  # start -> death
        X3 = Xall[idx3, , drop = FALSE],  # M1 -> M2
        X4 = Xall[idx4, , drop = FALSE],  # M1 -> death
        X5 = Xall[idx5, , drop = FALSE]   # M2 -> death
    )
}

# ============================================================
# Two sequential mediators: marginal interventional estimator

# ============================================================

# ============================================================
# Shared backend for sequential two time-to-event mediators
#
# This is a refactor of the computational setup in
# interventional_cmaverse_marginal_2med().  It performs the data
# preparation, five-transition Cox fit, baseline cumulative hazards,
# grid construction/thinning, design matrices, and tau-independent
# caches ONCE.  It does not itself define GE or PE contrasts.
# ============================================================
.fit_s2tem_backend <- function(
        dat,
        tau = 1:10,
        basec = character(0),
        a = 1,
        astar = 0,
        bh_method = "breslow",
        EMint = FALSE,
        verbose = FALSE,
        v_chunk_size = 4096L,
        max_grid_points = 200L
) {
    dat <- .im2_prepare_eff_vars(dat)
    
    tau <- sort(unique(as.numeric(tau)))
    basec <- intersect(basec, names(dat))
    v_chunk_size <- as.integer(v_chunk_size)
    max_grid_points <- as.integer(max_grid_points)
    if (!length(tau) || any(!is.finite(tau)) || any(tau <= 0)) {
        stop("tau must contain positive finite evaluation times.")
    }
    if (v_chunk_size <= 0L) stop("v_chunk_size must be a positive integer.")
    if (max_grid_points < 10L) stop("max_grid_points must be at least 10.")
    
    if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required.")
    if (!requireNamespace("mstate", quietly = TRUE)) stop("Package 'mstate' is required.")
    
    safe_exp <- function(x) exp(pmin(pmax(x, -50), 50))
    
    H_vec_at <- function(Hdf, tvec) {
        if (nrow(Hdf) == 0) return(rep(0, length(tvec)))
        idx <- findInterval(tvec, Hdf$time)
        out <- numeric(length(tvec))
        ok <- idx > 0
        out[ok] <- Hdf$Haz[idx[ok]]
        out
    }
    
    .get_coef_any <- function(b, cand_names) {
        hit <- cand_names[cand_names %in% names(b)]
        if (length(hit) == 0) return(0)
        unname(b[[hit[1]]])
    }
    
    .escape_regex <- function(x) {
        gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
    }
    
    .pick_base_terms <- function(v, cols, all_basec = character(0)) {
        # expand.covs() uses v.1 for numerical covariates and may use
        # vLevel.1 (e.g., sexMale.1) for factors. Match both formats.
        pat <- paste0("^", .escape_regex(v), ".*\\.[12345]$")
        out <- grep(pat, cols, value = TRUE)
        
        # Avoid prefix collisions, e.g. U1 must not select U10.1,
        # and age must not select ageGroupYes.3 when ageGroup is
        # another baseline covariate.
        longer_vars <- all_basec[nchar(all_basec) > nchar(v)]
        if (length(longer_vars) > 0L && length(out) > 0L) {
            for (w in longer_vars) {
                pat_w <- paste0("^", .escape_regex(w), ".*\\.[12345]$")
                out <- out[!grepl(pat_w, out)]
            }
        }
        
        out
    }
    
    .trim_H <- function(Hdf, tau_max_use) {
        Hdf <- Hdf[, c("time", "Haz"), drop = FALSE]
        Hdf <- Hdf[is.finite(Hdf$time), , drop = FALSE]
        Hdf <- Hdf[Hdf$time <= tau_max_use, , drop = FALSE]
        Hdf <- Hdf[!duplicated(Hdf$time), , drop = FALSE]
        Hdf <- Hdf[order(Hdf$time), , drop = FALSE]
        if (nrow(Hdf) == 0L || Hdf$time[1] > 0) {
            Hdf <- rbind(data.frame(time = 0, Haz = 0), Hdf)
        }
        Hdf
    }
    
    .thin_H_adaptive <- function(Hdf, max_pts = 200) {
        n <- nrow(Hdf)
        if (n <= max_pts) return(Hdf)
        
        if ("Haz" %in% names(Hdf)) {
            haz_col <- "Haz"
        } else if ("hazard" %in% names(Hdf)) {
            haz_col <- "hazard"
        } else {
            stop(".thin_H_adaptive: could not find a 'Haz' or 'hazard' column.")
        }
        
        Haz_vals <- Hdf[[haz_col]]
        Haz_min  <- Haz_vals[1]
        Haz_max  <- Haz_vals[n]
        
        if (Haz_max <= Haz_min) {
            idx <- unique(c(1L, round(seq(1, n, length.out = max_pts)), n))
            return(Hdf[sort(unique(idx)), , drop = FALSE])
        }
        
        target_Haz <- seq(Haz_min, Haz_max, length.out = max_pts)
        idx <- vapply(target_Haz, function(h) which.min(abs(Haz_vals - h)), integer(1))
        idx <- sort(unique(c(1L, idx, n)))
        Hdf[idx, , drop = FALSE]
    }
    
    # ----- 1) five-transition multistate Cox model -----
    trans <- mstate::transMat(
        x = list(c(2, 4), c(3, 4), c(4), c()),
        names = c("start", "M1", "M2", "death")
    )
    
    covs_df <- c("A", "t_med1_eff", "t_med2_eff", basec)
    
    mstate_dat <- mstate::msprep(
        time   = c(NA, "t_med1_eff", "t_med2_eff", "t_out"),
        status = c(NA, "mevent1_eff", "mevent2_eff", "yevent"),
        data   = dat,
        trans  = trans,
        keep   = covs_df
    )
    
    mstate_dat <- mstate::expand.covs(
        mstate_dat, covs_df,
        append = TRUE, longnames = FALSE
    )
    
    base_terms <- character(0)
    if (length(basec) > 0) {
        base_terms <- unique(unlist(lapply(
            basec, .pick_base_terms,
            cols = colnames(mstate_dat), all_basec = basec
        )))
    }
    
    rhs_terms <- c(
        paste0("A.", 1:5),
        "t_med1_eff.3", "t_med1_eff.4", "t_med2_eff.5",
        base_terms
    )
    if (EMint) {
        rhs_terms <- c(rhs_terms,
                       "A.3:t_med1_eff.3",
                       "A.4:t_med1_eff.4",
                       "A.5:t_med2_eff.5")
    }
    
    fml <- stats::as.formula(
        paste0("survival::Surv(Tstart, Tstop, status) ~ ",
               paste(rhs_terms, collapse = " + "),
               " + strata(trans)")
    )
    
    fit <- survival::coxph(
        fml, data = mstate_dat,
        method = bh_method, x = TRUE, model = TRUE
    )
    
    b <- stats::coef(fit)
    if (length(b) == 0) stop("coxph estimated no coefficients.")
    if (anyNA(b)) {
        if (verbose) warning(sprintf("coxph produced %d NA coefficients; they were set to 0.", sum(is.na(b))))
        b[is.na(b)] <- 0
    }
    
    # ----- 2) baseline cumulative hazards + thinning -----
    mf <- stats::model.frame(fit)
    need_cols <- setdiff(colnames(mf), c("(weights)", "(offset)"))
    need_cols <- setdiff(need_cols, "Surv(Tstart, Tstop, status)")
    need_cols <- intersect(need_cols, colnames(mstate_dat))
    
    newd0 <- mstate_dat[0, need_cols, drop = FALSE]
    newd0 <- newd0[rep(1, 5), , drop = FALSE]
    newd0[,] <- 0
    newd0$trans  <- 1:5
    newd0$strata <- 1:5
    attr(newd0, "trans") <- trans
    class(newd0) <- c("msdata", "data.frame")
    
    Haz0 <- mstate::msfit(fit, newdata = newd0, trans = trans)$Haz
    tau_max_use <- max(tau)
    
    H12_raw <- .trim_H(subset(Haz0, trans == 1), tau_max_use)
    H14_raw <- .trim_H(subset(Haz0, trans == 2), tau_max_use)
    H23_raw <- .trim_H(subset(Haz0, trans == 3), tau_max_use)
    H24_raw <- .trim_H(subset(Haz0, trans == 4), tau_max_use)
    H34_raw <- .trim_H(subset(Haz0, trans == 5), tau_max_use)
    
    n_outer_orig <- nrow(H12_raw)
    n_inner_orig <- nrow(H23_raw)
    
    H12 <- .thin_H_adaptive(H12_raw, max_grid_points)
    H14 <- .thin_H_adaptive(H14_raw, max_grid_points)
    H23 <- .thin_H_adaptive(H23_raw, max_grid_points)
    H24 <- .thin_H_adaptive(H24_raw, max_grid_points)
    H34 <- .thin_H_adaptive(H34_raw, max_grid_points)
    
    if (verbose) {
        cat(sprintf(
            "[grid thin] H12:%d->%d H14:%d->%d H23:%d->%d H24:%d->%d H34:%d->%d (max_grid_points=%d)\\n",
            n_outer_orig, nrow(H12),
            nrow(H14_raw), nrow(H14),
            n_inner_orig, nrow(H23),
            nrow(H24_raw), nrow(H24),
            nrow(H34_raw), nrow(H34),
            max_grid_points
        ))
    }
    
    if (nrow(H12) < 2) stop("H12 time grid is too short (<2).")
    if (nrow(H23) < 2) stop("H23 time grid is too short (<2).")
    
    time1    <- H12$time
    u1_left  <- time1[-length(time1)]
    u1_right <- time1[-1]
    du1_base <- diff(time1)
    if (any(du1_base <= 0)) stop("H12 time grid is not strictly increasing.")
    
    H12_u1 <- H12$Haz[-length(H12$Haz)]
    H14_u1 <- H_vec_at(H14, u1_left)
    H23_u1 <- H_vec_at(H23, u1_left)
    H24_u1 <- H_vec_at(H24, u1_left)
    
    h12_base <- diff(H12$Haz) / du1_base
    w1_full  <- diff(H12$Haz)
    
    time2    <- H23$time
    u2_left  <- time2[-length(time2)]
    u2_right <- time2[-1]
    du2_base <- diff(time2)
    if (any(du2_base <= 0)) stop("H23 time grid is not strictly increasing.")
    
    H23_u2 <- H23$Haz[-length(H23$Haz)]
    H24_u2 <- H_vec_at(H24, u2_left)
    H34_u2 <- H_vec_at(H34, u2_left)
    
    h23_base <- diff(H23$Haz) / du2_base
    w2_full  <- diff(H23$Haz)
    
    # ----- 3) design matrices and intervention-specific multipliers -----
    Xlist <- .build_all_X_2med(
        dat = dat, trans = trans,
        covs_df = covs_df, coef_cols = names(b), basec = basec
    )
    X1 <- Xlist$X1; X2 <- Xlist$X2; X3 <- Xlist$X3
    X4 <- Xlist$X4; X5 <- Xlist$X5
    
    colA1 <- "A.1"; colA2 <- "A.2"; colA3 <- "A.3"
    colA4 <- "A.4"; colA5 <- "A.5"
    colM1_3 <- "t_med1_eff.3"; colM1_4 <- "t_med1_eff.4"
    colM2_5 <- "t_med2_eff.5"
    
    set_A_5trans <- function(X, which_col, Aval) {
        for (nm in c(colA1, colA2, colA3, colA4, colA5)) {
            if (nm %in% colnames(X)) X[, nm] <- 0
        }
        if (which_col %in% colnames(X)) X[, which_col] <- Aval
        X
    }
    
    X3_const <- X3
    for (nm in c(colM1_3, "A.3:t_med1_eff.3", "t_med1_eff.3:A.3")) {
        if (nm %in% colnames(X3_const)) X3_const[, nm] <- 0
    }
    X4_const <- X4
    for (nm in c(colM1_4, "A.4:t_med1_eff.4", "t_med1_eff.4:A.4")) {
        if (nm %in% colnames(X4_const)) X4_const[, nm] <- 0
    }
    X5_const <- X5
    for (nm in c(colM2_5, "A.5:t_med2_eff.5", "t_med2_eff.5:A.5")) {
        if (nm %in% colnames(X5_const)) X5_const[, nm] <- 0
    }
    
    beta_m1_23 <- if (colM1_3 %in% names(b)) unname(b[[colM1_3]]) else 0
    beta_m1_24 <- if (colM1_4 %in% names(b)) unname(b[[colM1_4]]) else 0
    beta_m2_34 <- if (colM2_5 %in% names(b)) unname(b[[colM2_5]]) else 0
    
    beta_int_23 <- if (EMint) .get_coef_any(b, c("A.3:t_med1_eff.3", "t_med1_eff.3:A.3")) else 0
    beta_int_24 <- if (EMint) .get_coef_any(b, c("A.4:t_med1_eff.4", "t_med1_eff.4:A.4")) else 0
    beta_int_34 <- if (EMint) .get_coef_any(b, c("A.5:t_med2_eff.5", "t_med2_eff.5:A.5")) else 0
    
    eta12_astar <- as.numeric(set_A_5trans(X1, colA1, astar) %*% b)
    eta12_a     <- as.numeric(set_A_5trans(X1, colA1, a)     %*% b)
    eta14_astar <- as.numeric(set_A_5trans(X2, colA2, astar) %*% b)
    eta14_a     <- as.numeric(set_A_5trans(X2, colA2, a)     %*% b)
    eta23c_astar <- as.numeric(set_A_5trans(X3_const, colA3, astar) %*% b)
    eta23c_a     <- as.numeric(set_A_5trans(X3_const, colA3, a)     %*% b)
    eta24c_astar <- as.numeric(set_A_5trans(X4_const, colA4, astar) %*% b)
    eta24c_a     <- as.numeric(set_A_5trans(X4_const, colA4, a)     %*% b)
    eta34c_astar <- as.numeric(set_A_5trans(X5_const, colA5, astar) %*% b)
    eta34c_a     <- as.numeric(set_A_5trans(X5_const, colA5, a)     %*% b)
    
    e12_astar <- safe_exp(eta12_astar); e12_a <- safe_exp(eta12_a)
    e14_astar <- safe_exp(eta14_astar); e14_a <- safe_exp(eta14_a)
    e23c_astar <- safe_exp(eta23c_astar); e23c_a <- safe_exp(eta23c_a)
    e24c_astar <- safe_exp(eta24c_astar); e24c_a <- safe_exp(eta24c_a)
    e34c_astar <- safe_exp(eta34c_astar); e34c_a <- safe_exp(eta34c_a)
    
    N <- nrow(dat)
    
    E12  <- cbind(e12_astar, e12_astar, e12_a)
    E14  <- cbind(e14_astar, e14_a,     e14_a)
    E23c <- cbind(e23c_astar, e23c_astar, e23c_a)
    E24c <- cbind(e24c_astar, e24c_a,     e24c_a)
    E34c <- cbind(e34c_astar, e34c_a,     e34c_a)
    
    alpha23_vec <- c(beta_m1_23 + beta_int_23 * astar,
                     beta_m1_23 + beta_int_23 * astar,
                     beta_m1_23 + beta_int_23 * a)
    alpha24_vec <- c(beta_m1_24 + beta_int_24 * astar,
                     beta_m1_24 + beta_int_24 * a,
                     beta_m1_24 + beta_int_24 * a)
    alpha34_vec <- c(beta_m2_34 + beta_int_34 * astar,
                     beta_m2_34 + beta_int_34 * a,
                     beta_m2_34 + beta_int_34 * a)
    
    m1 <- length(u1_left)
    
    # ----- 4) tau-independent and tau-specific caches used by GE -----
    k23_u1_all <- safe_exp(outer(u1_left, alpha23_vec))
    k24_u1_all <- safe_exp(outer(u1_left, alpha24_vec))
    k34_u2_all <- safe_exp(outer(u2_left, alpha34_vec))
    
    H12_tau <- H_vec_at(H12, tau)
    H14_tau <- H_vec_at(H14, tau)
    H23_tau <- H_vec_at(H23, tau)
    H24_tau <- H_vec_at(H24, tau)
    H34_tau <- H_vec_at(H34, tau)
    
    dH23_tau_u1 <- outer(H23_tau, H23_u1, "-")
    dH24_tau_u1 <- outer(H24_tau, H24_u1, "-")
    dH23_tau_u1[dH23_tau_u1 < 0] <- 0
    dH24_tau_u1[dH24_tau_u1 < 0] <- 0
    
    q34_tau_u2 <- outer(H34_tau, H34_u2, "-")
    q34_tau_u2[q34_tau_u2 < 0] <- 0
    
    idx2_after_u1 <- vector("list", m1)
    q23_from_u1   <- vector("list", m1)
    q24_from_u1   <- vector("list", m1)
    
    for (ii in seq_len(m1)) {
        idx2 <- which(u2_left > u1_left[ii])
        idx2_after_u1[[ii]] <- idx2
        if (length(idx2) > 0L) {
            q23_from_u1[[ii]] <- pmax(H23_u2[idx2] - H23_u1[ii], 0)
            q24_from_u1[[ii]] <- pmax(H24_u2[idx2] - H24_u1[ii], 0)
        } else {
            q23_from_u1[[ii]] <- numeric(0)
            q24_from_u1[[ii]] <- numeric(0)
        }
    }
    
    build_w1_tau <- function(tt) {
        imax <- sum(u1_left < tt)
        if (imax <= 0L) return(numeric(0))
        w <- w1_full[seq_len(imax)]
        if (tt < u1_right[imax]) {
            w[imax] <- (tt - u1_left[imax]) * h12_base[imax]
        }
        w
    }
    build_w2_tau <- function(tt) {
        jmax <- sum(u2_left < tt)
        if (jmax <= 0L) return(numeric(0))
        w <- w2_full[seq_len(jmax)]
        if (tt < u2_right[jmax]) {
            w[jmax] <- (tt - u2_left[jmax]) * h23_base[jmax]
        }
        w
    }
    
    w1_tau_list <- lapply(tau, build_w1_tau)
    w2_tau_list <- lapply(tau, build_w2_tau)
    
    precomp_S0   <- vector("list", m1)
    precomp_e23  <- vector("list", m1)
    precomp_e24  <- vector("list", m1)
    
    for (ii in seq_len(m1)) {
        precomp_S0[[ii]] <- exp(-(H12_u1[ii] * E12 + H14_u1[ii] * E14))
        precomp_e23[[ii]] <- sweep(E23c, 2L, k23_u1_all[ii, ], "*")
        precomp_e24[[ii]] <- sweep(E24c, 2L, k24_u1_all[ii, ], "*")
    }
    
    settings <- list(
        v_chunk_size = v_chunk_size,
        max_grid_points = max_grid_points,
        n_outer_grid = length(time1),
        n_inner_grid = length(time2),
        n_outer_grid_original = n_outer_orig,
        n_inner_grid_original = n_inner_orig
    )
    
    context <- list(
        tau = tau,
        N = N,
        a = a,
        astar = astar,
        hazards = list(H12 = H12, H14 = H14, H23 = H23, H24 = H24, H34 = H34),
        mult = list(
            e12_astar = e12_astar, e12_a = e12_a,
            e14_astar = e14_astar, e14_a = e14_a,
            e23c_astar = e23c_astar, e23c_a = e23c_a,
            e24c_astar = e24c_astar, e24c_a = e24c_a,
            e34c_astar = e34c_astar, e34c_a = e34c_a
        ),
        beta = list(
            beta_m1_23 = beta_m1_23,
            beta_m1_24 = beta_m1_24,
            beta_m2_34 = beta_m2_34,
            beta_int_23 = beta_int_23,
            beta_int_24 = beta_int_24,
            beta_int_34 = beta_int_34
        ),
        grid = list(
            u1_left = u1_left,
            u2_left = u2_left,
            H12_u1 = H12_u1,
            H23_u1 = H23_u1,
            H23_u2 = H23_u2,
            w1_tau_list = w1_tau_list,
            w2_tau_list = w2_tau_list
        )
    )
    
    structure(list(
        dat = dat,
        tau = tau,
        basec = basec,
        a = a,
        astar = astar,
        bh_method = bh_method,
        EMint = EMint,
        fit = fit,
        mstate_dat = mstate_dat,
        trans = trans,
        settings = settings,
        context = context,
        ge_cache = list(
            N = N,
            E12 = E12, E14 = E14, E23c = E23c, E24c = E24c, E34c = E34c,
            k23_u1_all = k23_u1_all,
            k24_u1_all = k24_u1_all,
            k34_u2_all = k34_u2_all,
            H12_tau = H12_tau, H14_tau = H14_tau,
            dH23_tau_u1 = dH23_tau_u1,
            dH24_tau_u1 = dH24_tau_u1,
            q34_tau_u2 = q34_tau_u2,
            idx2_after_u1 = idx2_after_u1,
            q23_from_u1 = q23_from_u1,
            q24_from_u1 = q24_from_u1,
            w1_tau_list = w1_tau_list,
            w2_tau_list = w2_tau_list,
            precomp_S0 = precomp_S0,
            precomp_e23 = precomp_e23,
            precomp_e24 = precomp_e24,
            v_chunk_size = v_chunk_size
        )
    ), class = "mims_s2tem_backend")
}

# ============================================================
# GE prediction from the shared S2TEM backend
# ============================================================
.predict_s2tem_ge <- function(
        backend,
        scale = c("risk", "survival"),
        keep_fit_objects = FALSE
) {
    scale <- match.arg(scale)
    if (!inherits(backend, "mims_s2tem_backend")) {
        stop("backend must be created by .fit_s2tem_backend().")
    }
    
    tau <- backend$tau
    z <- backend$ge_cache
    N <- z$N
    
    calc_all_psi_tau <- function(k) {
        w1_tau <- z$w1_tau_list[[k]]
        if (!length(w1_tau)) return(c(0, 0, 0))
        
        jmax_global <- length(z$w2_tau_list[[k]])
        
        P0 <- exp(-(z$H12_tau[k] * z$E12 + z$H14_tau[k] * z$E14))
        P1 <- matrix(0, nrow = N, ncol = 3)
        P2 <- matrix(0, nrow = N, ncol = 3)
        
        n_outer <- length(w1_tau)
        
        for (ii in seq_len(n_outer)) {
            w1 <- w1_tau[ii]
            if (w1 <= 0) next
            
            S0_u1  <- z$precomp_S0[[ii]]
            e23_u1 <- z$precomp_e23[[ii]]
            e24_u1 <- z$precomp_e24[[ii]]
            haz12_w1 <- z$E12 * w1
            
            S1_ut <- exp(-(
                z$dH23_tau_u1[k, ii] * e23_u1 +
                    z$dH24_tau_u1[k, ii] * e24_u1
            ))
            P1 <- P1 + S0_u1 * haz12_w1 * S1_ut
            
            idx2_full <- z$idx2_after_u1[[ii]]
            if (!length(idx2_full) || jmax_global <= 0L) next
            
            keep_n <- sum(idx2_full <= jmax_global)
            if (keep_n <= 0L) next
            
            idx2    <- idx2_full[seq_len(keep_n)]
            q23_all <- z$q23_from_u1[[ii]][seq_len(keep_n)]
            q24_all <- z$q24_from_u1[[ii]][seq_len(keep_n)]
            q34_all <- z$q34_tau_u2[k, idx2]
            wc_all  <- z$w2_tau_list[[k]][idx2]
            
            J <- length(idx2)
            inner_sum <- matrix(0, nrow = N, ncol = 3)
            mem_cells <- as.numeric(N) * J
            
            if (mem_cells <= 5e7) {
                for (r in 1:3) {
                    term23 <- outer(e23_u1[, r], q23_all)
                    term24 <- outer(e24_u1[, r], q24_all)
                    term34 <- outer(z$E34c[, r], q34_all * z$k34_u2_all[idx2, r])
                    inner_sum[, r] <- exp(-(term23 + term24 + term34)) %*% wc_all
                }
            } else {
                chunk_starts <- seq.int(1L, J, by = z$v_chunk_size)
                for (st in chunk_starts) {
                    ed  <- min(st + z$v_chunk_size - 1L, J)
                    loc <- st:ed
                    q23c <- q23_all[loc]
                    q24c <- q24_all[loc]
                    q34c <- q34_all[loc]
                    wc   <- wc_all[loc]
                    
                    for (r in 1:3) {
                        k34c_r <- z$k34_u2_all[idx2[loc], r]
                        term23 <- outer(e23_u1[, r], q23c)
                        term24 <- outer(e24_u1[, r], q24c)
                        term34 <- outer(z$E34c[, r], q34c * k34c_r)
                        inner_sum[, r] <- inner_sum[, r] +
                            as.numeric(exp(-(term23 + term24 + term34)) %*% wc)
                    }
                }
            }
            
            P2 <- P2 + S0_u1 * haz12_w1 * e23_u1 * inner_sum
        }
        
        S <- P0 + P1 + P2
        S <- pmin(pmax(S, 0), 1)
        colMeans(1 - S)
    }
    
    psi_mat <- t(vapply(seq_along(tau), calc_all_psi_tau, numeric(3)))
    colnames(psi_mat) <- c("psi00", "psi10", "psi11")
    
    out <- data.frame(
        tau = tau,
        psi00 = psi_mat[, "psi00"],
        psi10 = psi_mat[, "psi10"],
        psi11 = psi_mat[, "psi11"]
    )
    out$S00 <- 1 - out$psi00
    out$S10 <- 1 - out$psi10
    out$S11 <- 1 - out$psi11
    
    out$DE <- out$psi10 - out$psi00
    out$IE <- out$psi11 - out$psi10
    out$TE <- out$psi11 - out$psi00
    
    if (scale == "survival") {
        out$DE <- -out$DE
        out$IE <- -out$IE
        out$TE <- -out$TE
    }
    
    # Restore the legacy column order before adding RMTL columns.
    out <- out[, c("tau", "DE", "IE", "TE",
                   "psi00", "psi10", "psi11",
                   "S00", "S10", "S11")]
    out <- add_rmtl_to_output(out, integration = "trapezoid")
    
    if (keep_fit_objects) {
        attr(out, "fit") <- backend$fit
        attr(out, "mstate_dat") <- backend$mstate_dat
    }
    attr(out, "exactfast_settings") <- backend$settings
    attr(out, "backend_class") <- class(backend)
    out
}

# Public-facing point-estimator candidate for chapter 2 / GE.
# For the first refactor stage this intentionally mirrors the legacy argument names.
mims_s2tem_ge_point <- function(
        dat,
        tau = 1:10,
        basec = character(0),
        a = 1,
        astar = 0,
        bh_method = "breslow",
        EMint = FALSE,
        scale = c("risk", "survival"),
        verbose = FALSE,
        keep_fit_objects = FALSE,
        v_chunk_size = 4096L,
        max_grid_points = 200L
) {
    scale <- match.arg(scale)
    backend <- .fit_s2tem_backend(
        dat = dat, tau = tau, basec = basec, a = a, astar = astar,
        bh_method = bh_method, EMint = EMint, verbose = verbose,
        v_chunk_size = v_chunk_size, max_grid_points = max_grid_points
    )
    .predict_s2tem_ge(
        backend = backend,
        scale = scale,
        keep_fit_objects = keep_fit_objects
    )
}

# Numerical regression test against the current chapter-2 function.
verify_s2tem_ge_equivalence <- function(
        dat,
        legacy_fun,
        tau = 1:10,
        basec = character(0),
        a = 1,
        astar = 0,
        bh_method = "breslow",
        EMint = FALSE,
        scale = "risk",
        v_chunk_size = 4096L,
        max_grid_points = 200L,
        tolerance = 1e-12
) {
    old <- legacy_fun(
        dat = dat, tau = tau, basec = basec, a = a, astar = astar,
        bh_method = bh_method, EMint = EMint, scale = scale,
        verbose = FALSE, keep_fit_objects = FALSE,
        v_chunk_size = v_chunk_size, max_grid_points = max_grid_points
    )
    new <- mims_s2tem_ge_point(
        dat = dat, tau = tau, basec = basec, a = a, astar = astar,
        bh_method = bh_method, EMint = EMint, scale = scale,
        verbose = FALSE, keep_fit_objects = FALSE,
        v_chunk_size = v_chunk_size, max_grid_points = max_grid_points
    )
    
    cols <- intersect(
        c("psi00", "psi10", "psi11", "DE", "IE", "TE",
          "rmtl00", "rmtl10", "rmtl11", "DE_RMTL", "IE_RMTL", "TE_RMTL"),
        intersect(names(old), names(new))
    )
    diffs <- vapply(cols, function(nm) max(abs(old[[nm]] - new[[nm]]), na.rm = TRUE), numeric(1))
    
    list(
        passed = all(diffs <= tolerance),
        tolerance = tolerance,
        max_error = if (length(diffs)) max(diffs) else NA_real_,
        by_metric = data.frame(metric = names(diffs), max_abs_error = unname(diffs)),
        old = old,
        new = new
    )
}
