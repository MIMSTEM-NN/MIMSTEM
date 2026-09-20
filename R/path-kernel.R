# ============================================================
# MIMS-S2TEM-PE refactor v2
# Shared S2TEM backend + independent PE path kernel.
#
# Required first:
#   source("s2tem_backend_refactor_v1.R", encoding = "UTF-8")
#
# This file no longer copies the chapter-2 estimator and no longer
# captures/modifies another function body.  It uses the verified
# .fit_s2tem_backend() directly and keeps the frozen PE kernel logic.
# ============================================================

.pe2_version <- "mims_s2tem_pe_shared_backend_v2_2026-09-17"
.pe2_risks <- c("psi0000", "psi1000", "psi1100", "psi1110", "psi1111")
.pe2_effects <- c(
    "P_A_Y", "P_A_T2_Y", "P_A_T1_Y", "P_A_T1_T2_Y",
    "DE", "IE_total", "TE_path", "IE", "TE"
)

.pe2_integer <- function(x, nm, minimum = 1L) {
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
        x != floor(x) || x < minimum || x > .Machine$integer.max) {
        stop(nm, " must be an integer >= ", minimum, ".")
    }
    as.integer(x)
}

.pe2_validate <- function(dat, tau, basec, a, astar) {
    if (!is.data.frame(dat) || nrow(dat) < 2L) {
        stop("dat must be a subject-level data.frame.")
    }
    if (!is.numeric(tau) || !length(tau) || any(!is.finite(tau)) || any(tau <= 0)) {
        stop("tau must contain positive finite times.")
    }
    if (length(a) != 1L || length(astar) != 1L || anyNA(c(a, astar)) ||
        any(!c(a, astar) %in% 0:1)) {
        stop("a and astar must be 0 or 1.")
    }
    if (is.null(basec)) basec <- character(0)
    if (!is.character(basec) || anyNA(basec) || anyDuplicated(basec)) {
        stop("Invalid basec.")
    }
    reserved <- c(
        "A", "id", "t_out", "d_out", "yevent",
        "t_med1", "t_med2", "d_med1", "d_med2",
        "t_med1_eff", "t_med2_eff", "mevent1_eff", "mevent2_eff"
    )
    if (any(basec %in% reserved)) {
        stop("basec must contain baseline covariates, not exposure or outcome variables.")
    }
    cols <- c("A", "t_out", "t_med1", "d_med1", "t_med2", "d_med2", basec)
    if (!any(c("d_out", "yevent") %in% names(dat))) {
        stop("Missing d_out or yevent.")
    }
    if (!all(cols %in% names(dat))) {
        stop("Missing columns: ", paste(setdiff(cols, names(dat)), collapse = ", "))
    }
    if (anyNA(dat[cols])) {
        stop("Missing data must be handled explicitly before fitting.")
    }
    if (any(vapply(dat[cols], function(x) is.numeric(x) && any(!is.finite(x)), logical(1)))) {
        stop("Non-finite observed input values are not supported.")
    }
    if (any(dat$t_out <= 0) || any(dat$t_med1 < 0) || any(dat$t_med2 < 0)) {
        stop("Invalid times.")
    }
    Avec <- as.numeric(as.character(dat$A))
    if (!all(c(0, 1) %in% Avec)) {
        stop("Both exposure groups are required.")
    }
    if ("id" %in% names(dat) && (anyNA(dat$id) || anyDuplicated(dat$id))) {
        stop("Supply one row per subject with unique id; bootstrap resets duplicated ids internally.")
    }
    for (nm in intersect(c("d_out", "yevent"), names(dat))) {
        if (anyNA(dat[[nm]]) || any(!dat[[nm]] %in% 0:1)) {
            stop("Invalid death indicator: ", nm)
        }
    }
    if (all(c("d_out", "yevent") %in% names(dat)) && any(dat$d_out != dat$yevent)) {
        stop("d_out and yevent disagree; remove or correct stale derived columns.")
    }
    invisible(TRUE)
}

# ------------------------------------------------------------
# Baseline cumulative-hazard utilities for the PE kernel.
# ------------------------------------------------------------
.pe2_hazard <- function(df) {
    if (!is.data.frame(df) || !all(c("time", "Haz") %in% names(df))) {
        stop("Each baseline hazard must contain time and Haz.")
    }
    if (any(!is.finite(df$time)) || any(!is.finite(df$Haz)) ||
        any(df$time < 0) || any(df$Haz < 0)) {
        stop("Invalid baseline hazard values.")
    }
    df <- df[order(df$time), , drop = FALSE]
    df <- df[!duplicated(df$time, fromLast = TRUE), , drop = FALSE]
    if (nrow(df) && any(df$time == 0 & df$Haz != 0)) {
        stop("The current continuous-time implementation does not support a hazard jump at time zero.")
    }
    df <- df[df$time > 0, , drop = FALSE]
    tt <- c(0, df$time)
    hh <- c(0, df$Haz)
    if (any(diff(hh) < 0)) stop("Cumulative baseline hazard is not monotone.")
    if (length(tt) == 1L) {
        tt <- c(0, 1)
        hh <- c(0, 0)
    }
    H <- stats::approxfun(tt, hh, method = "linear", rule = 2, ties = "ordered")
    inv <- function(target) {
        if (anyNA(target) || any(target < 0)) stop("Invalid hazard-inversion threshold.")
        ans <- rep(Inf, length(target))
        ok <- is.finite(target) & target <= utils::tail(hh, 1)
        if (!any(ok)) return(ans)
        y <- target[ok]
        exact <- match(y, hh, nomatch = 0L)
        z <- numeric(length(y))
        eq <- exact > 0L
        z[eq] <- tt[exact[eq]]
        if (any(!eq)) {
            j <- findInterval(y[!eq], hh)
            z[!eq] <- tt[j] + (y[!eq] - hh[j]) *
                (tt[j + 1L] - tt[j]) / (hh[j + 1L] - hh[j])
        }
        ans[ok] <- z
        ans
    }
    list(H = H, inv = inv)
}

.pe2_step_H <- function(df, t) {
    idx <- findInterval(t, df$time)
    x <- numeric(length(t))
    ok <- idx > 0
    x[ok] <- df$Haz[idx[ok]]
    x
}

# ------------------------------------------------------------
# Frozen structural PE kernel.
# switches are ordered (zY, z2, z1, z12).
# ------------------------------------------------------------
.pe2_path_kernel <- function(ctx, switches, rule = "t2_overrides") {
    if (length(switches) != 4L || anyNA(switches) || any(!switches %in% 0:1)) {
        stop("switches must contain four 0/1 values in the order (zY,z2,z1,z12).")
    }
    rule <- match.arg(rule, c("t2_overrides", "require_both"))
    zy <- switches[1]
    z2 <- switches[2]
    z1 <- switches[3]
    zm <- switches[4]
    N <- ctx$N
    gr <- ctx$grid

    if (N < 1 || any(!is.finite(unlist(ctx$mult))) || any(unlist(ctx$mult) <= 0)) {
        stop("Invalid transition-risk multipliers.")
    }

    ex <- function(x) exp(pmin(pmax(x, -50), 50))
    m <- function(st, z) ctx$mult[[paste0(st, if (z == 0) "_astar" else "_a")]]
    actual <- c(ctx$astar, ctx$a)
    al <- function(st, z) {
        stem <- if (st == 34) "beta_m2_" else "beta_m1_"
        ctx$beta[[paste0(stem, st)]] +
            ctx$beta[[paste0("beta_int_", st)]] * actual[z + 1L]
    }
    hm <- function(st, t) .pe2_step_H(ctx$hazards[[paste0("H", st)]], t)
    inverse12 <- .pe2_hazard(ctx$hazards$H12)$inv

    em <- m("e12", zm)
    ey <- m("e12", z1)
    d14 <- m("e14", zy)
    d24 <- m("e24c", zy)
    d34 <- m("e34c", zy)
    a23 <- al(23, z2)
    a24 <- al(24, zy)
    a34 <- al(34, zy)

    death_survival <- function(t, uy, v) {
        L <- max(length(uy), length(v))
        uy <- rep_len(uy, L)
        v <- rep_len(v, L)
        w <- if (rule == "t2_overrides") v else pmax(uy, v)

        D <- rep_len(d14, L) * hm(14, pmin(t, uy, w))
        end24 <- pmin(t, w)
        ok <- is.finite(uy) & uy < end24
        if (any(ok)) {
            D[ok] <- D[ok] + rep_len(d24, L)[ok] * ex(a24 * uy[ok]) *
                pmax(hm(24, end24[ok]) - hm(24, uy[ok]), 0)
        }
        ok <- is.finite(w) & w < t
        if (any(ok)) {
            D[ok] <- D[ok] + rep_len(d34, L)[ok] * ex(a34 * v[ok]) *
                pmax(hm(34, t) - hm(34, w[ok]), 0)
        }
        if (anyNA(D) || any(D < 0)) stop("Invalid cumulative death hazard.")
        exp(-D)
    }

    K <- length(ctx$tau)
    risk <- numeric(K)
    diagnostics <- vector("list", K)

    for (k in seq_len(K)) {
        t <- ctx$tau[k]
        w1 <- gr$w1_tau_list[[k]]
        w2 <- gr$w2_tau_list[[k]]
        if (!length(w1)) stop("No outer integration interval at an evaluation time.")

        qm <- em * hm(12, t)
        qy <- ey * hm(12, t)

        # E1 exceeds both thresholds: neither T1 copy occurs before t.
        S <- exp(-pmax(qm, qy)) * exp(-d14 * hm(14, t))

        # Region uM > t but uY < t, needed when the Y-copy is earlier.
        early <- ey > em
        tail_mass <- numeric(N)
        if (any(early)) {
            h12 <- ctx$hazards$H12
            for (j in seq_len(nrow(h12) - 1L)) {
                left <- pmax(qm, ey * h12$Haz[j])
                right <- pmin(qy, ey * h12$Haz[j + 1L])
                ok <- early & right > left
                if (!any(ok)) next

                uy <- rep(Inf, N)
                uy[ok] <- inverse12(left[ok] / ey[ok])
                exact_left <- ok & left == ey * h12$Haz[j]
                uy[exact_left] <- h12$time[j]

                wt <- numeric(N)
                wt[ok] <- exp(-left[ok]) * (right[ok] - left[ok])
                S <- S + wt * death_survival(t, uy, Inf)
                tail_mass <- tail_mass + wt
            }
        }

        mixed_mass <- numeric(N)
        for (ii in seq_along(w1)) {
            if (w1[ii] <= 0) next

            um <- gr$u1_left[ii]
            uy <- if (z1 == zm) rep(um, N) else inverse12(em * gr$H12_u1[ii] / ey)
            uy[em == ey] <- um

            wt1 <- exp(-em * gr$H12_u1[ii]) * em * w1[ii]
            e23 <- m("e23c", z2) * ex(a23 * um)
            Gt <- e23 * pmax(hm(23, t) - gr$H23_u1[ii], 0)

            S <- S + wt1 * exp(-Gt) * death_survival(t, uy, Inf)

            idx <- which(gr$u2_left > um & seq_along(gr$u2_left) <= length(w2))
            if (!length(idx)) next
            idx <- idx[w2[idx] > 0]
            if (!length(idx)) next

            vv <- gr$u2_left[idx]
            J <- length(idx)
            Gv <- outer(e23, pmax(gr$H23_u2[idx] - gr$H23_u1[ii], 0))
            weights <- sweep(exp(-Gv), 1, e23, "*")
            weights <- sweep(weights, 2, w2[idx], "*")
            sy <- matrix(
                death_survival(t, rep(uy, J), rep(vv, each = N)),
                nrow = N, ncol = J
            )
            S <- S + wt1 * rowSums(weights * sy)

            is_mixed <- outer(uy, vv, ">")
            mixed_mass <- mixed_mass + wt1 * rowSums(weights * is_mixed)
        }

        if (any(!is.finite(S))) stop("Integration produced non-finite survival probabilities.")

        diagnostics[[k]] <- data.frame(
            tau = t,
            min_unclipped_survival = min(S),
            max_unclipped_survival = max(S),
            fraction_clipped = mean(S < 0 | S > 1),
            mixed_history_quadrature_mass = mean(mixed_mass),
            earlier_Y_copy_tail_quadrature_mass = mean(tail_mass)
        )

        risk[k] <- mean(1 - pmin(pmax(S, 0), 1))
    }

    list(risk = risk, diagnostics = do.call(rbind, diagnostics))
}

.pe2_kernel <- function(ctx, switches, rule, subject_chunk_size = 512L) {
    size <- .pe2_integer(subject_chunk_size, "subject_chunk_size")
    if (ctx$N <= size) return(.pe2_path_kernel(ctx, switches, rule))

    blocks <- split(seq_len(ctx$N), ceiling(seq_len(ctx$N) / size))
    risk <- numeric(length(ctx$tau))
    dg <- NULL

    for (ids in blocks) {
        cc <- ctx
        cc$N <- length(ids)
        cc$mult <- lapply(ctx$mult, function(x) x[ids])
        z <- .pe2_path_kernel(cc, switches, rule)
        weight <- length(ids) / ctx$N
        risk <- risk + weight * z$risk
        dd <- z$diagnostics

        avg <- c(
            "fraction_clipped",
            "mixed_history_quadrature_mass",
            "earlier_Y_copy_tail_quadrature_mass"
        )
        if (is.null(dg)) {
            dg <- dd
            dg[avg] <- lapply(dd[avg], function(x) x * weight)
        } else {
            dg$min_unclipped_survival <- pmin(
                dg$min_unclipped_survival,
                dd$min_unclipped_survival
            )
            dg$max_unclipped_survival <- pmax(
                dg$max_unclipped_survival,
                dd$max_unclipped_survival
            )
            for (nm in avg) dg[[nm]] <- dg[[nm]] + weight * dd[[nm]]
        }
    }

    list(risk = risk, diagnostics = dg)
}

.pe2_format_risks <- function(q, scale = c("risk", "survival")) {
    scale <- match.arg(scale)
    if (!all(c("tau", .pe2_risks) %in% names(q))) {
        stop("Missing five PE risk curves.")
    }
    q <- q[c("tau", .pe2_risks)]
    if (any(!is.finite(as.matrix(q)))) stop("Non-finite PE risk output.")

    r <- list(
        P_A_Y = q$psi1000 - q$psi0000,
        P_A_T2_Y = q$psi1100 - q$psi1000,
        P_A_T1_Y = q$psi1110 - q$psi1100,
        P_A_T1_T2_Y = q$psi1111 - q$psi1110
    )
    r$DE <- r$P_A_Y
    r$IE_total <- q$psi1111 - q$psi1000
    r$TE_path <- q$psi1111 - q$psi0000
    # Legacy aliases retained during regression testing.
    r$IE <- r$IE_total
    r$TE <- r$TE_path

    sign <- if (scale == "risk") 1 else -1
    for (nm in names(r)) {
        q[[nm]] <- sign * r[[nm]]
        q[[paste0(nm, "_RMTL")]] <- .compute_rmtl_from_psi(
            q$tau, r[[nm]], method = "trapezoid"
        )
        q[[paste0(nm, "_RMST")]] <- -q[[paste0(nm, "_RMTL")]]
    }

    core <- c("tau", .pe2_risks, .pe2_effects)
    q <- q[c(core, setdiff(names(q), core))]
    for (nm in .pe2_risks) {
        q[[sub("psi", "S", nm, fixed = TRUE)]] <- 1 - q[[nm]]
    }
    attr(q, "effect_scale") <- scale
    attr(q, "rmtl_settings") <- list(
        integration = "trapezoid",
        origin = c(tau = 0, risk = 0),
        grid = q$tau,
        meaning = "RMTL risk differences; RMST is their negative, on either effect scale"
    )
    q
}

# ------------------------------------------------------------
# Fit the shared observational S2TEM model once for PE.
# ------------------------------------------------------------
.fit_s2tem_pe_base <- function(
        dat,
        basec = character(0),
        a = 1,
        astar = 0,
        bh_method = "breslow",
        EMint = FALSE,
        verbose = FALSE,
        tau = 1:10,
        max_grid_points = 200L,
        v_chunk_size = 4096L
) {
    if (!exists(".fit_s2tem_backend", mode = "function") ||
        !exists(".predict_s2tem_ge", mode = "function")) {
        stop("Source s2tem_backend_refactor_v1.R before this file.")
    }

    if (is.null(basec)) basec <- character(0)
    .pe2_validate(dat, tau, basec, a, astar)
    tau <- sort(unique(as.numeric(tau)))

    if (!bh_method %in% c("breslow", "efron")) {
        stop("bh_method must be breslow or efron.")
    }
    if (!is.logical(EMint) || length(EMint) != 1L || is.na(EMint)) {
        stop("EMint must be TRUE or FALSE.")
    }
    max_grid_points <- if (is.null(max_grid_points)) {
        .Machine$integer.max
    } else {
        .pe2_integer(max_grid_points, "max_grid_points", 10L)
    }
    v_chunk_size <- .pe2_integer(v_chunk_size, "v_chunk_size")

    dat_use <- as.data.frame(dat)
    for (nm in basec) {
        if (is.character(dat_use[[nm]])) dat_use[[nm]] <- factor(dat_use[[nm]])
    }

    warnings_seen <- character(0)
    backend <- withCallingHandlers(
        .fit_s2tem_backend(
            dat = dat_use,
            tau = tau,
            basec = basec,
            a = a,
            astar = astar,
            bh_method = bh_method,
            EMint = EMint,
            verbose = verbose,
            max_grid_points = max_grid_points,
            v_chunk_size = v_chunk_size
        ),
        warning = function(w) {
            warnings_seen <<- c(warnings_seen, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    fit <- backend$fit
    mstate_dat_fit <- backend$mstate_dat
    if (is.null(fit)) stop("Cox fit object is missing.")
    if (is.null(mstate_dat_fit) || !all(c("trans", "status") %in% names(mstate_dat_fit))) {
        stop("Cannot determine event counts by transition from the fitted multistate data.")
    }

    # Structural-zero handling retained from the frozen chapter-3 estimator.
    tr_id <- suppressWarnings(as.integer(as.character(mstate_dat_fit$trans)))
    if (anyNA(tr_id)) tr_id <- suppressWarnings(as.integer(mstate_dat_fit$trans))
    if (anyNA(tr_id) || any(!tr_id %in% 1:5)) stop("Invalid transition coding.")

    event_count <- vapply(
        1:5,
        function(j) sum(mstate_dat_fit$status == 1L & tr_id == j, na.rm = TRUE),
        numeric(1)
    )
    names(event_count) <- as.character(1:5)
    structural_none <- (event_count["2"] == 0 && event_count["4"] == 0)

    cf <- stats::coef(fit)
    if (length(cf) == 0L) stop("Cox fit contains no coefficients.")
    if (any(is.infinite(cf))) {
        stop("Cox fit contains infinite coefficients; structural-zero handling only permits NA coefficients.")
    }

    na_idx <- which(is.na(cf))
    allowed_na <- logical(length(na_idx))
    .coef_on_transition <- function(nm, j) {
        grepl(paste0("\\.", j, "(?:$|:)"), nm, perl = TRUE)
    }

    if (length(na_idx)) {
        na_names <- names(cf)[na_idx]
        allowed_na <- vapply(
            na_names,
            function(nm) {
                structural_none &&
                    (.coef_on_transition(nm, 2L) || .coef_on_transition(nm, 4L))
            },
            logical(1)
        )
        if (!all(allowed_na)) {
            stop(
                "Cox fit contains NA coefficients outside structurally zero transitions 2/4: ",
                paste(na_names[!allowed_na], collapse = ", ")
            )
        }
        cf[na_idx] <- 0
        fit$coefficients <- cf
        backend$fit <- fit
    }

    fatal_warn <- warnings_seen[
        grepl("converg|infinite", warnings_seen, ignore.case = TRUE)
    ]
    if (length(fatal_warn)) {
        stop("Cox fit failed convergence checks: ", paste(unique(fatal_warn), collapse = " | "))
    }

    singular_warn <- warnings_seen[
        grepl("singular", warnings_seen, ignore.case = TRUE)
    ]
    structural_singularity_ok <- structural_none &&
        length(na_idx) > 0L && all(allowed_na)
    if (length(singular_warn) && !structural_singularity_ok) {
        stop("Cox fit failed singularity checks: ", paste(unique(singular_warn), collapse = " | "))
    }

    chapter2_reference <- .predict_s2tem_ge(
        backend = backend,
        scale = "risk",
        keep_fit_objects = TRUE
    )
    if (any(!is.finite(as.matrix(chapter2_reference[c("psi00", "psi10", "psi11")])))) {
        stop("Non-finite chapter-2 endpoint result.")
    }

    report_warn <- warnings_seen[
        !grepl("singular", warnings_seen, ignore.case = TRUE)
    ]
    if (length(report_warn)) {
        warning(paste(unique(report_warn), collapse = " | "), call. = FALSE)
    }

    structure(list(
        version = .pe2_version,
        dat = backend$dat,
        basec = basec,
        a = a,
        astar = astar,
        bh_method = bh_method,
        EMint = EMint,
        tau = tau,
        fit = backend$fit,
        mstate_dat = backend$mstate_dat,
        context = backend$context,
        chapter2_reference = chapter2_reference,
        hazards = backend$context$hazards,
        mult = backend$context$mult,
        beta = backend$context$beta,
        fit_warnings = unique(warnings_seen),
        settings = backend$settings,
        backend = backend,
        event_count = event_count,
        structural_zero_14_24 = structural_none
    ), class = "mims_s2tem_pe_base")
}

.predict_s2tem_pe <- function(
        base,
        scale = c("risk", "survival"),
        primary_rule = c("t2_overrides", "require_both"),
        t1_coupling = "rank_preserving",
        check_controls = FALSE,
        rule_sensitivity = FALSE,
        keep_base_fit = FALSE,
        subject_chunk_size = 512L,
        verbose = FALSE
) {
    scale <- match.arg(scale)
    primary_rule <- match.arg(primary_rule)
    if (!inherits(base, "mims_s2tem_pe_base") || !identical(base$version, .pe2_version)) {
        stop("base must be refitted with .fit_s2tem_pe_base() from this file.")
    }
    if (!identical(t1_coupling, "rank_preserving")) {
        stop("Only rank_preserving coupling is implemented.")
    }

    ctx <- base$context
    ref <- base$chapter2_reference
    if (verbose) message("One shared Cox fit; computing psi1100 and psi1110 on the GE grid.")

    k1100 <- .pe2_kernel(ctx, c(1, 1, 0, 0), primary_rule, subject_chunk_size)
    k1110 <- .pe2_kernel(ctx, c(1, 1, 1, 0), primary_rule, subject_chunk_size)

    q <- data.frame(
        tau = ctx$tau,
        psi0000 = ref$psi00,
        psi1000 = ref$psi10,
        psi1100 = k1100$risk,
        psi1110 = k1110$risk,
        psi1111 = ref$psi11
    )

    controls <- NULL
    diagnostics <- rbind(
        data.frame(curve = "psi1100", k1100$diagnostics),
        data.frame(curve = "psi1110", k1110$diagnostics)
    )

    if (isTRUE(check_controls)) {
        regimes <- list(
            psi0000 = c(0, 0, 0, 0),
            psi1000 = c(1, 0, 0, 0),
            psi1111 = c(1, 1, 1, 1)
        )
        controls <- do.call(rbind, lapply(names(regimes), function(nm) {
            z <- .pe2_kernel(ctx, regimes[[nm]], primary_rule, subject_chunk_size)
            diagnostics <<- rbind(
                diagnostics,
                data.frame(curve = nm, z$diagnostics)
            )
            data.frame(
                tau = ctx$tau,
                curve = nm,
                new_kernel = z$risk,
                reference = q[[nm]],
                difference = z$risk - q[[nm]]
            )
        }))
        if (max(abs(controls$difference)) > 1e-8) {
            stop(
                "Generalized kernel failed the common-endpoint check; max error=",
                format(max(abs(controls$difference)))
            )
        }
    }

    out <- .pe2_format_risks(q, scale)
    bounded <- all(as.matrix(q[.pe2_risks]) >= 0 & as.matrix(q[.pe2_risks]) <= 1)
    monotone <- all(vapply(
        q[.pe2_risks],
        function(x) all(diff(x) >= -1e-10),
        logical(1)
    ))

    if (any(diagnostics$fraction_clipped > 0)) {
        warning(
            "Numerical survival probabilities were clipped; inspect path_diagnostics.",
            call. = FALSE
        )
    }
    if (!monotone) {
        warning("At least one risk curve decreases on the supplied grid.", call. = FALSE)
    }

    if (!isTRUE(keep_base_fit)) {
        attr(ref, "fit") <- NULL
        attr(ref, "mstate_dat") <- NULL
    }

    attr(out, "chapter2_reference") <- ref
    attr(out, "exactfast_settings") <- base$settings
    attr(out, "pathspecific_grid_settings") <- c(
        base$settings,
        list(
            scheme = "shared_ge_backend",
            subject_chunk_size = subject_chunk_size,
            tau_inserted = FALSE
        )
    )
    attr(out, "path_diagnostics") <- diagnostics
    attr(out, "control_checks") <- controls
    attr(out, "fit_warnings") <- base$fit_warnings
    attr(out, "estimand_definition") <- list(
        version = .pe2_version,
        switches = c("zY", "z2", "z1", "z12"),
        order = c("0000", "1000", "1100", "1110", "1111"),
        t1_coupling = t1_coupling,
        latent_innovations_independent = TRUE,
        primary_rule = primary_rule,
        rank_map_interpolation = "linear inverse of thinned H12",
        common_endpoints = "same shared GE backend",
        causal_identification_established = FALSE,
        numerical_accuracy_established = FALSE,
        old_path_oracle_compatible = FALSE
    )

    sg <- if (scale == "risk") 1 else -1
    attr(out, "summary") <- data.frame(
        chapter2_endpoints_preserved = TRUE,
        IE_matches_reference = max(abs(out$IE_total - sg * ref$IE)) <= 1e-12,
        decomposition_adds_up = max(abs(
            out$P_A_T2_Y + out$P_A_T1_Y + out$P_A_T1_T2_Y - out$IE_total
        )) <= 1e-12,
        controls_evaluated = isTRUE(check_controls),
        max_control_error = if (is.null(controls)) NA_real_ else max(abs(controls$difference)),
        all_risks_bounded = bounded,
        all_risks_nondecreasing = monotone,
        any_evaluated_survival_clipped = any(diagnostics$fraction_clipped > 0),
        causal_identification_established = FALSE
    )

    if (isTRUE(rule_sensitivity)) {
        alt_rule <- if (primary_rule == "t2_overrides") "require_both" else "t2_overrides"
        zz <- .pe2_kernel(ctx, c(1, 1, 1, 0), alt_rule, subject_chunk_size)
        qa <- q
        qa$psi1110 <- zz$risk
        alt <- .pe2_format_risks(qa, scale)
        attr(out, "rule_sensitivity") <- do.call(
            rbind,
            lapply(.pe2_effects, function(nm) {
                data.frame(
                    tau = ctx$tau,
                    metric = nm,
                    primary = out[[nm]],
                    alternative = alt[[nm]],
                    difference = alt[[nm]] - out[[nm]]
                )
            })
        )
        attr(out, "alternative_rule") <- alt_rule
        attr(out, "alternative_result") <- alt
        attr(out, "alternative_diagnostics") <- zz$diagnostics
    }

    if (isTRUE(keep_base_fit)) attr(out, "base_fit") <- base
    out
}

# ------------------------------------------------------------
# Stage-2 point estimator: no dependency on the legacy PE file.
# ------------------------------------------------------------
mims_s2tem_pe_point <- function(
        dat,
        tau = 1:10,
        basec = character(0),
        a = 1,
        astar = 0,
        bh_method = "breslow",
        EMint = FALSE,
        scale = c("risk", "survival"),
        verbose = FALSE,
        keep_base_fit = FALSE,
        max_grid_points = 200L,
        v_chunk_size = 4096L,
        primary_rule = c("t2_overrides", "require_both"),
        t1_coupling = "rank_preserving",
        check_controls = FALSE,
        rule_sensitivity = FALSE,
        subject_chunk_size = 512L
) {
    scale <- match.arg(scale)
    primary_rule <- match.arg(primary_rule)

    base <- .fit_s2tem_pe_base(
        dat = dat,
        basec = basec,
        a = a,
        astar = astar,
        bh_method = bh_method,
        EMint = EMint,
        verbose = verbose,
        tau = tau,
        max_grid_points = max_grid_points,
        v_chunk_size = v_chunk_size
    )

    .predict_s2tem_pe(
        base = base,
        scale = scale,
        primary_rule = primary_rule,
        t1_coupling = t1_coupling,
        check_controls = check_controls,
        rule_sensitivity = rule_sensitivity,
        keep_base_fit = keep_base_fit,
        subject_chunk_size = subject_chunk_size,
        verbose = verbose
    )
}

# ------------------------------------------------------------
# Audit PE common endpoints against a GE result.
# ------------------------------------------------------------
audit_s2tem_pe_ge <- function(path_out, ge_out, tolerance = 1e-6) {
    if (!is.data.frame(path_out) || !is.data.frame(ge_out)) {
        stop("Both inputs must be data.frames.")
    }
    if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
        stop("Invalid tolerance.")
    }
    if (anyDuplicated(path_out$tau) || anyDuplicated(ge_out$tau)) {
        stop("Duplicate tau values.")
    }

    ix <- match(path_out$tau, ge_out$tau)
    if (anyNA(ix)) stop("GE output lacks requested tau values.")
    j <- ge_out[ix, , drop = FALSE]
    p <- path_out
    pairs <- c(psi0000 = "psi00", psi1000 = "psi10", psi1111 = "psi11")

    comparison <- do.call(rbind, lapply(names(pairs), function(nm) {
        if (!nm %in% names(p) || !pairs[nm] %in% names(j)) stop("Missing risk columns.")
        d <- p[[nm]] - j[[pairs[nm]]]
        data.frame(
            tau = p$tau,
            metric = paste0(nm, "_vs_", pairs[nm]),
            PE = p[[nm]],
            GE = j[[pairs[nm]]],
            difference = d,
            abs_difference = abs(d),
            within_tolerance = abs(d) <= tolerance
        )
    }))

    pp <- list(
        DE = p$psi1000 - p$psi0000,
        IE = p$psi1111 - p$psi1000,
        TE = p$psi1111 - p$psi0000
    )
    jj <- list(
        DE = j$psi10 - j$psi00,
        IE = j$psi11 - j$psi10,
        TE = j$psi11 - j$psi00
    )
    for (nm in names(pp)) {
        d <- pp[[nm]] - jj[[nm]]
        comparison <- rbind(
            comparison,
            data.frame(
                tau = p$tau,
                metric = nm,
                PE = pp[[nm]],
                GE = jj[[nm]],
                difference = d,
                abs_difference = abs(d),
                within_tolerance = abs(d) <= tolerance
            )
        )
    }

    ep <- comparison$abs_difference[
        comparison$metric %in% paste0(names(pairs), "_vs_", pairs)
    ]
    ie <- comparison$abs_difference[comparison$metric == "IE"]

    list(
        summary = data.frame(
            endpoint_agreement = all(ep <= tolerance),
            IE_agreement = all(ie <= tolerance),
            max_endpoint_error = max(ep),
            max_IE_error = max(ie),
            tolerance = tolerance
        ),
        comparison = comparison
    )
}

# ------------------------------------------------------------
# Regression test: frozen legacy PE vs refactor-v2 PE.
# ------------------------------------------------------------
verify_s2tem_pe_v2_equivalence <- function(
        dat,
        legacy_pe_fun,
        tau = 1:10,
        basec = character(0),
        a = 1,
        astar = 0,
        EMint = FALSE,
        bh_method = "breslow",
        scale = "risk",
        max_grid_points = 200L,
        v_chunk_size = 4096L,
        primary_rule = "t2_overrides",
        t1_coupling = "rank_preserving",
        check_controls = TRUE,
        rule_sensitivity = FALSE,
        subject_chunk_size = 512L,
        tolerance = 1e-12
) {
    if (!is.function(legacy_pe_fun)) stop("legacy_pe_fun must be a function.")

    capture_run <- function(fun) {
        ww <- character(0)
        value <- withCallingHandlers(
            fun(),
            warning = function(w) {
                ww <<- c(ww, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        )
        list(value = value, warnings = unique(ww))
    }

    old_call <- capture_run(function() {
        legacy_pe_fun(
            dat = dat,
            tau = tau,
            basec = basec,
            a = a,
            astar = astar,
            bh_method = bh_method,
            EMint = EMint,
            scale = scale,
            verbose = FALSE,
            keep_base_fit = FALSE,
            max_grid_points = max_grid_points,
            v_chunk_size = v_chunk_size,
            primary_rule = primary_rule,
            t1_coupling = t1_coupling,
            check_controls = check_controls,
            rule_sensitivity = rule_sensitivity,
            subject_chunk_size = subject_chunk_size
        )
    })

    new_call <- capture_run(function() {
        mims_s2tem_pe_point(
            dat = dat,
            tau = tau,
            basec = basec,
            a = a,
            astar = astar,
            bh_method = bh_method,
            EMint = EMint,
            scale = scale,
            verbose = FALSE,
            keep_base_fit = FALSE,
            max_grid_points = max_grid_points,
            v_chunk_size = v_chunk_size,
            primary_rule = primary_rule,
            t1_coupling = t1_coupling,
            check_controls = check_controls,
            rule_sensitivity = rule_sensitivity,
            subject_chunk_size = subject_chunk_size
        )
    })

    old <- old_call$value
    new <- new_call$value
    common_num <- intersect(names(old), names(new))
    common_num <- common_num[vapply(old[common_num], is.numeric, logical(1)) &
                             vapply(new[common_num], is.numeric, logical(1))]

    by_metric <- do.call(rbind, lapply(common_num, function(nm) {
        data.frame(
            metric = nm,
            max_abs_error = max(abs(old[[nm]] - new[[nm]]), na.rm = TRUE)
        )
    }))

    main_metrics <- c(
        .pe2_risks, .pe2_effects,
        paste0(.pe2_effects, "_RMTL"),
        paste0(.pe2_effects, "_RMST"),
        paste0("S", sub("psi", "", .pe2_risks, fixed = TRUE))
    )
    main_metrics <- intersect(main_metrics, by_metric$metric)
    main_err <- by_metric$max_abs_error[match(main_metrics, by_metric$metric)]
    max_error <- if (length(main_err)) max(main_err, na.rm = TRUE) else NA_real_

    old_summary <- attr(old, "summary", exact = TRUE)
    new_summary <- attr(new, "summary", exact = TRUE)
    summary_comparison <- all.equal(old_summary, new_summary, check.attributes = FALSE)
    summary_equal <- isTRUE(summary_comparison)

    old_ctrl <- attr(old, "control_checks", exact = TRUE)
    new_ctrl <- attr(new, "control_checks", exact = TRUE)
    control_max_error <- NA_real_
    controls_equal <- TRUE
    if (!is.null(old_ctrl) || !is.null(new_ctrl)) {
        controls_equal <- !is.null(old_ctrl) && !is.null(new_ctrl) &&
            nrow(old_ctrl) == nrow(new_ctrl)
        if (controls_equal) {
            nums <- intersect(
                names(old_ctrl)[vapply(old_ctrl, is.numeric, logical(1))],
                names(new_ctrl)[vapply(new_ctrl, is.numeric, logical(1))]
            )
            if (length(nums)) {
                control_max_error <- max(vapply(
                    nums,
                    function(nm) max(abs(old_ctrl[[nm]] - new_ctrl[[nm]]), na.rm = TRUE),
                    numeric(1)
                ))
                controls_equal <- is.finite(control_max_error) && control_max_error <= tolerance
            }
        }
    }

    old_dg <- attr(old, "path_diagnostics", exact = TRUE)
    new_dg <- attr(new, "path_diagnostics", exact = TRUE)
    diagnostic_max_error <- NA_real_
    diagnostics_equal <- TRUE
    if (!is.null(old_dg) || !is.null(new_dg)) {
        diagnostics_equal <- !is.null(old_dg) && !is.null(new_dg) &&
            nrow(old_dg) == nrow(new_dg) && identical(old_dg$curve, new_dg$curve)
        if (diagnostics_equal) {
            nums <- intersect(
                names(old_dg)[vapply(old_dg, is.numeric, logical(1))],
                names(new_dg)[vapply(new_dg, is.numeric, logical(1))]
            )
            if (length(nums)) {
                diagnostic_max_error <- max(vapply(
                    nums,
                    function(nm) max(abs(old_dg[[nm]] - new_dg[[nm]]), na.rm = TRUE),
                    numeric(1)
                ))
                diagnostics_equal <- is.finite(diagnostic_max_error) &&
                    diagnostic_max_error <= tolerance
            }
        }
    }

    list(
        passed = is.finite(max_error) && max_error <= tolerance &&
            summary_equal && controls_equal && diagnostics_equal,
        max_error = max_error,
        by_metric = by_metric,
        summary_equal = summary_equal,
        summary_comparison = summary_comparison,
        controls_equal = controls_equal,
        control_max_error = control_max_error,
        diagnostics_equal = diagnostics_equal,
        diagnostic_max_error = diagnostic_max_error,
        legacy_warnings = old_call$warnings,
        refactor_warnings = new_call$warnings,
        legacy = old,
        refactor = new
    )
}
