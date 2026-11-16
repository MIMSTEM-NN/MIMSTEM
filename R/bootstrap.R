# ============================================================
# MIMS unified nonparametric bootstrap engine -- Wald CI
# Version: 2026-11-14 (PM Delta-Wald extension)
#
# Purpose
#   One bootstrap engine for:
#     MIMS-TEM
#     MIMS-S2TEM-GE
#     MIMS-S2TEM-PE
#
# Statistical convention
#   Individual-level nonparametric bootstrap.
#   SE_boot = SD(theta_hat^b)
#   Wald CI = theta_hat +/- z_(1-alpha/2) * SE_boot
#   PM = indirect / total; Delta-method SE uses PAIRED bootstrap covariance.
#
# Notes
#   * Default stratified = FALSE to match an ordinary subject-level
#     nonparametric bootstrap.  stratified = TRUE is available only as
#     an explicit sensitivity / implementation option.
#   * Every bootstrap sample refits the entire estimator once and all
#     requested metrics are extracted from that same fit.
#   * Failed replicates remain recorded and are never silently replaced.
# ============================================================


.mims_boot_integer <- function(x, nm, minimum = 1L) {
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
        x != floor(x) || x < minimum || x > .Machine$integer.max) {
        stop(nm, " must be an integer >= ", minimum, ".")
    }
    as.integer(x)
}


.mims_boot_validate_ci <- function(conf_level) {
    if (length(conf_level) != 1L || !is.numeric(conf_level) ||
        !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
        stop("conf_level must be a single number strictly between 0 and 1.")
    }
    invisible(TRUE)
}


.mims_boot_restore_seed <- function(old_seed_exists, old_seed) {
    if (old_seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = ".Random.seed", envir = .GlobalEnv)
    }
    invisible(NULL)
}


.mims_boot_sample_index <- function(dat, stratified = FALSE, exposure = "A") {
    n <- nrow(dat)
    
    if (!isTRUE(stratified)) {
        return(sample.int(n, size = n, replace = TRUE))
    }
    
    if (!exposure %in% names(dat)) {
        stop("stratified=TRUE requires exposure column: ", exposure)
    }
    
    grp <- split(seq_len(n), dat[[exposure]], drop = TRUE)
    if (length(grp) < 2L) {
        stop("stratified bootstrap requires at least two exposure groups.")
    }
    
    idx <- unlist(
        lapply(grp, function(ii) {
            ii[sample.int(length(ii), size = length(ii), replace = TRUE)]
        }),
        use.names = FALSE
    )
    
    # Randomize row order after within-stratum resampling.
    idx[sample.int(length(idx))]
}


.mims_boot_extract <- function(out, tau, metrics) {
    if (!is.data.frame(out)) {
        stop("The point estimator must return a data.frame.")
    }
    if (!"tau" %in% names(out)) {
        stop("Estimator output must contain a `tau` column.")
    }
    if (anyDuplicated(out$tau)) {
        stop("Estimator output contains duplicated tau values.")
    }
    
    ix <- match(tau, out$tau)
    if (anyNA(ix)) {
        stop("Estimator output does not contain every requested tau value.")
    }
    
    missing_metrics <- setdiff(metrics, names(out))
    if (length(missing_metrics)) {
        stop(
            "Estimator output is missing metrics: ",
            paste(missing_metrics, collapse = ", ")
        )
    }
    
    z <- out[ix, c("tau", metrics), drop = FALSE]
    
    mat <- as.matrix(z[metrics])
    storage.mode(mat) <- "double"
    if (any(!is.finite(mat))) {
        stop("Estimator returned non-finite values in requested bootstrap metrics.")
    }
    
    z
}


# Ratio specifications are derived from joint bootstrap replicates without
# refitting the point estimator or changing the sampling scheme.
.mims_boot_validate_ratios <- function(ratio_specs, metrics) {
    if (is.null(ratio_specs)) {
        return(data.frame(
            metric = character(), numerator = character(),
            denominator = character(), stringsAsFactors = FALSE
        ))
    }
    required <- c("metric", "numerator", "denominator")
    if (!is.data.frame(ratio_specs) ||
        !all(required %in% names(ratio_specs))) {
        stop("ratio_specs must be a data.frame with metric, numerator, denominator.")
    }
    specs <- ratio_specs[, required, drop = FALSE]
    if (anyNA(specs) || !all(vapply(specs, is.character, logical(1))) ||
        any(!nzchar(as.matrix(specs))) || anyDuplicated(specs$metric) ||
        any(specs$metric %in% metrics)) {
        stop("ratio_specs has missing, duplicate or invalid derived metric names.")
    }
    absent <- setdiff(unique(c(specs$numerator, specs$denominator)), metrics)
    if (length(absent)) {
        stop("ratio_specs needs jointly bootstrapped metrics: ",
             paste(absent, collapse = ", "))
    }
    specs
}

.mims_boot_ratio_specs_tem_ge <- function(scale = "risk") {
    data.frame(
        metric = c(paste0("PM_", scale), "PM_RMTL"),
        numerator = c("IE", "IE_RMTL"),
        denominator = c("TE", "TE_RMTL"),
        stringsAsFactors = FALSE
    )
}

.mims_boot_ratio_specs_pe <- function(scale = "risk") {
    data.frame(
        metric = c(
            paste0("PM_path_", scale), paste0("prop_T2_", scale),
            paste0("prop_T1_", scale), paste0("prop_T1_T2_", scale),
            "PM_path_RMTL", "prop_T2_RMTL", "prop_T1_RMTL", "prop_T1_T2_RMTL"
        ),
        numerator = c(
            "IE_total", "P_A_T2_Y", "P_A_T1_Y", "P_A_T1_T2_Y",
            "IE_total_RMTL", "P_A_T2_Y_RMTL", "P_A_T1_Y_RMTL",
            "P_A_T1_T2_Y_RMTL"
        ),
        denominator = c(
            rep("TE", 4L), rep("TE_RMTL", 4L)
        ),
        stringsAsFactors = FALSE
    )
}


.bootstrap_mims_wald <- function(
        dat,
        estimator,
        estimator_args = list(),
        tau,
        metrics,
        B = 200L,
        seed = 1L,
        conf_level = 0.95,
        stratified = FALSE,
        exposure = "A",
        reset_id = TRUE,
        min_success_fraction = 0.80,
        keep_replicates = FALSE,
        verbose = TRUE,
        ratio_specs = NULL,
        ratio_epsilon = 1e-6
) {
    if (!is.data.frame(dat) || nrow(dat) < 2L) {
        stop("dat must be a subject-level data.frame with at least two rows.")
    }
    if (!is.function(estimator)) {
        stop("estimator must be a function.")
    }
    if (!is.list(estimator_args)) {
        stop("estimator_args must be a list.")
    }
    
    tau <- sort(unique(as.numeric(tau)))
    if (!length(tau) || any(!is.finite(tau)) || any(tau <= 0)) {
        stop("tau must contain positive finite evaluation times.")
    }
    
    if (!is.character(metrics) || !length(metrics) ||
        anyNA(metrics) || any(!nzchar(metrics)) || anyDuplicated(metrics)) {
        stop("metrics must be a non-empty character vector with unique names.")
    }
    
    B <- .mims_boot_integer(B, "B", 2L)
    seed <- .mims_boot_integer(seed, "seed", 0L)
    .mims_boot_validate_ci(conf_level)
    ratio_specs <- .mims_boot_validate_ratios(ratio_specs, metrics)
    if (length(ratio_epsilon) != 1L || !is.numeric(ratio_epsilon) ||
        !is.finite(ratio_epsilon) || ratio_epsilon < 0) {
        stop("ratio_epsilon must be one finite nonnegative numeric value.")
    }
    
    if (length(stratified) != 1L || !is.logical(stratified) || is.na(stratified)) {
        stop("stratified must be TRUE or FALSE.")
    }
    if (length(reset_id) != 1L || !is.logical(reset_id) || is.na(reset_id)) {
        stop("reset_id must be TRUE or FALSE.")
    }
    if (length(min_success_fraction) != 1L ||
        !is.numeric(min_success_fraction) ||
        !is.finite(min_success_fraction) ||
        min_success_fraction <= 0 ||
        min_success_fraction > 1) {
        stop("min_success_fraction must lie in (0, 1].")
    }
    
    # Preserve the caller's RNG state.
    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (old_seed_exists) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    on.exit(.mims_boot_restore_seed(old_seed_exists, old_seed), add = TRUE)
    
    # Point estimate: use the same estimator and arguments as every replicate.
    point_args <- utils::modifyList(
        list(dat = dat, tau = tau),
        estimator_args,
        keep.null = TRUE
    )
    
    point <- do.call(estimator, point_args)
    point_extract <- .mims_boot_extract(point, tau = tau, metrics = metrics)
    
    K <- length(tau)
    M <- length(metrics)
    
    arr <- array(
        NA_real_,
        dim = c(K, M, B),
        dimnames = list(
            tau = as.character(tau),
            metric = metrics,
            replicate = paste0("b", seq_len(B))
        )
    )
    
    status <- vector("list", B)
    ok_index <- logical(B)
    
    set.seed(seed)
    
    for (b in seq_len(B)) {
        start_time <- proc.time()[[3]]
        
        idx <- .mims_boot_sample_index(
            dat = dat,
            stratified = stratified,
            exposure = exposure
        )
        
        dd <- dat[idx, , drop = FALSE]
        rownames(dd) <- NULL
        
        # Duplicated subject IDs after resampling are expected.  The estimator
        # works on one row per resampled subject, so reset the identifier.
        if (isTRUE(reset_id) && "id" %in% names(dd)) {
            dd$id <- seq_len(nrow(dd))
        }
        
        warnings_seen <- character(0)
        error_text <- ""
        
        val <- tryCatch(
            withCallingHandlers(
                {
                    args_b <- utils::modifyList(
                        list(dat = dd, tau = tau),
                        estimator_args,
                        keep.null = TRUE
                    )
                    z <- do.call(estimator, args_b)
                    .mims_boot_extract(z, tau = tau, metrics = metrics)
                },
                warning = function(w) {
                    warnings_seen <<- c(warnings_seen, conditionMessage(w))
                    invokeRestart("muffleWarning")
                }
            ),
            error = function(e) {
                error_text <<- conditionMessage(e)
                NULL
            }
        )
        
        ok <- !is.null(val)
        ok_index[b] <- ok
        
        if (ok) {
            arr[, , b] <- as.matrix(val[metrics])
        }
        
        status[[b]] <- data.frame(
            replicate = b,
            status = if (!ok) {
                "failed"
            } else if (length(warnings_seen)) {
                "ok_with_warning"
            } else {
                "ok"
            },
            seconds = proc.time()[[3]] - start_time,
            warnings = paste(unique(warnings_seen), collapse = " | "),
            error = error_text,
            stringsAsFactors = FALSE
        )
        
        if (isTRUE(verbose) && (b == 1L || b %% 10L == 0L || b == B)) {
            message(
                "Bootstrap ", b, "/", B, ": ",
                sum(ok_index[seq_len(b)]), " successful."
            )
        }
    }
    
    n_ok <- sum(ok_index)
    required_ok <- max(2L, ceiling(B * min_success_fraction))
    enough <- n_ok >= required_ok
    
    status_df <- do.call(rbind, status)
    
    alpha <- 1 - conf_level
    zcrit <- stats::qnorm(1 - alpha / 2)
    
    # Long-format interval table: one row per tau x estimand.
    interval_rows <- vector("list", K * M)
    rr <- 0L
    
    for (k in seq_len(K)) {
        for (m in seq_len(M)) {
            rr <- rr + 1L
            metric <- metrics[m]
            estimate <- point_extract[[metric]][k]
            
            se <- lower <- upper <- NA_real_
            
            if (enough) {
                x <- arr[k, m, ok_index]
                se <- stats::sd(x, na.rm = TRUE)
                lower <- estimate - zcrit * se
                upper <- estimate + zcrit * se
            }
            
            interval_rows[[rr]] <- data.frame(
                tau = tau[k],
                metric = metric,
                estimate = estimate,
                se_boot = se,
                lower = lower,
                upper = upper,
                conf_level = conf_level,
                stringsAsFactors = FALSE
            )
        }
    }
    
    intervals <- do.call(rbind, interval_rows)

    # Add ratios using the PAIRED bootstrap covariance of (numerator, total).
    # This is the first-order Delta SE, NOT sd(IE_b/TE_b) and NOT a percentile CI.
    # The original-sample ratio is the Wald interval center.
    ratio_diag_rows <- list()
    ratio_interval_rows <- list()
    if (nrow(ratio_specs) > 0L) {
        row_id <- 0L
        for (j in seq_len(nrow(ratio_specs))) {
            ratio_name <- ratio_specs$metric[j]
            num_name <- ratio_specs$numerator[j]
            den_name <- ratio_specs$denominator[j]
            mi <- match(num_name, metrics)
            md <- match(den_name, metrics)

            # Expose point estimates without changing the point estimator itself.
            pn <- as.numeric(point[[num_name]])
            pd <- as.numeric(point[[den_name]])
            good_point <- is.finite(pn) & is.finite(pd) &
                (abs(pd) > ratio_epsilon)
            ratio_points <- rep(NA_real_, nrow(point))
            ratio_points[good_point] <- pn[good_point] / pd[good_point]
            point[[ratio_name]] <- ratio_points

            for (k in seq_len(K)) {
                row_id <- row_id + 1L
                num0 <- point_extract[[num_name]][k]
                den0 <- point_extract[[den_name]][k]
                is_defined <- is.finite(num0) && is.finite(den0) &&
                    abs(den0) > ratio_epsilon
                ratio0 <- if (is_defined) num0 / den0 else NA_real_
                ratio_se <- ratio_lo <- ratio_hi <- NA_real_
                paired_cov <- NA_real_
                den_cross_zero <- NA
                b_pair <- n_ok

                if (n_ok >= 2L) {
                    boot_num <- as.numeric(arr[k, mi, ok_index])
                    boot_den <- as.numeric(arr[k, md, ok_index])
                    den_cross_zero <- min(boot_den) <= 0 && max(boot_den) >= 0
                    paired_cov <- stats::cov(boot_num, boot_den)

                    if (enough && is_defined) {
                        sigma <- stats::cov(cbind(boot_num, boot_den))
                        grad <- c(1 / den0, -num0 / den0^2)
                        delta_var <- as.numeric(crossprod(grad, sigma %*% grad))
                        # Tiny negative values can arise from floating point roundoff.
                        if (is.finite(delta_var) && delta_var >= -1e-12) {
                            ratio_se <- sqrt(max(0, delta_var))
                            ratio_lo <- ratio0 - zcrit * ratio_se
                            ratio_hi <- ratio0 + zcrit * ratio_se
                        }
                    }
                }
                ratio_interval_rows[[row_id]] <- data.frame(
                    tau = tau[k], metric = ratio_name, estimate = ratio0,
                    se_boot = ratio_se, lower = ratio_lo, upper = ratio_hi,
                    conf_level = conf_level, stringsAsFactors = FALSE
                )
                ratio_diag_rows[[row_id]] <- data.frame(
                    tau = tau[k], metric = ratio_name,
                    numerator = num_name, denominator = den_name,
                    denominator_estimate = den0,
                    denominator_near_zero = !is_defined,
                    denominator_bootstrap_crosses_zero = den_cross_zero,
                    B_valid_paired = b_pair,
                    paired_covariance = paired_cov,
                    stringsAsFactors = FALSE
                )
            }
        }
        intervals <- rbind(intervals, do.call(rbind, ratio_interval_rows))
    }
    ratio_diagnostics <- if (length(ratio_diag_rows)) {
        do.call(rbind, ratio_diag_rows)
    } else {
        data.frame()
    }
    
    diagnostics <- list(
        B_total = B,
        B_ok = n_ok,
        B_fail = B - n_ok,
        fail_rate = (B - n_ok) / B,
        success_fraction = n_ok / B,
        min_success_fraction = min_success_fraction,
        intervals_reported = enough,
        seed = seed,
        stratified = stratified,
        exposure = exposure,
        ci_method = "wald_bootstrap_se",
        conf_level = conf_level,
        ratio_ci_method = "paired_bootstrap_covariance_delta_wald",
        ratio_epsilon = ratio_epsilon,
        ratio_diagnostics = ratio_diagnostics,
        replicate_status = status_df,
        sample_errors = utils::head(status_df$error[status_df$status == "failed"], 10L),
        note = paste(
            "All requested metrics are estimated jointly within each resample.",
            "SEs are bootstrap standard deviations and intervals are Wald intervals",
            "centered on the original-sample point estimate."
        )
    )
    
    if (!enough) {
        warning(
            "Too few successful bootstrap replicates; SEs and Wald intervals are NA. ",
            "Inspect $diagnostics$replicate_status.",
            call. = FALSE
        )
    } else if (n_ok < B) {
        warning(
            "Some bootstrap replicates failed; intervals use successful replicates only. ",
            "Inspect $diagnostics$replicate_status.",
            call. = FALSE
        )
    }
    
    out <- list(
        point = point,
        intervals = intervals,
        diagnostics = diagnostics
    )
    
    if (isTRUE(keep_replicates)) {
        out$replicates <- arr
    }
    
    class(out) <- c("mims_bootstrap", "list")
    out
}


# ============================================================
# S2TEM-GE wrapper
# ============================================================

bootstrap_mims_s2tem_ge <- function(
        dat,
        tau = 1:10,
        basec = character(0),
        a = 1,
        astar = 0,
        bh_method = "breslow",
        EMint = FALSE,
        scale = c("risk", "survival"),
        B = 200L,
        seed = 1L,
        conf_level = 0.95,
        stratified = FALSE,
        min_success_fraction = 0.80,
        keep_replicates = FALSE,
        verbose = TRUE,
        v_chunk_size = 4096L,
        max_grid_points = 200L
) {
    if (!exists("mims_s2tem_ge_point", mode = "function")) {
        stop("Source the verified S2TEM backend before calling this function.")
    }
    
    scale <- match.arg(scale)
    
    metrics <- c(
        "DE", "IE", "TE",
        "DE_RMTL", "IE_RMTL", "TE_RMTL"
    )
    
    .bootstrap_mims_wald(
        dat = dat,
        estimator = mims_s2tem_ge_point,
        estimator_args = list(
            basec = basec,
            a = a,
            astar = astar,
            bh_method = bh_method,
            EMint = EMint,
            scale = scale,
            verbose = FALSE,
            keep_fit_objects = FALSE,
            v_chunk_size = v_chunk_size,
            max_grid_points = max_grid_points
        ),
        tau = tau,
        metrics = metrics,
        ratio_specs = .mims_boot_ratio_specs_tem_ge(scale),
        B = B,
        seed = seed,
        conf_level = conf_level,
        stratified = stratified,
        exposure = "A",
        reset_id = TRUE,
        min_success_fraction = min_success_fraction,
        keep_replicates = keep_replicates,
        verbose = verbose
    )
}


# ============================================================
# S2TEM-PE wrapper
# ============================================================

bootstrap_mims_s2tem_pe <- function(
        dat,
        tau = 1:10,
        basec = character(0),
        a = 1,
        astar = 0,
        bh_method = "breslow",
        EMint = FALSE,
        scale = c("risk", "survival"),
        B = 200L,
        seed = 1L,
        conf_level = 0.95,
        stratified = FALSE,
        min_success_fraction = 0.80,
        keep_replicates = FALSE,
        verbose = TRUE,
        max_grid_points = 200L,
        v_chunk_size = 4096L,
        primary_rule = c("t2_overrides", "require_both"),
        t1_coupling = "rank_preserving",
        subject_chunk_size = 512L
) {
    if (!exists("mims_s2tem_pe_point", mode = "function")) {
        stop("Source the verified S2TEM-PE refactor before calling this function.")
    }
    
    scale <- match.arg(scale)
    primary_rule <- match.arg(primary_rule)
    
    # PE bootstrap reports the total, direct, and total indirect effects,
    # alongside the three path-specific components and their RMTL analogues.
    # Each metric's SE is computed from the same joint bootstrap replicates;
    # the shared engine reports Wald confidence intervals.
    # Total/path proportions get Delta-Wald CIs from paired IE/path and TE covariance.
    metrics <- c(
        "TE",
        "DE",
        "IE_total",
        "P_A_T2_Y",
        "P_A_T1_Y",
        "P_A_T1_T2_Y",
        "TE_RMTL",
        "DE_RMTL",
        "IE_total_RMTL",
        "P_A_T2_Y_RMTL",
        "P_A_T1_Y_RMTL",
        "P_A_T1_T2_Y_RMTL"
    )
    
    .bootstrap_mims_wald(
        dat = dat,
        estimator = mims_s2tem_pe_point,
        estimator_args = list(
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
            check_controls = FALSE,
            rule_sensitivity = FALSE,
            subject_chunk_size = subject_chunk_size
        ),
        tau = tau,
        metrics = metrics,
        ratio_specs = .mims_boot_ratio_specs_pe(scale),
        B = B,
        seed = seed,
        conf_level = conf_level,
        stratified = stratified,
        exposure = "A",
        reset_id = TRUE,
        min_success_fraction = min_success_fraction,
        keep_replicates = keep_replicates,
        verbose = verbose
    )
}


# ============================================================
# Minimal print method for development-stage inspection
# ============================================================

print.mims_bootstrap <- function(x, ...) {
    cat("MIMS bootstrap result\n")
    cat("CI method: ", x$diagnostics$ci_method, "\n", sep = "")
    cat("Confidence level: ", x$diagnostics$conf_level, "\n", sep = "")
    cat(
        "Successful replicates: ",
        x$diagnostics$B_ok, "/", x$diagnostics$B_total, "\n",
        sep = ""
    )
    print(x$intervals, row.names = FALSE)
    invisible(x)
}
# ============================================================
# MIMS-TEM wrapper for the unified Wald bootstrap engine
# Requires:
#   mims_tem_refactor_v1.R
#   mims_bootstrap_wald_v1.R
# ============================================================

bootstrap_mims_tem <- function(
        dat,
        tau = 1:10,
        basec = character(0),
        a = 1,
        astar = 0,
        bh_method = "breslow",
        EMint = FALSE,
        scale = c("risk", "survival"),
        B = 200L,
        seed = 1L,
        conf_level = 0.95,
        stratified = FALSE,
        min_success_fraction = 0.80,
        keep_replicates = FALSE,
        verbose = TRUE
) {
    if (!exists("mims_tem_point", mode = "function")) {
        stop("Source mims_tem_refactor_v1.R before calling this function.")
    }
    if (!exists(".bootstrap_mims_wald", mode = "function")) {
        stop("Source mims_bootstrap_wald_v1.R before calling this function.")
    }
    
    scale <- match.arg(scale)
    
    metrics <- c(
        "DE", "IE", "TE",
        "DE_RMTL", "IE_RMTL", "TE_RMTL"
    )
    
    .bootstrap_mims_wald(
        dat = dat,
        estimator = mims_tem_point,
        estimator_args = list(
            basec = basec,
            a = a,
            astar = astar,
            bh_method = bh_method,
            EMint = EMint,
            scale = scale,
            verbose = FALSE,
            keep_fit_objects = FALSE
        ),
        tau = tau,
        metrics = metrics,
        ratio_specs = .mims_boot_ratio_specs_tem_ge(scale),
        B = B,
        seed = seed,
        conf_level = conf_level,
        stratified = stratified,
        exposure = "A",
        reset_id = TRUE,
        min_success_fraction = min_success_fraction,
        keep_replicates = keep_replicates,
        verbose = verbose
    )
}
