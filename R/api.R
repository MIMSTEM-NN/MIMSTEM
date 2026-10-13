# Public API ---------------------------------------------------------------

.mims_validate_public_input <- function(dat, tau, basec, a, astar) {
    if (!is.data.frame(dat)) stop("`dat` must be a data.frame.", call. = FALSE)
    if (nrow(dat) < 2L) stop("`dat` must contain at least two subjects.", call. = FALSE)
    
    tau <- as.numeric(tau)
    if (!length(tau) || any(!is.finite(tau)) || any(tau <= 0)) {
        stop("`tau` must contain positive finite evaluation times.", call. = FALSE)
    }
    
    if (is.null(basec)) basec <- character(0)
    if (!is.character(basec) || anyNA(basec) || anyDuplicated(basec)) {
        stop("`basec` must be a character vector of unique baseline covariate names.", call. = FALSE)
    }
    missing_basec <- setdiff(basec, names(dat))
    if (length(missing_basec)) {
        stop(
            "Requested baseline covariates are missing from `dat`: ",
            paste(missing_basec, collapse = ", "),
            call. = FALSE
        )
    }
    
    if (!"A" %in% names(dat)) stop("`dat` must contain exposure column `A`.", call. = FALSE)
    A_num <- suppressWarnings(as.numeric(as.character(dat$A)))
    if (anyNA(A_num) || !all(A_num %in% c(0, 1))) {
        stop("Exposure `A` must be coded 0/1.", call. = FALSE)
    }
    if (!all(c(0, 1) %in% unique(A_num))) {
        stop("Both exposure groups A=0 and A=1 are required.", call. = FALSE)
    }
    
    if (length(a) != 1L || length(astar) != 1L ||
        anyNA(c(a, astar)) || any(!c(a, astar) %in% c(0, 1))) {
        stop("`a` and `astar` must each be either 0 or 1.", call. = FALSE)
    }
    if (identical(as.numeric(a), as.numeric(astar))) {
        warning("`a` and `astar` are identical; all exposure contrasts should be zero.", call. = FALSE)
    }
    
    invisible(TRUE)
}

.mims_effect_metrics <- function(method) {
    switch(
        method,
        MIMS_TEM = c("DE", "IE", "TE", "DE_RMTL", "IE_RMTL", "TE_RMTL"),
        MIMS_S2TEM_GE = c("DE", "IE", "TE", "DE_RMTL", "IE_RMTL", "TE_RMTL"),
        MIMS_S2TEM_PE = c(
            "TE", "DE", "IE_total", "P_A_T2_Y", "P_A_T1_Y", "P_A_T1_T2_Y",
            "TE_RMTL", "DE_RMTL", "IE_total_RMTL", "P_A_T2_Y_RMTL",
            "P_A_T1_Y_RMTL", "P_A_T1_T2_Y_RMTL"
        ),
        stop("Unknown MIMS method.")
    )
}

.mims_make_effect_table <- function(point, method, boot = NULL) {
    metrics <- .mims_effect_metrics(method)
    miss <- setdiff(c("tau", metrics), names(point))
    if (length(miss)) {
        stop("Internal result is missing: ", paste(miss, collapse = ", "), call. = FALSE)
    }
    
    rows <- do.call(rbind, lapply(metrics, function(metric) {
        data.frame(
            tau = point$tau,
            estimand = metric,
            estimate = as.numeric(point[[metric]]),
            se = NA_real_,
            lower = NA_real_,
            upper = NA_real_,
            stringsAsFactors = FALSE
        )
    }))
    
    if (!is.null(boot)) {
        ints <- boot$intervals
        key_rows <- paste(rows$tau, rows$estimand, sep = "\r")
        key_ints <- paste(ints$tau, ints$metric, sep = "\r")
        ix <- match(key_rows, key_ints)
        ok <- !is.na(ix)
        rows$se[ok] <- ints$se_boot[ix[ok]]
        rows$lower[ok] <- ints$lower[ix[ok]]
        rows$upper[ok] <- ints$upper[ix[ok]]
    }
    
    rows$effect_type <- ifelse(grepl("_RMTL$", rows$estimand), "RMTL", "risk")
    rownames(rows) <- NULL
    rows
}

.mims_make_fit <- function(method, call, point, boot, settings) {
    effects <- .mims_make_effect_table(point, method, boot)
    
    endpoint_names <- switch(
        method,
        MIMS_TEM = intersect(c("tau", "psi00", "psi10", "psi11", "S00", "S10", "S11",
                               "rmtl00", "rmtl10", "rmtl11"), names(point)),
        MIMS_S2TEM_GE = intersect(c("tau", "psi00", "psi10", "psi11", "S00", "S10", "S11",
                                    "rmtl00", "rmtl10", "rmtl11"), names(point)),
        MIMS_S2TEM_PE = intersect(c("tau", "psi0000", "psi1000", "psi1100", "psi1110", "psi1111",
                                    "S0000", "S1000", "S1100", "S1110", "S1111"), names(point))
    )
    
    diagnostics <- list()
    if (!is.null(boot)) diagnostics$bootstrap <- boot$diagnostics
    pd <- attr(point, "path_diagnostics", exact = TRUE)
    if (!is.null(pd)) diagnostics$path <- pd
    controls <- attr(point, "controls", exact = TRUE)
    if (!is.null(controls)) diagnostics$controls <- controls
    
    out <- list(
        call = call,
        method = method,
        point = point,
        endpoints = point[endpoint_names],
        effects = effects,
        bootstrap = boot,
        diagnostics = diagnostics,
        settings = settings
    )
    class(out) <- c(
        switch(
            method,
            MIMS_TEM = "mims_tem",
            MIMS_S2TEM_GE = "mims_s2tem_ge",
            MIMS_S2TEM_PE = "mims_s2tem_pe"
        ),
        "mims_fit"
    )
    out
}

#' MIMS-TEM: single time-to-event mediator
#'
#' Estimate marginal interventional direct, indirect and total effects with a
#' three-state multi-state model for one time-to-event mediator.
#'
#' @param dat Subject-level data frame. Required canonical columns are `A`,
#'   `t_out`, `d_out` (or `yevent`), `t_med`, and `d_med`.
#' @param tau Positive evaluation times.
#' @param basec Character vector of baseline covariate names.
#' @param a,astar Exposure levels defining the contrast.
#' @param bootstrap Logical; if TRUE, compute individual-level nonparametric
#'   bootstrap standard errors and Wald confidence intervals.
#' @param B Number of bootstrap replicates.
#' @param seed Random seed for bootstrap resampling.
#' @param conf_level Confidence level for Wald intervals.
#' @param stratified Logical; use exposure-stratified resampling when TRUE.
#' @param bh_method Baseline Cox tie method, default `"breslow"`.
#' @param EMint Include the exposure-by-mediator-time interaction supported by
#'   the frozen estimator.
#' @param scale Effect display scale, `"risk"` or `"survival"`.
#' @param min_success_fraction Minimum successful bootstrap fraction needed to
#'   report standard errors and intervals.
#' @param keep_replicates Keep the bootstrap replicate array.
#' @param keep_fit_objects Keep fitted multi-state objects in the point result.
#' @param verbose Show bootstrap progress.
#' @return An object inheriting from `mims_fit`.
#' @export
mims_tem <- function(
        dat,
        tau = 1:10,
        basec = character(0),
        a = 1,
        astar = 0,
        bootstrap = TRUE,
        B = 200L,
        seed = 1L,
        conf_level = 0.95,
        stratified = FALSE,
        bh_method = "breslow",
        EMint = FALSE,
        scale = c("risk", "survival"),
        min_success_fraction = 0.80,
        keep_replicates = FALSE,
        keep_fit_objects = FALSE,
        verbose = TRUE
) {
    cl <- match.call()
    scale <- match.arg(scale)
    .mims_validate_public_input(dat, tau, basec, a, astar)
    
    if (isTRUE(bootstrap)) {
        boot <- bootstrap_mims_tem(
            dat = dat, tau = tau, basec = basec, a = a, astar = astar,
            bh_method = bh_method, EMint = EMint, scale = scale,
            B = B, seed = seed, conf_level = conf_level,
            stratified = stratified,
            min_success_fraction = min_success_fraction,
            keep_replicates = keep_replicates,
            verbose = verbose
        )
        point <- boot$point
    } else {
        boot <- NULL
        point <- mims_tem_point(
            dat = dat, tau = tau, basec = basec, a = a, astar = astar,
            bh_method = bh_method, EMint = EMint, scale = scale,
            verbose = verbose, keep_fit_objects = keep_fit_objects
        )
    }
    
    .mims_make_fit(
        method = "MIMS_TEM", call = cl, point = point, boot = boot,
        settings = list(
            tau = sort(unique(as.numeric(tau))), basec = basec, a = a, astar = astar,
            bootstrap = isTRUE(bootstrap), B = if (bootstrap) as.integer(B) else 0L,
            seed = seed, conf_level = conf_level, stratified = stratified,
            bh_method = bh_method, EMint = EMint, scale = scale
        )
    )
}

#' MIMS-S2TEM-GE: overall effects with two sequential time-to-event mediators
#'
#' Estimate total, direct and overall indirect marginal interventional effects
#' using the four-state/five-transition sequential-mediator model.
#'
#' @inheritParams mims_tem
#' @param v_chunk_size Internal chunk size used by the numerical integration.
#' @param max_grid_points Maximum adaptive hazard-grid size.
#' @return An object inheriting from `mims_fit`.
#' @export
mims_s2tem_ge <- function(
        dat,
        tau = 1:10,
        basec = character(0),
        a = 1,
        astar = 0,
        bootstrap = TRUE,
        B = 200L,
        seed = 1L,
        conf_level = 0.95,
        stratified = FALSE,
        bh_method = "breslow",
        EMint = FALSE,
        scale = c("risk", "survival"),
        min_success_fraction = 0.80,
        keep_replicates = FALSE,
        keep_fit_objects = FALSE,
        verbose = TRUE,
        v_chunk_size = 4096L,
        max_grid_points = 200L
) {
    cl <- match.call()
    scale <- match.arg(scale)
    .mims_validate_public_input(dat, tau, basec, a, astar)
    
    if (isTRUE(bootstrap)) {
        boot <- bootstrap_mims_s2tem_ge(
            dat = dat, tau = tau, basec = basec, a = a, astar = astar,
            bh_method = bh_method, EMint = EMint, scale = scale,
            B = B, seed = seed, conf_level = conf_level,
            stratified = stratified,
            min_success_fraction = min_success_fraction,
            keep_replicates = keep_replicates,
            verbose = verbose,
            v_chunk_size = v_chunk_size,
            max_grid_points = max_grid_points
        )
        point <- boot$point
    } else {
        boot <- NULL
        point <- mims_s2tem_ge_point(
            dat = dat, tau = tau, basec = basec, a = a, astar = astar,
            bh_method = bh_method, EMint = EMint, scale = scale,
            verbose = verbose, keep_fit_objects = keep_fit_objects,
            v_chunk_size = v_chunk_size, max_grid_points = max_grid_points
        )
    }
    
    .mims_make_fit(
        method = "MIMS_S2TEM_GE", call = cl, point = point, boot = boot,
        settings = list(
            tau = sort(unique(as.numeric(tau))), basec = basec, a = a, astar = astar,
            bootstrap = isTRUE(bootstrap), B = if (bootstrap) as.integer(B) else 0L,
            seed = seed, conf_level = conf_level, stratified = stratified,
            bh_method = bh_method, EMint = EMint, scale = scale,
            v_chunk_size = v_chunk_size, max_grid_points = max_grid_points
        )
    )
}

#' MIMS-S2TEM-PE: path-specific effects with two sequential time-to-event mediators
#'
#' Preserve the MIMS-S2TEM-GE total/direct/overall-indirect endpoints and
#' decompose the overall indirect effect into three mediator-related
#' path-specific intervention effects under the frozen switch order
#' `(zY, z2, z1, z12)`.
#'
#' @inheritParams mims_s2tem_ge
#' @param primary_rule Mixed-history rule. The primary rule is
#'   `"t2_overrides"`; `"require_both"` is available for sensitivity analysis.
#' @param t1_coupling Currently only `"rank_preserving"` is implemented.
#' @param subject_chunk_size Subject block size for path-kernel calculations.
#' @return An object inheriting from `mims_fit`. The main effect table reports
#'   TE, DE, `IE_total`, `P_A_T2_Y`, `P_A_T1_Y`, and `P_A_T1_T2_Y`
#'   (and their RMTL analogues). `IE_total` is also available in `$point`
#'   for decomposition checks, and its bootstrap SE and Wald confidence
#'   interval are included in the main effect table.
#' @export
mims_s2tem_pe <- function(
        dat,
        tau = 1:10,
        basec = character(0),
        a = 1,
        astar = 0,
        bootstrap = TRUE,
        B = 200L,
        seed = 1L,
        conf_level = 0.95,
        stratified = FALSE,
        bh_method = "breslow",
        EMint = FALSE,
        scale = c("risk", "survival"),
        min_success_fraction = 0.80,
        keep_replicates = FALSE,
        verbose = TRUE,
        v_chunk_size = 4096L,
        max_grid_points = 200L,
        primary_rule = c("t2_overrides", "require_both"),
        t1_coupling = "rank_preserving",
        subject_chunk_size = 512L
) {
    cl <- match.call()
    scale <- match.arg(scale)
    primary_rule <- match.arg(primary_rule)
    .mims_validate_public_input(dat, tau, basec, a, astar)
    
    if (!identical(t1_coupling, "rank_preserving")) {
        stop("Only `t1_coupling = \"rank_preserving\"` is currently implemented.", call. = FALSE)
    }
    
    if (isTRUE(bootstrap)) {
        boot <- bootstrap_mims_s2tem_pe(
            dat = dat, tau = tau, basec = basec, a = a, astar = astar,
            bh_method = bh_method, EMint = EMint, scale = scale,
            B = B, seed = seed, conf_level = conf_level,
            stratified = stratified,
            min_success_fraction = min_success_fraction,
            keep_replicates = keep_replicates,
            verbose = verbose,
            max_grid_points = max_grid_points,
            v_chunk_size = v_chunk_size,
            primary_rule = primary_rule,
            t1_coupling = t1_coupling,
            subject_chunk_size = subject_chunk_size
        )
        point <- boot$point
    } else {
        boot <- NULL
        point <- mims_s2tem_pe_point(
            dat = dat, tau = tau, basec = basec, a = a, astar = astar,
            bh_method = bh_method, EMint = EMint, scale = scale,
            verbose = verbose, keep_base_fit = FALSE,
            max_grid_points = max_grid_points, v_chunk_size = v_chunk_size,
            primary_rule = primary_rule, t1_coupling = t1_coupling,
            check_controls = FALSE, rule_sensitivity = FALSE,
            subject_chunk_size = subject_chunk_size
        )
    }
    
    .mims_make_fit(
        method = "MIMS_S2TEM_PE", call = cl, point = point, boot = boot,
        settings = list(
            tau = sort(unique(as.numeric(tau))), basec = basec, a = a, astar = astar,
            bootstrap = isTRUE(bootstrap), B = if (bootstrap) as.integer(B) else 0L,
            seed = seed, conf_level = conf_level, stratified = stratified,
            bh_method = bh_method, EMint = EMint, scale = scale,
            v_chunk_size = v_chunk_size, max_grid_points = max_grid_points,
            primary_rule = primary_rule, t1_coupling = t1_coupling,
            subject_chunk_size = subject_chunk_size
        )
    )
}
