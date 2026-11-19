# Shift-share (Bartik) instrumental variables: estimation and the exposure-robust standard errors of
# Adao, Kolesar and Morales (2019, QJE). Base R only. Checked against the R package ShiftShareSE and Stata ivreg_ss.
#
# Region l, sector k: exposure shares s_lk (matrix W, regions x sectors), sector-level shocks g_k, instrument z_l = sum_k s_lk g_k.
# The usual standard errors treat regions as independent; regions with similar shares share the same shocks, so the
# residuals are correlated. AKM variance is built from the sector-level sums of the share-weighted residuals.

# --- data in long form: one row per region and sector ------------------------------------------------------------
# data: data frame. Region-level variables (outcome, endogenous regressor, controls, weights, region cluster) must be
# constant within a region; the shock and the sector cluster constant within a sector.
iv_ss_prepare <- function(data, region, sector, share, shock, outcome, endog, controls = NULL, weights = NULL,
                          region_cluster = NULL, sector_cluster = NULL, sum_shares = FALSE) {
    vars <- unique(c(region, sector, share, shock, outcome, endog, controls, weights, region_cluster, sector_cluster))
    df <- data[, vars, drop = FALSE]
    num <- function(x, nm) {
        if (is.factor(x)) stop(sprintf("The variable '%s' must be numeric.", nm), call. = FALSE)
        v <- as.numeric(x)
        if (any(is.infinite(v))) stop(sprintf("The variable '%s' has infinite values.", nm), call. = FALSE)
        v
    }
    for (v in c(share, shock, outcome, endog, weights)) df[[v]] <- num(df[[v]], v)
    for (v in controls) if (!is.factor(df[[v]])) df[[v]] <- num(df[[v]], v)
    n_all <- nrow(df)
    # rows without region, sector, share or shock are dropped; regions with a missing region-level value are dropped whole
    ok <- stats::complete.cases(df[, c(region, sector, share, shock)])
    df <- df[ok, , drop = FALSE]
    rl_vars <- c(outcome, endog, controls, weights, region_cluster)
    bad_reg <- unique(df[[region]][!stats::complete.cases(df[, rl_vars, drop = FALSE])])
    if (length(bad_reg)) df <- df[!(df[[region]] %in% bad_reg), , drop = FALSE]
    if (nrow(df) < 10) stop("Too few complete observations.", call. = FALSE)
    if (anyDuplicated(df[, c(region, sector)])) stop("Each region may appear only once per sector.", call. = FALSE)
    if (any(df[[share]] < 0)) stop("Exposure shares must not be negative.", call. = FALSE)
    rid <- as.character(df[[region]]); sid <- as.character(df[[sector]])
    regs <- sort(unique(rid)); secs <- sort(unique(sid))
    L <- length(regs); K <- length(secs)
    if (L < 10) stop("Fewer than 10 regions.", call. = FALSE)
    if (K < 2) stop("At least two sectors are needed.", call. = FALSE)
    ri <- match(rid, regs); si <- match(sid, secs)
    const_by <- function(v, idx, n, what) {
        x <- df[[v]]
        if (is.factor(x)) x <- as.character(x)
        first <- x[match(seq_len(n), idx)]
        if (any(!identical(x, first[idx]) & (if (is.numeric(x)) abs(x - first[idx]) > 1e-9 * (1 + abs(first[idx])) else x != first[idx]), na.rm = TRUE))
            stop(sprintf("The variable '%s' varies within a %s. Region-level variables must repeat the same value on every row of a region; sector-level variables on every row of a sector.", v, what), call. = FALSE)
        first
    }
    reg_level <- function(v) const_by(v, ri, L, "region")
    sec_level <- function(v) const_by(v, si, K, "sector")
    y <- as.numeric(reg_level(outcome)); d <- as.numeric(reg_level(endog))
    g <- as.numeric(sec_level(shock))
    W <- matrix(0, L, K); W[cbind(ri, si)] <- df[[share]]
    z <- as.numeric(W %*% g)
    w <- if (!is.null(weights)) { ww <- as.numeric(reg_level(weights)); if (any(ww <= 0)) stop("Weights must be positive.", call. = FALSE); ww } else NULL
    Z <- matrix(1, L, 1, dimnames = list(NULL, "(Intercept)"))
    if (length(controls)) {
        cd <- data.frame(lapply(controls, function(v) reg_level(v)), stringsAsFactors = TRUE); names(cd) <- controls
        for (v in controls) if (is.character(cd[[v]])) cd[[v]] <- factor(cd[[v]])
        mm <- stats::model.matrix(stats::reformulate(paste0("`", controls, "`")), cd)
        Z <- mm
    }
    ssum <- rowSums(W)
    if (isTRUE(sum_shares) && stats::sd(ssum) > 1e-12) Z <- cbind(Z, `Sum of shares` = ssum)
    rc <- if (!is.null(region_cluster)) as.character(reg_level(region_cluster)) else NULL
    sc <- if (!is.null(sector_cluster)) as.character(sec_level(sector_cluster)) else NULL
    list(y = y, d = d, z = z, W = W, Z = Z, w = w, g = g, region_cvar = rc, sector_cvar = sc, L = L, K = K,
         regions = regs, sectors = secs, share_sum = ssum, n_dropped_rows = n_all - nrow(df), n_dropped_regions = length(bad_reg),
         sum_shares_added = isTRUE(sum_shares) && stats::sd(ssum) > 1e-12)
}

# Residual-maker with optional weights.
iv_ss_resid <- function(y, X, w = NULL) if (is.null(w)) stats::lm.fit(X, y)$residuals else stats::lm.wfit(X, y, w)$residuals

# Common AKM / AKM0 / robust inference. beta: estimate; resid: residuals of the estimating equation; ddX: the instrument
# (or regressor) after removing the controls; RX: the "denominator"; ddY / ddW: the pieces for the null-imposed variant.
iv_ss_infer <- function(beta, resid, ddX, ddY1, ddY2, W, w, RX, beta0, alpha, region_cvar, sector_cvar, hom_se, ehw_fun, clu_fun, iv = TRUE) {
    wgt <- if (is.null(w)) 1 else w
    hX <- if (is.null(w)) stats::lm.fit(W, ddX)$coefficients else stats::lm.wfit(W, ddX, w)$coefficients
    hX[is.na(hX)] <- 0
    cR <- hX * drop(crossprod(wgt * resid, W))
    cR0 <- hX * drop(crossprod(wgt * (ddY1 - ddY2 * beta0), W))
    cW <- hX * drop(crossprod(wgt * ddY2, W))
    if (!is.null(sector_cvar)) {
        cR <- tapply(cR, factor(sector_cvar), sum); cR0 <- tapply(cR0, factor(sector_cvar), sum); cW <- tapply(cW, factor(sector_cvar), sum)
    }
    se_akm <- sqrt(sum(cR^2)) / abs(RX)
    cv <- stats::qnorm(1 - alpha / 2)
    se0 <- sqrt(sum(cR0^2)) / abs(RX)
    Q <- RX^2 / cv^2 - sum(cW^2)
    Q2 <- sum(cR * cW) / Q
    mid <- beta - Q2
    dis <- Q2^2 + sum(cR^2) / Q
    if (Q > 0) { cir <- mid + sqrt(dis); cil <- mid - sqrt(dis); se_akm0 <- sqrt(dis) / cv }
    else if (dis > 0) { cir <- mid - sqrt(dis); cil <- mid + sqrt(dis); se_akm0 <- Inf }
    else { cir <- Inf; cil <- -Inf; se_akm0 <- Inf }
    list(akm = se_akm, akm0 = se_akm0, se0_akm0 = se0, akm0_ci = c(cil, cir))
}

# IV with a shift-share instrument. y1: outcome, y2: endogenous regressor, z: instrument, W: shares, Z: controls (with intercept).
iv_ss_ivreg <- function(y1, y2, z, W, Z, w = NULL, beta0 = 0, alpha = 0.05, region_cvar = NULL, sector_cvar = NULL) {
    ff <- function(y, X) if (is.null(w)) stats::lm.fit(X, y) else stats::lm.wfit(X, y, w)
    mm <- cbind(z, Z)
    ddX <- iv_ss_resid(z, Z, w); ddY1 <- iv_ss_resid(y1, Z, w); ddY2 <- iv_ss_resid(y2, Z, w)
    r1 <- ff(y1, mm); r2 <- ff(y2, mm)
    beta <- unname(r1$coefficients[1] / r2$coefficients[1])
    resid <- r1$residuals - r2$residuals * beta
    wgt <- if (is.null(w)) 1 else w
    RX <- sum(wgt * ddY2 * ddX)
    n <- NROW(Z)
    se_h <- unname(sqrt(sum(wgt * resid^2) / n) / (sqrt(sum(wgt * ddX^2)) * abs(r2$coefficients[1])))
    u <- wgt * resid * ddX
    se_r <- sqrt(drop(crossprod(u)) / RX^2)
    se_c <- if (!is.null(region_cvar)) sqrt(drop(crossprod(tapply(u, factor(region_cvar), sum))) / RX^2) else NA_real_
    inf <- iv_ss_infer(beta, resid, ddX, ddY1, ddY2, W, w, RX, beta0, alpha, region_cvar, sector_cvar)
    cv <- stats::qnorm(1 - alpha / 2)
    se <- c(Homoskedastic = se_h, `Robust (EHW)` = se_r, `Region cluster` = se_c, AKM = inf$akm, `AKM0` = inf$akm0)
    p <- 2 * (1 - stats::pnorm(abs(beta - beta0) / c(se[1:4], inf$se0_akm0)))
    lo <- c(beta - cv * se[1:4], inf$akm0_ci[1]); hi <- c(beta + cv * se[1:4], inf$akm0_ci[2])
    list(beta = beta, se = se, p = p, lo = lo, hi = hi, first = unname(r2$coefficients[1]), reduced = unname(r1$coefficients[1]))
}

# Regression of y on the shift-share regressor x (first stage: y = endogenous regressor; reduced form: y = outcome).
iv_ss_reg <- function(y, x, W, Z, w = NULL, beta0 = 0, alpha = 0.05, region_cvar = NULL, sector_cvar = NULL) {
    mm <- cbind(x, Z)
    r <- if (is.null(w)) stats::lm.fit(mm, y) else stats::lm.wfit(mm, y, w)
    ddX <- iv_ss_resid(x, Z, w); ddY <- iv_ss_resid(y, Z, w)
    wgt <- if (is.null(w)) 1 else w
    beta <- unname(r$coefficients[1]); n <- NROW(mm); p <- r$rank
    RX <- sum(wgt * ddX^2)
    u <- wgt * r$residuals * ddX
    se_h <- sqrt(sum(wgt * r$residuals^2) / r$df.residual / RX)
    se_r <- sqrt((n / (n - p)) * drop(crossprod(u))) / RX
    se_c <- if (!is.null(region_cvar)) { nc <- length(unique(region_cvar)); sqrt((nc / (nc - 1)) * (n - 1) / (n - p) * drop(crossprod(tapply(u, factor(region_cvar), sum)))) / RX } else NA_real_
    inf <- iv_ss_infer(beta, r$residuals, ddX, ddY, ddX, W, w, RX, beta0, alpha, region_cvar, sector_cvar)
    cv <- stats::qnorm(1 - alpha / 2)
    se <- c(Homoskedastic = se_h, `Robust (EHW)` = se_r, `Region cluster` = se_c, AKM = inf$akm, `AKM0` = inf$akm0)
    p <- 2 * (1 - stats::pnorm(abs(beta - beta0) / c(se[1:4], inf$se0_akm0)))
    list(beta = beta, se = se, p = p, lo = c(beta - cv * se[1:4], inf$akm0_ci[1]), hi = c(beta + cv * se[1:4], inf$akm0_ci[2]),
         partial_r2 = sum(wgt * ddX * ddY)^2 / (sum(wgt * ddX^2) * sum(wgt * ddY^2)))
}

# Effective number of sectors: inverse Herfindahl of the (weighted) exposure weights (Borusyak, Hull and Jaravel, 2022).
iv_ss_effective_sectors <- function(W, w = NULL) {
    ww <- if (is.null(w)) rep(1, nrow(W)) else w
    s <- colSums(ww * W)
    s <- s / sum(s)
    1 / sum(s^2)
}
