# Weak-instrument-robust inference for one endogenous regressor: Anderson-Rubin (AR), Kleibergen (K) and the
# conditional likelihood ratio (CLR) tests, with homoskedastic, heteroskedasticity-robust and cluster-robust
# moment covariances, and confidence sets obtained by inverting the tests.
# After partialling out the controls, the moments are g(b) = Z'(y - b x) / N.

iv_ar_setup <- function(pp, res, vce = res$vce) {
    if (length(res$endog) != 1)
        stop("The weak-instrument-robust tests are available for one endogenous regressor.", call. = FALSE)
    mm <- res$mm
    pw <- function(M) if (ncol(mm$W)) M - iv_qr_fit(mm$W, M) else M
    yt <- drop(pw(matrix(mm$y, ncol = 1))); xt <- drop(pw(mm$X1)); Zt <- pw(mm$Z1)
    N <- length(yt); q <- ncol(Zt)
    Q <- crossprod(Zt) / N
    A <- Zt * yt; B <- Zt * xt
    cl <- NULL
    if (vce == "classical") {
        # Moreira (2003) statistics: orthonormal instruments and the reduced-form covariance matrix
        Y <- cbind(yt, xt)
        Zo <- Zt %*% solve(chol(crossprod(Zt)))
        Om <- crossprod(Y - iv_qr_fit(Zt, Y)) / N
        cl <- list(ZoY = crossprod(Zo, Y), Om = Om, Omi = solve(Om), uu = c(sum(yt * yt), sum(yt * xt), sum(xt * xt)) / N)
        S <- NULL
    } else {
        grp <- if (vce == "robust" || !length(pp$cl)) NULL else pp$cl
        sm <- function(X, Y) {
            if (is.null(grp)) return(crossprod(X, Y) / N)
            out <- crossprod(rowsum(X, grp[[1]]), rowsum(Y, grp[[1]])) / N
            if (vce == "cluster2" && length(grp) > 1) {
                ii <- interaction(grp[[1]], grp[[2]], drop = TRUE)
                out <- out + crossprod(rowsum(X, grp[[2]]), rowsum(Y, grp[[2]])) / N - crossprod(rowsum(X, ii), rowsum(Y, ii)) / N
            }
            out
        }
        S <- list(aa = sm(A, A), ab = sm(A, B), bb = sm(B, B))
    }
    list(N = N, q = q, abar = colMeans(A), bbar = colMeans(B), S = S, cl = cl, yt = yt, xt = xt, Zt = Zt, vce = vce)
}

# Moreira conditional p-value of the CLR statistic given the rank statistic rk (df = q instruments).
iv_clr_pvalue <- function(lr, rk, q) {
    if (q == 1) return(stats::pchisq(lr, 1, lower.tail = FALSE))
    if (!is.finite(lr) || lr <= 0) return(1)
    f <- function(z) {
        thr <- lr + rk - z^2 * (rk / lr + 1)
        stats::dnorm(z) * ifelse(thr <= 0, 1, stats::pchisq(thr, q - 1, lower.tail = FALSE))
    }
    2 * stats::integrate(f, 0, 12, subdivisions = 400L, rel.tol = 1e-10)$value
}

# AR, K and CLR statistics and p-values at the hypothesised coefficient b0.
# Classical errors follow weakiv: AR and K use the moments scaled by u'u / N, the CLR uses the reduced-form covariance matrix
# (Moreira 2003). Robust and cluster errors use the moment covariances of Z u and Z x (Kleibergen 2005).
iv_weak_tests <- function(st, b0) {
    N <- st$N; q <- st$q
    if (st$vce == "classical") {
        cl <- st$cl
        u2 <- cl$uu[1] - 2 * b0 * cl$uu[2] + b0^2 * cl$uu[3]
        Sone <- (st$cl$ZoY %*% c(1, -b0)) / sqrt(u2)
        bv <- c(1, -b0); av <- c(b0, 1)
        Srf <- (cl$ZoY %*% bv) / sqrt(drop(t(bv) %*% cl$Om %*% bv))
        Tv <- (cl$ZoY %*% cl$Omi %*% av) / sqrt(drop(t(av) %*% cl$Omi %*% av))
        AR <- sum(Sone^2); Kst <- sum(Sone * Tv)^2 / sum(Tv^2); rk <- sum(Tv^2)
        QS <- sum(Srf^2); QST <- sum(Srf * Tv)
        LR <- 0.5 * (QS - rk + sqrt(max((QS + rk)^2 - 4 * (QS * rk - QST^2), 0)))
    } else {
        S <- st$S
        g <- st$abar - b0 * st$bbar
        Sg <- S$aa - b0 * (S$ab + t(S$ab)) + b0^2 * S$bb               # Var of the moments Z u
        Vdg <- t(S$ab) - b0 * S$bb                                      # Cov(Z x, Z u)
        Vdd <- S$bb
        Sinv <- solve(Sg)
        Dt <- st$bbar - Vdg %*% Sinv %*% g                              # part of Z'x / N uncorrelated with g
        AR <- drop(N * t(g) %*% Sinv %*% g)
        Kst <- drop(N * (t(g) %*% Sinv %*% Dt)^2 / (t(Dt) %*% Sinv %*% Dt))
        rk <- drop(N * t(Dt) %*% solve(Vdd - Vdg %*% Sinv %*% t(Vdg), Dt))
        LR <- 0.5 * (AR - rk + sqrt(max((AR + rk)^2 - 4 * rk * (AR - Kst), 0)))
    }
    list(b0 = b0, AR = AR, K = Kst, rk = rk, CLR = LR, df = q,
         p_ar = stats::pchisq(AR, q, lower.tail = FALSE), p_k = stats::pchisq(Kst, 1, lower.tail = FALSE),
         p_clr = iv_clr_pvalue(LR, rk, q))
}

# Confidence set {b : p(b) >= alpha} on a grid, with the end points refined by root finding. Returns a data frame of
# intervals; lo = -Inf or hi = Inf when the set reaches the edge of the grid.
iv_cset <- function(pfun, grid, alpha = 0.05) {
    p <- vapply(grid, pfun, numeric(1))
    inside <- p >= alpha
    if (!any(inside)) return(data.frame(lo = numeric(0), hi = numeric(0)))
    runs <- rle(inside)
    ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1
    out <- list()
    refine <- function(i0, i1) stats::uniroot(function(b) pfun(b) - alpha, c(grid[i0], grid[i1]), tol = 1e-10)$root
    for (r in which(runs$values)) {
        s <- starts[r]; e <- ends[r]
        lo <- if (s == 1) -Inf else refine(s - 1, s)
        hi <- if (e == length(grid)) Inf else refine(e, e + 1)
        out[[length(out) + 1]] <- data.frame(lo = lo, hi = hi)
    }
    do.call(rbind, out)
}

iv_cset_text <- function(cs, d = 3) {
    if (nrow(cs) == 0) return("empty")
    f <- function(x) if (is.infinite(x)) (if (x < 0) "-Inf" else "Inf") else formatC(x, format = "f", digits = d, drop0trailing = TRUE)
    paste(sprintf("[%s, %s]", vapply(cs$lo, f, ""), vapply(cs$hi, f, "")), collapse = " U ")
}
