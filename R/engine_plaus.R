# Sensitivity to violations of the exclusion restriction: plausibly exogenous instruments (Conley, Hansen and Rossi 2012).
# The structural equation is y = X b + Z g + u, where g is the direct effect of the excluded instruments Z on the outcome
# (zero under exclusion). With g known, 2SLS on y - Z g gives b(g) = b - A g, A = (X'P_Z X)^-1 X'Z with X and Z purified of the controls.

# Returns A (p x q): row j says how much the coefficient of endogenous regressor j moves per unit of direct effect.
iv_plaus_A <- function(res) {
    mm <- res$mm
    pw <- function(M) if (ncol(mm$W)) M - iv_qr_fit(mm$W, M) else M
    Xt <- pw(mm$X1); Zt <- pw(mm$Z1)
    Xh <- iv_qr_fit(Zt, Xt)
    solve(crossprod(Xh, Xt), crossprod(Xt, Zt))
}

# Union of confidence intervals for a box of direct effects g_k in [gmin_k, gmax_k] (all instruments).
iv_plaus_uci <- function(res, gmin, gmax, level = 95) {
    A <- iv_plaus_A(res)
    q <- ncol(A); gmin <- rep_len(gmin, q); gmax <- rep_len(gmax, q)
    z <- stats::qnorm(1 - (1 - level / 100) / 2)
    idx <- match(res$endog, names(res$beta))
    b <- res$beta[idx]; se <- res$se[idx]
    hi_shift <- vapply(seq_along(idx), function(j) sum(pmax(A[j, ] * gmin, A[j, ] * gmax)), 0)   # largest of A g
    lo_shift <- vapply(seq_along(idx), function(j) sum(pmin(A[j, ] * gmin, A[j, ] * gmax)), 0)
    data.frame(term = res$endog, estimate = unname(b), se = unname(se), lo = unname(b - z * se - hi_shift), hi = unname(b + z * se - lo_shift),
               stringsAsFactors = FALSE)
}

# Smallest delta for which the union of intervals with all |g_k| <= delta contains zero.
iv_plaus_tipping <- function(res, level = 95) {
    A <- iv_plaus_A(res)
    z <- stats::qnorm(1 - (1 - level / 100) / 2)
    idx <- match(res$endog, names(res$beta))
    b <- res$beta[idx]; se <- res$se[idx]
    pmax(0, (abs(b) - z * se) / rowSums(abs(A)))
}

# Local-to-zero approximation: g ~ N(mu, Omega); the estimate is shifted by A mu and its variance grows by A Omega A'.
iv_plaus_ltz <- function(res, mu, omega, level = 95) {
    A <- iv_plaus_A(res)
    q <- ncol(A); mu <- rep_len(mu, q)
    Om <- if (is.matrix(omega)) omega else diag(rep_len(omega, q), q)
    z <- stats::qnorm(1 - (1 - level / 100) / 2)
    idx <- match(res$endog, names(res$beta))
    b <- res$beta[idx] - drop(A %*% mu)
    v <- diag(res$V)[idx] + diag(A %*% Om %*% t(A))
    data.frame(term = res$endog, estimate = unname(b), se = sqrt(v), lo = unname(b - z * sqrt(v)), hi = unname(b + z * sqrt(v)), stringsAsFactors = FALSE)
}

# Union of confidence intervals evaluated as plausexog does it: 2SLS and its standard errors are recomputed at every corner of
# the box of direct effects (with robust or cluster errors the standard error changes with g). Falls back to the linear
# shift with a fixed standard error when there are more than 8 instruments (2^q corners).
iv_plaus_uci_exact <- function(pp, res, gmin, gmax, level = 95) {
    q <- ncol(pp$Z1)
    if (q > 8) return(iv_plaus_uci(res, gmin, gmax, level))
    gmin <- rep_len(gmin, q); gmax <- rep_len(gmax, q)
    z <- stats::qnorm(1 - (1 - level / 100) / 2)
    corners <- as.matrix(expand.grid(lapply(seq_len(q), function(k) unique(c(gmin[k], gmax[k])))))
    lo <- rep(Inf, length(res$endog)); hi <- rep(-Inf, length(res$endog))
    for (i in seq_len(nrow(corners))) {
        g <- corners[i, ]
        pp2 <- pp; pp2$y <- pp$y - drop(pp$Z1 %*% g)
        r <- iv_fit(pp2, "2sls", vce = res$vce, small = res$small)
        b <- r$beta[res$endog]; se <- r$se[res$endog]
        lo <- pmin(lo, b - z * se); hi <- pmax(hi, b + z * se)
    }
    est <- res$beta[res$endog]
    data.frame(term = res$endog, estimate = unname(est), se = unname(res$se[res$endog]), lo = lo, hi = hi, stringsAsFactors = FALSE)
}

# Tipping delta with the exact union: the smallest delta whose union interval reaches zero (root of the bound facing zero).
iv_plaus_tipping_exact <- function(pp, res, level = 95) {
    q <- ncol(pp$Z1)
    if (q > 8) return(iv_plaus_tipping(res, level))
    start <- iv_plaus_tipping(res, level)
    vapply(seq_along(res$endog), function(j) {
        est <- res$beta[[res$endog[j]]]
        f <- function(d) {
            u <- iv_plaus_uci_exact(pp, res, rep(-d, q), rep(d, q), level)[j, ]
            if (est >= 0) u$lo else -u$hi
        }
        if (f(0) <= 0) return(0)
        hi <- max(start[j] * 2, 1e-8)
        while (f(hi) > 0 && hi < 1e6) hi <- hi * 2
        stats::uniroot(f, c(0, hi), tol = 1e-10)$root
    }, 0)
}
