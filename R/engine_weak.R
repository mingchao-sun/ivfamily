# Weak-instrument diagnostics: Cragg-Donald Wald F, Kleibergen-Paap rk Wald F, Sanderson-Windmeijer conditional F
# and the Montiel Olea-Pflueger effective F. Everything works on the matrices built by iv_matrices().

# Partialled-out instruments and endogenous regressors (the included controls are removed).
iv_partial <- function(res) {
    mm <- res$mm
    pw <- function(M) if (ncol(mm$W)) M - iv_qr_fit(mm$W, M) else M
    list(Zt = pw(mm$Z1), Xt = pw(mm$X1))
}

# Variance matrix of vec(Pi-hat), the first-stage coefficients of the excluded instruments for all endogenous
# regressors (column-stacked, endogenous regressor by endogenous regressor). The small-sample factors of iv_vcov()
# are always applied (as ivreg2 does for its F statistics).
iv_pi_var <- function(pp, res, vce = res$vce) {
    pt <- iv_partial(res)
    Zt <- pt$Zt; Xt <- pt$Xt
    N <- nrow(Zt); q <- ncol(Zt); p <- ncol(Xt)
    A <- crossprod(Zt)
    Ai <- solve(A)
    Pi <- Ai %*% crossprod(Zt, Xt)
    E <- Xt - Zt %*% Pi
    K <- ncol(res$mm$Z) + res$extra_df
    Smat <- function(g) {              # sum over groups of (Z~ e_j)(Z~ e_k)' for every pair of endogenous regressors
        out <- matrix(0, p * q, p * q)
        gs <- lapply(seq_len(p), function(j) { M <- Zt * E[, j]; if (is.null(g)) M else rowsum(M, g) })
        for (j in seq_len(p)) for (k in seq_len(p))
            out[((j - 1) * q + 1):(j * q), ((k - 1) * q + 1):(k * q)] <- crossprod(gs[[j]], gs[[k]])
        out
    }
    bread <- kronecker(diag(p), Ai)
    if (vce == "classical") {
        Sig <- crossprod(E) / (N - K)
        V <- kronecker(Sig, Ai)
    } else if (vce == "robust") {
        V <- bread %*% Smat(NULL) %*% bread * N / (N - K)
    } else {
        cl <- pp$cl
        g1 <- cl[[1]]; G1 <- length(unique(g1))
        V <- bread %*% Smat(g1) %*% bread
        G <- G1
        if (vce == "cluster2" && length(cl) > 1) {
            V <- V + bread %*% Smat(cl[[2]]) %*% bread - bread %*% Smat(interaction(cl[[1]], cl[[2]], drop = TRUE)) %*% bread
            G <- min(G1, length(unique(cl[[2]])))
        }
        V <- V * (N - 1) / (N - K) * G / (G - 1)
    }
    list(Pi = Pi, V = (V + t(V)) / 2, E = E, N = N, q = q, p = p, K = K)
}

# Cragg-Donald Wald F: smallest eigenvalue of the (iid) first-stage F matrix. Equals the first-stage F when there is one regressor.
iv_cragg_donald <- function(res) {
    pt <- iv_partial(res)
    Zt <- pt$Zt; Xt <- pt$Xt
    N <- nrow(Zt); q <- ncol(Zt)
    K <- ncol(res$mm$Z) + res$extra_df
    Xh <- iv_qr_fit(Zt, Xt)
    Sres <- crossprod(Xt - Xh) / (N - K)
    Ma <- crossprod(Xh) / q
    min(Re(eigen(solve(Sres, Ma), only.values = TRUE)$values))
}

# Symmetric square root of a positive semi-definite matrix.
iv_sqrtm <- function(A) {
    e <- eigen((A + t(A)) / 2)
    e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow(A)) %*% t(e$vectors)
}

# Kleibergen-Paap rk Wald statistic for the rank of Pi (H0: rank = p - 1, i.e. underidentified), computed as in the
# Stata command ranktest (which ivreg2 uses): Pi is standardised with the Cholesky factors of Z'Z and X'X, decomposed by
# SVD, and the statistic is built from the blocks of the singular vectors. The F version follows ivreg2:
# F = chi2 * (N - L) / N / q with heteroskedasticity-robust errors, and chi2 * (N - L) * (G - 1) / (G * (N - 1)) / q with clusters.
iv_kleibergen_paap <- function(pp, res, vce = res$vce) {
    pt <- iv_partial(res)
    Zt <- pt$Zt; Yt <- pt$Xt
    N <- nrow(Zt); q <- ncol(Zt); p <- ncol(Yt)
    K <- ncol(res$mm$Z) + res$extra_df
    zhzh <- crossprod(Zt)
    Pi <- solve(zhzh, crossprod(Zt, Yt))
    E <- Yt - Zt %*% Pi
    Uz <- chol(zhzh); Uy <- chol(crossprod(Yt))
    that <- Uz %*% Pi %*% solve(Uy)
    sv <- svd(that, nu = q, nv = p)
    U <- sv$u; V <- sv$v
    u12 <- U[seq_len(p - 1), p:q, drop = FALSE]; u22 <- U[p:q, p:q, drop = FALSE]
    v12 <- V[seq_len(p - 1), p, drop = FALSE];   v22 <- V[p, p, drop = FALSE]
    aq <- rbind(u12, u22) %*% solve(u22) %*% iv_sqrtm(u22 %*% t(u22))
    bq <- iv_sqrtm(v22 %*% t(v22)) %*% solve(t(v22)) %*% t(rbind(v12, v22))
    L <- kronecker(bq, t(aq))
    lab <- L %*% as.vector(that)
    gi <- matrix(0, N, p * q)
    for (j in seq_len(p)) gi[, ((j - 1) * q + 1):(j * q)] <- Zt * E[, j]
    Sm <- function(g) if (is.null(g)) crossprod(gi) / N else crossprod(rowsum(gi, g)) / N
    cl <- pp$cl
    G <- NULL
    if (vce %in% c("cluster", "cluster2") && length(cl)) {
        S <- Sm(cl[[1]]); G <- length(unique(cl[[1]]))
        if (vce == "cluster2" && length(cl) > 1) {
            S <- S + Sm(cl[[2]]) - Sm(interaction(cl[[1]], cl[[2]], drop = TRUE))
            G <- min(G, length(unique(cl[[2]])))
        }
    } else S <- Sm(NULL)
    Kr <- kronecker(t(solve(Uy)), t(solve(Uz)))
    vhat <- Kr %*% S %*% t(Kr) * N
    chi2 <- drop(t(lab) %*% solve(L %*% vhat %*% t(L), lab))
    F <- if (is.null(G)) chi2 * (N - K) / N / q else chi2 * (N - K) * (G - 1) / (G * (N - 1)) / q
    list(chi2 = chi2, df = q - p + 1, F = F, q = q, p = p)
}

# Sanderson-Windmeijer conditional F for endogenous regressor j: first-stage F of the excluded instruments for x_j
# purged of the other endogenous regressors with their 2SLS coefficients.
iv_sw_f <- function(pp, res, j, vce = res$vce) {
    pt <- iv_partial(res)
    Zt <- pt$Zt; Xt <- pt$Xt
    N <- nrow(Zt); q <- ncol(Zt); p <- ncol(Xt)
    if (p < 2) return(NULL)
    xj <- Xt[, j]; Xo <- Xt[, -j, drop = FALSE]
    Xoh <- iv_qr_fit(Zt, Xo)
    gam <- solve(crossprod(Xoh, Xo), crossprod(Xoh, xj))
    xs <- drop(xj - Xo %*% gam)
    K <- ncol(res$mm$Z) + res$extra_df
    # OLS of xs on the partialled instruments (controls already removed)
    A <- crossprod(Zt); Ai <- solve(A)
    b <- Ai %*% crossprod(Zt, xs); e <- drop(xs - Zt %*% b)
    Sg <- function(g) { M <- Zt * e; if (is.null(g)) crossprod(M) else crossprod(rowsum(M, g)) }
    if (vce == "classical") V <- sum(e^2) / (N - K) * Ai
    else if (vce == "robust") V <- Ai %*% Sg(NULL) %*% Ai * N / (N - K)
    else {
        cl <- pp$cl; G <- length(unique(cl[[1]]))
        V <- Ai %*% Sg(cl[[1]]) %*% Ai
        if (vce == "cluster2" && length(cl) > 1) {
            V <- V + Ai %*% Sg(cl[[2]]) %*% Ai - Ai %*% Sg(interaction(cl[[1]], cl[[2]], drop = TRUE)) %*% Ai
            G <- min(G, length(unique(cl[[2]])))
        }
        V <- V * (N - 1) / (N - K) * G / (G - 1)
    }
    wald <- drop(t(b) %*% solve(V, as.vector(b)))
    df1 <- q - (p - 1)
    list(chi2 = wald, df1 = df1, F = wald / df1)
}

# Montiel Olea-Pflueger effective F for one endogenous regressor, with the chosen heteroskedasticity or cluster-robust variance.
iv_effective_f <- function(pp, res, vce = res$vce) {
    if (length(res$endog) != 1) return(NULL)
    pv <- iv_pi_var(pp, res, if (vce == "classical") "robust" else vce)
    Qzz <- crossprod(iv_partial(res)$Zt) / pv$N
    VQ <- pv$V %*% Qzz
    F <- drop(t(pv$Pi) %*% Qzz %*% pv$Pi) / sum(diag(VQ))
    # Simplified critical values (Montiel Olea and Pflueger 2013): a noncentral chi-square approximation that
    # holds for 2SLS and LIML; thresholds tau = 5, 10, 20 and 30 percent give x = 20, 10, 5 and 3.33 as in weakivtest.
    ev <- sort(Re(eigen(VQ, only.values = TRUE)$values) / sum(Re(eigen(VQ, only.values = TRUE)$values)))
    crit <- vapply(c(`5%` = 20, `10%` = 10, `20%` = 5, `30%` = 3.33), function(x) {
        EK <- 2 * (1 + 2 * x) / (2 * sum(ev^2) + 4 * x * max(ev))
        stats::qchisq(0.95, df = EK, ncp = EK * x) / EK
    }, numeric(1))
    list(F = F, crit = crit)
}

# Reduced form: the outcome regressed on the excluded instruments and controls, with a joint Wald F for the excluded instruments.
iv_reduced_form <- function(pp, res, vce = res$vce, small = res$small) {
    mm <- res$mm
    fm <- list(y = mm$y, X = mm$Z, Z = mm$Z, X1 = mm$Z1, W = mm$W)
    f <- iv_kclass(fm, 0)
    vc <- iv_vcov(f, fm, pp, vce, small, res$extra_df)
    q <- ncol(mm$Z1)
    b <- f$beta[seq_len(q)]
    Vq <- vc$V[seq_len(q), seq_len(q), drop = FALSE]
    se <- sqrt(diag(Vq))
    st <- b / se
    p <- if (is.finite(vc$df)) 2 * stats::pt(-abs(st), vc$df) else 2 * stats::pnorm(-abs(st))
    wald <- drop(t(b) %*% solve(Vq, b))
    list(term = colnames(mm$Z1), estimate = unname(b), se = unname(se), p = unname(p),
         F = wald / q, df1 = q, df2 = vc$df,
         pF = if (is.finite(vc$df)) stats::pf(wald / q, q, vc$df, lower.tail = FALSE) else stats::pchisq(wald, q, lower.tail = FALSE))
}
