# Overidentification, orthogonality (C / difference-in-Hansen) and endogeneity tests.
# All tests are differences of GMM criterion values J, so they are non-negative by construction.

# Moment covariance of the instruments Zs implied by residuals u: homoskedastic for classical errors, robust or cluster otherwise.
iv_moment_s <- function(Zs, u, pp, vce) {
    N <- nrow(Zs)
    if (vce == "classical") sum(u^2) / N * crossprod(Zs) / N else iv_moment_cov(Zs, u, pp, vce)
}

# GMM criterion minimised over beta with instruments Zs and weight matrix solve(S).
iv_jcrit <- function(y, X, Zs, S) {
    N <- length(y)
    Wg <- solve(S)
    XZ <- crossprod(X, Zs) / N
    b <- solve(XZ %*% Wg %*% t(XZ), XZ %*% Wg %*% crossprod(Zs, y) / N)
    g <- crossprod(Zs, y - X %*% b) / N
    drop(N * t(g) %*% Wg %*% g)
}

iv_check_clusters <- function(res, ncolZ) {
    if (res$vce %in% c("cluster", "cluster2") && res$G <= ncolZ)
        stop("Too few clusters for this test: it needs more clusters than instruments (including controls).", call. = FALSE)
}

# Difference-in-Hansen (C) test that a subset of the excluded instruments is orthogonal to the error.
iv_orthog_test <- function(pp, res, suspects) {
    mm <- res$mm; y <- mm$y; X <- mm$X
    sus <- match(suspects, colnames(mm$Z1))
    if (anyNA(sus)) stop("The instruments to test must be among the excluded instruments.", call. = FALSE)
    Zf <- mm$Z
    keep <- setdiff(seq_len(ncol(Zf)), sus)
    if (length(keep) < ncol(X))
        stop("Dropping these instruments would leave the model underidentified; test fewer instruments.", call. = FALSE)
    iv_check_clusters(res, ncol(Zf))
    u <- iv_kclass(mm, 1)$u
    S <- iv_moment_s(Zf, u, pp, res$vce)
    Jf <- iv_jcrit(y, X, Zf, S)
    Jr <- iv_jcrit(y, X, Zf[, keep, drop = FALSE], S[keep, keep, drop = FALSE])
    df <- length(sus)
    C <- Jf - Jr
    list(C = C, df = df, p = stats::pchisq(C, df, lower.tail = FALSE), J_full = Jf, J_restricted = Jr)
}

# Endogeneity test: treat the listed endogenous regressors as exogenous (add them to the instruments) and compare.
iv_endog_test <- function(pp, res, tested) {
    mm <- res$mm; y <- mm$y; X <- mm$X
    te <- match(tested, colnames(mm$X1))
    if (anyNA(te)) stop("The regressors to test must be among the endogenous regressors.", call. = FALSE)
    Zs <- mm$Z
    Zb <- cbind(mm$Z1, mm$X1[, te, drop = FALSE], mm$W)
    iv_check_clusters(res, ncol(Zb))
    if (qr(Zb, tol = 1e-10)$rank < ncol(Zb)) stop("The tested regressors are collinear with the instruments and controls.", call. = FALSE)
    ub <- iv_kclass(list(y = y, X = X, Z = Zb), 1)$u
    S <- iv_moment_s(Zb, ub, pp, res$vce)
    small <- seq_len(ncol(Zs))                 # Zs columns are the first ones of Zb only if the order matches; rebuild explicitly
    idx <- match(colnames(Zs), colnames(Zb))
    if (anyNA(idx)) idx <- seq_len(ncol(Zs))
    Jb <- iv_jcrit(y, X, Zb, S)
    Js <- iv_jcrit(y, X, Zs, S[idx, idx, drop = FALSE])
    df <- length(te)
    C <- Jb - Js
    list(C = C, df = df, p = stats::pchisq(C, df, lower.tail = FALSE))
}
