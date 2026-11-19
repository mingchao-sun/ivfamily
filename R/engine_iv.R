# Instrumental-variable estimation engine: OLS, 2SLS, LIML, Fuller, k-class and two-step GMM,
# with classical, robust (HC1), cluster and two-way cluster variance estimators.
# Base R only; nothing here depends on jamovi, so it can be tested from a plain R session.

iv_num <- function(x) if (is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)

# Dummy-coded design matrix (treatment contrasts, no intercept) for a set of control variables.
iv_design <- function(data, vars) {
    if (!length(vars)) return(NULL)
    df <- data[, vars, drop = FALSE]
    for (v in vars) {
        if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
        else if (is.factor(df[[v]])) df[[v]] <- droplevels(df[[v]])
    }
    mm <- stats::model.matrix(~ ., data = df)
    mm[, -1, drop = FALSE]
}

# Weighted within transformation for one or two sets of fixed effects (alternating projections).
iv_demean <- function(M, fes, w = NULL, tol = 1e-10, maxit = 5000) {
    M <- as.matrix(M)
    if (!length(fes)) return(M)
    if (is.null(w)) w <- rep(1, nrow(M))
    gi <- lapply(fes, function(f) as.integer(factor(f)))
    ws <- lapply(gi, function(g) as.numeric(rowsum(w, g)))
    proj <- function(M, g, s) M - (rowsum(M * w, g) / s)[g, , drop = FALSE]
    if (length(fes) == 1) return(proj(M, gi[[1]], ws[[1]]))
    for (it in seq_len(maxit)) {
        M0 <- M
        for (j in seq_along(gi)) M <- proj(M, gi[[j]], ws[[j]])
        if (max(abs(M - M0)) <= tol * max(1, max(abs(M0)))) break
    }
    M
}

# Build the estimation sample from a data frame. Rows with a missing value in any used variable are dropped.
iv_prepare <- function(data, y, endog, instr, exog = NULL, absorb = NULL, cluster = NULL, weights = NULL, extra = NULL) {
    if (length(endog) < 1) stop("Select at least one endogenous regressor.", call. = FALSE)
    if (length(instr) < length(endog))
        stop(sprintf("The model is underidentified: %d endogenous regressor(s) but only %d excluded instrument(s).",
                     length(endog), length(instr)), call. = FALSE)
    if (any(duplicated(c(y, endog, instr, exog))))
        stop("A variable appears in more than one role (outcome, endogenous, instrument, control).", call. = FALSE)
    used <- unique(c(y, endog, instr, exog, absorb, cluster, weights, extra))
    keep <- stats::complete.cases(data[, used, drop = FALSE])
    for (v in c(y, endog, instr, weights, extra)) keep <- keep & is.finite(iv_num(data[[v]]))
    n_drop <- sum(!keep)
    d <- data[keep, , drop = FALSE]
    if (nrow(d) < 10) stop("Fewer than 10 complete observations.", call. = FALSE)
    n <- nrow(d)
    yv <- iv_num(d[[y]])
    X1 <- matrix(vapply(endog, function(v) iv_num(d[[v]]), numeric(n)), n, dimnames = list(NULL, endog))
    Zx <- matrix(vapply(instr, function(v) iv_num(d[[v]]), numeric(n)), n, dimnames = list(NULL, instr))
    W <- iv_design(d, exog)
    w <- if (!is.null(weights)) iv_num(d[[weights]]) else rep(1, n)
    if (any(w <= 0)) stop("Weights must be positive.", call. = FALSE)
    fes <- lapply(absorb, function(a) as.character(d[[a]]))
    names(fes) <- absorb
    cl <- if (length(cluster)) lapply(cluster, function(a) as.character(d[[a]])) else NULL
    if (length(cl)) names(cl) <- cluster
    xtra <- if (length(extra)) matrix(vapply(extra, function(v) iv_num(d[[v]]), numeric(n)), n, dimnames = list(NULL, extra)) else NULL
    list(y = yv, X1 = X1, Z1 = Zx, W = W, w = w, fes = fes, cl = cl, n_drop = n_drop, extra = xtra,
         endog = endog, instr = instr, exog = colnames(W), has_const = !length(absorb))
}

# Full regressor and instrument matrices after fixed-effect removal and weighting.
iv_matrices <- function(pp) {
    cols <- list(y = matrix(pp$y, ncol = 1), X1 = pp$X1, Z1 = pp$Z1)
    if (!is.null(pp$W)) cols$W <- pp$W
    if (pp$has_const) cols$C <- matrix(1, length(pp$y), 1)
    all <- do.call(cbind, cols)
    nm <- c("y", colnames(pp$X1), colnames(pp$Z1), if (!is.null(pp$W)) colnames(pp$W), if (pp$has_const) "(Intercept)")
    if (length(pp$fes)) all <- iv_demean(all, pp$fes, pp$w)
    all <- all * sqrt(pp$w)
    colnames(all) <- nm
    wn <- c(if (!is.null(pp$W)) colnames(pp$W), if (pp$has_const) "(Intercept)")
    X1 <- all[, colnames(pp$X1), drop = FALSE]
    Z1 <- all[, colnames(pp$Z1), drop = FALSE]
    Wm <- all[, wn, drop = FALSE]
    list(y = all[, 1], X = cbind(X1, Wm), Z = cbind(Z1, Wm), X1 = X1, Z1 = Z1, W = Wm)
}

iv_qr_fit <- function(M, X) qr.fitted(qr(M, tol = 1e-10), X)

# Number of parameters absorbed by the fixed effects (levels nested in the first cluster variable are not counted).
iv_fe_df <- function(pp) {
    if (!length(pp$fes)) return(0L)
    tot <- 0L; nn <- 0L
    for (f in pp$fes) {
        lv <- length(unique(f))
        nested <- !is.null(pp$cl) && all(tapply(pp$cl[[1]], f, function(z) length(unique(z)) == 1L))
        if (!nested) { tot <- tot + lv; nn <- nn + 1L }
    }
    if (nn > 1L) tot <- tot - (nn - 1L)
    if (nn == 0L) tot <- 1L          # all fixed effects nested in clusters: only the constant is absorbed
    as.integer(tot)
}

# k-class estimation: k = 0 OLS, k = 1 2SLS, k = LIML eigenvalue, Fuller k = LIML - alpha / (N - L).
iv_kclass <- function(mm, k) {
    X <- mm$X; y <- mm$y
    PzX <- if (k == 0) X else iv_qr_fit(mm$Z, X)
    Xk <- (1 - k) * X + k * PzX
    A <- crossprod(Xk, X)
    beta <- drop(solve(A, crossprod(Xk, y)))
    list(beta = beta, u = drop(y - X %*% beta), Xk = Xk, A = A)
}

iv_liml_k <- function(mm) {
    Y <- cbind(mm$y, mm$X1)
    Wm <- mm$W
    MW <- if (ncol(Wm)) Y - iv_qr_fit(Wm, Y) else Y
    MZ <- Y - iv_qr_fit(mm$Z, Y)
    min(Re(eigen(solve(crossprod(Y, MZ), crossprod(Y, MW)), only.values = TRUE)$values))
}

# Variance estimators built from the score matrix Xk * u.
iv_vcov <- function(fit, mm, pp, vce, small = TRUE, extra_df = 0L) {
    N <- length(fit$u); K <- ncol(mm$X) + extra_df
    Ainv <- solve(fit$A)
    out <- list(N = N, K = K, G = NULL)
    if (vce == "classical") {
        V <- sum(fit$u^2) / (if (small) N - K else N) * Ainv
        out$df <- if (small) N - K else Inf
    } else if (vce == "robust") {
        V <- Ainv %*% crossprod(fit$Xk * fit$u) %*% Ainv * (if (small) N / (N - K) else 1)
        out$df <- if (small) N - K else Inf
    } else {
        ncl <- length(pp$cl)
        if (ncl < 1) stop("Select a cluster variable for cluster-robust standard errors.", call. = FALSE)
        S <- fit$Xk * fit$u
        raw <- function(g) {
            if (length(unique(g)) < 2) stop("A cluster variable needs at least two clusters.", call. = FALSE)
            Ainv %*% crossprod(rowsum(S, g)) %*% Ainv
        }
        G1 <- length(unique(pp$cl[[1]]))
        V <- raw(pp$cl[[1]]); G <- G1
        if (vce == "cluster2") {
            if (ncl < 2) stop("Two-way clustering needs two cluster variables.", call. = FALSE)
            V <- V + raw(pp$cl[[2]]) - raw(interaction(pp$cl[[1]], pp$cl[[2]], drop = TRUE))
            G <- min(G1, length(unique(pp$cl[[2]])))
        }
        if (small) V <- V * (N - 1) / (N - K) * G / (G - 1)
        out$G <- G
        out$df <- if (small) G - 1 else Inf
    }
    out$V <- (V + t(V)) / 2
    out
}

# Covariance matrix of the moments Z'u / sqrt(N): heteroskedasticity or cluster robust.
iv_moment_cov <- function(Z, u, pp, vce) {
    N <- nrow(Z); g <- Z * u
    if (vce == "cluster2" && length(pp$cl) > 1)
        return((crossprod(rowsum(g, pp$cl[[1]])) + crossprod(rowsum(g, pp$cl[[2]])) -
                crossprod(rowsum(g, interaction(pp$cl[[1]], pp$cl[[2]], drop = TRUE)))) / N)
    if (vce == "cluster" && length(pp$cl)) crossprod(rowsum(g, pp$cl[[1]])) / N
    else crossprod(g) / N
}

# Main entry: estimate and return the coefficients, variance matrix, tests and the matrices used.
iv_fit <- function(pp, method = c("2sls", "liml", "fuller", "kclass", "gmm", "ols"), alpha = 1, kvalue = 1,
                   vce = c("classical", "robust", "cluster", "cluster2"), small = TRUE) {
    method <- match.arg(method); vce <- match.arg(vce)
    mm <- iv_matrices(pp)
    N <- length(mm$y)
    X <- mm$X; Z <- mm$Z
    if (qr(Z, tol = 1e-10)$rank < ncol(Z)) stop("The instruments and controls are perfectly collinear.", call. = FALSE)
    if (qr(X, tol = 1e-10)$rank < ncol(X)) stop("The regressors are perfectly collinear.", call. = FALSE)
    if (qr(cbind(iv_qr_fit(Z, mm$X1), mm$W), tol = 1e-8)$rank < ncol(mm$X1) + ncol(mm$W))
        stop("The instruments do not identify the endogenous regressors (the first stage has deficient rank).", call. = FALSE)
    efd <- iv_fe_df(pp)
    k <- switch(method, ols = 0, `2sls` = 1, gmm = 1, kclass = kvalue, liml = iv_liml_k(mm),
                fuller = iv_liml_k(mm) - alpha / (N - ncol(Z)))
    fit <- iv_kclass(mm, k)
    vc <- iv_vcov(fit, mm, pp, vce, small, efd)
    beta <- fit$beta; V <- vc$V
    nov <- ncol(pp$Z1) - ncol(pp$X1)
    J <- NULL
    if (method == "gmm") {
        if (vce == "classical") stop("Two-step GMM needs a robust or cluster-robust variance estimator.", call. = FALSE)
        if (vce %in% c("cluster", "cluster2") && vc$G <= ncol(Z))
            stop("Two-step GMM needs more clusters than instruments (including controls).", call. = FALSE)
        Wg <- solve(iv_moment_cov(Z, fit$u, pp, vce))
        XZ <- crossprod(X, Z) / N
        Vb <- solve(XZ %*% Wg %*% t(XZ)) / N
        beta <- drop(Vb %*% XZ %*% Wg %*% crossprod(Z, mm$y))
        fit$u <- drop(mm$y - X %*% beta)
        g <- crossprod(Z, fit$u) / N
        J <- list(stat = drop(N * t(g) %*% Wg %*% g), df = nov)
        V <- Vb * (if (small) N / (N - vc$K) else 1)
    } else if (method != "ols" && nov > 0) {
        if (vce == "classical") {
            J <- list(stat = drop(sum(iv_qr_fit(Z, fit$u)^2) / (sum(fit$u^2) / N)), df = nov)      # Sargan
        } else if (method == "2sls" && !(vce %in% c("cluster", "cluster2") && vc$G <= ncol(Z))) {
            # Hansen J: the GMM criterion minimised with the weight matrix from the 2SLS residuals (as in ivreg2)
            Wg <- solve(iv_moment_cov(Z, fit$u, pp, vce))
            XZ <- crossprod(X, Z) / N
            b2 <- drop(solve(XZ %*% Wg %*% t(XZ), XZ %*% Wg %*% crossprod(Z, mm$y) / N))
            g <- crossprod(Z, mm$y - X %*% b2) / N
            J <- list(stat = drop(N * t(g) %*% Wg %*% g), df = nov)
        }
    }
    if (!is.null(J)) { J$p <- stats::pchisq(J$stat, J$df, lower.tail = FALSE); J$label <- if (vce == "classical" && method != "gmm") "Sargan" else "Hansen J" }
    se <- sqrt(diag(V))
    names(beta) <- colnames(X); names(se) <- colnames(X)
    list(beta = beta, se = se, V = V, u = fit$u, N = N, K = vc$K, df = vc$df, G = vc$G, J = J, k = k,
         method = method, vce = vce, small = small, mm = mm, extra_df = efd,
         endog = colnames(pp$X1), controls = colnames(mm$W))
}

iv_coef_table <- function(res, level = 95) {
    df <- res$df
    stat <- res$beta / res$se
    p <- if (is.finite(df)) 2 * stats::pt(-abs(stat), df) else 2 * stats::pnorm(-abs(stat))
    q <- if (is.finite(df)) stats::qt(1 - (1 - level / 100) / 2, df) else stats::qnorm(1 - (1 - level / 100) / 2)
    data.frame(term = names(res$beta), estimate = unname(res$beta), se = unname(res$se), stat = unname(stat),
               p = unname(p), lo = unname(res$beta - q * res$se), hi = unname(res$beta + q * res$se),
               type = ifelse(names(res$beta) %in% res$endog, "Endogenous", "Control"), stringsAsFactors = FALSE)
}

# First-stage regression of one endogenous regressor on the excluded instruments and controls,
# with the Wald (F) test that the excluded instruments are jointly zero.
iv_first_stage <- function(pp, res, j, vce = res$vce, small = res$small) {
    mm <- res$mm
    fm <- list(y = mm$X1[, j], X = mm$Z, Z = mm$Z, X1 = mm$Z1, W = mm$W)
    f <- iv_kclass(fm, 0)
    vc <- iv_vcov(f, fm, pp, vce, small, res$extra_df)
    q <- ncol(mm$Z1)
    b <- f$beta[seq_len(q)]
    Vq <- vc$V[seq_len(q), seq_len(q), drop = FALSE]
    wald <- drop(t(b) %*% solve(Vq, b))
    Fs <- wald / q
    df2 <- vc$df
    yw <- if (ncol(mm$W)) fm$y - iv_qr_fit(mm$W, fm$y) else fm$y
    list(F = Fs, df1 = q, df2 = df2,
         p = if (is.finite(df2)) stats::pf(Fs, q, df2, lower.tail = FALSE) else stats::pchisq(wald, q, lower.tail = FALSE),
         coef = b, se = sqrt(diag(Vq)),
         partial_r2 = 1 - sum(f$u^2) / sum(yw^2))
}
