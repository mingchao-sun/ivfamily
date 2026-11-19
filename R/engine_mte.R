# Marginal treatment effects (Heckman and Vytlacil, 2005) by local instrumental variables with a polynomial specification, following the Stata command mtefe
# (Andresen, 2018) in its default (non-separate, polynomial) form.
#
# Model (generalised Roy): Y = D Y1 + (1 - D) Y0,  Y_j = X b_j + U_j,  D = 1(P(Z) >= U_D),  U_D ~ uniform(0, 1).  With separability E[U_j | X, U_D] = E[U_j | U_D],
#   E[Y | X, P = p] = X b_0 + p X (b_1 - b_0) + K(p),   K(p) = sum_{j = 2}^{L + 1} psi_j p^j,
#   MTE(x, u) = x (b_1 - b_0) + sum_{j = 2}^{L + 1} j psi_j u^(j - 1)      (a polynomial of degree L in u).
# The propensity score P comes from a probit or logit of D on the instruments and the covariates. Treatment effect parameters are averages of the MTE over the
# 99 points u = 0.01, ..., 0.99 with the weights of mtefe: ATE equal weights at the mean of X; ATT weights Pr(P > u) / E[P] at the mean of X weighted by P;
# ATUT weights Pr(P <= u) / (1 - E[P]) at the mean of X weighted by 1 - P.

mte_grid <- seq(0.01, 0.99, by = 0.01)

mte_first <- function(d, W, link = c("probit", "logit")) {
    link <- match.arg(link)
    fit <- suppressWarnings(stats::glm.fit(W, d, family = stats::binomial(link), control = list(maxit = 100, epsilon = 1e-10)))
    if (anyNA(fit$coefficients)) stop("The first stage has a collinear variable: remove one of the instruments or covariates that are linear combinations of the others.", call. = FALSE)
    p <- as.numeric(fit$fitted.values)
    if (mean(p < 1e-6 | p > 1 - 1e-6) > 0.02) stop("The instruments and covariates predict the treatment almost perfectly; the marginal treatment effect cannot be estimated.", call. = FALSE)
    list(p = p, coef = fit$coefficients)
}

# 2SLS coefficient of d with X1 (incl. constant) as controls and W (instruments, covariates, constant) as instruments.
mte_2sls <- function(y, d, X1, W) {
    dh <- stats::lm.fit(W, d)$fitted.values
    R <- cbind(dh, X1)
    b <- stats::lm.fit(R, y)$coefficients
    b[[1]]
}

# One fit. y outcome, d 0/1, X covariates (matrix or NULL), Z instruments (matrix). L degree of the MTE in u.
mte_fit <- function(y, d, X = NULL, Z, L = 2, link = "probit") {
    n <- length(y)
    X1 <- if (is.null(X)) matrix(1, n, 1, dimnames = list(NULL, "Constant")) else cbind(Constant = 1, X)
    W <- cbind(X1, Z)
    fs <- mte_first(d, W, link)
    p <- fs$p; k <- ncol(X1)
    H <- cbind(X1, X1 * p, vapply(2:(L + 1), function(j) p^j, numeric(n)))
    colnames(H) <- c(paste0("b0:", colnames(X1)), paste0("delta:", colnames(X1)), paste0("psi", 2:(L + 1)))
    qh <- qr(H)
    if (qh$rank < ncol(H)) stop("The second stage is not identified: the propensity score varies too little (weak instruments) for this polynomial degree or these covariates. Use a lower degree, fewer covariates or stronger instruments.", call. = FALSE)
    cf <- stats::lm.fit(H, y)$coefficients
    delta <- cf[k + seq_len(k)]; psi <- cf[2 * k + seq_len(L)]
    kfun <- function(u) vapply(u, function(v) sum((2:(L + 1)) * psi * v^(1:L)), numeric(1))        # k(u) = sum_j j psi_j u^(j-1)
    ku <- kfun(mte_grid)
    mte_at <- function(x1) as.numeric(sum(x1 * delta)) + ku
    pbar <- mean(p)
    xbar <- colMeans(X1); xatt <- colSums(X1 * p) / sum(p); xatut <- colSums(X1 * (1 - p)) / sum(1 - p)
    prop <- vapply(mte_grid, function(u) mean(p > u), numeric(1))
    mte_ate <- mte_at(xbar)
    ate <- mean(mte_ate)
    att <- sum((1 / 99) * (prop / pbar) * mte_at(xatt))
    atut <- sum((1 / 99) * ((1 - prop) / (1 - pbar)) * mte_at(xatut))
    iv <- tryCatch(mte_2sls(y, d, X1, W), error = function(e) NA_real_)
    list(coef = cf, b0 = cf[seq_len(k)], delta = delta, psi = psi, p = p, first = fs$coef, mte = mte_ate, mte_x = xbar, ku = ku,
         ate = ate, att = att, atut = atut, iv = iv, L = L, link = link, k = k, n = n, names = colnames(X1), xbar = xbar, xatt = xatt, xatut = xatut,
         support = list(all = range(p), treated = range(p[d == 1]), untreated = range(p[d == 0]), share_outside = mean(p < 0.01 | p > 0.99)))
}

# Bootstrap (units, or clusters) of the whole procedure: first stage, second stage and parameters. Returns the point estimates, standard errors, intervals and the tests.
mte_boot <- function(y, d, X = NULL, Z, L = 2, link = "probit", cluster = NULL, nboot = 300, level = 0.95, seed = 1234, checkpoint = NULL) {
    est <- mte_fit(y, d, X, Z, L, link)
    n <- length(y)
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
    cl <- if (!is.null(cluster)) split(seq_len(n), factor(cluster)) else NULL
    idx_draw <- function() if (is.null(cl)) sample.int(n, replace = TRUE) else unlist(cl[sample.int(length(cl), replace = TRUE)], use.names = FALSE)
    pars <- function(f) c(ate = f$ate, att = f$att, atut = f$atut, iv = f$iv)
    B <- matrix(NA_real_, nboot, 4 + length(mte_grid) + length(est$psi) + est$k - 1)
    for (b in seq_len(nboot)) {
        if (!is.null(checkpoint)) checkpoint()
        i <- idx_draw()
        f <- tryCatch(mte_fit(y[i], d[i], if (!is.null(X)) X[i, , drop = FALSE], Z[i, , drop = FALSE], L, link), error = function(e) NULL)
        if (!is.null(f)) B[b, ] <- c(pars(f), f$mte, f$psi, f$delta[-1])
    }
    ok <- stats::complete.cases(B)
    if (sum(ok) < 0.5 * nboot) stop("Most bootstrap samples could not be estimated (the propensity score varies too little). Use stronger instruments or a lower polynomial degree.", call. = FALSE)
    B <- B[ok, , drop = FALSE]
    z <- stats::qnorm(1 - (1 - level) / 2)
    ip <- 1:4; im <- 4 + seq_along(mte_grid); ips <- 4 + length(mte_grid) + seq_along(est$psi); idl <- 4 + length(mte_grid) + length(est$psi) + seq_len(est$k - 1)
    se <- apply(B, 2, stats::sd)
    wald <- function(b, idx) { if (!length(idx)) return(c(chisq = NA, df = 0, p = NA)); V <- stats::cov(B[, idx, drop = FALSE]); w <- tryCatch(as.numeric(t(b) %*% solve(V, b)), error = function(e) NA_real_)
                               c(chisq = w, df = length(idx), p = if (is.finite(w)) stats::pchisq(w, length(idx), lower.tail = FALSE) else NA_real_) }
    parm <- cbind(estimate = unlist(pars(est)), se = se[ip]); parm <- cbind(parm, lo = parm[, 1] - z * parm[, 2], hi = parm[, 1] + z * parm[, 2], p = 2 * stats::pnorm(-abs(parm[, 1] / parm[, 2])))
    rownames(parm) <- c("ATE", "ATT", "ATUT", "2SLS")
    curve <- cbind(u = mte_grid, mte = est$mte, se = se[im]); curve <- cbind(curve, lo = curve[, "mte"] - z * curve[, "se"], hi = curve[, "mte"] + z * curve[, "se"])
    list(est = est, parameters = parm, curve = curve, nboot = sum(ok), test_unobserved = wald(est$psi, ips), test_observed = wald(est$delta[-1], idl), level = level)
}
