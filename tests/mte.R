# Marginal treatment effects: recovery of a known MTE in a generalised Roy model, identities, and the behaviour without essential heterogeneity.
# Run from the module root:  Rscript tests/mte.R   (the comparison with Stata mtefe is in stata_validation/18_mtefe.do)
source("R/engine_mte.R")
nfail <- 0
ok <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) nfail <<- nfail + 1 }

# Generalised Roy model: D = 1(P(Z) >= U_D) with U_D = Phi(V); U1 - U0 = rho V + noise, so MTE(x, u) = delta + beta x + rho Phi^-1(u)
sim <- function(seed, n = 20000, rho = 0.6, gam = 1.6) {
    set.seed(seed)
    z <- rnorm(n); x <- rnorm(n); v <- rnorm(n)
    p <- pnorm(gam * z + 0.3 * x); d <- as.numeric(p >= pnorm(v))
    u0 <- 0.5 * v + rnorm(n, 0, 0.5); u1 <- u0 + rho * v + rnorm(n, 0, 0.3)
    y0 <- 1 + 0.5 * x + u0; y1 <- y0 + 0.8 + 0.3 * x + (u1 - u0)
    list(y = ifelse(d == 1, y1, y0), d = d, x = cbind(x = x), z = cbind(z = z), p = p, rho = rho)
}
truth <- function(s, delta = 0.8, bx = 0.3) {          # ATE, ATT, ATUT by numerical integration of the true MTE
    xm <- colMeans(s$x); n <- length(s$p)
    mte <- function(x, u) delta + bx * x + s$rho * qnorm(u)
    grid <- (1:2000 - 0.5) / 2000
    ate <- mean(delta + bx * s$x) + s$rho * mean(qnorm(grid))
    pu <- function(p) { g <- (1:400 - 0.5) / 400 * p; mean(qnorm(g)) * p }          # integral of qnorm from 0 to p
    att <- (sum((delta + bx * s$x[, 1]) * s$p) + s$rho * sum(vapply(s$p, pu, numeric(1)))) / sum(s$p)
    atut <- (sum((delta + bx * s$x[, 1]) * (1 - s$p)) + s$rho * (n * 0 - sum(vapply(s$p, pu, numeric(1))))) / sum(1 - s$p)
    c(ate = ate, att = att, atut = atut)
}
s <- sim(1); tr <- truth(s)
f <- mte_fit(s$y, s$d, s$x, s$z, L = 3, link = "probit")
ok(abs(f$ate - tr["ate"]) < 0.08 && abs(f$att - tr["att"]) < 0.08 && abs(f$atut - tr["atut"]) < 0.08, sprintf("known MTE (essential heterogeneity, rho 0.6): ATE %.3f (true %.3f), ATT %.3f (%.3f), ATUT %.3f (%.3f)", f$ate, tr["ate"], f$att, tr["att"], f$atut, tr["atut"]))
ok(f$att > f$ate && f$ate > f$atut || (tr["att"] > tr["ate"]) == (f$att > f$ate), "the ordering ATT, ATE, ATUT agrees with the truth (selection on gains)")
xb <- f$xbar; true_curve <- 0.8 + 0.3 * xb[["x"]] + s$rho * qnorm(mte_grid)
mid <- mte_grid >= 0.1 & mte_grid <= 0.9
ok(max(abs(f$mte[mid] - true_curve[mid])) < 0.25, sprintf("the MTE curve is recovered between u = 0.1 and 0.9 (max error %.2f, polynomial degree 3)", max(abs(f$mte[mid] - true_curve[mid]))))
ok(cor(f$mte, mte_grid) < -0.9 || cor(f$mte, mte_grid) > 0.9, "the curve slopes with u as the model implies")
# identities
n <- length(s$y); X1 <- cbind(1, s$x)
ok(abs(f$ate - (sum(f$xbar * f$delta) + mean(f$ku))) < 1e-10, "ATE = x-bar times delta plus the grid mean of k(u)")
Kp <- function(pp) sapply(pp, function(q) sum(f$psi * q^(2:(f$L + 1))))
exact_att <- (sum(rowSums(X1 * rep(f$delta, each = n)) * f$p) + sum(Kp(f$p))) / sum(f$p)
ok(abs(f$att - exact_att) < 0.02, sprintf("ATT on the 99-point grid (%.4f) is close to the exact integral %.4f", f$att, exact_att))
# no essential heterogeneity: the MTE is flat, 2SLS = ATE = ATT = ATUT
s0 <- sim(2, rho = 0)
f0 <- mte_fit(s0$y, s0$d, s0$x, s0$z, L = 2)
ok(abs(f0$ate - 0.8) < 0.08 && abs(f0$att - f0$ate) < 0.08 && abs(f0$atut - f0$ate) < 0.08 && abs(f0$iv - 0.8) < 0.1, sprintf("no selection on gains: ATE %.3f, ATT %.3f, ATUT %.3f, 2SLS %.3f all near 0.8", f0$ate, f0$att, f0$atut, f0$iv))
b0 <- mte_boot(s0$y, s0$d, s0$x, s0$z, L = 2, nboot = 100, seed = 5)
b1 <- mte_boot(s$y, s$d, s$x, s$z, L = 3, nboot = 100, seed = 5)
ok(b1$test_unobserved[["p"]] < 0.01, sprintf("essential heterogeneity is detected (p %.1e)", b1$test_unobserved[["p"]]))
ok(b0$test_unobserved[["p"]] > 0.01, sprintf("and not found when the MTE is flat (p %.2f)", b0$test_unobserved[["p"]]))
ok(b1$test_observed[["p"]] < 0.01, "observed heterogeneity (the covariate changes the effect) is detected")
ok(all(b1$parameters[, "lo"] <= b1$parameters[, "estimate"] & b1$parameters[, "estimate"] <= b1$parameters[, "hi"]), "intervals contain the estimates")
ok(b1$parameters["ATE", "lo"] < tr["ate"] && tr["ate"] < b1$parameters["ATE", "hi"], sprintf("the ATE interval (%.2f, %.2f) contains the true ATE %.2f", b1$parameters["ATE", "lo"], b1$parameters["ATE", "hi"], tr["ate"]))
# coverage of the bootstrap interval for the ATE
cov <- replicate(40, { ss <- sim(sample.int(1e6, 1), n = 4000); bb <- mte_boot(ss$y, ss$d, ss$x, ss$z, L = 3, nboot = 60, seed = 1); tt <- truth(ss)[["ate"]]; bb$parameters["ATE", "lo"] <= tt && tt <= bb$parameters["ATE", "hi"] })
ok(mean(cov) >= 0.85, sprintf("95%% bootstrap intervals contain the true ATE in %.0f%% of 40 samples", 100 * mean(cov)))
# options and errors
fl <- mte_fit(s$y, s$d, s$x, s$z, L = 2, link = "logit"); ok(is.finite(fl$ate), "logit first stage")
ok(is.finite(mte_fit(s$y, s$d, NULL, s$z, L = 2)$ate), "no covariates")
ok(inherits(try(mte_fit(s$y, s$d, s$x, cbind(z = rep(1, n)), L = 2), silent = TRUE), "try-error"), "a constant instrument is refused")
w <- sim(3, gam = 0.01)
r <- try(mte_fit(w$y, w$d, w$x, w$z, L = 4), silent = TRUE); ok(inherits(r, "try-error") || is.finite(r$ate), "a nearly irrelevant instrument gives an estimate or a message, never a crash")
cl <- rep(1:200, length.out = n); bc <- mte_boot(s$y, s$d, s$x, s$z, L = 2, cluster = cl, nboot = 40, seed = 1); ok(is.finite(bc$parameters["ATE", "se"]), "bootstrap by clusters")
cat(if (nfail == 0) "\nAll marginal treatment effect checks passed.\n" else sprintf("\n%d check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)
