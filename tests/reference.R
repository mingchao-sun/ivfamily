# Engine vs reference packages (AER::ivreg + sandwich, fixest). Run from the module root:  Rscript tests/reference.R
source("R/engine_iv.R")
suppressMessages({library(AER); library(sandwich); library(fixest)})
nfail <- 0
chk <- function(name, got, want, tol) {
    d <- max(abs(unname(got) - unname(want)))
    ok <- is.finite(d) && d < tol
    cat(sprintf("%-58s %s  max |diff| = %.1e\n", name, if (ok) "ok  " else "FAIL", d))
    if (!ok) nfail <<- nfail + 1
}

# ---- cigarette demand (Stock and Watson): log quantity on log price, instrument = sales tax --------------------
data("CigarettesSW")
cg <- subset(CigarettesSW, year == "1995")
cg$rprice <- cg$price / cg$cpi; cg$rincome <- cg$income / cg$population / cg$cpi
cg$lq <- log(cg$packs); cg$lp <- log(cg$rprice); cg$li <- log(cg$rincome)
cg$tdiff <- (cg$taxs - cg$tax) / cg$cpi; cg$tax_r <- cg$tax / cg$cpi
cg$rp <- cg$price
m <- ivreg(lq ~ lp + li | li + tdiff + tax_r, data = cg)          # over-identified
pp <- iv_prepare(cg, "lq", "lp", c("tdiff", "tax_r"), "li")
r <- iv_fit(pp, "2sls", vce = "classical")
chk("2SLS coefficients (over-identified)", r$beta, coef(m)[c("lp", "li", "(Intercept)")], 1e-8)
chk("2SLS classical s.e.", r$se, sqrt(diag(vcov(m)))[c("lp", "li", "(Intercept)")], 1e-8)
r <- iv_fit(pp, "2sls", vce = "robust")
chk("2SLS HC1 s.e. vs sandwich", r$se, sqrt(diag(vcovHC(m, "HC1")))[c("lp", "li", "(Intercept)")], 1e-8)
pp2 <- iv_prepare(cg, "lq", "lp", c("tdiff", "tax_r"), "li", cluster = "state")
r <- iv_fit(pp2, "2sls", vce = "cluster")
chk("2SLS cluster s.e. vs sandwich", r$se, sqrt(diag(vcovCL(m, cluster = cg$state, type = "HC1")))[c("lp", "li", "(Intercept)")], 1e-8)
# Sargan statistic against the textbook formula
r <- iv_fit(pp, "2sls", vce = "classical")
u <- resid(m); Zm <- model.matrix(~ li + tdiff + tax_r, cg)
sar <- sum(fitted(lm(u ~ Zm - 1))^2) / (sum(u^2) / length(u))
chk("Sargan statistic", r$J$stat, sar, 1e-8)
# first-stage F vs a direct regression
fs <- iv_first_stage(pp, r, "lp")
a <- anova(lm(lp ~ li, cg), lm(lp ~ li + tdiff + tax_r, cg))
chk("first-stage F (classical)", fs$F, a$F[2], 1e-8)

# ---- college distance: controls with factors, cluster, weights, fixed effects -----------------------------
data("CollegeDistance")
cd <- CollegeDistance
cd$ed <- cd$education; cd$hi <- as.numeric(cd$distance)
m2 <- ivreg(wage ~ ed + gender + ethnicity + unemp | distance + tuition + gender + ethnicity + unemp, data = cd)
pp <- iv_prepare(cd, "wage", "ed", c("distance", "tuition"), c("gender", "ethnicity", "unemp"))
r <- iv_fit(pp, "2sls", vce = "classical")
cn <- names(coef(m2)); idx <- match(c("ed", "genderfemale", "ethnicityafam", "ethnicityhispanic", "unemp", "(Intercept)"), cn)
nm <- c("ed", "genderfemale", "ethnicityafam", "ethnicityhispanic", "unemp", "(Intercept)")
chk("college distance: 2SLS with factor controls", r$beta[nm], coef(m2)[nm], 1e-7)
chk("college distance: classical s.e.", r$se[nm], sqrt(diag(vcov(m2)))[nm], 1e-7)
set.seed(7); cd$grp <- sample(1:120, nrow(cd), TRUE)
pp <- iv_prepare(cd, "wage", "ed", c("distance", "tuition"), c("gender", "unemp"), cluster = "grp")
r <- iv_fit(pp, "2sls", vce = "cluster")
mf <- feols(wage ~ gender + unemp | ed ~ distance + tuition, data = cd, cluster = ~grp, ssc = ssc(adj = TRUE, cluster.adj = TRUE))
chk("cluster s.e. vs fixest (t(G-1) corrections)", r$se[["ed"]], se(mf)[["fit_ed"]], 1e-6)

# fixed effects: one and two sets, against fixest
set.seed(1)
n <- 3000; g1 <- sample(1:60, n, TRUE); g2 <- sample(1:12, n, TRUE); x2 <- rnorm(n)
z1 <- rnorm(n) + 0.3 * g1 / 60; z2 <- rnorm(n)
e <- rnorm(n); xe <- 0.5 * z1 + 0.5 * z2 + 0.6 * e + 0.2 * (g1 %% 7) + rnorm(n)
yy <- 1 + 2 * xe + 0.5 * x2 + 0.3 * (g1 %% 5) + 0.2 * g2 + e
d <- data.frame(y = yy, x = xe, z1, z2, x2, g1 = factor(g1), g2 = factor(g2))
pp <- iv_prepare(d, "y", "x", c("z1", "z2"), "x2", absorb = "g1")
r <- iv_fit(pp, "2sls", vce = "classical")
mf <- feols(y ~ x2 | g1 | x ~ z1 + z2, data = d, vcov = "iid", ssc = ssc(adj = TRUE, fixef.K = "full"))
chk("one set of fixed effects: coefficients", r$beta[c("x", "x2")], coef(mf)[c("fit_x", "x2")], 1e-8)
chk("one set of fixed effects: classical s.e.", r$se[c("x", "x2")], se(mf)[c("fit_x", "x2")], 1e-7)
pp <- iv_prepare(d, "y", "x", c("z1", "z2"), "x2", absorb = c("g1", "g2"))
r <- iv_fit(pp, "2sls", vce = "robust")
mf <- feols(y ~ x2 | g1 + g2 | x ~ z1 + z2, data = d, vcov = "hetero", ssc = ssc(adj = TRUE, fixef.K = "full"))
chk("two sets of fixed effects: coefficients", r$beta[c("x", "x2")], coef(mf)[c("fit_x", "x2")], 1e-7)
chk("two sets of fixed effects: HC1 s.e.", r$se[c("x", "x2")], se(mf)[c("fit_x", "x2")], 1e-6)

# weights: a weight of w equals w copies of the row
d$w <- sample(1:3, n, TRUE)
dd <- d[rep(seq_len(n), d$w), ]
pw <- iv_prepare(d, "y", "x", c("z1", "z2"), "x2", weights = "w"); rw <- iv_fit(pw, "2sls", vce = "classical", small = FALSE)
pc <- iv_prepare(dd, "y", "x", c("z1", "z2"), "x2");                rc <- iv_fit(pc, "2sls", vce = "classical", small = FALSE)
chk("weights = replicated rows: coefficients (frequency weights)", rw$beta, rc$beta, 1e-8)
# weights are analytic (aweight) weights as in ivreg2 [aw]: standard errors do not change when all weights are rescaled
d$w10 <- 10 * d$w
pw10 <- iv_prepare(d, "y", "x", c("z1", "z2"), "x2", weights = "w10"); rw10 <- iv_fit(pw10, "2sls", vce = "robust")
rw1 <- iv_fit(pw, "2sls", vce = "robust")
chk("weights rescaled by 10: coefficients and robust s.e.", c(rw10$beta, rw10$se), c(rw1$beta, rw1$se), 1e-9)

# LIML / Fuller / GMM identities
pp <- iv_prepare(d, "y", "x", c("z1", "z2"), "x2")
r1 <- iv_fit(pp, "2sls"); rl <- iv_fit(pp, "liml"); rf <- iv_fit(pp, "fuller", alpha = 1); rk <- iv_fit(pp, "kclass", kvalue = 1)
chk("k-class with k = 1 equals 2SLS", rk$beta, r1$beta, 1e-10)
cat(sprintf("%-58s %.6f (>= 1)\n", "LIML k", rl$k)); if (rl$k < 1 - 1e-9) nfail <- nfail + 1
chk("Fuller k = LIML k - 1/(N - L)", rf$k, rl$k - 1 / (rl$N - 4), 1e-10)
# just-identified: all estimators coincide
pj <- iv_prepare(d, "y", "x", "z1", "x2")
a <- iv_fit(pj, "2sls"); b <- iv_fit(pj, "liml"); cc <- iv_fit(pj, "gmm", vce = "robust")
chk("just-identified: LIML equals 2SLS", b$beta, a$beta, 1e-8)
chk("just-identified: GMM equals 2SLS", cc$beta, a$beta, 1e-8)
# homoskedastic GMM equals 2SLS in the over-identified case only asymptotically; check it is close and J is non-negative
g <- iv_fit(pp, "gmm", vce = "robust")
cat(sprintf("%-58s %.4f vs 2SLS %.4f, Hansen J %.3f\n", "GMM two-step coefficient on x", g$beta[["x"]], r1$beta[["x"]], g$J$stat))
if (!(g$J$stat >= 0)) nfail <- nfail + 1

cat(if (nfail == 0) "\nAll IV engine reference checks passed.\n" else sprintf("\n%d check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)

# ---- plausibly exogenous: b(g) = b - A g must equal 2SLS on y - Z g (linear-shift version) -------------------------------------------------
source("R/engine_plaus.R")
pp <- iv_prepare(d, "y", "x", c("z1", "z2"), "x2"); r <- iv_fit(pp, "2sls", vce = "robust", small = FALSE)
A <- iv_plaus_A(r); g <- c(0.3, -0.2)
d2 <- d; d2$ys <- d$y - g[1] * d$z1 - g[2] * d$z2
r2 <- iv_fit(iv_prepare(d2, "ys", "x", c("z1", "z2"), "x2"), "2sls", vce = "robust", small = FALSE)
chk("plausibly exogenous: b(g) = b - A g", r$beta[["x"]] - drop(A %*% g), r2$beta[["x"]], 1e-9)
u <- iv_plaus_uci(r, c(-0.2, -0.2), c(0.2, 0.2))
tip <- iv_plaus_tipping(r)
chk("tipping delta puts a union bound exactly at zero", abs(iv_plaus_uci(r, rep(-tip, 2), rep(tip, 2))[c("lo", "hi")][[if (r$beta[["x"]] > 0) "lo" else "hi"]]), 0, 1e-9)
if (nfail > 0) quit(status = 1)
