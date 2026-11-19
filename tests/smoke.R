# Headless smoke test of the jamovi analysis classes. Run from the module root:
#   Rscript tests/smoke.R                         (local R)
# or in jamovi's own R with the bundled libraries (see PLAN.md).
suppressMessages({library(jmvcore); library(R6); library(ggplot2)})
for (f in c("R/00jmv.R", "R/engine_iv.R", "R/engine_weak.R", "R/engine_tests.R", "R/engine_ar.R", "R/engine_late.R", "R/engine_plaus.R", "R/engine_mte.R", "R/common.R", "R/guide.R", list.files("R", "[.]h[.]R$", full.names = TRUE),
            list.files("R", "[.]b[.]R$", full.names = TRUE))) source(f)
cg <- read.csv("data/cigarettes_sw.csv"); cd <- read.csv("data/college_distance.csv")
cd$gender <- factor(cd$gender); cd$ethnicity <- factor(cd$ethnicity)
bad <- 0
run <- function(label, data, ...) {
    a <- list(...)
    fm <- formals(iv)
    for (n in setdiff(names(fm)[vapply(fm, function(z) identical(z, quote(expr = )), logical(1))], c("data", names(a)))) a[n] <- list(NULL)
    r <- tryCatch({ res <- do.call(iv, c(list(data = data), a)); as.data.frame(res$coef); "ok" },
                  error = function(e) paste("ERROR:", conditionMessage(e)))
    expect_err <- grepl("^EXPECT", label)
    cat(sprintf("%-52s %s\n", label, if (expect_err) sub("^ERROR: ", "message: ", r) else r))
    if (!expect_err && r != "ok") bad <<- bad + 1
    invisible(r)
}
ci <- list(outcome = "lq", endog = "lp", instruments = c("tdiff", "tax_r"), controls = "li")
do.call(run, c(list("cigarettes 2SLS classical", cg), ci))
do.call(run, c(list("cigarettes 2SLS robust", cg), ci, list(vce = "robust")))
do.call(run, c(list("cigarettes LIML", cg), ci, list(method = "liml")))
do.call(run, c(list("cigarettes Fuller", cg), ci, list(method = "fuller", alpha = 1)))
do.call(run, c(list("cigarettes k-class", cg), ci, list(method = "kclass", kvalue = 0.5)))
do.call(run, c(list("cigarettes GMM robust", cg), ci, list(method = "gmm", vce = "robust")))
do.call(run, c(list("cigarettes notes + guide", cg), ci, list(showHelp = TRUE, showGuide = TRUE)))
do.call(run, c(list("cigarettes large-sample", cg), ci, list(small = FALSE, showControls = FALSE)))
cd2 <- list(outcome = "wage", endog = "education", instruments = c("distance", "tuition"), controls = c("gender", "ethnicity", "unemp"))
do.call(run, c(list("college 2SLS robust, factor controls", cd), cd2, list(vce = "robust")))
set.seed(1); cd$grp <- factor(sample(1:100, nrow(cd), TRUE)); cd$grp2 <- factor(sample(1:50, nrow(cd), TRUE))
do.call(run, c(list("college cluster", cd), cd2, list(vce = "cluster", cluster = "grp")))
do.call(run, c(list("college two-way cluster", cd), cd2, list(vce = "cluster2", cluster = c("grp", "grp2"))))
do.call(run, c(list("college absorb fixed effects", cd), cd2[c("outcome", "endog", "instruments")], list(controls = "unemp", absorb = "grp", vce = "robust")))
do.call(run, c(list("college two sets of fixed effects", cd), cd2[c("outcome", "endog", "instruments")], list(controls = "unemp", absorb = c("grp", "grp2"))))
cd$w <- runif(nrow(cd), 0.5, 2)
do.call(run, c(list("college weights", cd), cd2, list(weights = "w", vce = "robust")))
# expected messages for bad input
do.call(run, c(list("EXPECT underidentified", cg), list(outcome = "lq", endog = c("lp", "li"), instruments = "tdiff")))
do.call(run, c(list("EXPECT cluster without variable", cg), ci, list(vce = "cluster")))
do.call(run, c(list("EXPECT GMM with classical errors", cg), ci, list(method = "gmm")))
do.call(run, c(list("EXPECT variable in two roles", cg), list(outcome = "lq", endog = "lp", instruments = c("lp", "tdiff"))))
cg2 <- cg; cg2$tdiff <- cg2$tax_r * 2
do.call(run, c(list("EXPECT collinear instruments", cg2), list(outcome = "lq", endog = "lp", instruments = c("tdiff", "tax_r"), controls = "li")))
cg3 <- cg; cg3$lq[1:3] <- NA
do.call(run, c(list("missing values are dropped", cg3), ci))

# ---- diagnostics: ivfirst and ivover on the simulated two-endogenous data -------------------------------------
sm <- read.csv("data/sim_two_endogenous.csv"); sm$g <- factor(sm$g); sm$g2 <- factor(sample(1:40, nrow(sm), TRUE))
run2 <- function(label, fn, data, tab, ...) {
    a <- list(...)
    fm <- formals(fn)
    for (n in setdiff(names(fm)[vapply(fm, function(z) identical(z, quote(expr = )), logical(1))], c("data", names(a)))) a[n] <- list(NULL)
    r <- tryCatch({ res <- do.call(fn, c(list(data = data), a)); n <- nrow(as.data.frame(res[[tab]])); if (n < 1) stop("empty table"); "ok" },
                  error = function(e) paste("ERROR:", conditionMessage(e)))
    expect_err <- grepl("^EXPECT", label)
    cat(sprintf("%-52s %s
", label, if (expect_err) sub("^ERROR: ", "message: ", r) else r))
    if (!expect_err && r != "ok") bad <<- bad + 1
}
b2 <- list(outcome = "y", endog = c("x1", "x2"), instruments = c("z1", "z2", "z3", "z4"), controls = "w")
for (v in c("classical", "robust")) {
    run2(paste("ivfirst", v), ivfirst, sm, "joint", outcome = b2$outcome, endog = b2$endog, instruments = b2$instruments, controls = b2$controls, vce = v)
    run2(paste("ivover", v), ivover, sm, "overid", outcome = b2$outcome, endog = b2$endog, instruments = c("z1", "z2", "z4"), controls = b2$controls, vce = v, suspects = "z3")
}
run2("ivfirst cluster + notes", ivfirst, sm, "each", outcome = "y", endog = b2$endog, instruments = b2$instruments, controls = "w", vce = "cluster", cluster = "g", showHelp = TRUE, showGuide = TRUE)
run2("ivfirst two-way cluster", ivfirst, sm, "each", outcome = "y", endog = b2$endog, instruments = b2$instruments, controls = "w", vce = "cluster2", cluster = c("g", "g2"))
run2("ivfirst single endogenous", ivfirst, sm, "joint", outcome = "y", endog = "x1", instruments = c("z1", "z2"), controls = "w", vce = "robust")
run2("ivfirst with fixed effects", ivfirst, sm, "each", outcome = "y", endog = b2$endog, instruments = b2$instruments, absorb = "g", vce = "robust")
run2("ivover cluster", ivover, sm, "overid", outcome = "y", endog = "x2", endogTest = "x1", instruments = b2$instruments, controls = "w", vce = "cluster", cluster = "g", showHelp = TRUE, showGuide = TRUE)
run2("ivover exactly identified", ivover, sm, "endog", outcome = "y", endog = "x1", instruments = "z1", controls = "w", vce = "robust")
run2("EXPECT ivover suspects too many", ivover, sm, "ctest", outcome = "y", endog = b2$endog, instruments = "z4", controls = "w", suspects = c("z1", "z2", "z3"))
run2("ivover suspects only, separate box", ivover, sm, "ctest", outcome = "y", endog = b2$endog, instruments = c("z1", "z2", "z4"), controls = "w", suspects = "z3", vce = "robust")
run2("EXPECT ivfirst cluster without variable", ivfirst, sm, "joint", outcome = "y", endog = b2$endog, instruments = b2$instruments, vce = "cluster")
sw <- read.csv("data/sim_weak.csv"); sw$g <- factor(sw$g)
for (v in c("classical", "robust")) run2(paste("ivweak", v), ivweak, sw, "tests", outcome = "y", endog = "x", instruments = c("z1", "z2", "z3"), controls = "w", vce = v, nullValue = 0.5)
run2("ivweak cluster + notes", ivweak, sw, "csets", outcome = "y", endog = "x", instruments = c("z1", "z2", "z3"), controls = "w", vce = "cluster", cluster = "g", showHelp = TRUE, showGuide = TRUE)
run2("ivweak just identified", ivweak, sw, "csets", outcome = "y", endog = "x", instruments = "z1", controls = "w", vce = "robust")
run2("ivweak strong instruments", ivweak, sm, "csets", outcome = "y", endog = "x1", instruments = c("z1", "z2", "z3", "z4"), controls = "w", vce = "robust")
run2("EXPECT ivweak two endogenous", ivweak, sm, "tests", outcome = "y", endog = c("x1", "x2"), instruments = c("z1", "z2", "z3"), controls = "w")

# ---- ivcmp, ivlate, ivplaus ------------------------------------------------------------------------------------
run3 <- function(label, fn, data, tab, ...) {
    a <- list(...)
    fm <- formals(fn)
    for (n in setdiff(names(fm)[vapply(fm, function(z) identical(z, quote(expr = )), logical(1))], c("data", names(a)))) a[n] <- list(NULL)
    r <- tryCatch({ res <- do.call(fn, c(list(data = data), a)); n <- nrow(as.data.frame(res[[tab]])); if (n < 1) stop("empty table"); "ok" },
                  error = function(e) paste("ERROR:", conditionMessage(e)))
    expect_err <- grepl("^EXPECT", label)
    cat(sprintf("%-52s %s", label, if (expect_err) sub("^ERROR: ", "message: ", r) else r), "\n", sep = "")
    if (!expect_err && r != "ok") bad <<- bad + 1
}
ft <- read.csv("data/fertility.csv"); ft <- ft[seq(1, nrow(ft), by = 10), ]
for (v in c("classical", "robust"))
    run3(paste("ivcmp", v), ivcmp, sw, "est", outcome = "y", endog = "x", instruments = c("z1", "z2", "z3"), controls = "w", vce = v, showGmm = TRUE, showFuller = TRUE)
run3("ivcmp two endogenous + notes", ivcmp, sm, "est", outcome = "y", endog = c("x1", "x2"), instruments = c("z1", "z2", "z3", "z4"), controls = "w", vce = "robust", showHelp = TRUE, showGuide = TRUE)
run3("ivlate", ivlate, ft, "effects", outcome = "weeks", treatment = "morekids", instrument = "samesex", profile = c("age", "boy1", "afam"), vce = "robust", showHelp = TRUE, showGuide = TRUE)
run3("ivlate with controls", ivlate, ft, "shares", outcome = "weeks", treatment = "morekids", instrument = "samesex", controls = c("age", "boy1"), vce = "classical")
ft2 <- ft; ft2$morekids <- factor(ifelse(ft2$morekids == 1, "yes", "no"))
run3("ivlate with a two-level factor treatment", ivlate, ft2, "effects", outcome = "weeks", treatment = "morekids", instrument = "samesex", vce = "robust")
run3("EXPECT ivlate non-binary instrument", ivlate, ft, "effects", outcome = "weeks", treatment = "morekids", instrument = "age", vce = "robust")
for (v in c("classical", "robust"))
    run3(paste("ivplaus", v), ivplaus, cg, "uci", outcome = "lq", endog = "lp", instruments = c("tdiff", "tax_r"), controls = "li", vce = v, delta = 0.05, ltzSd = 0.02)
run3("ivplaus two endogenous + notes", ivplaus, sm, "tip", outcome = "y", endog = c("x1", "x2"), instruments = c("z1", "z2", "z3", "z4"), controls = "w", vce = "robust", showHelp = TRUE, showGuide = TRUE)
mt <- read.csv("data/sim_mte.csv", stringsAsFactors = TRUE)
for (dg in 1:3) run3(paste("ivmte degree", dg), ivmte, mt, "pars", outcome = "logwage", treat = "college", instruments = c("distance", "tuition"), covs = c("score", "female"), degree = dg, nboot = 60)
run3("ivmte logit, no covariates, no plot", ivmte, mt, "pars", outcome = "logwage", treat = "college", instruments = c("distance", "tuition"), link = "logit", plot = FALSE, nboot = 60)
run3("ivmte clusters, notes and guide", ivmte, transform(mt, grp = id %% 100), "curve", outcome = "logwage", treat = "college", instruments = "distance", covs = "score", cluster = "grp", showHelp = TRUE, showGuide = TRUE, nboot = 60)
run3("ivmte factor treatment", ivmte, mt, "tests", outcome = "logwage", treat = "college", instruments = c("distance", "tuition"), covs = "female", nboot = 60)
run3("EXPECT ivmte treatment not binary", ivmte, transform(mt, college = score), "pars", outcome = "logwage", treat = "college", instruments = c("distance", "tuition"), nboot = 60)
run3("EXPECT ivmte variable in two roles", ivmte, mt, "pars", outcome = "logwage", treat = "college", instruments = c("distance", "tuition"), covs = "distance", nboot = 60)
run3("ivmte irrelevant instrument (wide bands or a message)", ivmte, transform(mt, junk = rnorm(nrow(mt))), "pars", outcome = "logwage", treat = "college", instruments = "junk", nboot = 60)
cat(if (bad == 0) "\nAll IV analyses ran.\n" else sprintf("\n%d run(s) failed.\n", bad))
if (bad > 0) quit(status = 1)
